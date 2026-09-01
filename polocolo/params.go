// Package polocolo implements the Polocolo permutation and hash (Ha, Hwang,
// Lee, Park, Son, Eurocrypt 2025, https://eprint.iacr.org/2025/926) as a gnark
// circuit, mirroring the Python reference (ref/polocolo/).
//
// Polocolo is a lookup-based ("type-4") SPN over F_p^t:
//
//	Polocolo_pi = LinLayer^(R) o SBoxLayer o ... o LinLayer^(1) o SBoxLayer o LinLayer^(0)
//
// where LinLayer^(i)(x) = M*x + c^(i) (c^(R) = 0) with a low-addition MDS
// matrix M, and SBoxLayer applies the power-residue S-box
// S(x) = x^{-1} * T[x^((p-1)/m)] to every element. In-circuit, the S-box is
// the shared algebra.PowerResidueSbox gadget (the paper's own Plonk constraint
// system: witness the residue class, never evaluate the ~254-bit annihilator
// exponent), and the linear layer is applied through the published
// straight-line programs (Appendix A.1) via algebra.SLP, realizing the
// low-addition PLONK gate counts a dense product cannot.
package polocolo

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded parameter set of one Polocolo instance,
// mirroring ref polocolo/params.py PolocoloParams. Structural numbers the
// reference fixes per instance (t, m, R, g, the field label, r/c/d) are passed
// in; the derived tables (sigma, the MDS matrix + its straight-line program,
// the round constants) are rebuilt from them.
type Parameters struct {
	P *big.Int // field modulus (a curve scalar field)
	T int      // state size
	M int      // power-residue order, a power of two dividing p-1
	R int      // number of rounds

	G     *big.Int // generator of F_p^*
	Ann   *big.Int // the "annihilator" exponent (p-1)/m; documentation/tests only — never evaluated in-circuit
	Sigma []int    // S-box permutation of {0,...,m-1} (pinned per m, see sigma.go)

	Mat   [][]*big.Int // t x t low-addition MDS matrix (Appendix A.1)
	Slp   algebra.SLP  // the published addition program computing Mat (5/8/13/17/24/31 gates for t=3..8)
	Rcons [][]*big.Int // R x t round constants c^(0)..c^(R-1); c^(R) = 0 is implicit

	// GrMode picks how the S-box gadget obtains g^r: a second lookup table
	// (once-per-circuit cost ~m) or a bit recomposition (~2*log2(m) per S-box).
	// NewParameters selects whichever estimate is cheaper for (m, t, R);
	// overwrite before building circuits to force a mode.
	GrMode algebra.GrMode

	Rate     int
	Capacity int
	Digest   int
}

// mdsMatrices pins the published MDS matrices (Appendix A.1; module constant
// MDS in ref polocolo/params.py). t = 3 is circ(2,1,1) and t = 4 the
// Duval-Leurent M_{4,4}^{8,4} (alpha = 2); t = 5..8 come from the paper's
// randomized low-addition search and exist only as pinned constants (the
// search is unseeded). All entries are small integers, so the same matrices
// serve every field.
var mdsMatrices = map[int][][]int64{
	3: {
		{2, 1, 1},
		{1, 2, 1},
		{1, 1, 2}},
	4: {
		{5, 7, 1, 3},
		{4, 6, 1, 1},
		{1, 3, 5, 7},
		{1, 1, 4, 6}},
	5: {
		{39, 6, 10, 28, 8},
		{174, 28, 32, 80, 16},
		{348, 58, 42, 84, 2},
		{39, 4, 54, 100, 44},
		{204, 20, 300, 560, 244}},
	6: {
		{1011, 1470, 42, 140, 508, 1700},
		{232, 70, 48, 48, 264, 1280},
		{4227, 7371, 3, 490, 1420, 2900},
		{6744, 11760, 60, 844, 2272, 4670},
		{13281, 23163, 9, 1540, 4460, 9100},
		{48, 84, 12, 35, 40, 200}},
	7: {
		{3538, 3090, 768, 480, 720, 96, 336},
		{470862, 470750, 1120, 16380, 94284, 136, 924},
		{10112885, 10113269, 24960, 352496, 2023200, 768, 18048},
		{3799380, 3799524, 9024, 132256, 760128, 288, 6783},
		{94120, 94080, 232, 3276, 18816, 5, 198},
		{1357780, 1357788, 3240, 47268, 271632, 101, 2454},
		{270260, 270260, 640, 9402, 54108, 64, 480}},
	8: {
		{3840, 24, 4728, 2952, 258912, 99840, 94222, 74400},
		{1386, 78, 280, 1218, 32256, 13044, 8120, 6496},
		{6180, 743, 10416, 4428, 508032, 194858, 193984, 153056},
		{432, 400, 1920, 144, 73728, 27776, 30400, 23936},
		{10122, 1246, 5320, 8526, 346752, 136724, 108570, 86184},
		{950, 1052, 5424, 240, 202944, 76333, 84683, 66656},
		{2564, 16, 3072, 1920, 172128, 66380, 62528, 49408},
		{661, 35, 908, 585, 43008, 16448, 14512, 11456}},
}

// mdsSLPs transcribes the published straight-line programs Constraints_3..8
// (Appendix A.1) computing y = M*x: each step is one two-term addition gate
// w = c1*a + c2*b. Index convention: 0..t-1 are the inputs x_1..x_t and t-1+k
// is w_k. NewParameters asserts each program's symbolic matrix equals the
// pinned mdsMatrices entry, so a transcription error cannot ship.
var mdsSLPs = map[int]algebra.SLP{
	3: algebra.NewSLP(3, [][4]int64{
		{1, 0, 1, 1}, // w1 = x1 + x2
		{1, 3, 1, 2}, // w2 = w1 + x3
		{1, 4, 1, 0}, // w3 = w2 + x1
		{1, 4, 1, 1}, // w4 = w2 + x2
		{1, 4, 1, 2}, // w5 = w2 + x3
	}, []int{5, 6, 7}),
	4: algebra.NewSLP(4, [][4]int64{
		{1, 0, 1, 1}, // w1 = x1 + x2
		{1, 2, 1, 3}, // w2 = x3 + x4
		{2, 1, 1, 5}, // w3 = 2*x2 + w2
		{2, 3, 1, 4}, // w4 = 2*x4 + w1
		{4, 5, 1, 7}, // w5 = 4*w2 + w4
		{4, 4, 1, 6}, // w6 = 4*w1 + w3
		{1, 7, 1, 9}, // w7 = w4 + w6
		{1, 6, 1, 8}, // w8 = w3 + w5
	}, []int{10, 9, 11, 8}),
	5: algebra.NewSLP(5, [][4]int64{
		{6, 0, 1, 1},   // w1 = 6*x1 + x2
		{2, 2, 4, 3},   // w2 = 2*x3 + 4*x4
		{3, 0, 8, 4},   // w3 = 3*x1 + 8*x5
		{8, 5, 3, 6},   // w4 = 8*w1 + 3*w2
		{5, 6, 1, 7},   // w5 = 5*w2 + w3
		{4, 2, 5, 9},   // w6 = 4*x3 + 5*w5
		{8, 3, 6, 5},   // w7 = 8*x4 + 6*w1
		{2, 5, 2, 4},   // w8 = 2*w1 + 2*x5
		{1, 11, 1, 9},  // w9 = w7 + w5
		{2, 8, 2, 13},  // w10 = 2*w4 + 2*w9
		{1, 12, 7, 8},  // w11 = w8 + 7*w4
		{1, 10, 2, 12}, // w12 = w6 + 2*w8
		{3, 9, 5, 16},  // w13 = 3*w5 + 5*w12
	}, []int{13, 14, 15, 16, 17}),
	6: algebra.NewSLP(6, [][4]int64{
		{4, 0, 7, 1},   // w1 = 4*x1 + 7*x2
		{2, 2, 2, 3},   // w2 = 2*x3 + 2*x4
		{1, 4, 5, 5},   // w3 = x5 + 5*x6
		{3, 0, 4, 8},   // w4 = 3*x1 + 4*w3
		{5, 6, 4, 4},   // w5 = 5*w1 + 4*x5
		{5, 7, 5, 5},   // w6 = 5*w2 + 5*x6
		{3, 6, 3, 2},   // w7 = 3*w1 + 3*x3
		{8, 9, 3, 7},   // w8 = 8*w4 + 3*w2
		{8, 8, 7, 3},   // w9 = 8*w3 + 7*x4
		{2, 14, 6, 10}, // w10 = 2*w9 + 6*w5
		{1, 9, 7, 15},  // w11 = w4 + 7*w10
		{7, 13, 1, 16}, // w12 = 7*w8 + w11
		{2, 10, 8, 13}, // w13 = 2*w5 + 8*w8
		{5, 16, 1, 12}, // w14 = 5*w11 + w7
		{8, 16, 6, 11}, // w15 = 8*w11 + 6*w6
		{3, 19, 5, 15}, // w16 = 3*w14 + 5*w10
		{5, 14, 4, 12}, // w17 = 5*w9 + 4*w7
	}, []int{17, 18, 19, 20, 21, 22}),
	7: algebra.NewSLP(7, [][4]int64{
		{5, 0, 5, 1},   // w1 = 5*x1 + 5*x2
		{8, 2, 5, 3},   // w2 = 8*x3 + 5*x4
		{3, 4, 2, 5},   // w3 = 3*x5 + 2*x6
		{8, 0, 6, 6},   // w4 = 8*x1 + 6*x7
		{6, 4, 6, 7},   // w5 = 6*x5 + 6*w1
		{2, 8, 2, 11},  // w6 = 2*w2 + 2*w5
		{8, 12, 7, 7},  // w7 = 8*w6 + 7*w1
		{7, 11, 7, 3},  // w8 = 7*w5 + 7*x4
		{7, 10, 6, 9},  // w9 = 7*w4 + 6*w3
		{1, 5, 1, 10},  // w10 = x6 + w4
		{8, 14, 3, 6},  // w11 = 8*w8 + 3*x7
		{8, 17, 4, 8},  // w12 = 8*w11 + 4*w2
		{8, 2, 7, 18},  // w13 = 8*x3 + 7*w12
		{8, 18, 6, 1},  // w14 = 8*w12 + 6*x2
		{1, 18, 2, 7},  // w15 = w12 + 2*w1
		{8, 9, 5, 21},  // w16 = 8*w3 + 5*w15
		{6, 22, 8, 20}, // w17 = 6*w16 + 8*w14
		{8, 15, 6, 13}, // w18 = 8*w9 + 6*w7
		{7, 22, 2, 15}, // w19 = 7*w16 + 2*w9
		{8, 23, 7, 13}, // w20 = 8*w17 + 7*w7
		{5, 17, 3, 23}, // w21 = 5*w11 + 3*w17
		{1, 19, 5, 16}, // w22 = w13 + 5*w10
		{1, 28, 1, 23}, // w23 = w22 + w17
		{4, 22, 6, 14}, // w24 = 4*w16 + 6*w8
	}, []int{24, 25, 26, 27, 28, 29, 30}),
	8: algebra.NewSLP(8, [][4]int64{
		{3, 0, 5, 1},   // w1 = 3*x1 + 5*x2
		{4, 2, 3, 3},   // w2 = 4*x3 + 3*x4
		{8, 4, 3, 5},   // w3 = 8*x5 + 3*x6
		{5, 6, 4, 7},   // w4 = 5*x7 + 4*x8
		{4, 5, 7, 11},  // w5 = 4*x6 + 7*w4
		{5, 9, 1, 8},   // w6 = 5*w2 + w1
		{3, 10, 2, 11}, // w7 = 3*w3 + 2*w4
		{2, 0, 6, 10},  // w8 = 2*x1 + 6*w3
		{4, 2, 6, 14},  // w9 = 4*x3 + 6*w7
		{8, 12, 1, 1},  // w10 = 8*w5 + x2
		{3, 3, 2, 15},  // w11 = 3*x4 + 2*w8
		{6, 16, 6, 6},  // w12 = 6*w9 + 6*x7
		{5, 8, 5, 19},  // w13 = 5*w1 + 5*w12
		{4, 20, 2, 12}, // w14 = 4*w13 + 2*w5
		{1, 6, 6, 20},  // w15 = x7 + 6*w13
		{8, 18, 4, 12}, // w16 = 8*w11 + 4*w5
		{6, 23, 2, 13}, // w17 = 6*w16 + 2*w6
		{2, 17, 5, 5},  // w18 = 2*w10 + 5*x6
		{5, 23, 8, 19}, // w19 = 5*w16 + 8*w12
		{2, 26, 1, 25}, // w20 = 2*w19 + w18
		{3, 9, 8, 6},   // w21 = 3*w2 + 8*x7
		{1, 17, 4, 20}, // w22 = w10 + 4*w13
		{6, 27, 4, 28}, // w23 = 6*w20 + 4*w21
		{2, 30, 1, 19}, // w24 = 2*w23 + w12
		{4, 25, 7, 24}, // w25 = 4*w18 + 7*w17
		{7, 29, 3, 30}, // w26 = 7*w22 + 3*w23
		{6, 23, 4, 21}, // w27 = 6*w16 + 4*w14
		{7, 32, 7, 21}, // w28 = 7*w25 + 7*w14
		{1, 27, 7, 22}, // w29 = w20 + 7*w15
		{2, 15, 8, 27}, // w30 = 2*w8 + 8*w20
		{4, 26, 7, 13}, // w31 = 4*w19 + 7*w6
	}, []int{31, 32, 33, 34, 35, 36, 37, 38}),
}

// NewParameters expands the structural instance numbers into a full parameter
// set, mirroring ref polocolo/params.py PolocoloParams for the official
// instances: sigma comes from the pinned per-m tables, the MDS matrix and its
// addition program from the pinned Appendix A.1 constants, and the round
// constants from a SHAKE128 "bitshift" big-endian stream seeded with
// "Polocolo-{m}-{R}-{t}-{fieldLabel}". fieldLabel is "BLS12" or "BN254" for
// the official fields.
func NewParameters(p *big.Int, t, m, rounds int, g int64, fieldLabel string, rate, capacity, digest int) (*Parameters, error) {
	if t < 3 || t > 8 {
		return nil, fmt.Errorf("polocolo: no published matrix/program for t=%d (official instances use 3 <= t <= 8)", t)
	}
	if m < 2 || m&(m-1) != 0 {
		return nil, fmt.Errorf("polocolo: m=%d must be a power of two >= 2", m)
	}
	pMinus1 := new(big.Int).Sub(p, big.NewInt(1))
	if new(big.Int).Mod(pMinus1, big.NewInt(int64(m))).Sign() != 0 {
		return nil, fmt.Errorf("polocolo: m=%d does not divide p-1", m)
	}
	if rounds < 1 {
		return nil, fmt.Errorf("polocolo: R=%d must be >= 1", rounds)
	}
	if rate+capacity != t || digest > rate {
		return nil, fmt.Errorf("polocolo: need r+c == t and d <= r, got r=%d c=%d d=%d t=%d", rate, capacity, digest, t)
	}
	sigma, ok := sigmaTables[m]
	if !ok {
		return nil, fmt.Errorf("polocolo: no pinned sigma for m=%d (official orders: 32, 64, 128, 512, 1024)", m)
	}

	params := &Parameters{
		P: new(big.Int).Set(p), T: t, M: m, R: rounds,
		G:     big.NewInt(g),
		Ann:   new(big.Int).Div(pMinus1, big.NewInt(int64(m))),
		Sigma: sigma,
		Slp:   mdsSLPs[t],
		Rate:  rate, Capacity: capacity, Digest: digest,
	}

	// The pinned matrix, and the assertion that the transcribed straight-line
	// program computes exactly it.
	params.Mat = make([][]*big.Int, t)
	for i, row := range mdsMatrices[t] {
		params.Mat[i] = make([]*big.Int, t)
		for j, v := range row {
			params.Mat[i][j] = big.NewInt(v)
		}
	}
	slpMat := params.Slp.Matrix(nil) // the pinned matrix is over the integers
	for i := range params.Mat {
		for j := range params.Mat[i] {
			if slpMat[i][j].Cmp(params.Mat[i][j]) != 0 {
				return nil, fmt.Errorf("polocolo: t=%d addition program disagrees with the pinned matrix at [%d][%d]: %s != %s",
					t, i, j, slpMat[i][j], params.Mat[i][j])
			}
		}
	}

	// Round constants: R x t grid off one sequential SHAKE128 stream (bitshift
	// trim, big-endian draws), exactly ref params.py _init_cons. The final
	// linear layer has no constant (c^(R) = 0), so no row is drawn for it.
	seed := fmt.Sprintf("Polocolo-%d-%d-%d-%s", m, rounds, t, fieldLabel)
	xof := sampler.NewXOFSamplerBigEndian([]byte(seed), p, sampler.SHAKE128, sampler.Bitshift)
	params.Rcons = make([][]*big.Int, rounds)
	for r := range params.Rcons {
		params.Rcons[r] = make([]*big.Int, t)
		for i := range params.Rcons[r] {
			params.Rcons[r][i] = xof.Next()
		}
	}

	// g^r mode: measured under gnark's logderiv emulation (see
	// TestGrModeConstraintCost), the bit recomposition wins whenever the G
	// table's once-per-circuit cost (~m) exceeds the S-box count t*R — each
	// bit path costs ~2*log2(m) constraints, but each lookup query carries
	// its own overhead, so the break-even sits near m ≈ t*R rather than the
	// naive m ≈ 2*log2(m)*t*R. Instances can overwrite GrMode.
	if m >= t*rounds {
		params.GrMode = algebra.GrBits
	} else {
		params.GrMode = algebra.GrLookup
	}

	return params, nil
}
