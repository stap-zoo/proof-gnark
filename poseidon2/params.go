// Package poseidon2 implements the Poseidon2 permutation
// (https://eprint.iacr.org/2023/323) as a gnark circuit. Poseidon2 is a
// Hades-strategy design like Poseidon, but with improved linear layers: a
// leading external matrix, a distinct external matrix M_ext (full rounds) and
// internal matrix M_int = J + diag(mat_diag) (partial rounds), and internal-round
// constants placed on the first u branches only. The S-box is the same power map.
// See ref/hades for the reference (Poseidon2Params / the Poseidon2 class).
package poseidon2

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific instantiation of Poseidon2 —
// the single source of truth the permutation reads. Everything below the chosen
// numbers (field, t, alpha, R_ext, R_int, r/c/d, version, mat_diag) is derived by
// NewParameters via the Grain LFSR (round constants) and fixed matrix builders,
// matching hades/params.py Poseidon2Params. A SNARK only ever evaluates the
// permutation forward, so everything here is forward data.
type Parameters struct {
	Field   *big.Int // field characteristic p
	T       int      // state size (branches)
	Alpha   int      // power-map S-box exponent
	RExt    int      // external (full) rounds
	RInt    int      // internal (partial) rounds
	RExtBeg int      // external rounds before the internal block (R_ext/2)
	RExtEnd int      // external rounds after
	R       int      // RExt + RInt
	U       int      // branches the internal S-box (and internal round constants) hit (1)
	Version string   // Grain seeding version: "isec" or "circom"

	Rate     int
	Capacity int
	Digest   int

	MExt           [][]*big.Int // external-round matrix (also the leading matrix)
	MInt           [][]*big.Int // internal-round matrix M_I = J + diag(mat_diag)
	RoundConstants [][]*big.Int // R rows: external rows full width, internal rows on branch 0 only

	// Addition programs for the two linear layers, each verified against its
	// matrix by algebra.SelectSLP and nil where the dense product is already
	// cheaper (t = 2, whose 2x2 layers need two additions either way). They cut
	// only the PLONK column — R1CS folds a constant linear layer into wire
	// coefficients for free — and M_int's matters most: it runs in every one of
	// the ~56 internal rounds.
	SlpExt *algebra.SLP // circ(2,1,...,1) = J+I for t in {2,3}, else the M4 block-circulant program
	SlpInt *algebra.SLP // shared-sum program for M_int = J + diag(mat_diag)
}

// NewParameters builds and validates Poseidon2 parameters for the given field.
// The round constants are drawn from the Grain LFSR (rejection sampling only —
// unlike Poseidon, Poseidon2's matrices are fixed and draw nothing off the
// stream): external rounds take t constants, internal rounds take u constants on
// the first u branches (the rest zero). The matrices are the fixed Poseidon2
// linear layers: M_ext a small circulant (t in {2,3}) or the block-circulant
// built from the Duval-Leurent M4 (t a multiple of 4), and M_int = J +
// diag(matDiag). matDiag (the MAT_DIAG_M_1 vector) is field-specific and supplied
// per instance. version selects the Grain seed layout: "isec" (the reference
// implementation) or "circom".
func NewParameters(field *big.Int, t, alpha, rExt, rInt, rate, capacity, digest int, matDiag []*big.Int, version string) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("poseidon2: field modulus must be a positive integer")
	}
	if t != 2 && t != 3 && t%4 != 0 {
		return nil, fmt.Errorf("poseidon2: state size t must be 2, 3, or a multiple of 4, got %d", t)
	}
	if alpha < 3 {
		return nil, fmt.Errorf("poseidon2: alpha must be >= 3, got %d", alpha)
	}
	if rExt < 2 || rExt%2 != 0 {
		return nil, fmt.Errorf("poseidon2: R_ext must be even and >= 2, got %d", rExt)
	}
	if rInt < 1 {
		return nil, fmt.Errorf("poseidon2: R_int must be >= 1, got %d", rInt)
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("poseidon2: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("poseidon2: digest size must be >= 1, got %d", digest)
	}
	if len(matDiag) != t {
		return nil, fmt.Errorf("poseidon2: mat_diag must have t=%d entries, got %d", t, len(matDiag))
	}
	if version != "isec" && version != "circom" {
		return nil, fmt.Errorf("poseidon2: unknown version %q (want isec or circom)", version)
	}

	p := &Parameters{
		Field: new(big.Int).Set(field),
		T:     t, Alpha: alpha,
		RExt: rExt, RInt: rInt,
		RExtBeg: rExt / 2, RExtEnd: rExt - rExt/2,
		R: rExt + rInt, U: 1, Version: version,
		Rate: rate, Capacity: capacity, Digest: digest,
	}

	// Only the round constants come off the Grain stream; the matrices are fixed.
	g := sampler.NewPoseidonGrain(field, t, rExt, rInt, version)
	p.RoundConstants = sampledRoundConstants(g, p)
	p.MExt = externalMatrix(field, t)
	p.MInt = onesPlusDiag(field, matDiag)
	if err := p.initSLPs(); err != nil {
		return nil, err
	}

	return p, nil
}

// initSLPs derives the addition program for each linear layer and keeps it only
// where it beats the dense product. Both layers have exploitable structure:
// M_int = J + diag(mat_diag) is the shared-sum family, and so is M_ext for
// t in {2,3} since circ(2,1,...,1) = J + I; for t a multiple of 4, M_ext is the
// Duval-Leurent M4 block-circulant, whose published sequence generalizes to any
// number of blocks. algebra.SelectSLP checks each program against the matrix it
// replaces before accepting it.
func (p *Parameters) initSLPs() error {
	var ext algebra.SLP
	if p.T%4 == 0 {
		ext = algebra.BlockCirculantM4SLP(p.T, 2)
	} else {
		var ok bool
		if ext, ok = algebra.JPlusDiagSLP(p.MExt, p.Field); !ok {
			return fmt.Errorf("poseidon2: external matrix for t=%d is not of the form J+diag", p.T)
		}
	}
	slpExt, err := algebra.SelectSLP(ext, p.MExt, p.Field)
	if err != nil {
		return fmt.Errorf("poseidon2: external matrix program: %w", err)
	}

	int_, ok := algebra.JPlusDiagSLP(p.MInt, p.Field)
	if !ok {
		return fmt.Errorf("poseidon2: internal matrix is not of the form J+diag")
	}
	slpInt, err := algebra.SelectSLP(int_, p.MInt, p.Field)
	if err != nil {
		return fmt.Errorf("poseidon2: internal matrix program: %w", err)
	}

	p.SlpExt, p.SlpInt = slpExt, slpInt
	return nil
}

func (p *Parameters) isInternal(r int) bool {
	return p.RExtBeg <= r && r < p.RExtBeg+p.RInt
}

// ---------------------------------------------------------------------------
// Round constants (matches Poseidon2Params._init_rcons)
// ---------------------------------------------------------------------------

// sampledRoundConstants draws the R x t round-constant grid: external rounds draw
// t constants (full width), internal rounds draw u constants placed on the first
// u branches with the rest zero. The draws come off the Grain stream in round
// order (external, internal, external), matching the HorizenLabs reference for
// u=1.
func sampledRoundConstants(g *sampler.Grain, p *Parameters) [][]*big.Int {
	rc := make([][]*big.Int, p.R)
	for i := 0; i < p.R; i++ {
		row := make([]*big.Int, p.T)
		if p.isInternal(i) {
			for j := 0; j < p.T; j++ {
				if j < p.U {
					row[j] = g.Next()
				} else {
					row[j] = big.NewInt(0)
				}
			}
		} else {
			for j := 0; j < p.T; j++ {
				row[j] = g.Next()
			}
		}
		rc[i] = row
	}
	return rc
}

// ---------------------------------------------------------------------------
// Fixed linear layers (matches Poseidon2Params._init_M_ext / _init_M_int)
// ---------------------------------------------------------------------------

// externalMatrix returns the fixed Poseidon2 external matrix M_ext (which is also
// the leading matrix). For t=2 it is circ(2,1); for t=3, circ(2,1,1); for t a
// multiple of 4 it is built from the Duval-Leurent M4 = dl_m44_84(alpha=2):
// t=4 gives 2*M4, and larger multiples give the block-circulant circ(2,1,...,1)
// (x) M4. Entries are reduced mod field. Matches _init_M_ext in ref/hades.
func externalMatrix(field *big.Int, t int) [][]*big.Int {
	switch {
	case t == 2:
		return reduceMatrix(circulant([]int64{2, 1}), field)
	case t == 3:
		return reduceMatrix(circulant([]int64{2, 1, 1}), field)
	case t%4 == 0:
		m4 := dlM44(2)
		if t == 4 {
			scaled := make([][]int64, 4)
			for i := range m4 {
				scaled[i] = make([]int64, 4)
				for j := range m4[i] {
					scaled[i][j] = 2 * m4[i][j]
				}
			}
			return reduceMatrix(scaled, field)
		}
		return reduceMatrix(blockCirculantM4(t, m4), field)
	default:
		panic(fmt.Sprintf("poseidon2: external matrix undefined for t=%d", t))
	}
}

// circulant builds the right-circulant matrix from its first row:
// M[i][j] = row[(j-i) mod n]. Matches utils.matrix.circulant.
func circulant(row []int64) [][]int64 {
	n := len(row)
	m := make([][]int64, n)
	for i := 0; i < n; i++ {
		m[i] = make([]int64, n)
		for j := 0; j < n; j++ {
			m[i][j] = row[((j-i)%n+n)%n]
		}
	}
	return m
}

// dlM44 is the Duval-Leurent M^{8,4}_{4,4} matrix (DL18 Fig. 13), the 4x4 diffusion
// block used by Poseidon2 (and Griffin) for state sizes that are multiples of 4.
// At alpha=2 it is [[5,7,1,3],[4,6,1,1],[1,3,5,7],[1,1,4,6]]. Matches
// utils.matrix.dl_m44_84_matrix.
func dlM44(alpha int64) [][]int64 {
	a, a2 := alpha, alpha*alpha
	return [][]int64{
		{a2 + 1, a2 + a + 1, 1, a + 1},
		{a2, a2 + a, 1, 1},
		{1, a + 1, a2 + 1, a2 + a + 1},
		{1, 1, a2, a2 + a},
	}
}

// blockCirculantM4 builds the t x t block-circulant matrix circ(2,1,...,1) (x) M4
// for t a multiple of 4 (> 4): 4x4 block (I,J) is 2*M4 on the block diagonal and
// M4 off it. Matches utils.matrix.m4_to_block_circulant_matrix.
func blockCirculantM4(t int, m4 [][]int64) [][]int64 {
	m := make([][]int64, t)
	for row := 0; row < t; row++ {
		m[row] = make([]int64, t)
		for col := 0; col < t; col++ {
			val := m4[row%4][col%4]
			if row/4 == col/4 {
				val *= 2
			}
			m[row][col] = val
		}
	}
	return m
}

// onesPlusDiag returns the internal matrix M_int = J + diag(diag): the all-ones
// matrix with diag[i] added on the diagonal, reduced mod field. Matches
// utils.matrix.ones_plus_diag_matrix (Poseidon2's M_I from MAT_DIAG_M_1).
func onesPlusDiag(field *big.Int, diag []*big.Int) [][]*big.Int {
	t := len(diag)
	m := make([][]*big.Int, t)
	for i := 0; i < t; i++ {
		m[i] = make([]*big.Int, t)
		for j := 0; j < t; j++ {
			if i == j {
				v := new(big.Int).Add(big.NewInt(1), diag[i])
				m[i][j] = v.Mod(v, field)
			} else {
				m[i][j] = big.NewInt(1)
			}
		}
	}
	return m
}

// reduceMatrix maps an int64 matrix to field elements reduced mod field.
func reduceMatrix(m [][]int64, field *big.Int) [][]*big.Int {
	out := make([][]*big.Int, len(m))
	for i, row := range m {
		out[i] = make([]*big.Int, len(row))
		for j, v := range row {
			e := big.NewInt(v)
			out[i][j] = e.Mod(e, field)
		}
	}
	return out
}
