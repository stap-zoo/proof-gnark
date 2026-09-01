package arion

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific Arion instance the
// permutation reads, mirroring ArionParams in the reference arion/params.py.
// The chosen numbers (field, t, R, alpha1, alpha2, r/c/d) are passed to
// NewParameters; the derived tables — the circulant matrix M, the affine round
// constants, and the GTDS map coefficients — are rebuilt from them.
//
// No inverse data (M_inv, alpha1_inv, alpha2_inv) is stored: the in-circuit
// permutation only evaluates forward, and the single high-degree inverse S-box
// (x^{1/alpha2} on the last branch) is handled by algebra.InvPow, whose hint
// derives the inverse exponent itself.
type Parameters struct {
	Field  *big.Int // field characteristic p
	T      int      // state size (branches)
	R      int      // number of rounds
	Alpha1 int      // small forward power on branches 0..t-2 in the GTDS
	Alpha2 int      // high-degree power on branch t-1 (inverted via a hint)

	M       [][]*big.Int   // t x t circulant MDS matrix
	Rcons   [][]*big.Int   // R rows of t affine round constants
	CoeffsG [][][]*big.Int // R rows of (t-1) [a, b] pairs for the g_i quadratics
	CoeffsH [][]*big.Int   // R rows of (t-1) coefficients c for the h_i quadratics

	// Slp is the addition program for M, verified against it by
	// algebra.SelectSLP, and nil where the dense product is at least as cheap.
	// At t=3 the circulant of [1,2,3] has a shared linear form
	// V = 6*x0 + 2*x1 + 3*x2: every row is a multiple of V off a single
	// coordinate, so evaluating V once (2 gates) and reading each output off it
	// (1 gate each) is 5 gates against 6 dense (see algebra.SharedFormSLP).
	// Only the PLONK column changes.
	Slp *algebra.SLP

	Rate     int
	Capacity int
	Digest   int
}

// NewParameters builds and validates Arion parameters for the given field. R is
// the round count, alpha1/alpha2 the two GTDS power-map exponents, and r/c/d the
// sponge parameters — all taken explicitly per instance, as in the reference
// instances.py. The circulant matrix, round constants, and GTDS coefficients are
// derived deterministically from (p, t, R).
func NewParameters(field *big.Int, t, R, alpha1, alpha2, rate, capacity, digest int) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("arion: field modulus must be a positive integer")
	}
	if field.Bit(0) == 0 {
		return nil, fmt.Errorf("arion: field modulus must be odd, got an even value")
	}
	if t < 2 {
		return nil, fmt.Errorf("arion: t must be >= 2, got %d", t)
	}
	if R < 1 {
		return nil, fmt.Errorf("arion: R must be >= 1, got %d", R)
	}
	pMinus1 := new(big.Int).Sub(field, big.NewInt(1))
	if alpha1 < 2 || !coprime(alpha1, pMinus1) {
		return nil, fmt.Errorf("arion: x^%d is not a permutation (gcd(alpha1, p-1) != 1)", alpha1)
	}
	if alpha2 < 2 || !coprime(alpha2, pMinus1) {
		return nil, fmt.Errorf("arion: x^%d is not a permutation (gcd(alpha2, p-1) != 1)", alpha2)
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("arion: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("arion: digest size must be >= 1, got %d", digest)
	}

	coeffsG, coeffsH, rcons, err := deriveConstants(field, t, R)
	if err != nil {
		return nil, fmt.Errorf("arion: constants: %w", err)
	}

	p := &Parameters{
		Field:    new(big.Int).Set(field),
		T:        t,
		R:        R,
		Alpha1:   alpha1,
		Alpha2:   alpha2,
		M:        circulantMatrix(t),
		Rcons:    rcons,
		CoeffsG:  coeffsG,
		CoeffsH:  coeffsH,
		Rate:     rate,
		Capacity: capacity,
		Digest:   digest,
	}

	// The circulant of [1,2,...,t] has a shared linear form at t=3 and, as the
	// search reports, not at t>=4 (each row would have to agree with a multiple
	// of v on t-1 coordinates, which t rows over-determine). A miss is a cost
	// question, not an error, so the layer simply stays dense.
	if slp, ok := algebra.SharedFormSLP(p.M, p.Field); ok {
		if p.Slp, err = algebra.SelectSLP(slp, p.M, p.Field); err != nil {
			return nil, fmt.Errorf("arion: linear layer program: %w", err)
		}
	}
	return p, nil
}

func coprime(a int, n *big.Int) bool {
	return new(big.Int).GCD(nil, nil, big.NewInt(int64(a)), n).Cmp(big.NewInt(1)) == 0
}

// ---------------------------------------------------------------------------
// Circulant MDS matrix (ref utils.matrix.simple_circulant_matrix): the
// right-circulant matrix of the row [1, 2, ..., t]. Row i is a cyclic right-shift
// of the row, i.e. row[(t-i)%t:] ++ row[:(t-i)%t]. Guaranteed MDS for t in {2,3,4}
// (Arion paper, Remark 4). No MDS search or generator is involved.
// ---------------------------------------------------------------------------

func circulantMatrix(t int) [][]*big.Int {
	row := make([]*big.Int, t)
	for j := 0; j < t; j++ {
		row[j] = big.NewInt(int64(j + 1))
	}
	m := make([][]*big.Int, t)
	for i := 0; i < t; i++ {
		s := (t - i) % t
		m[i] = make([]*big.Int, t)
		for j := 0; j < t; j++ {
			m[i][j] = new(big.Int).Set(row[(s+j)%t])
		}
	}
	return m
}

// ---------------------------------------------------------------------------
// Constants (ref params.py ArionParams._init_constants). rcons, coeffs_h and
// coeffs_g are three independent SHAKE256 "mod" streams, all seeded by the
// instance label "Arion(<p-decimal>,<t>,<R>)" with a per-table suffix. kappa is
// fixed at 128 (the reference default; it does not affect these tables).
// ---------------------------------------------------------------------------

func deriveConstants(field *big.Int, t, R int) (coeffsG [][][]*big.Int, coeffsH, rcons [][]*big.Int, err error) {
	base := fmt.Sprintf("Arion(%s,%d,%d)", field.String(), t, R)
	seed := func(suffix string) []byte { return []byte(base + suffix) }

	rcons = sampler.Grid(seed("aff"), field, R, t)
	coeffsH = sampler.Grid(seed("h"), field, R, t-1)

	coeffsG, err = deriveCoeffsG(seed("g"), field, t, R)
	if err != nil {
		return nil, nil, nil, err
	}
	return coeffsG, coeffsH, rcons, nil
}

// deriveCoeffsG rejection-samples the R*(t-1) [a, b] pairs for the g_i
// quadratics. Candidates are drawn flat (round-major) from one SHAKE256 stream
// as 4*n rows of 2 elements (n = R*(t-1)); a pair is kept iff a^2 - 4*b is a
// quadratic non-residue mod p, so g_i(x) = x^2 + a*x + b is root-free (never
// invertible to a field element). The rejection breaks the row alignment of the
// stream, so the kept pairs are reshaped to R x (t-1) afterwards.
func deriveCoeffsG(seed []byte, field *big.Int, t, R int) ([][][]*big.Int, error) {
	n := R * (t - 1)
	candidates := sampler.Grid(seed, field, 4*n, 2)

	kept := make([][]*big.Int, 0, n)
	for _, ab := range candidates {
		a, b := ab[0], ab[1]
		// disc = a^2 - 4*b mod p
		disc := new(big.Int).Mul(a, a)
		disc.Sub(disc, new(big.Int).Lsh(b, 2)) // -4*b
		disc.Mod(disc, field)
		if sampler.IsQuadraticNonResidue(disc, field) {
			kept = append(kept, []*big.Int{a, b})
			if len(kept) == n {
				break
			}
		}
	}
	if len(kept) < n {
		return nil, fmt.Errorf("not enough valid coeffs_g candidates (%d of %d); widen the candidate pool", len(kept), n)
	}

	// Reshape row-major to [round][branch].
	w := t - 1
	out := make([][][]*big.Int, R)
	for r := 0; r < R; r++ {
		out[r] = kept[r*w : (r+1)*w]
	}
	return out, nil
}
