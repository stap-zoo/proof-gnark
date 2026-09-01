package griffin

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific Griffin instance the
// permutation reads, mirroring GriffinParams in the reference griffin/params.py.
// The chosen numbers (field, t, R, alpha, r/c/d) are passed to NewParameters; the
// derived tables — the mixing matrix M, the affine round constants, and the
// quadratic-map coefficients — are rebuilt from them.
//
// No inverse data (M_inv, alpha_inv) is stored: the in-circuit permutation only
// evaluates forward, and the single high-degree inverse S-box (x^{1/alpha} on
// branch 0) is handled by algebra.InvPow, whose hint derives the inverse exponent
// itself.
type Parameters struct {
	Field *big.Int // field characteristic p
	T     int      // state size (branches)
	R     int      // number of rounds
	Alpha int      // power-map exponent d: forward x^d on branch 1, inverse x^{1/d} on branch 0

	M       [][]*big.Int // t x t mixing matrix
	Rcons   [][]*big.Int // R rows of t affine round constants (the last row is all zero)
	CoeffsG [][]*big.Int // t-2 [a, b] pairs for the root-free quadratics G_i (round-independent)

	// Slp is the addition program for M, verified against it by
	// algebra.SelectSLP. Which program applies follows the matrix family (initM):
	// at t=3 the matrix is circ(2,1,1) = J + I, so the shared sum s = x0+x1+x2
	// gives every output in 5 additions instead of 6; at t a multiple of 4 it is
	// built from the Duval-Leurent M4, and its published sequence takes 8 additions
	// against 12 dense at t=4 (algebra.BlockCirculantM4SLP's doc comment has the
	// count for wider blocks). Either way a PLONK-only saving, applied in each of
	// the R rounds plus the leading multiplication.
	Slp *algebra.SLP

	Rate     int
	Capacity int
	Digest   int
}

// NewParameters builds and validates Griffin parameters for the given field. R is
// the round count, alpha the power-map exponent, and r/c/d the sponge parameters —
// all taken explicitly per instance, as in the reference instances.py. The mixing
// matrix, round constants, and quadratic-map coefficients are derived
// deterministically from (p, t, R).
func NewParameters(field *big.Int, t, R, alpha, rate, capacity, digest int) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("griffin: field modulus must be a positive integer")
	}
	if field.Bit(0) == 0 {
		return nil, fmt.Errorf("griffin: field modulus must be odd, got an even value")
	}
	if !(t == 3 || t%4 == 0) {
		return nil, fmt.Errorf("griffin: state size t must be 3 or a multiple of 4, got %d", t)
	}
	if R < 1 {
		return nil, fmt.Errorf("griffin: R must be >= 1, got %d", R)
	}
	if alpha != 3 && alpha != 5 && alpha != 7 {
		return nil, fmt.Errorf("griffin: alpha must be 3, 5, or 7, got %d", alpha)
	}
	pMinus1 := new(big.Int).Sub(field, big.NewInt(1))
	if !coprime(alpha, pMinus1) {
		return nil, fmt.Errorf("griffin: x^%d is not a permutation (gcd(alpha, p-1) != 1)", alpha)
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("griffin: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("griffin: digest size must be >= 1, got %d", digest)
	}

	m, err := initM(t)
	if err != nil {
		return nil, err
	}
	rcons, coeffsG, err := deriveConstants(field, t, R)
	if err != nil {
		return nil, fmt.Errorf("griffin: constants: %w", err)
	}

	p := &Parameters{
		Field:    new(big.Int).Set(field),
		T:        t,
		R:        R,
		Alpha:    alpha,
		M:        m,
		Rcons:    rcons,
		CoeffsG:  coeffsG,
		Rate:     rate,
		Capacity: capacity,
		Digest:   digest,
	}

	// One program per matrix family, the same split initM makes: circ(2,1,1) at
	// t=3 is J + I, so the shared-sum program applies; a multiple of 4 is the M4
	// block-circulant, so the Duval-Leurent sequence does (one block at t=4, the
	// blockwise generalisation above it). SelectSLP checks whichever we pick
	// against the matrix it stands in for, so a mismatched pairing fails here
	// rather than silently computing the wrong layer.
	var slp algebra.SLP
	switch {
	case t == 4:
		slp = algebra.DLM44SLP(2, 1) // M itself is M4: scale 1, not Poseidon2's 2*M4
	case t%4 == 0:
		slp = algebra.BlockCirculantM4SLP(t, 2)
	default:
		var ok bool
		if slp, ok = algebra.JPlusDiagSLP(p.M, p.Field); !ok {
			return nil, fmt.Errorf("griffin: mixing matrix for t=%d is not of the form J+diag", t)
		}
	}
	if p.Slp, err = algebra.SelectSLP(slp, p.M, p.Field); err != nil {
		return nil, fmt.Errorf("griffin: mixing matrix program: %w", err)
	}
	return p, nil
}

func coprime(a int, n *big.Int) bool {
	return new(big.Int).GCD(nil, nil, big.NewInt(int64(a)), n).Cmp(big.NewInt(1)) == 0
}

// ---------------------------------------------------------------------------
// Mixing matrix (ref params.py GriffinParams._init_mat), one per admissible width
// family: the pinned circulant of [2, 1, 1] at t=3, and the M4 block-circulant
// construction at t a multiple of 4. Those two families are all Griffin defines
// (eprint 2022/403), which is what NewParameters validates.
// ---------------------------------------------------------------------------

func initM(t int) ([][]*big.Int, error) {
	switch {
	case t == 3:
		// circulant([2, 1, 1]): each row is a cyclic right-shift of [2, 1, 1].
		return circulant([]int64{2, 1, 1}), nil
	case t%4 == 0:
		return m4BlockCirculant(t), nil
	}
	return nil, fmt.Errorf("griffin: no mixing matrix for t=%d, want 3 or a multiple of 4", t)
}

// m4BlockCirculant is ref utils.matrix.m4_to_block_circulant_matrix: the t x t
// matrix circ(2,1,...,1) (x) M4 for t a multiple of 4, where M4 is the
// Duval-Leurent M^{8,4}_{4,4} at alpha=2 (MDS for every p > 2^31). Block (I,J) is
// 2*M4 on the block diagonal and M4 off it — full diffusion in one round at O(t)
// additions, at the price of the whole matrix not being MDS, which is the
// trade Griffin and Poseidon2 both take at these widths.
//
// t=4 is the one block, and there the reference returns M4 *unscaled*. Poseidon2
// doubles it at that width (its own spec says so), so the two constructions do not
// share this builder even though they share the family; what they do share is the
// addition program in algebra, and algebra.SelectSLP re-checks it against whichever
// matrix it is handed.
func m4BlockCirculant(t int) [][]*big.Int {
	// dl_m44_84_matrix(alpha=2), DL18 Fig. 13.
	const a, a2 = 2, 4
	m4 := [4][4]int64{
		{a2 + 1, a2 + a + 1, 1, a + 1},
		{a2, a2 + a, 1, 1},
		{1, a + 1, a2 + 1, a2 + a + 1},
		{1, 1, a2, a2 + a},
	}
	m := make([][]*big.Int, t)
	for row := 0; row < t; row++ {
		m[row] = make([]*big.Int, t)
		for col := 0; col < t; col++ {
			v := m4[row%4][col%4]
			if t > 4 && row/4 == col/4 {
				v *= 2
			}
			m[row][col] = big.NewInt(v)
		}
	}
	return m
}

// circulant builds the right-circulant matrix whose first row is `row`: row i is
// row cyclically right-shifted by i (ref utils.matrix.circulant).
func circulant(row []int64) [][]*big.Int {
	n := len(row)
	m := make([][]*big.Int, n)
	for i := 0; i < n; i++ {
		s := (n - i) % n
		m[i] = make([]*big.Int, n)
		for j := 0; j < n; j++ {
			m[i][j] = big.NewInt(row[(s+j)%n])
		}
	}
	return m
}

// ---------------------------------------------------------------------------
// Constants (ref params.py GriffinParams._init_constants). The (R-1)xt round
// constants and the (t-2) quadratic-map [a, b] pairs are drawn from ONE SHAKE128
// "bitmask" stream, seeded with "Griffin" || p (p serialized little-endian, width
// rounded up to whole 64-bit limbs). The round constants are drawn first, then the
// coefficients, so both share and advance the same stream.
// ---------------------------------------------------------------------------

func deriveConstants(field *big.Int, t, R int) (rcons, coeffsG [][]*big.Int, err error) {
	// seed = b"Griffin" + p.to_bytes(ceil(bitlen/64)*8, "little")
	nSeed := ((field.BitLen() + 63) / 64) * 8
	seed := append([]byte("Griffin"), leBytes(field, nSeed)...)
	r := sampler.NewXOFSampler(seed, field, sampler.SHAKE128, sampler.Bitmask)

	// Round constants: (R-1) rows of t, then a zero row so the round loop can index
	// rcons[r] uniformly for r in 0..R-1 (the final round has no constants).
	rcons = make([][]*big.Int, 0, R)
	for row := 0; row < R-1; row++ {
		rr := make([]*big.Int, t)
		for c := 0; c < t; c++ {
			rr[c] = r.Next()
		}
		rcons = append(rcons, rr)
	}
	zeroRow := make([]*big.Int, t)
	for c := range zeroRow {
		zeroRow[c] = big.NewInt(0)
	}
	rcons = append(rcons, zeroRow)

	coeffsG = deriveCoeffsG(r, field, t)
	return rcons, coeffsG, nil
}

// deriveCoeffsG draws the (t-2) [a, b] pairs for the quadratics G_i(l) = l^2 +
// a*l + b off the shared stream, exactly as the reference _init_constants:
//
//   - the first pair: draw two distinct nonzero elements a, b and keep them iff
//     a^2 - 4*b is a quadratic non-residue (so G_2 is root-free); otherwise redraw
//     the whole pair.
//   - the remaining pairs (i = 2 .. t-2, i.e. for G_3 .. G_{t-1}): a_i = i*a and
//     b_i = i^2*b derived from the first pair, resampling b_i off the stream only
//     if it collides with a_i.
func deriveCoeffsG(r *sampler.XOFSampler, field *big.Int, t int) [][]*big.Int {
	coeffs := make([][]*big.Int, 0, t-2)

	var a0, b0 *big.Int
	for {
		a := r.NextNonzero()
		b := r.NextNonzero()
		for a.Cmp(b) == 0 {
			b = r.NextNonzero()
		}
		disc := new(big.Int).Mul(a, a)         // a^2
		disc.Sub(disc, new(big.Int).Lsh(b, 2)) // - 4*b
		disc.Mod(disc, field)
		if sampler.IsQuadraticNonResidue(disc, field) {
			a0, b0 = a, b
			coeffs = append(coeffs, []*big.Int{a, b})
			break
		}
	}

	for i := 2; i <= t-2; i++ { // ref range(2, t-1) => i = 2 .. t-2
		bi := new(big.Int).Mul(b0, big.NewInt(int64(i*i)))
		ai := new(big.Int).Mul(a0, big.NewInt(int64(i)))
		ai.Mod(ai, field)
		bi.Mod(bi, field)
		for ai.Cmp(bi) == 0 {
			bi = r.NextNonzero()
		}
		coeffs = append(coeffs, []*big.Int{ai, bi})
	}
	return coeffs
}

// leBytes returns x as an n-byte little-endian slice, zero-padded in the high
// bytes (matching Python's int.to_bytes(n, "little")). Panics if x does not fit.
func leBytes(x *big.Int, n int) []byte {
	be := x.Bytes() // big-endian, minimal length
	if len(be) > n {
		panic(fmt.Sprintf("griffin: value needs %d bytes, seed width is %d", len(be), n))
	}
	b := make([]byte, n)
	for i := 0; i < len(be); i++ {
		b[i] = be[len(be)-1-i]
	}
	return b
}
