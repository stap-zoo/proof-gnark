// Package anemoi implements the Anemoi permutation (eprint 2022/840) as a gnark
// circuit. Anemoi is a "type-2" design: its S-box (the open Flystel) is defined
// through a high-degree inverse power map but verified with the cheap forward
// map, which is what makes it constraint-efficient to prove.
//
// File layout mirrors the Python reference framework:
//   - anemoi.go     — the permutation round function                          ≈ hash.py
//   - params.go     — parameter generation and sanity checks                  ≈ params.py
//   - instances.go  — concrete instances + harness/mode wiring                ≈ instances.py
//
// Modes of operation (Jive compression, Hirose sponge) come from package mode.
package anemoi

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/algebra"
)

// Digits of pi used to derive the round constants via an open butterfly
// (ref params.py PI_0 / PI_1).
var (
	pi0, _ = new(big.Int).SetString("1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679", 10)
	pi1, _ = new(big.Int).SetString("8214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196", 10)
)

// Parameters is the fully-expanded, field-specific Anemoi instance the
// permutation reads. Everything below the chosen numbers (field, l, alpha, g, R,
// r/c/d) is derived by NewParameters, matching anemoi/params.py.
//
// No inverse data is stored: the forward permutation's only inverse power map
// lives inside algebra.InvPow, whose hint derives the inverse exponent itself.
type Parameters struct {
	Field *big.Int // field characteristic p
	L     int      // number of Flystel columns
	T     int      // state size t = 2*l
	R     int      // rounds
	Alpha int      // Flystel power-map exponent
	Quad  int      // Q-map exponent (2 in odd characteristic)

	Beta  *big.Int // Q-map quadratic coefficient (= g)
	Gamma *big.Int // Q_gamma constant (= 0)
	Delta *big.Int // Q_delta constant (= g^{-1})

	Mx [][]*big.Int // l x l MDS matrix (x-lane)
	My [][]*big.Int // l x l matrix (y-lane): Mx with each row rotated right by one

	// MxAlpha is the power of g the MDS search settled on — the `a` the DL18
	// shape is instantiated at (g itself for every in-scope instance).
	MxAlpha *big.Int

	// SlpMx / SlpMy are the addition programs for the two lanes, verified against
	// their matrices by algebra.SelectSLP. Only l=3 has anything to gain: its
	// M^{5,2}_{3,3} shape takes 5 additions against 6 dense. At l=1 M_x is the
	// identity (already free, all diffusion comes from the PHT) and at l=2 the
	// published pht_apply needs the same two additions the dense product does —
	// the DL18 sequences save *multiplications*, which PLONK gives away free, so
	// both are left dense and both fields are nil there.
	SlpMx *algebra.SLP
	SlpMy *algebra.SLP

	C [][]*big.Int // round constants, x-lane: R rows of l
	D [][]*big.Int // round constants, y-lane: R rows of l

	// LinearConstants[r] is L(C[r] || D[r]) — the round constants pushed through
	// the linear layer L, indexed like the state (x-lane then y-lane). The circuit
	// adds these *after* the layer instead of adding C/D before it, which is the
	// same map (L(s + c) = L(s) + L(c)) but lets every constant ride on the gate
	// that produces its slot instead of costing a PLONK gate. Derivation only —
	// C and D remain the reference's constants.
	LinearConstants [][]*big.Int

	Rate     int
	Capacity int
	Digest   int
}

// NewParameters builds and validates Anemoi parameters for the given field.
// g is a multiplicative generator (beta = g, delta = g^{-1}); R is the round
// count and r/c/d the sponge parameters (all taken explicitly per instance, as
// in the reference instances.py). The MDS matrix and round constants are derived.
func NewParameters(field *big.Int, l, alpha int, g *big.Int, R, rate, capacity, digest int) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("anemoi: field modulus must be a positive integer")
	}
	if l < 1 {
		return nil, fmt.Errorf("anemoi: l must be >= 1, got %d", l)
	}
	if alpha < 3 {
		return nil, fmt.Errorf("anemoi: alpha must be >= 3, got %d", alpha)
	}
	if new(big.Int).GCD(nil, nil, big.NewInt(int64(alpha)), new(big.Int).Sub(field, big.NewInt(1))).Cmp(big.NewInt(1)) != 0 {
		return nil, fmt.Errorf("anemoi: x^%d is not a permutation (gcd(alpha, p-1) != 1)", alpha)
	}
	if g == nil || g.Sign() <= 0 {
		return nil, fmt.Errorf("anemoi: generator g must be a positive integer")
	}
	t := 2 * l
	if rate+capacity != t {
		return nil, fmt.Errorf("anemoi: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if R < 1 {
		return nil, fmt.Errorf("anemoi: R must be >= 1, got %d", R)
	}
	if digest < 1 {
		return nil, fmt.Errorf("anemoi: digest size must be >= 1, got %d", digest)
	}

	gModP := new(big.Int).Mod(g, field)
	delta := new(big.Int).ModInverse(gModP, field)
	if delta == nil {
		return nil, fmt.Errorf("anemoi: g=%s is not invertible mod p", g)
	}

	mx, mxAlpha, err := buildMx(l, gModP, field)
	if err != nil {
		return nil, err
	}

	p := &Parameters{
		Field: new(big.Int).Set(field),
		L:     l, T: t, R: R, Alpha: alpha, Quad: 2,
		Beta:    new(big.Int).Set(gModP),
		Gamma:   big.NewInt(0),
		Delta:   delta,
		Mx:      mx,
		My:      rotateRowsRight(mx),
		MxAlpha: mxAlpha,
		Rate:    rate, Capacity: capacity, Digest: digest,
	}
	if err := p.initSLPs(); err != nil {
		return nil, err
	}
	p.C, p.D = buildRoundConstants(field, l, R, alpha, gModP, delta)
	p.LinearConstants = p.pushConstantsThroughLinearLayer()
	return p, nil
}

// pushConstantsThroughLinearLayer evaluates L(C[r] || D[r]) natively for every
// round, L being the circuit's linear layer: M_x on the x-lane, M_y on the y-lane,
// then the pseudo-Hadamard transform x' = 2*xm + ym, y' = xm + ym. Mirrors
// Permutation.linearLayer exactly — if one changes, so must the other, and the
// reference vectors are what catches a mismatch.
func (p *Parameters) pushConstantsThroughLinearLayer() [][]*big.Int {
	out := make([][]*big.Int, p.R)
	for r := 0; r < p.R; r++ {
		xm := matVecMod(p.Mx, p.C[r], p.Field)
		ym := matVecMod(p.My, p.D[r], p.Field)
		row := make([]*big.Int, p.T)
		for i := 0; i < p.L; i++ {
			x2 := new(big.Int).Lsh(xm[i], 1)
			row[i] = new(big.Int).Mod(new(big.Int).Add(x2, ym[i]), p.Field)
			row[p.L+i] = new(big.Int).Mod(new(big.Int).Add(xm[i], ym[i]), p.Field)
		}
		out[r] = row
	}
	return out
}

// matVecMod is the native m·v mod p.
func matVecMod(m [][]*big.Int, v []*big.Int, p *big.Int) []*big.Int {
	out := make([]*big.Int, len(m))
	for i, row := range m {
		acc := new(big.Int)
		for j, c := range row {
			acc.Add(acc, new(big.Int).Mul(c, v[j]))
		}
		out[i] = acc.Mod(acc, p)
	}
	return out
}

// initSLPs derives the lane programs for l=3, the only width where the published
// DL18 sequence costs fewer PLONK additions than the dense product. M_y is M_x
// with every row rotated right, i.e. M_x applied to the left-rotated input, so
// the same program serves both lanes with its inputs permuted — and both are
// still checked against their own matrix.
func (p *Parameters) initSLPs() error {
	if p.L != 3 {
		return nil
	}
	mx := algebra.DLM3352SLP(p.MxAlpha, p.Field)
	slpMx, err := algebra.SelectSLP(mx, p.Mx, p.Field)
	if err != nil {
		return fmt.Errorf("anemoi: x-lane matrix program: %w", err)
	}
	rotl := make([]int, p.L) // input i of the program reads state slot i+1
	for i := range rotl {
		rotl[i] = (i + 1) % p.L
	}
	slpMy, err := algebra.SelectSLP(mx.PermuteInputs(rotl), p.My, p.Field)
	if err != nil {
		return fmt.Errorf("anemoi: y-lane matrix program: %w", err)
	}
	p.SlpMx, p.SlpMy = slpMx, slpMy
	return nil
}

// ---------------------------------------------------------------------------
// MDS matrix (ref params.py _init_M + utils.matrix)
// ---------------------------------------------------------------------------

// buildMx returns Anemoi's l x l matrix M_x and the power of g it is
// instantiated at: the identity for l=1 (diffusion then comes entirely from the
// pseudo-Hadamard transform), otherwise the DL18 low-addition shape at the
// smallest power of g that is MDS. The exponent is returned because the addition
// program for the shape is parameterized by it.
func buildMx(l int, g, p *big.Int) ([][]*big.Int, *big.Int, error) {
	if l == 1 {
		return [][]*big.Int{{big.NewInt(1)}}, big.NewInt(1), nil
	}
	var shape func(a, p *big.Int) [][]*big.Int
	switch l {
	case 2:
		shape = phtMatrix
	case 3:
		shape = dlM3352Matrix
	default:
		return nil, nil, fmt.Errorf("anemoi: no MDS shape for l=%d, want l in {1, 2, 3}", l)
	}
	const maxTries = 1000
	gi := new(big.Int).Set(g)
	for try := 0; try < maxTries; try++ {
		m := shape(gi, p)
		if isMDS(m, p) {
			return m, new(big.Int).Set(gi), nil
		}
		gi.Mod(gi.Mul(gi, g), p)
	}
	return nil, nil, fmt.Errorf("anemoi: no MDS instance of the l=%d shape within %d powers of g", l, maxTries)
}

// phtMatrix is the generalized pseudo-Hadamard transform [[1, a], [a, 1+a^2]]
// (ref utils.matrix.pht_matrix), Anemoi's M_2 with a = g^i.
func phtMatrix(a, p *big.Int) [][]*big.Int {
	a2 := new(big.Int).Mod(new(big.Int).Mul(a, a), p)
	return [][]*big.Int{
		{big.NewInt(1), new(big.Int).Set(a)},
		{new(big.Int).Set(a), addMod(big.NewInt(1), a2, p)},
	}
}

// dlM3352Matrix is the DL18 M^{5,2}_{3,3} shape (ref utils.matrix.dl_m33_52_matrix),
// Anemoi's M_3 with a = g^i:
//
//	[[1+a, 1, 1+a], [1, 1, a], [a, 1, 1]]
func dlM3352Matrix(a, p *big.Int) [][]*big.Int {
	one := big.NewInt(1)
	onePlusA := addMod(one, a, p)
	return [][]*big.Int{
		{new(big.Int).Set(onePlusA), big.NewInt(1), new(big.Int).Set(onePlusA)},
		{big.NewInt(1), big.NewInt(1), new(big.Int).Set(a)},
		{new(big.Int).Set(a), big.NewInt(1), big.NewInt(1)},
	}
}

// rotateRowsRight returns M_y = M_x with each row rotated right by one, i.e.
// M_x @ P_rho (ref params.py _init_My_from_Mx).
func rotateRowsRight(m [][]*big.Int) [][]*big.Int {
	out := make([][]*big.Int, len(m))
	for i, row := range m {
		n := len(row)
		nr := make([]*big.Int, n)
		nr[0] = new(big.Int).Set(row[n-1])
		for j := 1; j < n; j++ {
			nr[j] = new(big.Int).Set(row[j-1])
		}
		out[i] = nr
	}
	return out
}

// isMDS reports whether every square submatrix of m has a nonzero determinant
// mod p, i.e. all minors of every order are nonzero (ref utils.matrix.is_mds).
func isMDS(m [][]*big.Int, p *big.Int) bool {
	n := len(m)
	for k := 1; k <= n; k++ {
		combos := combinations(n, k)
		for _, rows := range combos {
			for _, cols := range combos {
				if detMod(submatrix(m, rows, cols), p).Sign() == 0 {
					return false
				}
			}
		}
	}
	return true
}

// detMod is the determinant of a small matrix mod p, by cofactor expansion.
func detMod(m [][]*big.Int, p *big.Int) *big.Int {
	n := len(m)
	if n == 1 {
		return new(big.Int).Mod(m[0][0], p)
	}
	det := big.NewInt(0)
	for j := 0; j < n; j++ {
		minor := make([][]*big.Int, n-1)
		for r := 1; r < n; r++ {
			row := make([]*big.Int, 0, n-1)
			for c := 0; c < n; c++ {
				if c != j {
					row = append(row, m[r][c])
				}
			}
			minor[r-1] = row
		}
		term := new(big.Int).Mul(m[0][j], detMod(minor, p))
		if j%2 == 1 {
			term.Neg(term)
		}
		det.Add(det, term)
	}
	return det.Mod(det, p)
}

func submatrix(m [][]*big.Int, rows, cols []int) [][]*big.Int {
	out := make([][]*big.Int, len(rows))
	for i, r := range rows {
		row := make([]*big.Int, len(cols))
		for j, c := range cols {
			row[j] = m[r][c]
		}
		out[i] = row
	}
	return out
}

// combinations returns all k-subsets of {0,..,n-1} in lexicographic order.
func combinations(n, k int) [][]int {
	var out [][]int
	idx := make([]int, k)
	var rec func(start, depth int)
	rec = func(start, depth int) {
		if depth == k {
			c := make([]int, k)
			copy(c, idx)
			out = append(out, c)
			return
		}
		for i := start; i < n; i++ {
			idx[depth] = i
			rec(i+1, depth+1)
		}
	}
	rec(0, 0)
	return out
}

// ---------------------------------------------------------------------------
// Round constants (ref params.py _init_rcons): built from the digits of pi via
// an open butterfly. gamma = 0 is baked in (it does not appear in C).
// ---------------------------------------------------------------------------

func buildRoundConstants(p *big.Int, l, R, alpha int, beta, delta *big.Int) (C, D [][]*big.Int) {
	piF0 := new(big.Int).Mod(pi0, p)
	piF1 := new(big.Int).Mod(pi1, p)
	e := big.NewInt(int64(alpha))
	two := big.NewInt(2)

	C = make([][]*big.Int, R)
	D = make([][]*big.Int, R)
	for r := 0; r < R; r++ {
		pi0r := new(big.Int).Exp(piF0, big.NewInt(int64(r)), p)
		C[r] = make([]*big.Int, l)
		D[r] = make([]*big.Int, l)
		for i := 0; i < l; i++ {
			pi1i := new(big.Int).Exp(piF1, big.NewInt(int64(i)), p)
			powAlpha := new(big.Int).Exp(addMod(pi0r, pi1i, p), e, p)

			// C = beta*pi0r^2 + pow_alpha
			cr := new(big.Int).Mul(beta, new(big.Int).Exp(pi0r, two, p))
			cr.Add(cr, powAlpha)
			C[r][i] = cr.Mod(cr, p)

			// D = beta*pi1i^2 + pow_alpha + delta
			dr := new(big.Int).Mul(beta, new(big.Int).Exp(pi1i, two, p))
			dr.Add(dr, powAlpha)
			dr.Add(dr, delta)
			D[r][i] = dr.Mod(dr, p)
		}
	}
	return C, D
}

func addMod(a, b, p *big.Int) *big.Int {
	return new(big.Int).Mod(new(big.Int).Add(a, b), p)
}
