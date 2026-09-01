package algebra

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
)

// SF is the Lai-Massey lifting S_F(x, y) = (a*x + s, a*y + s) with s = b*(x - y)^2,
// the quadratic round function shared by the open-Flystel / Lai-Massey S-boxes.
// It costs one multiplication constraint (the square); a and b are constants and
// fold in. In PLONK it is four rows: the difference, the square, and one per
// output.
func SF(api frontend.API, x, y frontend.Variable, a, b *big.Int) (frontend.Variable, frontend.Variable) {
	return sfShift(api, x, y, a, b, nil)
}

// sfShift is SF with the constant kx added to its FIRST output. The constant rides
// the qC selector of the row that already sums a*x + s, so it is free — which is
// what lets LaiMassey below spend no gate on either of its two constant offsets.
func sfShift(api frontend.API, x, y frontend.Variable, a, b, kx *big.Int) (frontend.Variable, frontend.Variable) {
	d := api.Sub(x, y)
	s := api.Mul(b, api.Mul(d, d)) // b*(x-y)^2
	// One row for a*x + s, with kx in its qC when there is one. Computing the
	// unshifted sum first and then rebuilding it with kx would emit a row for a
	// value nothing reads — gnark's frontend does not eliminate it — so the branch
	// has to pick the spelling before anything is emitted.
	ax := api.Mul(a, x)
	var u frontend.Variable
	if kx != nil && kx.Sign() != 0 {
		u = api.Add(ax, s, kx)
	} else {
		u = api.Add(ax, s)
	}
	return u, api.Add(api.Mul(a, y), s)
}

// LaiMassey is Neptune's external quadratic pair-wise S-box (open-Flystel map),
// Eq. (22) of https://eprint.iacr.org/2021/1695.pdf:
//
//	(x, y) = SF(x, y)
//	(x, y) = M · (x, y)          // 2x2 constant matrix
//	(x, y) = (x + gamma, y)
//	(x, y) = SF(x, y)
//	return   (x - a*gamma, y)
//
// It costs two constraints (two SF squares); the 2x2 matrix and the constant
// offsets fold in. a, b, gamma are constants and M is the 2x2 constant matrix.
//
// Neither constant offset costs a PLONK gate either, because each sits directly
// after a row that has a free constant selector and so rides it: +gamma is folded
// into the matrix outputs, and -a*gamma follows the second SF, whose first output
// row absorbs it (sfShift). Written literally — an Add and a Sub of their own —
// they cost two gates per S-box, R times per permutation.
//
// In PLONK it is nine rows: five for the first SF fused with the matrix
// (sfMatrix), four for the second SF. R1CS is two constraints, the two squares.
//
// Nine is the floor for this decomposition, and reaching it is worth −24 gates per
// Neptune permutation at t=4 (782 → 758) and −16 at t=2 (386 → 370), R1CS untouched
// at 232 and 186. Those are 12 and 8 S-box applications at 2 gates each; the
// per-application saving is what nl_laimassey_test.go pins, because the totals move
// whenever anything else in the permutation does — the t=4 pair read 844 → 820
// before the internal rounds went through algebra.PartialRounds.
func LaiMassey(api frontend.API, x, y frontend.Variable, a, b, gamma *big.Int, m [][]*big.Int) (frontend.Variable, frontend.Variable) {
	p, q := sfMatrix(api, x, y, a, b, gamma, m)
	aGamma := new(big.Int).Neg(new(big.Int).Mul(a, orZero(gamma)))
	return sfShift(api, p, q, a, b, aGamma)
}

// sfMatrix computes the first SF and the 2x2 matrix TOGETHER — M·SF(x, y) +
// (gamma, 0) — in five PLONK rows where the layers written separately cost six.
// R1CS is one constraint either way, the square.
//
// # Why the obvious two spellings both cost six
//
// Materialising SF's outputs costs a row each and leaves the matrix two 2-term
// rows: with the difference and the square, 1+1+2+2 = 6. Substituting them out
// does no better. Writing s = b*(x-y)^2, SF is (a*x + s, a*y + s), so
//
//	p = (m00+m01)*s + a*m00*x + a*m01*y + gamma
//	q = (m10+m11)*s + a*m10*x + a*m11*y
//
// which never materialises u or v, but pays for it with a 3-term row per output:
// 1+1+2+2 = 6 again. Every arrangement of these terms is six rows, because both
// outputs genuinely depend on all three of s, x and y.
//
// # The square's row has a spare linear selector
//
// The row that computes the square can carry a linear term in its own wire for
// free (MulAffine's qL — see its doc comment). So computing
//
//	s' = b*d^2 + sigma*d       instead of      s = b*d^2,      d = x - y
//
// costs the same single row for any constant sigma. Since s = s' - sigma*d =
// s' - sigma*x + sigma*y, feeding s' to the outputs above shifts their x and y
// coefficients in step:
//
//	row_i = (m_i0+m_i1)*s' + (a*m_i0 - sum_i*sigma)*x + (a*m_i1 + sum_i*sigma)*y
//
// with sum_i = m_i0 + m_i1. The shift is one degree of freedom over both rows, so
// it can zero exactly one of the four x/y coefficients — turning that output into
// a 2-term row, which is one gate instead of two. Five rows: 1+1+1+2.
//
// sigma is picked by matrixShift, which evaluates the candidates that zero each
// coefficient and keeps the cheapest; it cannot lose, since sigma = 0 is the
// unfused count. What it buys is fixed by M alone, so the choice is a compile-time
// constant and the pinned saving in nl_laimassey_test.go is exact.
//
// Nothing here is PLONK-only in the sense MulAffine's doc warns about: sigma
// changes only constant coefficients, and in R1CS a shifted square costs the one
// constraint an unshifted one does, with every linear term free. Both backends run
// this spelling, and the reference vectors pass under it because it is an
// identity, not an approximation.
func sfMatrix(api frontend.API, x, y frontend.Variable, a, b, gamma *big.Int, m [][]*big.Int) (frontend.Variable, frontend.Variable) {
	field := api.Compiler().Field()
	sigma, ok := matrixShift(a, b, m, orZero(gamma), field)
	if !ok {
		// b = 0 leaves no square to shift; spell the layers out plainly.
		u, v := SF(api, x, y, a, b)
		out := MatVecMulConst(api, m, []frontend.Variable{u, v}, []*big.Int{gamma, nil})
		return out[0], out[1]
	}

	d := api.Sub(x, y)
	// w = d^2 + (sigma/b)*d in one row, so that b*w is the shifted square s'. b
	// stays in the output rows' coefficients rather than costing a row of its own.
	w := MulAffine(api, d, d, mulMod(sigma, invMod(b, field), field))

	row := func(i int, konst *big.Int) frontend.Variable {
		cw, cx, cy := shiftedCoeffs(a, b, m, sigma, i, field)
		var terms []frontend.Variable
		for _, t := range []struct {
			coeff *big.Int
			v     frontend.Variable
		}{{cw, w}, {cx, x}, {cy, y}} {
			if t.coeff.Sign() != 0 {
				terms = append(terms, api.Mul(t.coeff, t.v))
			}
		}
		return sumWithConstant(api, terms, konst)
	}
	return row(0, gamma), row(1, nil)
}

// shiftedCoeffs returns row i's coefficients on (w, x, y) once the square row
// carries sigma*d, where w = d^2 + (sigma/b)*d so that b*w = b*d^2 + sigma*d:
//
//	cw = (m_i0 + m_i1) * b        cx = a*m_i0 - sum_i*sigma
//	                              cy = a*m_i1 + sum_i*sigma
func shiftedCoeffs(a, b *big.Int, m [][]*big.Int, sigma *big.Int, i int, field *big.Int) (cw, cx, cy *big.Int) {
	sum := addMod(m[i][0], m[i][1], field)
	shift := mulMod(sum, sigma, field)
	cw = mulMod(sum, b, field)
	cx = subMod(mulMod(a, m[i][0], field), shift, field)
	cy = addMod(mulMod(a, m[i][1], field), shift, field)
	return cw, cx, cy
}

// matrixShift chooses the square row's linear coefficient sigma. Each candidate
// zeroes one x or y coefficient of one output row; the winner is whichever
// minimises the two rows' gate count, with sigma = 0 (the unfused spelling) as the
// floor, so the choice is never worse than not shifting at all. ok is false when
// b vanishes in the field, in which case there is no square to shift — tested on
// the reduced value, since b = p is as much a zero here as b = 0.
func matrixShift(a, b *big.Int, m [][]*big.Int, gamma, field *big.Int) (*big.Int, bool) {
	if b == nil || new(big.Int).Mod(b, field).Sign() == 0 {
		return nil, false
	}
	candidates := []*big.Int{big.NewInt(0)}
	for i := 0; i < 2; i++ {
		sum := addMod(m[i][0], m[i][1], field)
		if sum.Sign() == 0 {
			continue // this row does not read the square, so no shift moves it
		}
		inv := invMod(sum, field)
		// zero cx: sigma = a*m_i0 / sum_i;  zero cy: sigma = -a*m_i1 / sum_i
		candidates = append(candidates,
			mulMod(mulMod(a, m[i][0], field), inv, field),
			mulMod(negMod(mulMod(a, m[i][1], field), field), inv, field))
	}

	best, bestCost := candidates[0], -1
	for _, sigma := range candidates {
		cost := 0
		for i := 0; i < 2; i++ {
			cw, cx, cy := shiftedCoeffs(a, b, m, sigma, i, field)
			konst := gamma
			if i == 1 {
				konst = nil
			}
			cost += rowGates([]*big.Int{cw, cx, cy}, konst)
		}
		if bestCost < 0 || cost < bestCost {
			best, bestCost = sigma, cost
		}
	}
	return best, true
}

// rowGates is what sumWithConstant spends on a row whose terms have these
// coefficients: one gate per term past the first, and one for a lone term that
// still has a constant to carry (there is no addition for it to ride).
func rowGates(coeffs []*big.Int, konst *big.Int) int {
	terms := 0
	for _, c := range coeffs {
		if c != nil && c.Sign() != 0 {
			terms++
		}
	}
	switch {
	case terms >= 2:
		return terms - 1
	case terms == 1 && konst != nil && konst.Sign() != 0:
		return 1
	default:
		return 0
	}
}

// ---------------------------------------------------------------------------
// constant arithmetic, all mod the field
// ---------------------------------------------------------------------------

func orZero(z *big.Int) *big.Int {
	if z == nil {
		return new(big.Int)
	}
	return z
}

func addMod(u, v, field *big.Int) *big.Int {
	return new(big.Int).Mod(new(big.Int).Add(u, v), field)
}

func subMod(u, v, field *big.Int) *big.Int {
	return new(big.Int).Mod(new(big.Int).Sub(u, v), field)
}

func mulMod(u, v, field *big.Int) *big.Int {
	return new(big.Int).Mod(new(big.Int).Mul(u, v), field)
}

func negMod(u, field *big.Int) *big.Int {
	return new(big.Int).Mod(new(big.Int).Neg(u), field)
}

func invMod(u, field *big.Int) *big.Int {
	return new(big.Int).ModInverse(u, field)
}
