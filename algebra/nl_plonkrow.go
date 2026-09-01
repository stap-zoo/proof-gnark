// This file holds the gadgets whose only reason to exist is that gnark's frontend
// leaves PLONK row selectors empty. A row of the sparse constraint system is
//
//	qL*a + qR*b + qM*a*b + qO*o + qC = 0
//
// over three wires, but api.Mul emits qM alone and api.Add emits qL/qR/qC alone, so
// anything mixing a product with a linear term — or a product with the assertion
// that consumes it — pays two rows for what one row can hold. Each gadget here is
// that same row, written once, with the cost pinned in nl_plonkrow_test.go.
//
// None of them changes what a circuit computes, and none of them changes an R1CS
// count: there a constant shift is a free linear combination either way. They are
// PLONK-only, which is why each tests for frontend.PlonkAPI — the interface only the
// sparse-R1CS builder implements — and falls back to the plain spelling otherwise.
//
// What the pass over the type-2 designs was worth, per construction, against the
// plain spelling: −16…12% Anemoi, −14% Arion, −13% Griffin, −9% Rescue-Prime,
// −10…2% Polocolo, −2…1% Neptune, −0.3% Skyscraper. Nothing in Poseidon, Poseidon2
// or GMiMC, whose rows are all either a bare power map or a pure linear combination,
// and nothing in R1CS on any of the 48 permutation targets. The absolute counts
// these produced are pinned in registry/cost_test.go, which is what a regression
// would trip; the percentages are here because they are the claim about the
// arithmetization rather than about any one instance.

package algebra

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
)

// MulAffine returns u * (v + k) for a constant k — equivalently u*v + k*u — in a
// SINGLE PLONK gate, where the obvious spelling costs two. In R1CS it is one
// constraint either way.
//
// A PLONK row is qL*a + qR*b + qM*a*b + qO*o + qC = 0 over three wires, so one row
// can carry a product AND a linear term in the same wires: with a = u, b = v,
// o = result, qM = 1, qL = k, qO = -1 it evaluates u*v + k*u. gnark's scs builder
// never fills qL/qR on a multiplication row of its own accord — api.Mul emits qM
// alone — so api.Mul(u, api.Add(v, k)) pays a second row just for the shift.
// api.MulAcc's fast path does fill it, whenever its accumulator shares a wire with
// one of its factors; handing it k*u (free, since a constant only scales a term's
// coefficient) is what lands the whole thing on one row.
//
// Setting v = u gives the quadratic form u^2 + k*u, and v = q*u gives q*u^2 + k*u.
// That is the shape behind every quadratic-plus-linear in this repository: Griffin's
// and Arion's root-free G/g/h on a running sum, and Polocolo's g^r bit
// recomposition.
//
// # The constant term does not fit, and must be pushed onto the consumer
//
// A quadratic's CONSTANT term cannot ride this row, so x^2 + k*x + d is two gates,
// not one. The reason is worth stating exactly, because it is a limit of gnark's
// frontend and not of PLONK — the underlying sparseR1C holds qC as a full field
// element, and api.Add(t1, t2, d) fills it with one:
//
//	api.Mul(u, v)            qM only               (qL, qR, qC forced to 0)
//	api.Add(t1, t2, d)       qL, qR, qC            field-sized, but qM is 0
//	api.MulAcc fast path     qM, qL                field-sized, qC HARDCODED 0
//	EvaluatePlonkExpression  qM, qL, qR, qC        but every selector typed int
//
// No exported path sets qM together with a field-sized qC. The asymmetry between qL
// and qC in the MulAcc row is the crux: that qL *is* a wire's term coefficient
// (qL: aVar.Coeff), so scaling u by k puts a 254-bit value there for free, whereas
// qC belongs to no wire and there is nothing to scale. A round constant does not fit
// in an int.
//
// So d has to move onto whatever CONSUMES the result. Both cases occur here, and
// they work differently:
//
//   - A product consumes it: distribute, and d stops being a constant at all. It
//     becomes a wire coefficient, which is field-sized. Griffin's
//     x_i * (l^2 + a*l + b) is MulAffine(x_i, MulAffine(l, l, a), b) — the outer row
//     computes x_i*core + b*x_i, so b lands in the OUTER row's qL. Two rows.
//   - A sum consumes it: just merge, since addition rows take a field-sized qC
//     already. Anemoi's t = x - (beta*y^2 + gamma) is
//     api.Sub(x, api.Mul(beta, Pow(api, y, 2)), gamma) — one row that was going to
//     exist anyway, with gamma in its empty qC (the same trick MatVecMulConst plays
//     with round constants). Nothing is converted.
//
// What is NOT lost either way: at both shapes the pushed-constant form ties the row
// count an unrestricted qC would buy — Griffin is {l^2+a*l, x_i*(core+b)} = 2 rows
// against {l^2+a*l+b, x_i*g} = 2, and Arion is 3 either way. gnark's int-typed qC
// therefore costs this repository nothing; it would only bite a quadratic consumed by
// neither a product nor a sum, which none of the eleven constructions has.
//
// What IS a mistake is api.Add(MulAffine(...), d). That is the second gate, paid for
// nothing. Do not write it.
//
// # Why this is gated on the builder
//
// The fused form is used only on the sparse-R1CS builder — which is also the only
// one implementing frontend.PlonkAPI, so that interface is the test. R1CS gets the
// plain form, where a constant shift is a free linear combination and the count is
// identical. That split is not only cosmetic: r1cs's MulAcc mutates its accumulator
// argument in place when capacity allows, and keeping R1CS out of this path keeps
// that aliasing hazard out of reach.
func MulAffine(api frontend.API, u, v frontend.Variable, k *big.Int) frontend.Variable {
	if k == nil || k.Sign() == 0 {
		return api.Mul(u, v)
	}
	if _, sparse := api.(frontend.PlonkAPI); sparse {
		return api.MulAcc(api.Mul(k, u), u, v)
	}
	return api.Mul(u, api.Add(v, k))
}

// AssertProductIsEqual asserts u*v == x in a SINGLE PLONK gate, where materialising
// the product and then comparing it costs two. In R1CS it is one constraint either
// way.
//
// A verification row is the cheapest thing in a type-2 design and the easiest to
// overpay for. api.Mul(u, v) spends a row whose qO points at a fresh internal wire,
// and api.AssertIsEqual then spends a second row equating that wire with x — but the
// first row's qO can point at x directly, which is qM*u*v + qO*x = 0 with qO = -1.
// Nothing here needs a field-sized selector, so unlike MulAffine's constant term
// this fits frontend.PlonkAPI's int-typed interface exactly.
//
// This is what pins every hinted value in the repository: InvPow's forward check
// w^alpha == x (see AssertPowIsEqual) and the power-residue S-box's membership check
// wm*g^r == x.
func AssertProductIsEqual(api frontend.API, u, v, x frontend.Variable) {
	if sparse, ok := api.(frontend.PlonkAPI); ok {
		sparse.AddPlonkConstraint(u, v, x, 0, 0, -1, 1, 0) // u*v - x == 0
		return
	}
	api.AssertIsEqual(api.Mul(u, v), x)
}

// AssertPowIsEqual asserts w^degree == x, stopping the square-and-multiply chain one
// multiplication short so that the last multiplication and the equality share a row
// (AssertProductIsEqual). One PLONK gate cheaper than Pow followed by
// AssertIsEqual, for every degree above 1; R1CS is unchanged.
func AssertPowIsEqual(api frontend.API, w frontend.Variable, degree int, x frontend.Variable) {
	if degree < 0 {
		panic("algebra: negative degree")
	}
	if degree < 2 {
		// No multiplication row to share: w^0 is the constant 1, w^1 is w itself.
		api.AssertIsEqual(Pow(api, w, degree), x)
		return
	}
	if degree%2 == 0 {
		u := Pow(api, w, degree/2)
		AssertProductIsEqual(api, u, u, x) // u^2 == x
		return
	}
	AssertProductIsEqual(api, Pow(api, w, degree-1), w, x) // w^(degree-1) * w == x
}
