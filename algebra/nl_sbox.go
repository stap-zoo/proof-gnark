// Package algebra holds the in-circuit arithmetic building blocks that hash
// permutations are assembled from — power-map S-boxes and linear (matrix) layers.
// Almost every construction in this repository reuses these, so they live here
// once rather than being re-derived per hash.
package algebra

import "github.com/consensys/gnark/frontend"

// Pow returns x^degree using binary square-and-multiply, which minimises the
// number of multiplication constraints. Multiplications by the constant 1 fold
// away in gnark, so e.g. degree 5 costs exactly 3 constraints (two squarings and
// one multiply). degree 0 returns the constant 1; a negative degree panics
// (inversion is not a plain power map).
func Pow(api frontend.API, x frontend.Variable, degree int) frontend.Variable {
	if degree < 0 {
		panic("algebra: negative degree")
	}
	var result frontend.Variable = 1
	base := x
	for d := degree; d > 0; d >>= 1 {
		if d&1 == 1 {
			result = api.Mul(result, base)
		}
		if d > 1 {
			base = api.Mul(base, base)
		}
	}
	return result
}

// PowPlusConst returns x^degree + k for a SMALL INTEGER k, at exactly the cost of
// x^degree alone: the constant rides the qC selector of the row that performs the
// last multiplication of the square-and-multiply chain, which api.Mul leaves empty.
// Adding k afterwards would cost a PLONK gate (it is free in R1CS, where a shift is
// just a linear combination — that path is taken verbatim).
//
// k must be an int, not a field element, because gnark exposes qC only through
// frontend.PlonkAPI's int-typed selectors — see algebra.MulAffine, which documents
// that limitation and the way around it for field-sized constants. So this is for
// the small shifts: Arion's GTDS needs x^alpha1 + 1.
func PowPlusConst(api frontend.API, x frontend.Variable, degree, k int) frontend.Variable {
	if degree < 0 {
		panic("algebra: negative degree")
	}
	sparse, isSparse := api.(frontend.PlonkAPI)
	if degree == 0 || degree == 1 || !isSparse {
		// Nothing to fuse into (or nothing to gain): x^0 and x^1 have no
		// multiplication row, and in R1CS the addition is free.
		return api.Add(Pow(api, x, degree), k)
	}
	// Stop the chain one multiplication short and let that row carry qC = k.
	if degree%2 == 0 {
		u := Pow(api, x, degree/2)
		return sparse.EvaluatePlonkExpression(u, u, 0, 0, 1, k) // u^2 + k
	}
	u := Pow(api, x, degree-1)
	return sparse.EvaluatePlonkExpression(u, x, 0, 0, 1, k) // u*x + k
}

// PowerMap applies the power map x^degree elementwise to state, returning a new
// slice (a full S-box layer, as used by Poseidon-like constructions). Callers
// that power only one branch — GMiMC — call Pow directly.
func PowerMap(api frontend.API, degree int, state []frontend.Variable) []frontend.Variable {
	out := make([]frontend.Variable, len(state))
	for i, x := range state {
		out[i] = Pow(api, x, degree)
	}
	return out
}

// PowerMapPartial applies x^degree to the first u branches and leaves the rest
// unchanged — the partial S-box layer of Hades-strategy internal rounds
// (Poseidon2, Neptune). Returns a new slice.
func PowerMapPartial(api frontend.API, degree int, state []frontend.Variable, u int) []frontend.Variable {
	out := make([]frontend.Variable, len(state))
	for i, x := range state {
		if i < u {
			out[i] = Pow(api, x, degree)
		} else {
			out[i] = x
		}
	}
	return out
}
