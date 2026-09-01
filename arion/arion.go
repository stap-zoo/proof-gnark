// Package arion implements the Arion permutation (arxiv 2303.04639) as a gnark
// circuit. Arion is a "type-2" design: its nonlinear layer is a Generalized
// Triangular Dynamical System (GTDS) whose last branch applies the high-degree
// power map x^{1/alpha2}. That single inverse is never evaluated directly
// in-circuit — algebra.InvPow supplies the value with a hint and pins it with the
// cheap forward check w^alpha2 == x. Every other branch uses only the small
// forward power x^alpha1, so the whole permutation is forward/low-degree.
//
// File layout mirrors the Python reference framework (arion/):
//   - arion.go     — the permutation round function (GTDS + affine layer)     ≈ hash.py
//   - params.go    — parameter generation and sanity checks                    ≈ params.py
//   - instances.go — concrete instances + harness/mode wiring                  ≈ instances.py
//
// The sponge mode of operation comes from package mode.
package arion

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit Arion permutation.
type Permutation struct {
	params *Parameters
}

var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns an Arion permutation for the given parameters.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the Arion permutation: R rounds, each of the shape
//
//	nonlinear (GTDS) -> linear layer (M) -> constant addition (rcons[r])
//
// There is no pre- or post-round step (Arion does no work outside the loop).
// Mirrors Arion.permutation in the reference arion/hash.py.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	out := make([]frontend.Variable, len(state))
	copy(out, state)

	for r := 0; r < pr.R; r++ {
		out = p.nonlinearLayer(api, out, r)
		// The round constants follow the matrix, so the layer absorbs them: each
		// rides on the gate that produced its slot instead of costing one of its
		// own. The layer itself runs as an addition program where one is cheaper
		// than the dense product (algebra.ApplyMatrixConst).
		out = algebra.ApplyMatrixConst(api, pr.M, pr.Slp, out, pr.Rcons[r])
	}
	return out
}

// nonlinearLayer is the GTDS (arxiv 2303.04639, Def. 1) for round r. Working from
// the last branch inward, it maintains a running sum sigma:
//
//	y[t-1] = x[t-1]^{1/alpha2}                       # the single inverse S-box
//	sigma  = x[t-1] + y[t-1]
//	for i = t-2 .. 0:
//	    y[i]  = x[i]^alpha1 * g_i(sigma) + h_i(sigma)
//	    sigma = sigma + x[i] + y[i]
//
// where g_i(s) = s^2 + a*s + b and h_i(s) = s^2 + c*s share the same s^2 term.
// Only branch t-1 uses the high-degree inverse; every other branch is forward.
//
// Because g_i and h_i are both monic, h_i = g_i - b + (c-a)*s, and the whole branch
// collapses onto ONE quadratic instead of two:
//
//	y[i] = (x[i]^alpha1 + 1) * g_i(sigma) + (c-a)*sigma - b
//
// which is what is emitted below — three PLONK gates where the literal form costs
// five. See the per-line notes; each of the three fills selectors that gnark's
// api.Mul/api.Add leave empty, and none of it changes the arithmetic (the reference
// vectors are the proof). R1CS is unaffected: there the constant-folding is free
// either way and both forms cost the same two multiplications.
//
// What the layer costs as emitted, writing c(a) = floor(log2 a) + hw(a) - 1 for the
// square-and-multiply cost of x^a: (c(alpha1) + 5)(t-1) + c(alpha2) - 1 PLONK gates,
// which is 8t at alpha1 = 5 and alpha2 = 257 — the exponents both scalar fields
// select — so 24 a round at t=3 and 32 at t=4. Of the eight a branch, c(alpha1) = 3
// are the power map and three are the running sum — one for the quadratic core
// sigma^2 + a*sigma, two for extending sigma, which is a 3-term sum and so cannot be
// one row. R1CS is 5t + 5, additions being free there.
func (p *Permutation) nonlinearLayer(api frontend.API, x []frontend.Variable, r int) []frontend.Variable {
	pr := p.params
	t := pr.T
	y := make([]frontend.Variable, t)

	// Last branch: the sole inverse S-box, verified forward by algebra.InvPow.
	y[t-1] = algebra.InvPow(api, x[t-1], pr.Alpha2)
	sigma := api.Add(x[t-1], y[t-1])

	for i := t - 2; i >= 0; i-- {
		a, b := pr.CoeffsG[r][i][0], pr.CoeffsG[r][i][1]
		c := pr.CoeffsH[r][i]
		cMinusA := new(big.Int).Sub(c, a)

		// The quadratic core shared by g_i and h_i, as one gate.
		core := algebra.MulAffine(api, sigma, sigma, a) // sigma^2 + a*sigma
		// x[i]^alpha1 + 1, at the cost of the power map alone (the +1 rides the
		// last row of the square-and-multiply chain).
		z := algebra.PowPlusConst(api, x[i], pr.Alpha1, 1)
		// z * g_i(sigma) — g_i's constant b rides this row's linear selector
		// rather than costing a gate to add to core first.
		zg := algebra.MulAffine(api, z, core, b)
		// The linear remainder, both terms free on one row.
		y[i] = api.Add(zg, api.Mul(cMinusA, sigma), new(big.Int).Neg(b))
		if i > 0 {
			// Extend sigma_{i+1,n} to sigma_{i,n} — two rows, since three wires
			// do not fit one. Not on the last branch: nothing reads sigma_{0,n},
			// and gnark's SCS builder emits dead rows rather than eliminating them.
			sigma = api.Add(sigma, x[i], y[i])
		}
	}
	return y
}
