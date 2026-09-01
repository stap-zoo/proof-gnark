// Package griffin implements the Griffin-pi permutation (eprint 2022/403) as a
// gnark circuit. Griffin is a "type-2" design: its nonlinear layer applies the
// high-degree inverse power map x^{1/d} on the first branch, the small forward
// power x^d on the second, and a Horst-style x_i * G_i(L_i(...)) on the rest,
// where G_i is a root-free quadratic. The single inverse is never evaluated
// directly in-circuit — algebra.InvPow supplies the value with a hint and pins
// it with the cheap forward check w^d == x — so the whole permutation is
// forward/low-degree.
//
// File layout mirrors the Python reference framework (griffin/):
//   - griffin.go   — the permutation round function (nonlinear + affine layer)  ≈ hash.py
//   - params.go    — parameter generation and sanity checks                     ≈ params.py
//   - instances.go — concrete instances + harness/mode wiring                   ≈ instances.py
//
// The sponge mode of operation comes from package mode.
package griffin

import (
	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit Griffin permutation.
type Permutation struct {
	params *Parameters
}

var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns a Griffin permutation for the given parameters.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the Griffin permutation: an initial matrix multiplication,
// then R rounds each of the shape
//
//	nonlinear -> linear layer (M) -> constant addition (rcons[r])
//
// The final round's constants are all zero (rcons has R rows, the last a zero
// row), so it is a plain nonlinear+linear step. There is no post-round step.
// Mirrors Griffin.permutation in the reference griffin/hash.py.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params

	out := algebra.ApplyMatrix(api, pr.M, pr.Slp, state) // initial matrix multiplication (_pre_rounds)
	for r := 0; r < pr.R; r++ {
		out = p.nonlinearLayer(api, out)
		// The round constants follow the matrix, so the layer absorbs them: each
		// rides on the gate that produced its slot instead of costing one of its
		// own (algebra.ApplyMatrixConst).
		out = algebra.ApplyMatrixConst(api, pr.M, pr.Slp, out, pr.Rcons[r])
	}
	return out
}

// nonlinearLayer is Griffin's S-box layer (eprint 2022/403, Eq. 6). It is
// round-independent (the same map every round):
//
//	y0 = x0^{1/d}                          # the single inverse S-box
//	y1 = x1^d
//	y_i = x_i * G_i(L_i(y0, y1, z))        for i = 2 .. t-1   (Horst step)
//
// where L_i(y0, y1, z) = (i-1)*y0 + y1 + z with z = 0 for i = 2 and the *input*
// word x_{i-1} for i > 2, and G_i(l) = l^2 + a_i*l + b_i is root-free (its
// coefficients were sampled so a_i^2 - 4*b_i is a non-residue). Only branch 0
// uses the high-degree inverse; every other branch is forward, low degree.
func (p *Permutation) nonlinearLayer(api frontend.API, x []frontend.Variable) []frontend.Variable {
	pr := p.params
	t := pr.T
	y := make([]frontend.Variable, t)

	y[0] = algebra.InvPow(api, x[0], pr.Alpha) // x0^{1/d}, verified forward by InvPow
	y[1] = algebra.Pow(api, x[1], pr.Alpha)    // x1^d (forward, low degree)

	for i := 2; i < t; i++ {
		a, b := pr.CoeffsG[i-2][0], pr.CoeffsG[i-2][1]

		// L_i = (i-1)*y0 + y1 + z; the constant*var terms fold away (no constraint).
		l := api.Add(api.Mul(i-1, y[0]), y[1])
		if i > 2 {
			l = api.Add(l, x[i-1]) // feedback word (input x_{i-1}); none for i = 2
		}

		// Both remaining steps are a product plus a linear term in the same wires,
		// which is one PLONK gate each (algebra.MulAffine): the inner row holds
		// l^2 + a*l, and G_i's constant b rides the outer row's linear selector as
		// x_i*(core + b) rather than costing a gate to add to core first.
		core := algebra.MulAffine(api, l, l, a)      // l^2 + a*l
		y[i] = algebra.MulAffine(api, x[i], core, b) // x_i*(l^2 + a*l + b) = x_i*G_i(l)
	}
	return y
}
