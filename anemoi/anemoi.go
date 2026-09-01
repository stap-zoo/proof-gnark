package anemoi

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit Anemoi permutation.
type Permutation struct {
	params *Parameters
}

var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns an Anemoi permutation for the given parameters.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t = 2*l.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the Anemoi permutation. The state is laid out as x || y with
// x = state[:l], y = state[l:]. Each of the R rounds is
//
//	constant addition -> linear layer -> non-linear layer
//
// followed by one trailing linear layer (the reference's _post_rounds). The
// non-linear layer is the open Flystel per column, evaluated in-circuit through
// its closed form (algebra.ClosedFlystel); everything else folds into constant
// additions and constant matrix products.
//
// The constant addition is applied *after* the linear layer rather than before it,
// with the constants premultiplied by that layer: L(state + c) == L(state) + L(c),
// and LinearConstants holds the L(c) rows. That is what makes the round constants
// free — sitting after the layer, each rides on the gate that produced its slot
// instead of costing a PLONK gate of its own. The arithmetic is unchanged, as the
// reference vectors confirm.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	out := make([]frontend.Variable, len(state))
	copy(out, state)

	for r := 0; r < p.params.R; r++ {
		out = p.linearLayer(api, out, p.params.LinearConstants[r])
		out = p.nonlinearLayer(api, out)
	}
	return p.linearLayer(api, out, nil)
}

// linearLayer applies M_x to the x-lane and M_y to the y-lane, then the
// pseudo-Hadamard transform, adding consts (indexed like the state, x-lane then
// y-lane) to the result. Round-independent apart from those constants.
//
// The PHT is emitted as the two independent combinations x' = 2*xm + ym and
// y' = xm + ym rather than the sequential "y += x; x += y". Both forms cost two
// gates per column and compute the same thing, but in the sequential form x' reads
// the already-updated y', so folding y's constant into that gate would corrupt x'.
// Written independently, both constants fold.
func (p *Permutation) linearLayer(api frontend.API, state []frontend.Variable, consts []*big.Int) []frontend.Variable {
	l := p.params.L
	xm := algebra.ApplyMatrix(api, p.params.Mx, p.params.SlpMx, state[:l])
	ym := algebra.ApplyMatrix(api, p.params.My, p.params.SlpMy, state[l:])

	out := make([]frontend.Variable, p.params.T)
	for i := 0; i < l; i++ {
		x2 := api.Mul(2, xm[i])
		if consts == nil {
			out[i], out[l+i] = api.Add(x2, ym[i]), api.Add(xm[i], ym[i])
			continue
		}
		out[i] = api.Add(x2, ym[i], consts[i])
		out[l+i] = api.Add(xm[i], ym[i], consts[l+i])
	}
	return out
}

// nonlinearLayer applies the closed Flystel S-box to each column (x_i, y_i).
func (p *Permutation) nonlinearLayer(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	out := make([]frontend.Variable, pr.T)
	for i := 0; i < pr.L; i++ {
		out[i], out[pr.L+i] = algebra.ClosedFlystel(api, state[i], state[pr.L+i], pr.Alpha, pr.Quad, pr.Beta, pr.Gamma, pr.Delta)
	}
	return out
}
