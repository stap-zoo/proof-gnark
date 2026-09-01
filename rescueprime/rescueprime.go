// Package rescueprime implements the Rescue-Prime permutation (the RescuePrime
// class of the Marvellous family, Rescue-XLIX) as a gnark circuit. Rescue-Prime
// is a "type-2" design: its inverse S-box x -> x^{1/alpha} is a high-degree
// power map, so it is never evaluated directly in-circuit — algebra.InvPow
// supplies the value with a hint and pins it with the cheap forward check
// w^alpha == x. The forward S-box x -> x^alpha is a plain low-degree power map.
//
// File layout mirrors the Python reference framework (marvellous/):
//   - rescueprime.go — the permutation round function                        ≈ hash.py
//   - params.go      — parameter generation and sanity checks                 ≈ params.py
//   - instances.go   — concrete instances + harness/mode wiring              ≈ instances.py
//
// The sponge mode of operation comes from package mode.
package rescueprime

import (
	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit Rescue-Prime permutation.
type Permutation struct {
	params *Parameters
}

var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns a Rescue-Prime permutation for the given parameters.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the Rescue-Prime permutation: 2*R half-rounds forming R double
// rounds, each of the shape (F)(B). A half-round is
//
//	S-box -> linear layer (M) -> constant addition (rcons[r])
//
// where the S-box is the forward power map x^alpha on even half-rounds (F) and
// the inverse power map x^{1/alpha} on odd half-rounds (B). There is no pre- or
// post-round step (Rescue-Prime folds the constant into every round). Mirrors
// RescuePrime.permutation in the reference marvellous/hash.py.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	out := make([]frontend.Variable, len(state))
	copy(out, state)

	for r := 0; r < 2*pr.R; r++ {
		if r%2 == 0 {
			out = p.forwardSbox(api, out)
		} else {
			out = p.inverseSbox(api, out)
		}
		// The round constants follow the matrix, so the layer absorbs them: each
		// rides on the gate that produced its slot instead of costing one of its
		// own (algebra.MatVecMulConst).
		out = algebra.MatVecMulConst(api, pr.M, out, pr.RoundConstants[r])
	}
	return out
}

// forwardSbox applies x -> x^alpha element-wise (the cheap direction).
func (p *Permutation) forwardSbox(api frontend.API, state []frontend.Variable) []frontend.Variable {
	return algebra.PowerMap(api, p.params.Alpha, state)
}

// inverseSbox applies x -> x^{1/alpha} element-wise via algebra.InvPow: a hint
// supplies each value and the low-degree constraint w^alpha == x pins it, so the
// high-degree inverse exponent is never evaluated in-circuit.
func (p *Permutation) inverseSbox(api frontend.API, state []frontend.Variable) []frontend.Variable {
	out := make([]frontend.Variable, len(state))
	for i, x := range state {
		out[i] = algebra.InvPow(api, x, p.params.Alpha)
	}
	return out
}
