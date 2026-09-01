package poseidon2

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit Poseidon2 permutation.
type Permutation struct {
	params *Parameters
}

var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns a Poseidon2 permutation for the given parameters.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the Poseidon2 permutation: a leading external matrix, then R
// rounds of (ARK -> S-box -> matrix). Round constants are added to every branch
// in external rounds and to the first u branches in internal rounds (the internal
// rows are zero elsewhere, so a uniform ARK is equivalent). The S-box is the
// full-width power map in external rounds and the single-branch power map in
// internal rounds; the linear layer is M_ext in external rounds and M_int in
// internal rounds.
//
// The loop is written rotated — each layer absorbs the *next* round's constants
// instead of the round adding its own — because a constant added right after a
// linear layer rides on the gate that produced its slot for free, while a separate
// addition costs a PLONK gate (see algebra.ApplyMatrixConst). The leading matrix
// takes round 0's constants and the last round's layer takes none; the arithmetic
// is unchanged, as the reference vectors confirm.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	out := make([]frontend.Variable, len(state))
	copy(out, state)

	out = algebra.ApplyMatrixConst(api, pr.MExt, pr.SlpExt, out, pr.RoundConstants[0])

	for r := 0; r < pr.R; r++ {
		next := nextConstants(pr, r)
		if pr.isInternal(r) {
			out = algebra.PowerMapPartial(api, pr.Alpha, out, pr.U)
			out = algebra.ApplyMatrixConst(api, pr.MInt, pr.SlpInt, out, next)
		} else {
			out = algebra.PowerMap(api, pr.Alpha, out)
			out = algebra.ApplyMatrixConst(api, pr.MExt, pr.SlpExt, out, next)
		}
	}
	return out
}

// nextConstants returns the ARK row for round r+1 — the one the layer closing
// round r absorbs — and nil after the final round.
func nextConstants(pr *Parameters, r int) []*big.Int {
	if r+1 < pr.R {
		return pr.RoundConstants[r+1]
	}
	return nil
}
