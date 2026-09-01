package poseidon

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit Poseidon permutation.
type Permutation struct {
	params *Parameters
}

var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns a Poseidon permutation for the given parameters.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the Poseidon permutation: R rounds of (ARK -> S-box -> MDS),
// with no leading matrix or output whitening. The round constants are added to
// every branch each round; the S-box is the full-width power map in external
// rounds and the single-branch power map in internal rounds; the MDS is one
// matrix shared by all rounds.
//
// The loop is written rotated — each MDS multiplication absorbs the next round's
// constants rather than each round adding its own — because a constant added right
// after a linear layer rides on the gate that produced its slot for free, while a
// separate addition costs a PLONK gate (algebra.MatVecMulConst). Unlike Poseidon2
// there is no leading matrix, so round 0's ARK has nothing to fold into and stays
// explicit; the last round's MDS absorbs nothing. The arithmetic is unchanged, as
// the reference vectors confirm.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	out := make([]frontend.Variable, len(state))
	copy(out, state)

	out = addConstants(api, out, pr.RoundConstants[0]) // ARK for round 0: no layer precedes it
	for r := 0; r < pr.R; r++ {
		if pr.isInternal(r) {
			out = algebra.PowerMapPartial(api, pr.Alpha, out, pr.U)
		} else {
			out = algebra.PowerMap(api, pr.Alpha, out)
		}
		out = pr.applyLayer(api, r, out)
	}
	return out
}

// applyLayer applies round r's linear layer, carrying the constants of the round
// that follows it. That is the MDS matrix and RoundConstants[r+1] — except in the
// internal-round block, where the sparse factorisation supplies both, and in the
// round that feeds the block, whose matrix and constants carry the residual
// (Parameters.Internal). Which layer runs where is the whole of the difference;
// the state it computes is the same either way.
func (p *Parameters) applyLayer(api frontend.API, r int, state []frontend.Variable) []frontend.Variable {
	if f := p.Internal; f != nil {
		if p.isInternal(r) {
			return f.Apply(api, r-p.RExtBeg, state)
		}
		if r == p.RExtBeg-1 {
			return f.ApplyPre(api, state)
		}
	}
	var next []*big.Int
	if r+1 < p.R {
		next = p.RoundConstants[r+1]
	}
	return algebra.MatVecMulConst(api, p.M, state, next)
}

// addConstants returns state[i] + consts[i] for each branch.
func addConstants(api frontend.API, state []frontend.Variable, consts []*big.Int) []frontend.Variable {
	out := make([]frontend.Variable, len(state))
	for i := range state {
		out[i] = api.Add(state[i], consts[i])
	}
	return out
}
