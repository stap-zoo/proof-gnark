package neptune

import (
	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit Neptune permutation.
type Permutation struct {
	params *Parameters
}

var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns a Neptune permutation for the given parameters.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the Neptune permutation. The Hades structure — a leading
// external matrix, then R rounds of (ARK -> S-box -> matrix), then output
// whitening — is wired here from the algebra gadgets; the round constants and
// matrices come from the parameters.
//
// The S-box differs by round: external rounds apply the quadratic Lai-Massey
// S-box pair-wise across the state; internal rounds apply the degree-alpha power
// map to the first U branches. Likewise the matrix is M_ext (external) or M_int
// (internal).
//
// The loop is written rotated — each layer absorbs the constants of the round that
// follows it, rather than each round adding its own — because a constant added
// right after a linear layer rides on the gate that produced its slot for free,
// while a separate addition costs a PLONK gate (algebra.ApplyMatrixConst). The
// leading matrix takes row 0 (all zero here) and the final round's layer takes the
// whitening row, so no separate constant addition is left. The arithmetic is
// unchanged, as the reference vectors confirm.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	out := make([]frontend.Variable, len(state))
	copy(out, state)

	// Leading external matrix (_pre_rounds).
	out = algebra.MatVecMulConst(api, pr.MExt, out, pr.RoundConstants[0])

	for r := 0; r < pr.R; r++ {
		if pr.isInternal(r) {
			out = algebra.PowerMapPartial(api, pr.Alpha, out, pr.U)
		} else {
			out = p.laiMasseyLayer(api, out)
		}
		out = pr.applyLayer(api, r, out)
	}
	return out
}

// applyLayer applies round r's linear layer, carrying the constants of the round
// that follows it (rows 1..R of RoundConstants, the last being the whitening
// row). That is M_ext or M_int as the round dictates — except in the
// internal-round block, where the sparse factorisation supplies the matrix and
// the constants, and in the round that feeds the block, whose matrix and
// constants carry the residual (Parameters.Internal). Which layer runs where is
// the whole of the difference; the state it computes is the same either way.
func (p *Parameters) applyLayer(api frontend.API, r int, state []frontend.Variable) []frontend.Variable {
	if f := p.Internal; f != nil {
		if p.isInternal(r) {
			return f.Apply(api, r-p.RExtBeg, state)
		}
		if r == p.RExtBeg-1 {
			return f.ApplyPre(api, state)
		}
	}
	next := p.RoundConstants[r+1]
	if p.isInternal(r) {
		return algebra.ApplyMatrixConst(api, p.MInt, p.SlpInt, state, next)
	}
	return algebra.MatVecMulConst(api, p.MExt, state, next)
}

// laiMasseyLayer applies the Neptune external S-box to consecutive branch pairs.
func (p *Permutation) laiMasseyLayer(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	out := make([]frontend.Variable, len(state))
	for i := 0; i < len(state); i += 2 {
		out[i], out[i+1] = algebra.LaiMassey(api, state[i], state[i+1], pr.LMAlpha, pr.LMBeta, pr.LMGamma, pr.LMMatrix)
	}
	return out
}

func (p *Parameters) isInternal(r int) bool {
	return p.RExtBeg <= r && r < p.RExtBeg+p.RInt
}
