package polocolo

import (
	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
)

// Permutation is the in-circuit Polocolo permutation over F_p^t, mirroring
// Polocolo.permutation in the reference polocolo/hash.py.
type Permutation struct {
	params *Parameters
}

// NewPermutation wraps a parameter set as a permutation.Permutation.
func NewPermutation(params *Parameters) *Permutation { return &Permutation{params: params} }

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies R rounds of
//
//	linear layer (M*x, via the published addition program) ->
//	constant addition (rcons[r]) -> S-box layer (power-residue S-box)
//
// followed by the bare final linear layer (LinLayer^(R) carries no constant,
// c^(R) = 0).
//
// The S-box gadget is built once per *compilation*, not per call: its tables cost
// one constraint per entry (m of them, up to 1024 here), so a circuit that calls
// the permutation k times must pay for those entries once and route all k*t*R
// S-box evaluations through the same lookup argument. algebra.Shared caches it;
// see the note there on why that is both sound and what a real circuit does.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	pr := p.params
	sbox := algebra.Shared(api, sboxKey{pr}, func() *algebra.PowerResidueSbox {
		return algebra.NewPowerResidueSbox(api, pr.M, pr.G, pr.Sigma, pr.GrMode)
	})

	out := state
	for r := 0; r < pr.R; r++ {
		// The round constants follow the linear layer, so the program absorbs
		// them: each rides on the gate that produced its slot instead of costing
		// one of its own (algebra.SLP.ApplyConst).
		out = pr.Slp.ApplyConst(api, out, pr.Rcons[r])
		for i := range out {
			out[i] = sbox.Apply(out[i])
		}
	}
	return pr.Slp.Apply(api, out) // _post_rounds: the bare matrix product
}

// sboxKey identifies one S-box gadget per (compilation, parameter set). Instances
// differ in m, the generator and sigma, so their tables differ and must not be
// shared; the parameter pointer distinguishes them.
type sboxKey struct{ params *Parameters }
