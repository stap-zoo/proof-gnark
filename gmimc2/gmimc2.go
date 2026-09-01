// Package gmimc2 implements the GMiMC2 (GMiMC-erf2) permutation as a gnark circuit.
//
// GMiMC2 is the updated GMiMC of the SoK: same expanding-round-function Feistel,
// three changes. The S-box exponent is alpha = 2^k (optimal degree growth per
// multiplication), the round constant is added into the state before the power
// map, and an input/output matrix M_IO is applied at the start and the end of the
// permutation. Round numbers are re-derived against the current cryptanalysis;
// see instances.go for the table and params.go for what each parameter means.
//
// Unlike every other construction here, GMiMC2 has no counterpart in the Python
// reference framework (../ref implements GMiMC only). Its oracle is
// gnark-hashes/gmimc2_ref.py, a standalone pure-Python implementation in this
// repository's root, which produced testdata/vectors.json and follows the
// reference's conventions for the sampler, the sponge and the KAT shape.
//
// File layout mirrors the rest of the repo:
//   - gmimc2.go    — the permutation round function                          ≈ hash.py
//   - params.go    — parameter derivation and sanity checks                  ≈ params.py
//   - instances.go — concrete instances + harness/mode wiring                ≈ instances.py
package gmimc2

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Permutation is the in-circuit GMiMC2 permutation.
//
// It has two arithmetizations of the same function, and picks between them by
// width — see Permute. Direct forces the wide one, for the test that measures
// what the other saves.
type Permutation struct {
	params *Parameters
	Direct bool
}

// Compile-time check that Permutation satisfies the shared interface.
var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns a GMiMC2 permutation for the given parameters, in
// whichever arithmetization is cheaper at that width.
func NewPermutation(params *Parameters) *Permutation {
	return &Permutation{params: params}
}

// NewPermutationDirect returns the same permutation forced into the direct,
// round-by-round arithmetization. Only gmimc2_test.go builds one: it is the
// baseline the accumulator circuit is measured against.
func NewPermutationDirect(params *Parameters) *Permutation {
	return &Permutation{params: params, Direct: true}
}

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the GMiMC2 permutation:
//
//	pi(x) = M_IO . R^(R) . ... . R^(1) (M_IO . x)
//
// Each round R^(i) adds the round constant into the first branch, powers it —
// x_1 <- x_1 + c^(i), then y <- x_1^alpha — adds y into every other branch, and
// rotates the state one position left. Note where the constant lands: in the
// STATE, not merely in the S-box's input. Branch 1 keeps it and carries it away
// on the rotation. That is what separates this round from GMiMC's, where the
// constant feeds (x_0 + rc)^alpha and is then discarded, and it is what makes the
// constant free in PLONK — it merges into a branch addition the round already
// performs, rather than needing a gate of its own.
//
// M_IO is free at every registered instance: all have capacity 1, where the
// specification sets M_IO to the identity, and MatVecMul emits nothing for a row
// with a single unit coefficient. It is applied unconditionally rather than
// special-cased, so a non-identity M_IO would need no change here.
//
// Two arithmetizations compute this, and the width picks between them, matching
// the specification's own threshold (the improved circuit is proposed for t > 4):
//
//   - roundsDirect, t <= 4: one branch addition per branch past the first, so
//     R*(k + t - 1) gates for alpha = 2^k.
//   - roundsAccumulated, t > 4: the specification's alg:gmimc2-efficient-circuit,
//     R*(k + 3) gates plus a closing loop, independent of t.
//
// They cross over exactly at t = 4 (t-1 = 3), which is why t = 4 keeps the direct
// form: same cost, fewer moving parts. Both satisfy testdata/vectors.json, and
// gmimc2_test.go pins what the second saves at every width.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	// Leading M_IO, with the first round's constant folded into it:
	// MatVecMulConst rides it on row 0's last addition, so it costs nothing that
	// a plain M_IO.x did not already cost. Both forms start here.
	lead := make([]*big.Int, p.params.T)
	lead[0] = p.params.RoundConstants[0]
	out := algebra.MatVecMulConst(api, p.params.MIO, state, matVec(p.params.MIO, lead, p.params.Field))

	if p.Direct || p.params.T <= 4 {
		out = p.roundsDirect(api, out)
	} else {
		out = p.roundsAccumulated(api, out)
	}
	return algebra.MatVecMul(api, p.params.MIO, out)
}

// roundsDirect runs the R rounds as written: power the first branch, add the
// result into the other t-1, rotate. state arrives with the first round's
// constant already on branch 0, and each round pre-loads the next round's
// constant onto the branch the rotation will bring to the front — a PLONK gate
// carries a constant selector, so folding it into that branch's addition is free
// where a separate api.Add would cost a gate.
//
// Cost: R*(k + t - 1) gates for alpha = 2^k, k of them multiplications.
func (p *Permutation) roundsDirect(api frontend.API, state []frontend.Variable) []frontend.Variable {
	t := p.params.T
	out := make([]frontend.Variable, t)
	copy(out, state)

	for r := 0; r < p.params.Rounds; r++ {
		// Branch 0 already carries c^(r), so it is the S-box input as it stands.
		y := algebra.Pow(api, out[0], p.params.Alpha)
		for i := 1; i < t; i++ {
			// Branch 1 is the one the rotation moves to the front, so it takes the
			// next round's constant along with y.
			if i == 1 && r+1 < p.params.Rounds {
				out[i] = api.Add(out[i], y, p.params.RoundConstants[r+1])
			} else {
				out[i] = api.Add(out[i], y)
			}
		}
		// Linear layer: out = M . out, a rotation that costs nothing.
		out = algebra.MatVecMul(api, p.params.M, out)
	}
	return out
}

// roundsAccumulated is alg:gmimc2-efficient-circuit, the specification's improved
// circuit for t > 4. It computes exactly what roundsDirect does — the KAT vectors
// and gmimc2_ref.py's self-test both assert it — for R*(k + 3) gates instead of
// R*(k + t - 1), so the round cost stops growing with the width.
//
// The idea is that only one branch is ever read: the one at position 0, which the
// S-box consumes. Every other branch is just accumulating y's, so it need not
// exist as a wire until it comes back around. A branch takes t rounds to return to
// position 0 and receives y in t-1 of them, so what it missed is always the sum of
// the last t-1 S-box outputs. The circuit therefore keeps
//
//	a     — the last t-1 S-box outputs, a rotating window, and
//	acc   — their running sum, one eviction and one insertion per round,
//
// and settles a branch in one addition at the moment it arrives. Rotating a and
// the state are re-indexings of Go slices and cost nothing; only acc's two updates
// and the branch's own addition do, hence 3 per round at any width.
//
// The closing loop pays that back for the t-2 branches still in flight when the
// rounds end, shrinking the window by one each time. The last branch needs
// nothing: it has just left position 0, where it received no y at all. And nothing
// removes the constant a branch picked up at position 0 — under this round
// function the constant belongs to the state, so carrying it away is the point.
//
// Early rounds are cheaper than the count suggests: a and acc start as constants,
// so gnark folds the evictions and the first insertion away rather than emitting
// gates for them.
func (p *Permutation) roundsAccumulated(api frontend.API, state []frontend.Variable) []frontend.Variable {
	t := p.params.T
	xi := make([]frontend.Variable, t)
	copy(xi, state)

	var acc frontend.Variable = 0
	a := make([]frontend.Variable, t-1)
	for i := range a {
		a[i] = frontend.Variable(0)
	}

	for r := 0; r < p.params.Rounds; r++ {
		xAlpha := algebra.Pow(api, xi[0], p.params.Alpha)

		a = rotateLeft(a)
		acc = api.Sub(acc, a[0]) // evict the output that has fallen out of the window
		a[0] = xAlpha
		acc = api.Add(acc, xAlpha) // ... and take the new one in

		xi = rotateLeft(xi)
		// The branch arriving at position 0 has missed exactly acc, and takes the
		// next round's constant in the same gate.
		if r+1 < p.params.Rounds {
			xi[0] = api.Add(xi[0], acc, p.params.RoundConstants[r+1])
		} else {
			xi[0] = api.Add(xi[0], acc)
		}
	}

	// Settle the branches still in flight, oldest first; the window shrinks by one
	// per branch, and the last branch (index t-1) has nothing to collect.
	for j := 1; j < t-1; j++ {
		a = rotateLeft(a)
		acc = api.Sub(acc, a[0])
		xi[j] = api.Add(xi[j], acc)
	}
	return xi
}

// rotateLeft returns v rotated one position left (out[i] = v[i+1], out[n-1] =
// v[0]) as a new slice — the same permutation as the linear layer M, and the
// "rotate one position to the left" of the specification's algorithm. It is a
// re-indexing of wires, so it costs nothing in either constraint system.
func rotateLeft(v []frontend.Variable) []frontend.Variable {
	n := len(v)
	out := make([]frontend.Variable, n)
	copy(out, v[1:])
	out[n-1] = v[0]
	return out
}

// matVec is the native (out-of-circuit) product m.v over the field, for a vector
// of constants where a nil entry means zero and a zero result stays nil — the form
// algebra.MatVecMulConst takes. It transports the first round constant through M_IO.
func matVec(m [][]*big.Int, v []*big.Int, field *big.Int) []*big.Int {
	out := make([]*big.Int, len(m))
	for i, row := range m {
		acc := new(big.Int)
		for j, coeff := range row {
			if v[j] == nil || coeff.Sign() == 0 {
				continue
			}
			acc.Add(acc, new(big.Int).Mul(coeff, v[j]))
		}
		if acc.Sign() != 0 {
			out[i] = acc.Mod(acc, field)
		}
	}
	return out
}
