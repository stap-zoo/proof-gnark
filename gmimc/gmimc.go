// Package gmimc implements the S-GMiMC permutation as a gnark circuit.
//
// File layout mirrors the Python reference framework:
//   - gmimc.go     — the permutation round function                          ≈ hash.py
//   - params.go    — parameter generation and sanity checks                  ≈ params.py
//   - instances.go — concrete instances + harness/mode wiring (Sponges, Targets) ≈ instances.py
//
// Modes of operation come from package mode (a construction only declares which
// ones it is eligible for, in instances.go).
package gmimc

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// form selects which of Permute's two arithmetizations to build. The zero value
// picks by width, which is what every caller outside circuit_test.go wants; the
// other two force one, so the test can measure them against each other at widths
// where Permute would not choose them.
type form uint8

const (
	formByWidth form = iota
	formDirect
	formAccumulated
)

// Permutation is the in-circuit S-GMiMC permutation.
//
// It has two arithmetizations of the same function and picks between them by
// width — see Permute.
type Permutation struct {
	params *Parameters
	form   form
}

// Compile-time check that Permutation satisfies the shared interface.
var _ permutation.Permutation = (*Permutation)(nil)

// NewPermutation returns an S-GMiMC permutation for the given parameters, in
// whichever arithmetization is cheaper at that width.
func NewPermutation(params *Parameters) *Permutation {
	return &Permutation{params: params}
}

// NewPermutationDirect returns the same permutation forced into the direct,
// round-by-round arithmetization at any width. Only circuit_test.go builds one:
// it is the baseline the accumulator circuit is measured against.
func NewPermutationDirect(params *Parameters) *Permutation {
	return &Permutation{params: params, form: formDirect}
}

// NewPermutationAccumulated returns the same permutation forced into the
// sliding-window arithmetization at any width, including the widths where it is
// the more expensive of the two. Only circuit_test.go builds one, to pin where
// the crossover is.
func NewPermutationAccumulated(params *Parameters) *Permutation {
	return &Permutation{params: params, form: formAccumulated}
}

// Width returns the state size t.
func (p *Permutation) Width() int { return p.params.T }

// Permute applies the S-GMiMC permutation to state (an unbalanced Feistel: each
// round powers the first branch through the round constant and adds the result
// into every other branch, then applies the linear layer). The two generic
// building blocks — the power-map S-box and the matrix-vector product — come
// from the algebra package; this function only wires them per the GMiMC spec.
//
// # Round constants are free, in both arithmetizations
//
// A constant cannot fold into a multiplication, so `(x0 + rc)^alpha` would need a
// gate for the addition — but a PLONK gate carries a constant selector, and every
// branch gets an addition of its own before it next reaches position 0. Folding
// the constant into that addition costs nothing, where a separate api.Add costs a
// row. So each branch's wire is allowed to carry a known constant above the
// reference state, `shift` tracks what it carries, and every addition folds in
// (what it should carry) - (what it does carry). Both are constants, so the
// correction is free too.
//
// Two additions survive that, whichever arithmetization runs: round 0's constant,
// which has no preceding gate to ride on, and one residual correction at the end,
// for the branch that left position 0 in the final round and so was never touched
// again. The arithmetic is unchanged, as the reference vectors confirm.
//
// # Two arithmetizations
//
//   - roundsDirect, t <= 4: the round as written, one addition per branch past
//     the first, R*(ell(alpha) + t - 1) gates.
//   - roundsAccumulated, t > 4: the sliding-window circuit, R*(ell(alpha) + 3)
//     gates plus a closing loop — the round cost stops growing with the width.
//
// They cross over exactly at t = 4 (t - 1 = 3): the difference is (t-4)*(R-1)
// gates in the accumulated form's favour, so t = 4 keeps the direct form for the
// same cost and fewer moving parts, and t = 3 needs it. Both satisfy
// testdata/vectors.json, and circuit_test.go pins what the second saves.
func (p *Permutation) Permute(api frontend.API, state []frontend.Variable) []frontend.Variable {
	t := p.params.T
	out := make([]frontend.Variable, t)
	copy(out, state)

	// shift[i] is the constant branch i's wire carries relative to the reference
	// state; nil means none.
	shift := make([]*big.Int, t)
	out[0] = api.Add(out[0], p.params.RoundConstants[0])
	shift[0] = p.params.RoundConstants[0]

	if p.form == formDirect || (p.form == formByWidth && t <= 4) {
		out, shift = p.roundsDirect(api, out, shift)
	} else {
		out, shift = p.roundsAccumulated(api, out, shift)
	}

	// One wire still carries the constant it picked up at position 0 in the last
	// round, and nothing came after to absorb it.
	for i := range out {
		if fold := diffMod(nil, shift[i], p.params.Field); fold != nil {
			out[i] = api.Add(out[i], fold)
			shift[i] = nil
		}
	}
	return out
}

// roundsDirect runs the R rounds as written: power the first branch, add the
// result into the other t-1, shift. state arrives with the first round's constant
// already on branch 0, and each round pre-loads the next round's constant onto the
// branch the shift will bring to the front, on the addition that branch receives
// anyway.
//
// Cost: R*(ell(alpha) + t - 1) gates.
func (p *Permutation) roundsDirect(api frontend.API, state []frontend.Variable, shift []*big.Int) ([]frontend.Variable, []*big.Int) {
	t, field := p.params.T, p.params.Field
	out := make([]frontend.Variable, t)
	copy(out, state)

	for r := 0; r < p.params.Rounds; r++ {
		// Nonlinear layer: branch 0 already carries rc_r, so it is the S-box input.
		y := algebra.Pow(api, out[0], p.params.Alpha)
		for i := 1; i < t; i++ {
			// The wire at branch 1 is the one the shift moves into position 0, so it
			// pre-loads the next round's constant; every other branch carries none.
			var want *big.Int
			if i == 1 && r+1 < p.params.Rounds {
				want = p.params.RoundConstants[r+1]
			}
			if fold := diffMod(want, shift[i], field); fold != nil {
				out[i] = api.Add(out[i], y, fold)
			} else {
				out[i] = api.Add(out[i], y)
			}
			shift[i] = want
		}
		// Linear layer: out = M · out, a cyclic shift that costs nothing and moves
		// each wire's carried constant with it.
		out = algebra.MatVecMul(api, p.params.M, out)
		shift = shiftConstants(shift)
	}
	return out, shift
}

// roundsAccumulated computes exactly what roundsDirect does — the KAT vectors
// assert it, and circuit_test.go asserts the two agree in-circuit — for
// R*(ell(alpha) + 3) gates instead of R*(ell(alpha) + t - 1), so the round cost
// stops growing with the width.
//
// The idea is that only one branch is ever read: the one at position 0, which the
// S-box consumes. Every other branch is just accumulating y's, so it need not
// exist as a wire until it comes back around. A branch takes t rounds to return to
// position 0 and receives y in t-1 of them, so what it missed is always the sum of
// the last t-1 S-box outputs. The circuit therefore keeps
//
//	win   — the last t-1 S-box outputs, a rotating window, and
//	acc   — their running sum, one eviction and one insertion per round,
//
// and settles a branch in one addition at the moment it arrives. Rotating win and
// the state are re-indexings of Go slices and cost nothing; only acc's two updates
// and the branch's own addition do, hence 3 per round at any width.
//
// That one addition per branch per pass is also the one the round constant folds
// into (see Permute), which is why GMiMC's constant — which feeds the S-box rather
// than staying in the state, so it has to come back off — still costs nothing
// here. The window never sees it: it holds S-box outputs, which carry no shift.
//
// The closing loop pays the width back once, for the t-2 branches still in flight
// when the rounds end, shrinking the window by one each time. The last branch
// needs nothing from the window: it has just left position 0, where it received no
// y at all — only its constant comes off, in Permute.
//
// Early rounds are cheaper than the count suggests: win and acc start as
// constants, so gnark folds the first t-1 evictions and the first insertion away
// rather than emitting gates for them. Nothing here requires the round count to be
// a multiple of t: the window rotates in step with the state, so it stays aligned
// wherever the rounds stop.
func (p *Permutation) roundsAccumulated(api frontend.API, state []frontend.Variable, shift []*big.Int) ([]frontend.Variable, []*big.Int) {
	t, field := p.params.T, p.params.Field
	out := make([]frontend.Variable, t)
	copy(out, state)

	var acc frontend.Variable = 0
	win := make([]frontend.Variable, t-1)
	for i := range win {
		win[i] = frontend.Variable(0)
	}

	for r := 0; r < p.params.Rounds; r++ {
		y := algebra.Pow(api, out[0], p.params.Alpha) // out[0] carries exactly rc_r

		win = rotateLeft(win)
		acc = api.Sub(acc, win[0]) // evict the output that has fallen out of the window
		win[0] = y
		acc = api.Add(acc, y) // ... and take the new one in

		out, shift = rotateLeft(out), shiftConstants(shift)
		// The branch arriving at position 0 has missed exactly acc, and takes the
		// next round's constant — net of whatever it still carries — in the same gate.
		var want *big.Int
		if r+1 < p.params.Rounds {
			want = p.params.RoundConstants[r+1]
		}
		if fold := diffMod(want, shift[0], field); fold != nil {
			out[0] = api.Add(out[0], acc, fold)
		} else {
			out[0] = api.Add(out[0], acc)
		}
		shift[0] = want
	}

	// Settle the branches still in flight, oldest first; the window shrinks by one
	// per branch, and each branch's carried constant comes off in its own gate.
	for j := 1; j < t-1; j++ {
		win = rotateLeft(win)
		acc = api.Sub(acc, win[0])
		if fold := diffMod(nil, shift[j], field); fold != nil {
			out[j] = api.Add(out[j], acc, fold)
		} else {
			out[j] = api.Add(out[j], acc)
		}
		shift[j] = nil
	}
	return out, shift
}

// rotateLeft returns v rotated one position left (out[i] = v[i+1], out[n-1] =
// v[0]) as a new slice — the same permutation as the linear layer M. It is a
// re-indexing of wires, so it costs nothing in either constraint system.
func rotateLeft(v []frontend.Variable) []frontend.Variable {
	n := len(v)
	out := make([]frontend.Variable, n)
	copy(out, v[1:])
	out[n-1] = v[0]
	return out
}

// diffMod returns (want - have) mod field, or nil when that is zero (nil operands
// count as zero).
func diffMod(want, have, field *big.Int) *big.Int {
	d := new(big.Int)
	if want != nil {
		d.Set(want)
	}
	if have != nil {
		d.Sub(d, have)
	}
	if d.Sign() == 0 {
		return nil
	}
	return d.Mod(d, field)
}

// shiftConstants applies the same cyclic shift as the linear layer (M[i][i+1]=1,
// M[t-1][0]=1) to the per-branch carried constants.
func shiftConstants(shift []*big.Int) []*big.Int {
	t := len(shift)
	out := make([]*big.Int, t)
	for i := 0; i < t-1; i++ {
		out[i] = shift[i+1]
	}
	out[t-1] = shift[0]
	return out
}
