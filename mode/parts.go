package mode

import "github.com/consensys/gnark/frontend"

// The three pluggable pieces a sponge is configured with, ported from
// ref/utils/mode.py. Padding and IV are purely structural — a circuit's input
// length is fixed at compile time, so they emit constants, not wires — while
// Absorb mixes a block into the state and needs the api.

// Padding pads data to a whole number of rate-sized blocks, returning a new
// slice. Ported from the pad_* rules in mode.py.
type Padding func(data []frontend.Variable, rate int) []frontend.Variable

// Absorb mixes a rate-sized block into the state, returning a new state. Ported
// from add_to_start / replace_start in matrix.py.
type Absorb func(api frontend.API, state, block []frontend.Variable) []frontend.Variable

// IV builds the capacity-sized initial value from the (unpadded) input length,
// capacity and digest size — the domain-separation constants a sponge starts
// from. Returns capacity elements.
type IV func(inputLen, capacity, digest int) []frontend.Variable

// ---------------------------------------------------------------------------
// Padding rules
// ---------------------------------------------------------------------------

// PadZero zero-pads to the next multiple of rate (no-op if already aligned).
func PadZero(data []frontend.Variable, rate int) []frontend.Variable {
	n := (rate - len(data)%rate) % rate
	return appendZeros(clone(data), n)
}

// PadOne appends a single 1, then zero-pads to the next multiple of rate; always
// applied, even when already aligned.
func PadOne(data []frontend.Variable, rate int) []frontend.Variable {
	n := (rate - (len(data)+1)%rate) % rate
	return appendZeros(append(clone(data), frontend.Variable(1)), n)
}

// ---------------------------------------------------------------------------
// Absorb rules
// ---------------------------------------------------------------------------

// AddToStart adds block element-wise into the first len(block) state positions.
func AddToStart(api frontend.API, state, block []frontend.Variable) []frontend.Variable {
	out := clone(state)
	for i := range block {
		out[i] = api.Add(out[i], block[i])
	}
	return out
}

// ---------------------------------------------------------------------------
// IV builders
// ---------------------------------------------------------------------------

// ZeroIV is the all-zero initial value.
func ZeroIV(inputLen, capacity, digest int) []frontend.Variable {
	return zeros(capacity)
}

// LengthIV encodes the input length in the first capacity element, zeros after —
// the domain separation GMiMC's sponge uses.
func LengthIV(inputLen, capacity, digest int) []frontend.Variable {
	iv := zeros(capacity)
	iv[0] = inputLen
	return iv
}

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

func clone(v []frontend.Variable) []frontend.Variable {
	return append([]frontend.Variable{}, v...)
}

func zeros(n int) []frontend.Variable {
	out := make([]frontend.Variable, n)
	for i := range out {
		out[i] = frontend.Variable(0)
	}
	return out
}

func appendZeros(v []frontend.Variable, n int) []frontend.Variable {
	for i := 0; i < n; i++ {
		v = append(v, frontend.Variable(0))
	}
	return v
}
