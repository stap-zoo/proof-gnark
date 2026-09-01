package sampler

import "math/big"

// Grain is the Grain-LFSR field-element sampler that derives Poseidon and
// Poseidon2 parameters deterministically from an 80-bit seed. It mirrors
// ref/utils/sampler.py LFSRFieldElementSampler: an LFSR register clocked with
// feedback taps, from which candidate field elements are read MSB-first and
// mapped to [0, p) either by rejection ("bitmask") or by reduction ("mod").
//
// Poseidon draws its round constants under rejection, then switches to
// reduction (SwitchToMod) for the Cauchy MDS material, all off this one stream;
// Poseidon2 uses rejection only. Candidate width is the field bit length in both
// modes (the reference's "naive" width is unused by the Hades family and omitted).
type Grain struct {
	field     *big.Int
	candBits  int    // bits collected per candidate (= field bit length)
	state     []byte // one bit (0/1) per register cell
	taps      []int  // feedback tap indices
	shrink    bool   // self-shrinking output extraction
	reduceMod bool   // "mod" (reduce) vs "bitmask" (reject)
}

// NewGrain builds a Grain LFSR seeded with seedBits (one bit per entry), clocks
// it through warmup discarded steps, and starts in rejection ("bitmask") mode.
// shrink selects the self-shrinking output extraction (Poseidon's hadeshash
// layout uses shrink=true).
func NewGrain(seedBits []byte, field *big.Int, taps []int, warmup int, shrink bool) *Grain {
	g := &Grain{
		field:    field,
		candBits: field.BitLen(),
		state:    append([]byte(nil), seedBits...),
		taps:     taps,
		shrink:   shrink,
	}
	for i := 0; i < warmup; i++ {
		g.nextBit()
	}
	return g
}

// SwitchToMod switches from rejection to mod-reduction sampling, as Poseidon
// does between its round constants and its MDS material.
func (g *Grain) SwitchToMod() { g.reduceMod = true }

// Next returns the next field element in [0, p): under rejection, resample until
// the candidate is < p; under reduction, reduce mod p.
func (g *Grain) Next() *big.Int {
	for {
		v := g.draw()
		if g.reduceMod {
			return v.Mod(v, g.field)
		}
		if v.Cmp(g.field) < 0 {
			return v
		}
	}
}

// NextNonzero returns the next nonzero field element.
func (g *Grain) NextNonzero() *big.Int {
	for {
		if v := g.Next(); v.Sign() != 0 {
			return v
		}
	}
}

// Grid returns a rows x cols grid (row-major) of field elements.
func (g *Grain) Grid(rows, cols int) [][]*big.Int {
	out := make([][]*big.Int, rows)
	for i := range out {
		out[i] = make([]*big.Int, cols)
		for j := range out[i] {
			out[i][j] = g.Next()
		}
	}
	return out
}

// draw collects candBits output bits MSB-first into one candidate integer.
func (g *Grain) draw() *big.Int {
	v := new(big.Int)
	for i := 0; i < g.candBits; i++ {
		v.Lsh(v, 1)
		if g.nextOutputBit() == 1 {
			v.SetBit(v, 0, 1)
		}
	}
	return v
}

// nextOutputBit is one output bit: self-shrinking (read pairs (b0, b1), keep b1
// iff b0 == 1) when shrink is set, otherwise the raw clocked bit.
func (g *Grain) nextOutputBit() byte {
	if !g.shrink {
		return g.nextBit()
	}
	for {
		b0 := g.nextBit()
		b1 := g.nextBit()
		if b0 == 1 {
			return b1
		}
	}
}

// nextBit clocks the LFSR once: XOR the tapped cells, shift the register left,
// and append the feedback bit. Matches s = s[1:] + [new] in the reference.
func (g *Grain) nextBit() byte {
	var nb byte
	for _, i := range g.taps {
		nb ^= g.state[i]
	}
	copy(g.state, g.state[1:])
	g.state[len(g.state)-1] = nb
	return nb
}
