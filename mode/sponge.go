package mode

import (
	"fmt"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// PlainSponge is the standard sponge (hash_sponge): initialise the state as
// [0]*rate || IV, absorb rate-sized blocks (each followed by a permutation),
// then squeeze. It is configured with the pluggable IV, Absorb and Padding
// pieces from parts.go.
type PlainSponge struct {
	perm     permutation.Permutation
	rate     int
	capacity int
	digest   int
	iv       IV
	absorb   Absorb
	pad      Padding
}

var _ Sponge = (*PlainSponge)(nil)

// NewPlainSponge builds a standard sponge. It requires rate+capacity == width.
func NewPlainSponge(perm permutation.Permutation, rate, capacity, digest int, iv IV, absorb Absorb, pad Padding) (*PlainSponge, error) {
	if err := checkSpongeDims(perm, rate, capacity, digest); err != nil {
		return nil, err
	}
	if iv == nil || absorb == nil || pad == nil {
		return nil, fmt.Errorf("mode: sponge iv/absorb/pad must be non-nil")
	}
	return &PlainSponge{perm, rate, capacity, digest, iv, absorb, pad}, nil
}

func (s *PlainSponge) Name() string    { return "sponge" }
func (s *PlainSponge) Rate() int       { return s.rate }
func (s *PlainSponge) DigestSize() int { return s.digest }

func (s *PlainSponge) Hash(api frontend.API, input []frontend.Variable) []frontend.Variable {
	state := initState(s.rate, s.iv(len(input), s.capacity, s.digest))
	for _, blk := range blocksOf(s.pad(input, s.rate), s.rate) {
		state = s.absorb(api, state, blk)
		state = s.perm.Permute(api, state)
	}
	return squeeze(api, s.perm, state, s.rate, s.digest)
}

// HiroseSponge is the Hirose sponge variant (hash_sponge_hirose) used by Anemoi:
// a zero IV, and a domain separator sigma added to the last state element after
// the final absorb-permutation, before squeezing.
//
// The padding and sigma are not fixed at construction: matching the reference
// Anemoi.hash_sponge, they are chosen from the input length (a compile-time
// quantity — the number of input wires, not a witness value). A rate-aligned,
// non-empty input is absorbed unpadded with sigma=1; every other input (including
// the empty one) is PadOne-padded with sigma=0.
type HiroseSponge struct {
	perm     permutation.Permutation
	rate     int
	capacity int
	digest   int
	absorb   Absorb
}

var _ Sponge = (*HiroseSponge)(nil)

// NewHiroseSponge builds a Hirose sponge. Padding and the domain separator are
// derived per input length inside Hash (see the type doc), so unlike the other
// sponges it takes neither a Padding nor a sigma.
func NewHiroseSponge(perm permutation.Permutation, rate, capacity, digest int, absorb Absorb) (*HiroseSponge, error) {
	if err := checkSpongeDims(perm, rate, capacity, digest); err != nil {
		return nil, err
	}
	if absorb == nil {
		return nil, fmt.Errorf("mode: sponge absorb must be non-nil")
	}
	return &HiroseSponge{perm, rate, capacity, digest, absorb}, nil
}

func (s *HiroseSponge) Name() string    { return "sponge-hirose" }
func (s *HiroseSponge) Rate() int       { return s.rate }
func (s *HiroseSponge) DigestSize() int { return s.digest }

func (s *HiroseSponge) Hash(api frontend.API, input []frontend.Variable) []frontend.Variable {
	// Domain separation: sigma=1 for a rate-aligned non-empty message (absorbed
	// as-is), sigma=0 otherwise (PadOne applied). Matches Anemoi.hash_sponge.
	var padded []frontend.Variable
	var sigma int
	if len(input) != 0 && len(input)%s.rate == 0 {
		padded, sigma = clone(input), 1
	} else {
		padded, sigma = PadOne(input, s.rate), 0
	}

	state := zeros(s.perm.Width())
	for _, blk := range blocksOf(padded, s.rate) {
		state = s.absorb(api, state, blk)
		state = s.perm.Permute(api, state)
	}
	last := len(state) - 1
	state[last] = api.Add(state[last], sigma)
	return squeeze(api, s.perm, state, s.rate, s.digest)
}

// ---------------------------------------------------------------------------
// shared sponge helpers
// ---------------------------------------------------------------------------

func checkSpongeDims(perm permutation.Permutation, rate, capacity, digest int) error {
	if rate < 1 {
		return fmt.Errorf("mode: sponge rate must be >= 1, got %d", rate)
	}
	if capacity < 1 {
		return fmt.Errorf("mode: sponge capacity must be >= 1, got %d", capacity)
	}
	if digest < 1 {
		return fmt.Errorf("mode: sponge digest must be >= 1, got %d", digest)
	}
	if rate+capacity != perm.Width() {
		return fmt.Errorf("mode: sponge requires rate+capacity == width, got r=%d c=%d width=%d", rate, capacity, perm.Width())
	}
	return nil
}

// initState builds [0]*rate || iv.
func initState(rate int, iv []frontend.Variable) []frontend.Variable {
	return append(zeros(rate), iv...)
}

// blocksOf splits rate-aligned data into rate-sized blocks.
func blocksOf(data []frontend.Variable, rate int) [][]frontend.Variable {
	blocks := make([][]frontend.Variable, 0, len(data)/rate)
	for i := 0; i < len(data); i += rate {
		blocks = append(blocks, data[i:i+rate])
	}
	return blocks
}

// squeeze reads rate elements at a time, re-permuting between reads, until digest
// elements have been produced.
func squeeze(api frontend.API, perm permutation.Permutation, state []frontend.Variable, rate, digest int) []frontend.Variable {
	out := make([]frontend.Variable, 0, digest)
	for len(out) < digest {
		out = append(out, state[:rate]...)
		if len(out) < digest {
			state = perm.Permute(api, state)
		}
	}
	return out[:digest]
}
