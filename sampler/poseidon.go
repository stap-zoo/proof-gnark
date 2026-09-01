package sampler

import "math/big"

// Grain LFSR settings shared by the Poseidon family (Poseidon and Poseidon2):
// the 80-bit register and feedback taps of Appendix E of
// https://eprint.iacr.org/2019/458, warmed up by two full passes (2*80).
// Mirrors GRAIN_STATE_SIZE / GRAIN_TAPS in ref/hades/params.py.
const (
	poseidonGrainStateSize = 80
	poseidonGrainWarmup    = 2 * poseidonGrainStateSize
)

var poseidonGrainTaps = []int{0, 13, 23, 38, 51, 62}

// NewPoseidonGrain builds the Grain LFSR for a Poseidon-family instance. The
// 80-bit seed is the nothing-up-my-sleeve encoding of the instance parameters
// (field code, S-box marker, field bit length, t, R_ext, R_int), MSB-first,
// padded with ones — the _ORIGINAL layout of ref/hades/params.py. The only
// version-dependent field is the S-box marker: the "isec"/horizenlabs reference
// forces it to 1 even for a power map, while "circom" uses the plain rule (0 for
// a power map). Both Poseidon and Poseidon2 seed identically; they diverge only
// in what they draw off the resulting stream.
//
// The returned generator is self-shrinking (the hadeshash layout) and starts in
// rejection ("bitmask") sampling; callers switch it to reduction with SwitchToMod
// where the reference does (Poseidon, for its MDS material).
func NewPoseidonGrain(field *big.Int, t, rExt, rInt int, version string) *Grain {
	sbox := 0 // power-map marker (reference _sbox rule for alpha != -1)
	if version == "isec" || version == "horizenlabs" {
		sbox = 1 // isec/horizenlabs deviation: force the marker to 1
	}
	fields := []struct{ value, width int }{
		{1, 2},               // field code: 1 = prime field
		{sbox, 4},            // S-box marker
		{field.BitLen(), 12}, // n = field bit length
		{t, 12},              // state size
		{rExt, 10},           // R_F (external rounds)
		{rInt, 10},           // R_P (internal rounds)
	}
	seed := make([]byte, 0, poseidonGrainStateSize)
	for _, f := range fields {
		for i := 0; i < f.width; i++ {
			seed = append(seed, byte((f.value>>(f.width-1-i))&1)) // MSB-first
		}
	}
	for len(seed) < poseidonGrainStateSize {
		seed = append(seed, 1) // pad with ones
	}
	return NewGrain(seed, field, poseidonGrainTaps, poseidonGrainWarmup, true /*shrink*/)
}
