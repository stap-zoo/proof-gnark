// Package mode provides the modes of operation that turn a fixed-width
// permutation into a hash or compression function: sponges and compressions.
//
// Every mode is generic over permutation.Permutation (the only thing it needs is
// Permute + Width), so the concrete constructions (GMiMC, Anemoi, ...) never
// reimplement a mode. A construction instead declares which modes it is designed
// and proven for by having its instance implement SpongeHash and/or Compressor.
//
// Layout mirrors ref/utils/mode.py:
//   - compression.go — Davies-Meyer, Jive
//   - sponge.go       — plain and Hirose sponges
//   - parts.go        — the pluggable padding, absorb and IV pieces sponges use
//   - circuit.go      — generic SpongeCircuit / CompressionCircuit for test+bench
package mode

import "github.com/consensys/gnark/frontend"

// Sponge hashes a variable-length input to a fixed-size digest.
type Sponge interface {
	// Name identifies the mode variant (e.g. "sponge", "sponge-hirose").
	Name() string
	// Hash absorbs input and squeezes DigestSize() field elements.
	Hash(api frontend.API, input []frontend.Variable) []frontend.Variable
	Rate() int
	DigestSize() int
}

// Compression maps Arity() input blocks of BlockSize() elements to one block.
type Compression interface {
	Name() string
	Compress(api frontend.API, blocks [][]frontend.Variable) []frontend.Variable
	Arity() int     // number of input blocks
	BlockSize() int // elements per block (and per output)
}

// SpongeHash is implemented by a construction's instance when it is eligible for
// (i.e. designed and proven for) one or more sponge modes.
type SpongeHash interface {
	Sponges() []Sponge
}

// Compressor is implemented by a construction's instance when it is eligible for
// one or more compression modes.
type Compressor interface {
	Compressions() []Compression
}

// Sponges collects the sponge modes of every instance that is a SpongeHash,
// skipping those that are not. This is how the comparison harness gathers all
// sponges across constructions to benchmark them together.
func Sponges(instances ...any) []Sponge {
	var out []Sponge
	for _, inst := range instances {
		if sh, ok := inst.(SpongeHash); ok {
			out = append(out, sh.Sponges()...)
		}
	}
	return out
}

// Compressions collects the compression modes of every instance that is a
// Compressor, skipping those that are not.
func Compressions(instances ...any) []Compression {
	var out []Compression
	for _, inst := range instances {
		if c, ok := inst.(Compressor); ok {
			out = append(out, c.Compressions()...)
		}
	}
	return out
}
