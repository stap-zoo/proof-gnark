package mode

import "github.com/consensys/gnark/frontend"

// SpongeCircuit hashes Input with a Sponge and asserts the digest equals Output.
// It is the mode-level analog of permutation.PermutationCircuit, letting any
// sponge be tested against vectors and benchmarked through the shared harness.
// Non-witness config (the Sponge) is excluded with `gnark:"-"`.
type SpongeCircuit struct {
	Input  []frontend.Variable
	Output []frontend.Variable `gnark:",public"`

	Sponge Sponge `gnark:"-"`
}

// NewSpongeCircuit returns a SpongeCircuit sized for nbInputs input elements and
// the sponge's digest size, as required for compilation.
func NewSpongeCircuit(s Sponge, nbInputs int) *SpongeCircuit {
	return &SpongeCircuit{
		Input:  make([]frontend.Variable, nbInputs),
		Output: make([]frontend.Variable, s.DigestSize()),
		Sponge: s,
	}
}

// Define implements frontend.Circuit.
func (c *SpongeCircuit) Define(api frontend.API) error {
	got := c.Sponge.Hash(api, c.Input)
	for i := range c.Output {
		api.AssertIsEqual(got[i], c.Output[i])
	}
	return nil
}

// CompressionCircuit compresses Input with a Compression and asserts the result
// equals Output. Input holds the Arity input blocks concatenated row-major
// (Arity*BlockSize elements); Output is one block (BlockSize elements).
type CompressionCircuit struct {
	Input  []frontend.Variable
	Output []frontend.Variable `gnark:",public"`

	Comp Compression `gnark:"-"`
}

// NewCompressionCircuit returns a CompressionCircuit sized for the compression's
// arity and block size, as required for compilation.
func NewCompressionCircuit(c Compression) *CompressionCircuit {
	return &CompressionCircuit{
		Input:  make([]frontend.Variable, c.Arity()*c.BlockSize()),
		Output: make([]frontend.Variable, c.BlockSize()),
		Comp:   c,
	}
}

// Define implements frontend.Circuit.
func (c *CompressionCircuit) Define(api frontend.API) error {
	m := c.Comp.BlockSize()
	blocks := make([][]frontend.Variable, c.Comp.Arity())
	for j := range blocks {
		blocks[j] = c.Input[j*m : (j+1)*m]
	}
	got := c.Comp.Compress(api, blocks)
	for i := range c.Output {
		api.AssertIsEqual(got[i], c.Output[i])
	}
	return nil
}
