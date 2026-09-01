package permutation

import "github.com/consensys/gnark/frontend"

// PermutationCircuit applies a bare Permutation to Input and asserts the full
// output state equals Output. It is what permutation-level test vectors
// (state -> state) run against, and lets the permutation be benchmarked
// independently of any mode of operation.
//
// gnark builds a circuit by reflecting over the struct: every frontend.Variable
// field (including inside slices) becomes a wire. Configuration that is NOT part
// of the witness — here the Permutation — must be excluded with the `gnark:"-"`
// tag, and slice lengths must be set before frontend.Compile (see the
// constructor).
type PermutationCircuit struct {
	Input  []frontend.Variable
	Output []frontend.Variable `gnark:",public"`

	Perm Permutation `gnark:"-"`
}

// NewPermutationCircuit returns a PermutationCircuit with Input and Output
// pre-allocated to the permutation width, as required for compilation.
func NewPermutationCircuit(p Permutation) *PermutationCircuit {
	return &PermutationCircuit{
		Input:  make([]frontend.Variable, p.Width()),
		Output: make([]frontend.Variable, p.Width()),
		Perm:   p,
	}
}

// Define implements frontend.Circuit.
func (c *PermutationCircuit) Define(api frontend.API) error {
	got := c.Perm.Permute(api, c.Input)
	for i := range c.Output {
		api.AssertIsEqual(got[i], c.Output[i])
	}
	return nil
}
