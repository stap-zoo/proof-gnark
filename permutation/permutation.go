// Package permutation defines the core contract every hash construction plugs
// into: the [Permutation] interface — an in-circuit, fixed-width permutation
// over field elements. Every construction (S-GMiMC, Rescue, Anemoi, ...)
// implements it, which lets the modes of operation (package mode) and the
// comparison harness treat all of them uniformly. [PermutationCircuit] wraps a
// Permutation for permutation-level testing and benchmarking.
package permutation

import "github.com/consensys/gnark/frontend"

// Permutation is a fixed-width permutation evaluated inside a gnark circuit.
//
// Implementations operate on a state of Width() field elements, expressed as
// []frontend.Variable, using the constraint-emitting api. They must NOT read or
// branch on concrete values (those are unknown at circuit-definition time);
// they only describe the arithmetic relations between wires.
//
// The concrete field (BN254, BLS12-381, ...) is fixed when the circuit is
// compiled, not here. Anything that depends on the field — round constants, MDS
// matrices — belongs in each construction's params, expressed as *big.Int so it
// can be assigned into a frontend.Variable for any field.
type Permutation interface {
	// Permute applies the permutation to state and returns the new state.
	// len(state) must equal Width(). Implementations may return a new slice
	// or mutate and return state; callers must use the returned slice.
	Permute(api frontend.API, state []frontend.Variable) []frontend.Variable

	// Width is the number of field elements in the state.
	Width() int
}
