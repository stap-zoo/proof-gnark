// The registry layer: instances become tagged, runnable targets that the tests
// and the benchmark command both iterate and filter. Package doc: harness.go.
package harness

import (
	"fmt"
)

// Kind is the category of a target: the bare permutation, a sponge hash, or a
// compression function.
type Kind int

const (
	Permutation Kind = iota
	Sponge
	Compression
)

func (k Kind) String() string {
	switch k {
	case Permutation:
		return "permutation"
	case Sponge:
		return "sponge"
	case Compression:
		return "compression"
	default:
		return "unknown"
	}
}

// ParseKind maps a kind name (as used on the command line) to a Kind.
func ParseKind(s string) (Kind, error) {
	switch s {
	case "permutation":
		return Permutation, nil
	case "sponge":
		return Sponge, nil
	case "compression":
		return Compression, nil
	default:
		return 0, fmt.Errorf("harness: unknown kind %q (want permutation|sponge|compression)", s)
	}
}

// Target is one runnable (instance × mode): metadata labels plus the Case
// that knows how to build its circuit.
type Target struct {
	Construction string // e.g. "gmimc"
	Instance     string // e.g. "gmimc-bn254-t3"
	FieldName    string // e.g. "bn254"
	Kind         Kind
	Mode         string // e.g. "permutation", "sponge", "jive-2"
	Case         Case

	// Width is the permutation's state size in field elements — t, the number
	// every comparison of these designs is indexed by.
	//
	// It comes from the permutation itself and never from the instance name,
	// because the references do not agree on how to name it: Skyscraper's own
	// label for its state of two elements is "n1", so `skyscraper-bn254-n1` and
	// `poseidon-bn254-t2` are the same width under two spellings. Reporting the
	// number the permutation actually has is the only way a width column can be
	// read across constructions.
	Width int

	// VectorInstance is the instance name this target's reference vectors are
	// labeled with. It equals Instance for every instance that generated its own
	// vectors, and differs only for an alternative *arithmetization* of a hash
	// already registered — a circuit computing the identical function, which
	// therefore satisfies the original instance's KATs under the original's name
	// (Skyscraper's 16-bit-word Bar; see VectorSource). Vector matching goes
	// through this field, never through Instance, so the label a target carries
	// in the CSV is free to differ from the label its KATs carry.
	VectorInstance string
}

// MatchesVector reports whether v is a reference vector for this target: the
// vector's instance is the one this target's KATs are labeled with, and its mode
// is either the target's specific mode name ("jive-2") or the target's generic
// Kind ("compression").
//
// The kind fallback is there because the reference's known-answer generator labels
// vectors by kind only, which is how a specifically-named mode (Hirose sponge,
// Jive) picks up its vectors; for a target whose Mode already equals its Kind
// (permutation, the plain "sponge") the two coincide and matching is unchanged.
//
// This is the *one* definition of the match. RunTests runs what it selects, and
// the coverage guards in registry/ assert every target selects something and every
// in-scope vector is selected by someone — all three off this method, so they
// cannot drift apart.
func (t Target) MatchesVector(v Vector) bool {
	if v.Instance != t.VectorInstance {
		return false
	}
	return v.Mode == t.Mode || v.Mode == t.Kind.String()
}
