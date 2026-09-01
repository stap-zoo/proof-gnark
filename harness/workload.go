// The job layer: a whole circuit of many hash calls, rather than the single call
// the cost table measures. Package doc: harness.go.
//
// The per-call numbers answer "what does this permutation cost". They do not answer
// "what does it cost to prove something", and the difference is not just
// multiplication:
//
//   - **Setup cost amortizes.** A lookup-based design pays one constraint per table
//     entry, once per circuit (see [algebra.Shared]). Skyscraper's first permutation
//     costs 1025 R1CS constraints and every later one ~340, so a job with thousands of
//     calls approaches a third of the per-call figure, while designs with no table are
//     unaffected — a single-call table ranks the lookup designs pessimistically, by a
//     factor only a real job reveals.
//   - **The mode changes the call count.** Absorbing n elements at rate r is
//     ceil(n/r) permutations, so a wider instance does more per call.
package harness

import (
	"fmt"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/mode"
)

// Message is the case for hashing nbElements field elements with the instance's
// sponge — the long-data workload. The element count is the circuit's input size,
// so the sponge absorbs ceil(nbElements/rate) blocks and calls the permutation
// that many times.
//
// An instance with several sponges uses the first it declares (only Anemoi and
// Skyscraper declare a non-plain one, and they declare exactly one).
func Message(inst Instance, nbElements int) (Case, error) {
	if nbElements < 1 {
		return Case{}, fmt.Errorf("workload: nbElements must be >= 1, got %d", nbElements)
	}
	sh, ok := inst.(mode.SpongeHash)
	if !ok {
		return Case{}, fmt.Errorf("workload: %s has no sponge", inst.InstanceName())
	}
	sponges := sh.Sponges()
	if len(sponges) == 0 {
		return Case{}, fmt.Errorf("workload: %s declares no sponge", inst.InstanceName())
	}
	s := sponges[0]
	return Case{
		Name:       fmt.Sprintf("%s/%s/message-%d", inst.InstanceName(), s.Name(), nbElements),
		Field:      inst.Field(),
		InputSize:  nbElements,
		OutputSize: s.DigestSize(),
		Build: func(in, out []frontend.Variable) frontend.Circuit {
			c := mode.NewSpongeCircuit(s, len(in))
			copy(c.Input, in)
			copy(c.Output, out)
			return c
		},
		Eval: s.Hash,
	}, nil
}

// Merkle is the case for computing the root of a 2:1 Merkle tree over nbLeaves
// leaves, which is nbLeaves-1 hash calls. nbLeaves must be a power of two.
//
// A "node" is [Node.Size] field elements, and the tree is 2:1 in *nodes*: two
// children make one parent. That is the only formulation that composes for every
// construction here, because a digest is not always one element — Rescue-Prime and
// Arion squeeze two — so a parent absorbs 2*Size elements and emits Size.
//
// The node comes from the caller ([NewNode] for the mode the construction is proven
// for, [NewSpongeNode] to force the sponge) because which one it is has to be
// reported alongside the number.
func Merkle(inst Instance, node Node, nbLeaves int) (Case, error) {
	if nbLeaves < 2 || nbLeaves&(nbLeaves-1) != 0 {
		return Case{}, fmt.Errorf("workload: nbLeaves must be a power of two >= 2, got %d", nbLeaves)
	}
	sz := node.Size()
	return Case{
		Name:       fmt.Sprintf("%s/%s/merkle-%d", inst.InstanceName(), node.Name(), nbLeaves),
		Field:      inst.Field(),
		InputSize:  nbLeaves * sz,
		OutputSize: sz,
		Build: func(in, out []frontend.Variable) frontend.Circuit {
			c := newMerkleCircuit(node, nbLeaves)
			copy(c.Leaves, in)
			copy(c.Root, out)
			return c
		},
		Eval: func(api frontend.API, in []frontend.Variable) []frontend.Variable {
			return foldToRoot(api, node, in)
		},
	}, nil
}

// ---------------------------------------------------------------------------
// The tree node
// ---------------------------------------------------------------------------

// Node is one 2:1 step of a Merkle tree: it maps two child nodes to their parent.
// It exists so the tree is written once and works over either a compression
// function or a sponge, whichever the construction is actually proven for.
type Node interface {
	// Combine returns the parent of two Size()-element children.
	Combine(api frontend.API, left, right []frontend.Variable) []frontend.Variable
	// Size is the number of field elements in a node.
	Size() int
	// Name identifies the mode used, for the case name and CSV.
	Name() string
}

// NewNode picks the tree node for an instance: its arity-2 compression if it
// declares one, otherwise its sponge.
//
// A compression is the right primitive for a Merkle tree and is what the
// compression-mode designs are proven for (Anemoi's Jive-2, Skyscraper's
// Davies-Meyer), so it wins when available. The sponge is the fallback, and for
// the six constructions that declare no compression it is the only honest choice —
// the alternative would be to invent a compression the design was never analyzed
// in, which is exactly what mode eligibility exists to prevent.
func NewNode(inst Instance) (Node, error) {
	if c, ok := inst.(mode.Compressor); ok {
		for _, comp := range c.Compressions() {
			if comp.Arity() == 2 {
				return &compressionNode{comp: comp}, nil
			}
		}
	}
	return NewSpongeNode(inst)
}

// NewSpongeNode is the sponge tree node, even for an instance that also declares a
// compression. It is what gives a same-mode column across every construction: six
// of the eleven have no compression, so a comparison that uses the better primitive
// where it exists is comparing two different things, and this is the other side of
// that trade.
func NewSpongeNode(inst Instance) (Node, error) {
	sh, ok := inst.(mode.SpongeHash)
	if !ok {
		return nil, fmt.Errorf("workload: %s has no sponge", inst.InstanceName())
	}
	sponges := sh.Sponges()
	if len(sponges) == 0 {
		return nil, fmt.Errorf("workload: %s declares no sponge", inst.InstanceName())
	}
	return &spongeNode{s: sponges[0]}, nil
}

// compressionNode is a tree node backed by an arity-2 compression: the natural
// 2:1 primitive, one call per node.
type compressionNode struct{ comp mode.Compression }

func (n *compressionNode) Size() int    { return n.comp.BlockSize() }
func (n *compressionNode) Name() string { return n.comp.Name() }

func (n *compressionNode) Combine(api frontend.API, left, right []frontend.Variable) []frontend.Variable {
	return n.comp.Compress(api, [][]frontend.Variable{left, right})
}

// spongeNode is a tree node backed by a sponge: absorb both children, squeeze the
// parent. Costlier than a compression at the same width (padding and a full
// absorption per node), which is part of what the comparison shows.
type spongeNode struct{ s mode.Sponge }

func (n *spongeNode) Size() int    { return n.s.DigestSize() }
func (n *spongeNode) Name() string { return n.s.Name() }

func (n *spongeNode) Combine(api frontend.API, left, right []frontend.Variable) []frontend.Variable {
	in := make([]frontend.Variable, 0, len(left)+len(right))
	in = append(in, left...)
	in = append(in, right...)
	return n.s.Hash(api, in)
}

// ---------------------------------------------------------------------------
// The circuit
// ---------------------------------------------------------------------------

// merkleCircuit computes the root of a 2:1 tree over Leaves and asserts it equals
// Root. Leaves holds nbLeaves nodes concatenated, Size() elements each.
//
// The node is excluded from the witness with `gnark:"-"`, like every other
// non-witness config in this repo.
type merkleCircuit struct {
	Leaves []frontend.Variable
	Root   []frontend.Variable `gnark:",public"`

	Node Node `gnark:"-"`
}

// newMerkleCircuit returns a merkleCircuit sized for nbLeaves leaves, as required
// for compilation.
func newMerkleCircuit(node Node, nbLeaves int) *merkleCircuit {
	sz := node.Size()
	return &merkleCircuit{
		Leaves: make([]frontend.Variable, nbLeaves*sz),
		Root:   make([]frontend.Variable, sz),
		Node:   node,
	}
}

// Define implements frontend.Circuit.
func (c *merkleCircuit) Define(api frontend.API) error {
	got := foldToRoot(api, c.Node, c.Leaves)
	for i := range c.Root {
		api.AssertIsEqual(got[i], c.Root[i])
	}
	return nil
}

// foldToRoot folds the leaves up the tree, level by level, and returns the root. It is
// shared by Define and the case's Eval so the circuit and the native evaluation
// can never disagree about the tree's shape.
func foldToRoot(api frontend.API, node Node, leaves []frontend.Variable) []frontend.Variable {
	sz := node.Size()
	level := make([][]frontend.Variable, len(leaves)/sz)
	for i := range level {
		level[i] = leaves[i*sz : (i+1)*sz]
	}
	for len(level) > 1 {
		next := make([][]frontend.Variable, len(level)/2)
		for i := range next {
			next[i] = node.Combine(api, level[2*i], level[2*i+1])
		}
		level = next
	}
	return level[0]
}
