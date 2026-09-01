package gmimc2

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete, named parameter set pinned to a specific curve's
// scalar field. Each entry supplies only the numbers the specification chooses
// (t, alpha, R, and the sponge split r/c/d); Parameters derives the rest.
//
// Instance implements harness.Instance (metadata) and mode.SpongeHash. It is
// deliberately not a mode.Compressor: over ~256-bit fields the specification
// proposes sponge mode only, because compression with d=1 admits a univariate
// Meet-in-the-Middle attack and would need roughly twice the rounds — to digest
// t elements in one call it is cheaper to run a sponge of width t+1. Compression
// mode is proposed only for the 32/64-bit instances, which are out of scope here
// (see Instances).
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance: metadata for tagging targets ---

func (i Instance) Construction() string { return "gmimc2" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash: the sponge GMiMC2 is eligible for ---

// Sponges returns GMiMC2's sponge mode: a plain sponge with zero padding,
// add-to-start absorption, and a length-encoding IV — [len(input), 0, ..., 0]
// over the capacity. This is the reference framework's SpongeLE over a
// fixed-length input, the same mode GMiMC uses, and what gmimc2_ref.py's
// hash_sponge implements.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.LengthIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("gmimc2: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// rounds is tab:gmimc2-rounds subtable (a) — log2(q) ~ 256, c = d = 1, 128-bit
// security — indexed [t][alpha]. Every entry is the smallest multiple of t
// exceeding the specification's bound for that configuration; over 256-bit fields
// the dominating attack is univariate polynomial system solving.
//
// The other two subtables are not here because their fields are not: they give
// round numbers for log2(q) ~ 64 (Goldilocks, t in {8,12}, c = d = 4) and ~32
// (Mersenne31, t in {16,24}, c = d = 8), and this repo implements only over the
// BN254 and BLS12-381 scalar fields — see README, "Fields". Those are also the
// instances that use compression mode and a non-identity M_IO, so both are absent
// with them.
//
// alpha = 8 is the specification's bold row at every width: the best tradeoff
// between plain and proof-system performance. Note how flat the round count is in
// t compared with GMiMC's — the erf criterion there grew steeply with the branch
// count (t=3 needs R=228 at alpha=5), which is what makes wide GMiMC2 practical.
var rounds = map[int]map[int]int{
	3: {2: 168, 4: 90, 8: 63},
	4: {2: 172, 4: 92, 8: 64},
	5: {2: 175, 4: 95, 8: 75},
	8: {2: 184, 4: 104, 8: 80},
}

// Instances is the registry of concrete instances.
// Naming convention: gmimc2-<field>-t<width>-a<alpha>.
//
// All twelve in-scope configurations over BLS12-381, plus the four bold alpha = 8
// ones over BN254: the repo's usual split, where BLS12-381 carries every width and
// BN254 twins the recommended row. Every instance has c = d = 1, hence r = t-1 and
// M_IO = identity.
//
// Both the round numbers and the vectors that pin them come from outside the
// Python reference framework, which has no GMiMC2: gnark-hashes/gmimc2_ref.py is
// the oracle, and `python3 gmimc2_ref.py --check gmimc2/testdata/vectors.json`
// re-derives every vector from the table above. Change a round number here and
// there, then regenerate — the KATs are what catch a mismatch.
var Instances = buildInstances()

func buildInstances() []Instance {
	var out []Instance
	for _, f := range []struct {
		label string
		curve ecc.ID
		alpha []int // which alpha columns this field carries
	}{
		{"bls12", ecc.BLS12_381, []int{2, 4, 8}},
		{"bn254", ecc.BN254, []int{8}},
	} {
		for _, t := range []int{3, 4, 5, 8} {
			for _, alpha := range f.alpha {
				name := fmt.Sprintf("gmimc2-%s-t%d-a%d", f.label, t, alpha)
				out = append(out, newInstance(name, f.curve, t, rounds[t][alpha], alpha,
					t-1 /*r*/, 1 /*c*/, 1 /*d*/))
			}
		}
	}
	return out
}

// newInstance constructs an Instance, panicking if the parameters are invalid —
// Instances is static configuration validated at process start, so a bad entry is
// a programming error, not a recoverable runtime condition.
func newInstance(name string, curve ecc.ID, t, rounds, alpha, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, rounds, alpha, rate, capacity, digest)
	if err != nil {
		panic(fmt.Sprintf("gmimc2: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
