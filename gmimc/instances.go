package gmimc

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete, named parameter set pinned to a specific curve's
// scalar field, mirroring instances.py. Each entry supplies only the chosen
// numbers (t, R, r, c, d); Parameters derives the rest.
//
// Instance implements harness.Instance (metadata) and mode.SpongeHash (it is
// eligible for a sponge). It is not a mode.Compressor, so it declares no
// compression targets.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance: metadata for tagging targets ---

func (i Instance) Construction() string { return "gmimc" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash: the sponge(s) GMiMC is eligible for ---

// Sponges returns GMiMC's sponge mode: a plain sponge with zero padding,
// add-to-start absorption, and a length-encoding IV — [len(input), 0, ..., 0]
// over the capacity — matching the Python reference (hash.py hash_sponge, whose
// IV is [to_field(len(data))] + [0]*(c-1)).
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.LengthIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("gmimc: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// Instances is the registry of concrete instances.
// Naming convention: gmimc-<field>-t<width>.
//
// t=4 needs nothing derived that t=3 did not: the linear layer is the same cyclic
// shift at any width and the round constants are one per round. Only R changes, and
// it is not transferable between widths — the erf construction's round criterion
// grows with the branch count, so 228 is a t=3 number. R=231 is the updated t=4
// count at alpha=5, which is the exponent defaultAlpha derives for both scalar
// fields (the smallest >= 2 coprime to p-1: 3 divides p-1 in both). At alpha=3 it
// would be 334. The reference has no t=4 instance and its own round derivation is
// unimplemented, so these vectors come from tools/export_kat_t4.py.
var Instances = []Instance{
	newInstance("gmimc-bn254-t3", ecc.BN254, 3 /*t*/, 228 /*R*/, 2 /*r*/, 1 /*c*/, 1 /*d*/),
	newInstance("gmimc-bls12-t3", ecc.BLS12_381, 3, 228, 2, 1, 1),
	newInstance("gmimc-bls12-t4", ecc.BLS12_381, 4, 231, 3, 1, 1),
}

// newInstance constructs an Instance, panicking if the parameters are invalid —
// Instances is static configuration validated at process start, so a bad entry
// is a programming error, not a recoverable runtime condition.
func newInstance(name string, curve ecc.ID, t, rounds, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, rounds, rate, capacity, digest)
	if err != nil {
		panic(fmt.Sprintf("gmimc: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
