package polocolo

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete, named parameter set pinned to a curve's scalar
// field, mirroring ref polocolo/instances.py. m and R are spelled out from
// Table 1 (Table 7 for the tight variants, which drop the security margin:
// best attack >= 2^128 instead of 2^160); everything else is derived in
// NewParameters. The sponge uses one capacity element and a one-element digest
// (r = t-1, c = 1, d = 1), giving 128-bit security on these ~255-bit fields.
//
// Instance implements harness.Instance (metadata) and mode.SpongeHash. The
// reference defines no compression mode for Polocolo.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance: metadata for tagging targets ---

func (i Instance) Construction() string { return "polocolo" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Polocolo's sponge mode: a plain sponge with zero padding,
// add-to-start absorption, and an all-zero IV, matching the reference
// hash.py hash_sponge (IV = [0]*c).
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.ZeroIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("polocolo: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// generators of F_p^* for the two official fields (ref utils/field.py).
const (
	genBLS12 = 7
	genBN254 = 5
)

// Instances is the registry of all 20 official instances: recommended
// (Table 1) and tight (Table 7) parameters over both scalar fields. Names
// follow the reference instance names (and hence the KAT labels):
// polocolo-<field>-scalar-t<width>[-tight].
var Instances = []Instance{
	// BLS12-381 scalar field, recommended (Table 1)
	newInstance("polocolo-bls12-381-scalar-t3", ecc.BLS12_381, 3, 1024, 6, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t4", ecc.BLS12_381, 4, 512, 5, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t5", ecc.BLS12_381, 5, 128, 5, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t6", ecc.BLS12_381, 6, 64, 5, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t7", ecc.BLS12_381, 7, 32, 5, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t8", ecc.BLS12_381, 8, 32, 5, genBLS12),
	// BLS12-381 scalar field, tight (Table 7)
	newInstance("polocolo-bls12-381-scalar-t3-tight", ecc.BLS12_381, 3, 1024, 5, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t4-tight", ecc.BLS12_381, 4, 1024, 4, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t6-tight", ecc.BLS12_381, 6, 64, 4, genBLS12),
	newInstance("polocolo-bls12-381-scalar-t8-tight", ecc.BLS12_381, 8, 32, 4, genBLS12),
	// BN254 scalar field, recommended (Table 1)
	newInstance("polocolo-bn254-scalar-t3", ecc.BN254, 3, 1024, 6, genBN254),
	newInstance("polocolo-bn254-scalar-t4", ecc.BN254, 4, 512, 5, genBN254),
	newInstance("polocolo-bn254-scalar-t5", ecc.BN254, 5, 128, 5, genBN254),
	newInstance("polocolo-bn254-scalar-t6", ecc.BN254, 6, 64, 5, genBN254),
	newInstance("polocolo-bn254-scalar-t7", ecc.BN254, 7, 32, 5, genBN254),
	newInstance("polocolo-bn254-scalar-t8", ecc.BN254, 8, 32, 5, genBN254),
	// BN254 scalar field, tight (Table 7)
	newInstance("polocolo-bn254-scalar-t3-tight", ecc.BN254, 3, 1024, 5, genBN254),
	newInstance("polocolo-bn254-scalar-t4-tight", ecc.BN254, 4, 1024, 4, genBN254),
	newInstance("polocolo-bn254-scalar-t6-tight", ecc.BN254, 6, 64, 4, genBN254),
	newInstance("polocolo-bn254-scalar-t8-tight", ecc.BN254, 8, 32, 4, genBN254),
}

// newInstance constructs an Instance, panicking on invalid parameters —
// Instances is static configuration, so a bad entry is a programming error.
// The field label of the round-constant seed is fixed per curve ("BLS12" /
// "BN254"); the sponge shape is always r = t-1, c = 1, d = 1.
func newInstance(name string, curve ecc.ID, t, m, rounds int, g int64) Instance {
	label := "BN254"
	if curve == ecc.BLS12_381 {
		label = "BLS12"
	}
	params, err := NewParameters(curve.ScalarField(), t, m, rounds, g, label, t-1, 1, 1)
	if err != nil {
		panic(fmt.Sprintf("polocolo: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
