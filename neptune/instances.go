package neptune

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Neptune parameter set pinned to a curve's scalar
// field. Mirrors hades/instances.py over the BN254 and BLS12-381 scalar fields.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "neptune" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Neptune's sponge mode: a plain sponge with zero padding,
// add-to-start absorption and a length-encoding IV, matching the reference
// hades hash_sponge (IV = [len(data)] + [0]*(c-1)).
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.LengthIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("neptune: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// Instances is the registry of concrete instances (bn254 + bls12-381).
// Naming convention: neptune-<field>-t<width>.
var Instances = []Instance{
	newInstance("neptune-bn254-t4", ecc.BN254, 4 /*t*/, 5 /*alpha*/, 6 /*R_ext*/, 68 /*R_int*/, 3 /*r*/, 1 /*c*/, 1 /*d*/),
	newInstance("neptune-bls12-t4", ecc.BLS12_381, 4, 5, 6, 68, 3, 1, 1),
	newInstance("neptune-bls12-t2", ecc.BLS12_381, 2, 5, 8, 56, 1, 1, 1),
}

func newInstance(name string, curve ecc.ID, t, alpha, rExt, rInt, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, alpha, rExt, rInt, rate, capacity, digest)
	if err != nil {
		panic(fmt.Sprintf("neptune: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
