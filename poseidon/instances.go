package poseidon

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Poseidon parameter set pinned to a curve's scalar
// field. Mirrors hades/instances.py over the BN254 and BLS12-381 scalar fields.
type Instance struct {
	Name    string
	Curve   ecc.ID
	Version string
	Params  *Parameters
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "poseidon" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Poseidon's sponge mode: a plain sponge with zero padding,
// add-to-start absorption and a length-encoding IV, matching the reference
// hades hash_sponge (IV = [len(data)] + [0]*(c-1)). Poseidon defines no
// compression mode of its own, so none is declared.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.LengthIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("poseidon: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// Instances is the registry of concrete instances (bn254 + bls12-381 scalar).
// Naming convention: poseidon-<field>-t<width>, matching the reference KAT labels.
//
// The t=4 instance runs Poseidon2's round numbers for that width (R_ext=8,
// R_int=56, the same pair poseidon2-bls12-t4 carries): the two share the Hades
// round structure and the same S-box, so the round criterion does not distinguish
// them, and giving Poseidon a different count would make the head-to-head a
// comparison of round budgets instead of arithmetizations. Everything else derives:
// the Cauchy MDS comes off the Grain stream at any width and so does the constant
// grid, which is why this instance needs no code beyond its own line (its vectors
// come from tools/export_kat_t4.py, the reference having no t=4 instance).
var Instances = []Instance{
	newInstance("poseidon-bn254-t3", ecc.BN254, "isec", 3 /*t*/, 5 /*alpha*/, 8 /*R_ext*/, 57 /*R_int*/, 2 /*r*/, 1 /*c*/, 1 /*d*/),
	newInstance("poseidon-bls12-t3", ecc.BLS12_381, "isec", 3, 5, 8, 57, 2, 1, 1),
	newInstance("poseidon-bls12-t2", ecc.BLS12_381, "circom", 2, 5, 8, 56, 1, 1, 1),
	newInstance("poseidon-bls12-t4", ecc.BLS12_381, "isec", 4, 5, 8, 56, 3, 1, 1),
}

func newInstance(name string, curve ecc.ID, version string, t, alpha, rExt, rInt, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, alpha, rExt, rInt, rate, capacity, digest, version)
	if err != nil {
		panic(fmt.Sprintf("poseidon: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Version: version, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
