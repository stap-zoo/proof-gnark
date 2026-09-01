package anemoi

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Anemoi parameter set pinned to a curve's scalar field.
// Mirrors anemoi/instances.py over the BN254 and BLS12-381 scalar fields.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "anemoi" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Anemoi's sponge mode: the Hirose variant with add-to-start
// absorption and the length-dependent padding/sigma of the reference.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewHiroseSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.AddToStart)
	if err != nil {
		panic(fmt.Sprintf("anemoi: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// --- mode.Compressor ---

// Compressions returns Anemoi's Jive_2 compression (the two lanes x and y of the
// width-t state compressed to l elements).
func (i Instance) Compressions() []mode.Compression {
	j, err := mode.NewJive(i.Permutation(), 2)
	if err != nil {
		panic(fmt.Sprintf("anemoi: instance %q jive: %v", i.Name, err))
	}
	return []mode.Compression{j}
}

// Instances is the registry of concrete instances (BN254 + BLS12-381 scalar).
// Naming (matching the reference KAT labels): anemoi-<field>-scalar-t<width>.
// r/c/d are taken from anemoi/instances.py; alpha=5 and g are the field's
// power-map exponent and generator (utils/field.py).
//
// R at t=4 is 20, not the reference's 14: the newer analysis of that width raises
// it, and at alpha=3 it would be 25 (no instance here uses alpha=3 — 3 divides p-1
// in both fields, so x^3 is not a permutation). The reference is left alone, so
// these two instances no longer match ANEMOI_*_SCALAR_T4 and their vectors come
// from tools/export_kat_t4.py instead of export_kat.py. t=2 and t=6 are unchanged
// and still the reference's own.
var Instances = []Instance{
	newInstance("anemoi-bn254-scalar-t2", ecc.BN254, 1 /*l*/, 5 /*alpha*/, 5 /*g*/, 21 /*R*/, 1 /*r*/, 1 /*c*/, 1 /*d*/),
	newInstance("anemoi-bn254-scalar-t4", ecc.BN254, 2, 5, 5, 20, 3, 1, 1),
	newInstance("anemoi-bn254-scalar-t6", ecc.BN254, 3, 5, 5, 12, 5, 1, 1),
	newInstance("anemoi-bls12-381-scalar-t2", ecc.BLS12_381, 1, 5, 7, 21, 1, 1, 1),
	newInstance("anemoi-bls12-381-scalar-t4", ecc.BLS12_381, 2, 5, 7, 20, 3, 1, 1),
	newInstance("anemoi-bls12-381-scalar-t6", ecc.BLS12_381, 3, 5, 7, 12, 5, 1, 1),
}

func newInstance(name string, curve ecc.ID, l, alpha int, g int64, R, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), l, alpha, big.NewInt(g), R, rate, capacity, digest)
	if err != nil {
		panic(fmt.Sprintf("anemoi: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
	_ mode.Compressor  = Instance{}
)
