package griffin

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Griffin parameter set pinned to a curve's scalar
// field. Mirrors the GRIFFIN_* instances in griffin/instances.py, restricted to
// the BN254 and BLS12-381 scalar fields per the project's prime-field policy. The
// reference's other instances (ST, Goldilocks) are over excluded fields and have
// no Go target.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "griffin" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Griffin's sponge: zero padding, add-to-start absorption, and a
// length-encoding IV — [len(input), 0, ..., 0] over the capacity — matching the
// reference Griffin.hash_sponge (IV = [to_field(len(data))] + [0]*(c-1)). Every
// instance here has 2*d != t, so none declares a compression mode.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.LengthIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("griffin: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// Instances is the registry of concrete instances. The reference defines t=3
// instances over both scalar fields (GRIFFIN_BN254_T3, GRIFFIN_BLS12_T3): t=3,
// alpha=5, R=12, r/c/d = 2/1/1 (kappa=128). The instance-name strings match the
// reference KAT labels griffin-bn254-t3 / griffin-bls12-t3.
//
// The t=4 instance is the second of Griffin's two admissible width families
// (t = 3 or a multiple of 4), so it needs no new design decision — only the M4
// mixing matrix initM already builds. Its round count is R=20, not the 12 the
// reference's own criterion returns for the width (eprint 2022/403 Section 5.2, the
// Groebner bound with a 20% margin): the newer analysis of t=4 supersedes that. Note
// that R is not only a loop bound here — the round constants and the quadratic
// coefficients come off one SHAKE stream in that order (deriveConstants), so R also
// shifts the coefficients. Every vector for this instance is therefore specific to
// R=20 and comes from tools/export_kat_t4.py, not the reference's export_kat.py; t=3
// is unchanged and still the reference's own.
var Instances = []Instance{
	newInstance("griffin-bn254-t3", ecc.BN254, 3 /*t*/, 12 /*R*/, 5 /*alpha*/, 2 /*r*/, 1 /*c*/, 1 /*d*/),
	newInstance("griffin-bls12-t3", ecc.BLS12_381, 3, 12, 5, 2, 1, 1),
	newInstance("griffin-bls12-t4", ecc.BLS12_381, 4, 20, 5, 3, 1, 1),
}

func newInstance(name string, curve ecc.ID, t, R, alpha, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, R, alpha, rate, capacity, digest)
	if err != nil {
		panic(fmt.Sprintf("griffin: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
