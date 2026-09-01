package arion

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Arion parameter set pinned to a curve's scalar field.
// Mirrors arion/instances.py, whose scalar-field instance is ARION_BLS12_T3 —
// every parameter here is the reference's, none invented.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "arion" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Arion's sponge: zero-padding, add-to-start absorption, and a
// conditional IV — an all-zero IV when the input is rate-aligned, otherwise the
// (unpadded) input length encoded in the first capacity slot. Matches the
// reference Arion.hash_sponge. Arion has no compression mode.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, arionIV(p.Rate), mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("arion: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// arionIV builds Arion's rate-dependent IV. mode.IV receives only
// (inputLen, capacity, digest) — not the rate — so the "was the input
// rate-aligned?" test that the reference makes cannot be expressed with the stock
// ZeroIV/LengthIV builders alone. This closure captures the rate and dispatches
// to them: ZeroIV on an aligned input, LengthIV otherwise.
func arionIV(rate int) mode.IV {
	return func(inputLen, capacity, digest int) []frontend.Variable {
		if inputLen%rate == 0 {
			return mode.ZeroIV(inputLen, capacity, digest)
		}
		return mode.LengthIV(inputLen, capacity, digest)
	}
}

// Instances is the registry of concrete instances over the BLS12-381 scalar field:
// the t=3 instance at the Arion.sage defaults — t=3, R=6, alpha1=5, alpha2=257,
// r/c/d = 2/1/2 (kappa=128) — and a t=4 one at R=10 with the same exponents. The
// t=3 instance-name string matches the reference KAT label arion-bls12-t3; the
// reference has no t=4 instance (and no round derivation of its own), so t=4's
// vectors come from tools/export_kat_t4.py.
//
// t=4 needs no new derivation: circ(1,2,3,4) is still MDS (Arion paper Remark 4
// covers t in {2,3,4}) and the GTDS coefficients are drawn per branch. What does
// change is cost — the shared-form addition program exists only at t=3, so the
// linear layer runs dense here; see slp_test.go.
var Instances = []Instance{
	newInstance("arion-bls12-t3", ecc.BLS12_381, 3 /*t*/, 6 /*R*/, 5 /*alpha1*/, 257 /*alpha2*/, 2 /*r*/, 1 /*c*/, 2 /*d*/),
	newInstance("arion-bls12-t4", ecc.BLS12_381, 4, 10, 5, 257, 3, 1, 2),
}

func newInstance(name string, curve ecc.ID, t, R, alpha1, alpha2, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, R, alpha1, alpha2, rate, capacity, digest)
	if err != nil {
		panic(fmt.Sprintf("arion: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
