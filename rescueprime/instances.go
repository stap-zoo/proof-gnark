package rescueprime

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Rescue-Prime parameter set pinned to a curve's scalar
// field. Mirrors the RESCUE_PRIME_* instances in marvellous/instances.py over the
// BN254 and BLS12-381 scalar fields.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "rescueprime" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Rescue-Prime's sponge: variable-length padding (pad-one), a
// zero IV of capacity length, and add-to-start absorption — matching the
// reference RescuePrime.hash_sponge.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.ZeroIV, mode.AddToStart, mode.PadOne)
	if err != nil {
		panic(fmt.Sprintf("rescueprime: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// Instances is the registry of concrete instances (BN254 + BLS12-381 scalar).
// Naming matches the reference KAT labels: rescue-prime-<field>-t<width>. The t=3
// pair uses t=3, alpha=5, R=14, r/c/d = 2/1/2 (kappa=128); g is the field generator
// (utils/field.py): 5 for BN254, 7 for BLS12-381.
//
// t=4 needs no new derivation — the Vandermonde MDS and the SHAKE round constants
// are built for any width — only a round count. R=13 is the updated t=4 number at
// alpha=5 (the reference's own criterion, marvellous/params.py _init_rounds, still
// returns 11 there; the newer analysis supersedes it, which is why this instance
// does not exist upstream and its vectors come from tools/export_kat_t4.py). The
// companion numbers, unused here because alpha=5 is the exponent every instance in
// this repo uses: alpha=3 -> R=18, alpha=7 -> R=11.
var Instances = []Instance{
	newInstance("rescue-prime-bn254-t3", ecc.BN254, 3 /*t*/, 5 /*alpha*/, 5 /*g*/, 14 /*R*/, 2 /*r*/, 1 /*c*/, 2 /*d*/),
	newInstance("rescue-prime-bls12-t3", ecc.BLS12_381, 3, 5, 7, 14, 2, 1, 2),
	newInstance("rescue-prime-bls12-t4", ecc.BLS12_381, 4, 5, 7, 13, 3, 1, 2),
}

func newInstance(name string, curve ecc.ID, t, alpha int, g int64, R, rate, capacity, digest int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, alpha, big.NewInt(g), R, rate, capacity, digest)
	if err != nil {
		panic(fmt.Sprintf("rescueprime: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
