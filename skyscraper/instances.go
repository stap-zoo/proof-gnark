package skyscraper

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Skyscraper parameter set pinned to a curve's scalar
// field. Mirrors the SKYSCRAPER_<FIELD>_N<n> instances in skyscraper/instances.py
// over the BN254 and BLS12-381 scalar fields, where n is the extension degree of the
// branch field GF(p^n) and the state is t = 2n base-field elements.
type Instance struct {
	Name   string
	Curve  ecc.ID
	Params *Parameters

	// Vectors is the instance name whose reference vectors this instance must
	// satisfy. It is Name for the byte-word instances, which is what the reference
	// generated its KATs under; the 16-bit-word variants set it to their byte-word
	// twin's name, because they compute the identical function under a different
	// arithmetization (see Instances16 and harness.VectorSource).
	Vectors string
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "skyscraper" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// VectorInstance implements harness.VectorSource: it is the name this instance's
// KATs are labeled with, which is its own except for the 16-bit-word variants.
func (i Instance) VectorInstance() string { return i.Vectors }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.Compressor ---

// Compressions returns Skyscraper's 2-to-1 Davies-Meyer compression over the two
// Feistel branches (perm(l||r)_L + l), matching the reference Skyscraper.compress.
// The block is one branch, i.e. half the state.
func (i Instance) Compressions() []mode.Compression {
	dm, err := mode.NewDaviesMeyer(i.Permutation(), i.Params.T/2)
	if err != nil {
		panic(fmt.Sprintf("skyscraper: instance %q davies-meyer: %v", i.Name, err))
	}
	return []mode.Compression{dm}
}

// --- mode.SpongeHash ---

// Sponges returns Skyscraper's sponge: zero padding, add-to-start absorption, and
// a zero IV over the capacity, matching the reference Skyscraper.hash_sponge
// (IV = [0]*c), which is what its vectors pin.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.ZeroIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("skyscraper: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// Instances is the registry of concrete instances with the reference's byte-word
// Bar. Every one has the reference's R=18 Feistel rounds and Bar at rounds
// {6,7,10,11}; what varies is the extension degree n of the branch field, which is
// Skyscraper's only width parameter (state t = 2n) and the sponge shape that follows
// from it (r=c=d=n).
//
//	skyscraper-bn254-n1        n=1  t=2  the prime field
//	skyscraper-bls12-381-n1    n=1  t=2
//	skyscraper-bls12-381-n2    n=2  t=4  GF(p^2) = F_p[X]/(X^2 + 5)
//
// The instance-name strings match the reference KAT labels (variable name lowercased,
// "_" -> "-"), and the "n<k>" is the reference's own label rather than a parameter of
// this implementation — it is what ties an instance to its vectors. beta = 5 is the
// reference's own choice for both curve scalar fields (its fmod = [5, 0, 1]);
// NewParameters re-derives nothing about it but does check that X^2 + 5 is irreducible
// over the field it is handed.
//
// n=2 is registered over BLS12-381 only, which is where the cross-construction plan
// compares (cmd/bench/plan.go: t=4 over the BLS12-381 scalar field). The reference
// also defines SKYSCRAPER_BN254_N2 and both N3 instances; adding one is a line here
// plus its four reference vectors, except N3, which params.go rejects — see
// checkExtension.
//
// The two-byte-word arithmetization of the same instances is [Instances16].
var Instances = []Instance{
	newInstance("skyscraper-bn254-n1", ecc.BN254, 1 /*n*/, nil /*beta*/, 18 /*R*/, 1 /*r*/, 1 /*c*/, 1 /*d*/, 8 /*wordBits*/),
	newInstance("skyscraper-bls12-381-n1", ecc.BLS12_381, 1, nil, 18, 1, 1, 1, 8),
	newInstance("skyscraper-bls12-381-n2", ecc.BLS12_381, 2, big.NewInt(5), 18, 2, 2, 2, 8),
}

// Instances16 mirrors Instances but with the two-byte (16-bit) word Bar S-box: the
// decomposition looks its words up in a 2^16-entry table (entry i = the byte S-box
// applied to each byte of i), which halves the word count. The output is
// byte-identical to the byte-word Bar, so these satisfy the SAME reference vectors
// (see TestVectors16) — which is why each keeps its twin's KAT label in Vectors
// while taking a name of its own, the "-w16" suffix.
//
// It is an alternative *arithmetization* of the same hash, not a distinct
// construction, and it is a SoK datum in its own right because the two rank in
// opposite orders depending on the unit: at n=1 the 2^16 table costs 65,280 R1CS
// constraints more than the 2^8 one for a single call (66,305 against 1,025), and then
// cuts the marginal call by nearly two fifths (212 against 340 per Merkle node), so it
// overtakes after ~510 calls. Measure it with `go run ./cmd/bench -workload merkle`;
// the amortization itself is guarded by harness.TestLookupTablesStayShared and the two
// registries are compared by TestBarWordSizeConstraintCost.
//
// The extension degree is a second axis on the same fact and pushes the same way. The
// table is one per circuit and the Bars are 2n per Bar round, so going to n=2 adds
// 1,724 PLONK gates to the byte-word permutation (3,566 -> 5,290) and 1,148 to the
// two-byte-word one (329,390 -> 330,538): the wider word is the wrong choice for one
// call by a smaller factor at every step up in width, for exactly the reason it is the
// right choice for a tree.
//
// The separate slice is what lets a plan name one arithmetization or the other
// (harness.Pick over this registry or over Instances).
//
// Derived from Instances rather than written out, so the two sets cannot drift in
// anything but word size.
var Instances16 = wordBits16(Instances)

// wordBits16 rebuilds a registry with the two-byte-word Bar: same field, rounds and
// sponge shape, wordBits 16, a "-w16" name, and the byte-word instance's KAT label.
func wordBits16(base []Instance) []Instance {
	out := make([]Instance, len(base))
	for i, b := range base {
		p := b.Params
		params, err := NewParameters(b.Curve.ScalarField(), p.N, p.Beta, p.R, p.Rate, p.Capacity, p.Digest, 16)
		if err != nil {
			panic(fmt.Sprintf("skyscraper: invalid 16-bit-word instance for %q: %v", b.Name, err))
		}
		out[i] = Instance{Name: b.Name + "-w16", Curve: b.Curve, Params: params, Vectors: b.Name}
	}
	return out
}

func newInstance(name string, curve ecc.ID, n int, beta *big.Int, R, rate, capacity, digest, wordBits int) Instance {
	params, err := NewParameters(curve.ScalarField(), n, beta, R, rate, capacity, digest, wordBits)
	if err != nil {
		panic(fmt.Sprintf("skyscraper: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Params: params, Vectors: name}
}

// Ensure Instance implements the harness metadata, the vector-label hook, and both
// mode-eligibility interfaces.
var (
	_ harness.Instance     = Instance{}
	_ harness.VectorSource = Instance{}
	_ mode.Compressor      = Instance{}
	_ mode.SpongeHash      = Instance{}
)
