package poseidon2

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is one concrete Poseidon2 parameter set pinned to a curve's scalar
// field. Mirrors hades/instances.py over the BN254 and BLS12-381 scalar fields.
type Instance struct {
	Name    string
	Curve   ecc.ID
	Version string
	Params  *Parameters
}

// --- harness.Instance ---

func (i Instance) Construction() string { return "poseidon2" }
func (i Instance) InstanceName() string { return i.Name }
func (i Instance) FieldName() string    { return i.Curve.String() }
func (i Instance) Field() *big.Int      { return i.Curve.ScalarField() }

// Permutation builds the in-circuit permutation for this instance.
func (i Instance) Permutation() permutation.Permutation { return NewPermutation(i.Params) }

// --- mode.SpongeHash ---

// Sponges returns Poseidon2's sponge mode: a plain sponge with zero padding,
// add-to-start absorption and a length-encoding IV, matching the reference hades
// hash_sponge (IV = [len(data)] + [0]*(c-1)), inherited unchanged from the Hades
// base. Poseidon2 defines no compression mode of its own, so none is declared.
func (i Instance) Sponges() []mode.Sponge {
	p := i.Params
	s, err := mode.NewPlainSponge(i.Permutation(), p.Rate, p.Capacity, p.Digest, mode.LengthIV, mode.AddToStart, mode.PadZero)
	if err != nil {
		panic(fmt.Sprintf("poseidon2: instance %q sponge: %v", i.Name, err))
	}
	return []mode.Sponge{s}
}

// matDiag parses the MAT_DIAG_M_1 vector (decimal strings) into big.Ints.
func matDiag(decs ...string) []*big.Int {
	out := make([]*big.Int, len(decs))
	for i, s := range decs {
		v, ok := new(big.Int).SetString(s, 10)
		if !ok {
			panic(fmt.Sprintf("poseidon2: invalid mat_diag entry %q", s))
		}
		out[i] = v
	}
	return out
}

// Instances is the registry of concrete instances (bn254 + bls12-381 scalar).
// Naming convention: poseidon2-<field>-t<width>, matching the reference KAT
// labels. All use version "isec" (the reference default for Poseidon2). The
// mat_diag (MAT_DIAG_M_1) vectors are field-specific spec data, taken verbatim
// from hades/instances.py; t in {2,3} use the small fixed diagonals.
var Instances = []Instance{
	newInstance("poseidon2-bls12-t2", ecc.BLS12_381, "isec", 2 /*t*/, 5 /*alpha*/, 8 /*R_ext*/, 56 /*R_int*/, 1 /*r*/, 1 /*c*/, 1, /*d*/
		matDiag("1", "2")),
	newInstance("poseidon2-bls12-t3", ecc.BLS12_381, "isec", 3, 5, 8, 56, 2, 1, 1,
		matDiag("1", "1", "2")),
	newInstance("poseidon2-bls12-t4", ecc.BLS12_381, "isec", 4, 5, 8, 56, 3, 1, 1,
		matDiag(
			"1655454839116271620575749293234574274222993456067245397034476501783047702986",
			"50438009192746386491573599403504532375212735259047687062053082186860235946862",
			"2333904567228711547445735898465837165392237851337731850779276628672040817654",
			"50596178259748710644710488613201561091625758003284944221322227344704022331809")),
	newInstance("poseidon2-bls12-t8", ecc.BLS12_381, "isec", 8, 5, 8, 57, 7, 1, 1,
		matDiag(
			"30096150855626815013648098211149190726024462037711422151903740709493260300796",
			"1670165440712639538351788804048794320334719999221164772595439843145725878342",
			"46650283196255499715295956296451823074186598793282238922518499397091394504987",
			"35540998704038434561306263855447791582438082003486383019176839889144643735624",
			"5729993341210689688259865343737097561899201100829809250544373980467298239994",
			"50448381164981534770780393450590297246005426700091703154028953045831153438287",
			"17212607781899146579730953203574521410127157101393873991200310682506890698706",
			"34706517679585386573351886038051072526558536116546116247060202675936040458533")),
	newInstance("poseidon2-bn254-t3", ecc.BN254, "isec", 3, 5, 8, 56, 2, 1, 1,
		matDiag("1", "1", "2")),
}

func newInstance(name string, curve ecc.ID, version string, t, alpha, rExt, rInt, rate, capacity, digest int, mDiag []*big.Int) Instance {
	params, err := NewParameters(curve.ScalarField(), t, alpha, rExt, rInt, rate, capacity, digest, mDiag, version)
	if err != nil {
		panic(fmt.Sprintf("poseidon2: invalid instance %q: %v", name, err))
	}
	return Instance{Name: name, Curve: curve, Version: version, Params: params}
}

// Ensure Instance implements the harness metadata and the mode-eligibility
// interfaces it declares.
var (
	_ harness.Instance = Instance{}
	_ mode.SpongeHash  = Instance{}
)
