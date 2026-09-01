// Package registry is the list of every construction: the one place that knows the
// full set. The benchmark command and the cross-construction tests both build on it.
//
// It is a separate package from harness for a structural reason rather than a
// stylistic one: the constructions import harness (to describe their instances and
// embed their vectors), so harness cannot import the constructions. This is where
// that cycle is broken, and it is why this package holds nothing but the list.
//
// Add a construction by appending its registry to [Instances] and its vectors to
// [Vectors]. Everything else — targets, filtering, coverage — follows from those.
package registry

import (
	"github.com/zkhash-sok/gnark-hashes/anemoi"
	"github.com/zkhash-sok/gnark-hashes/arion"
	"github.com/zkhash-sok/gnark-hashes/gmimc"
	"github.com/zkhash-sok/gnark-hashes/gmimc2"
	"github.com/zkhash-sok/gnark-hashes/griffin"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/neptune"
	"github.com/zkhash-sok/gnark-hashes/polocolo"
	"github.com/zkhash-sok/gnark-hashes/poseidon"
	"github.com/zkhash-sok/gnark-hashes/poseidon2"
	"github.com/zkhash-sok/gnark-hashes/rescueprime"
	"github.com/zkhash-sok/gnark-hashes/skyscraper"
)

// Instances is every instance of every construction, in a fixed order.
//
// Skyscraper appears twice: its byte-word Bar and the two-byte-word arithmetization
// of the same instances (named "…-n1-w16"). The two compute the identical hash and
// share its reference vectors — the label indirection is harness.VectorSource — but
// they are separate rows wherever cost is reported, because what a 2^16 lookup table
// costs and when it pays is one of the SoK's own questions.
func Instances() []harness.Instance {
	var out []harness.Instance
	out = append(out, harness.Lift(gmimc.Instances)...)
	out = append(out, harness.Lift(gmimc2.Instances)...)
	out = append(out, harness.Lift(neptune.Instances)...)
	out = append(out, harness.Lift(poseidon.Instances)...)
	out = append(out, harness.Lift(poseidon2.Instances)...)
	out = append(out, harness.Lift(anemoi.Instances)...)
	out = append(out, harness.Lift(rescueprime.Instances)...)
	out = append(out, harness.Lift(arion.Instances)...)
	out = append(out, harness.Lift(griffin.Instances)...)
	out = append(out, harness.Lift(skyscraper.Instances)...)
	out = append(out, harness.Lift(skyscraper.Instances16)...)
	out = append(out, harness.Lift(polocolo.Instances)...)
	return out
}

// Targets is every runnable (instance × mode) of every construction — each
// instance's permutation plus whichever sponge and compression it declares.
func Targets() []harness.Target { return harness.Targets(Instances()...) }

// Vectors is every construction's reference vectors.
func Vectors() []harness.Vector {
	var v []harness.Vector
	v = append(v, gmimc.Vectors()...)
	v = append(v, gmimc2.Vectors()...)
	v = append(v, neptune.Vectors()...)
	v = append(v, poseidon.Vectors()...)
	v = append(v, poseidon2.Vectors()...)
	v = append(v, anemoi.Vectors()...)
	v = append(v, rescueprime.Vectors()...)
	v = append(v, arion.Vectors()...)
	v = append(v, griffin.Vectors()...)
	v = append(v, skyscraper.Vectors()...)
	v = append(v, polocolo.Vectors()...)
	return v
}
