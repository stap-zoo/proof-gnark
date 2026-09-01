package main

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

// plan is what gets measured by default, for every workload: one instance per
// construction, all at width t=4 over the BLS12-381 scalar field, each carrying the
// round count the current analysis gives it (see the construction's instances.go).
// One width and one field is the smallest set in which every construction is
// comparable to every other, which is what makes it the default.
//
// EVERY LINE HERE IS t=4 OVER BLS12-381, the commented-out ones included: they are
// alternatives at that same width — a different exponent, a different arithmetization —
// so uncommenting one cannot quietly break the comparison. Other widths and the BN254
// twins stay in the registry and are reached with `-all` or `-instance` instead of from
// here. Comment a line out to drop an instance, a block to drop a construction.
//
// An entry is a typed reference into the construction's own registry: the registry
// fixes the construction *and* the arithmetization and is a symbol to jump to (its
// parameters are right there), while the name selects within it and panics at startup
// if that registry does not contain it. So a plan cannot be quietly wrong, and it can
// say which of two arithmetizations of one hash it means — Skyscraper's two Bar word
// sizes share their KAT labels and are told apart only by their registry.
//
// Skyscraper gets to t=4 by a different route from the rest, having no width parameter:
// its branch becomes an element of GF(p^2), so t = 2n. That makes its row a genuinely
// different permutation from the t=2 one rather than a wider version of it — the reason
// is in skyscraper.go's "Width" section and is not repeated here.
var plan = []harness.Instance{
	// ---- GMiMC ------------------------------------------------------- t=4
	harness.Pick(gmimc.Instances, "gmimc-bls12-t4"),

	// ---- GMiMC2 ------------------------------------------------------ t=4
	// alpha = 8 is the specification's recommended column at every width — the best
	// plain/proof-system tradeoff — so it is the row that stands for GMiMC2 here.
	// The other two exponents are one uncommented line away, and the alpha sweep is
	// the interesting one for this construction: fewer rounds buys more
	// multiplications per round, and the two cross over differently in R1CS (which
	// charges only the multiplications) than in PLONK (which at t=4 also charges the
	// t-1 branch additions a round makes; gmimc2.Permute's accumulator circuit, which
	// caps that at 3, only pays above this width and so is not exercised here).
	harness.Pick(gmimc2.Instances, "gmimc2-bls12-t4-a8"),
	// harness.Pick(gmimc2.Instances, "gmimc2-bls12-t4-a4"),
	// harness.Pick(gmimc2.Instances, "gmimc2-bls12-t4-a2"),

	// ---- Neptune ----------------------------------------------------- t=4
	harness.Pick(neptune.Instances, "neptune-bls12-t4"),

	// ---- Poseidon ---------------------------------------------------- t=4
	harness.Pick(poseidon.Instances, "poseidon-bls12-t4"),

	// ---- Poseidon2 --------------------------------------------------- t=4
	harness.Pick(poseidon2.Instances, "poseidon2-bls12-t4"),

	// ---- Anemoi ------------------------------------------------------ t=4
	harness.Pick(anemoi.Instances, "anemoi-bls12-381-scalar-t4"),

	// ---- Rescue-Prime ------------------------------------------------ t=4
	harness.Pick(rescueprime.Instances, "rescue-prime-bls12-t4"),

	// ---- Arion ------------------------------------------------------- t=4
	harness.Pick(arion.Instances, "arion-bls12-t4"),

	// ---- Griffin ----------------------------------------------------- t=4
	harness.Pick(griffin.Instances, "griffin-bls12-t4"),

	// ---- Skyscraper -------------------------------------------------- t=4
	// n=2: two Feistel branches over GF(p^2), which is how this construction gets a
	// four-element state (see the note above). Both Bar word sizes: the reference's
	// byte-word S-box (a 2^8 table) and the two-byte-word one (a 2^16 table, "-w16").
	// They satisfy the same KATs and differ only in cost, in opposite directions per
	// workload — the 2^16 table costs 65k R1CS more for one call and roughly half as
	// much per call once amortized — which is the whole reason both workloads exist.
	// It is also the expensive half of a -prove run (~40 s per instance, a 33-50 MB
	// proving key); drop the -w16 line for a quick one.
	harness.Pick(skyscraper.Instances, "skyscraper-bls12-381-n2"),
	harness.Pick(skyscraper.Instances16, "skyscraper-bls12-381-n2-w16"),

	// ---- Polocolo ---------------------------------------------------- t=4
	harness.Pick(polocolo.Instances, "polocolo-bls12-381-scalar-t4"),
	// harness.Pick(polocolo.Instances, "polocolo-bls12-381-scalar-t4-tight"),
}
