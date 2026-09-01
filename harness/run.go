package harness

import (
	"fmt"
	"testing"
)

// This file holds the *testing.T / *testing.B wrappers. Importing harness pulls
// in "testing"; in modern Go that does not register -test.* flags on a normal
// binary (that happens only when a test's M.Run executes), so cmd/bench can use
// Filter/Constraints without side effects.

// RunTests runs every reference vector against the target it belongs to — see
// [Target.MatchesVector] for the match. Targets with no matching vectors are
// skipped; if nothing matches at all the test is skipped.
func RunTests(t *testing.T, targets []Target, vectors []Vector) {
	t.Helper()

	ran := false
	for _, tg := range targets {
		i := 0
		for _, v := range vectors {
			if !tg.MatchesVector(v) {
				continue
			}
			ran = true
			v := v
			name := fmt.Sprintf("%s/%s/%d", tg.Instance, tg.Mode, i)
			i++
			t.Run(name, func(t *testing.T) {
				if err := Solve(tg.Case, v); err != nil {
					t.Fatalf("vector not satisfied: %v", err)
				}
			})
		}
	}
	if !ran {
		t.Skip("harness: no vectors matched any target")
	}
}

// RunBench reports each target's cost under both backends — R1CS constraints and
// PLONK gates (plus compile time). Bench subtests are named construction/field/mode
// for readability. For filtered comparisons prefer the CSV tool
// `go run ./cmd/bench -field ... -mode ... -backend ...`; `go test -bench` matches
// name segments positionally, which is fiddlier.
func RunBench(b *testing.B, targets []Target) {
	for _, tg := range targets {
		tg := tg
		name := fmt.Sprintf("%s/%s/%s", tg.Construction, tg.FieldName, tg.Mode)
		b.Run(name, func(b *testing.B) {
			var constraints, gates int
			for i := 0; i < b.N; i++ {
				var err error
				if constraints, err = Count(R1CS, tg.Case); err != nil {
					b.Fatal(err)
				}
				if gates, err = Count(PLONK, tg.Case); err != nil {
					b.Fatal(err)
				}
			}
			b.ReportMetric(float64(constraints), "constraints")
			b.ReportMetric(float64(gates), "gates")
		})
	}
}
