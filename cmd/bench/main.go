// Command bench measures what an arithmetization-oriented hash costs inside a
// circuit and prints one CSV row per measurement. It is the only benchmark command:
// three flags cover everything the SoK measures.
//
//	-workload   which circuit. "call" (the default) is one hash call — the unit the
//	            cost table compares. "message" is a sponge absorbing many elements
//	            and "merkle" a 2:1 tree: whole jobs, where a lookup table's setup
//	            amortizes over the calls and the mode's rate decides how many calls
//	            there are. Neither is derivable from the other.
//	-backend    which constraint system: r1cs (Groth16 proves it), plonk (its sparse
//	            system), or both. R1CS constraints and PLONK gates are different
//	            units — compare each against itself, never across.
//	-prove      also run a real setup, proof and verification and report their
//	            timings plus proof and key sizes. Off by default: compiling and
//	            counting is seconds where proving needs a proportionally large SRS
//	            and proving key, and the count is the primary metric.
//
// What gets measured is the plan in plan.go, or the whole registry with -all, or the
// instances named by -instance — narrowed by the label filters (-construction,
// -field, -kind, -mode).
//
//	go run ./cmd/bench                                   # the plan, one call each, R1CS + PLONK
//	go run ./cmd/bench -all                              # every target in the registry
//	go run ./cmd/bench -all -kind sponge -backend plonk  # PLONK gates, sponges only
//	go run ./cmd/bench -field bn254 -mode jive-2         # jive on bn254
//	go run ./cmd/bench -prove -runs 10                   # + setup/prove/verify timings
//	go run ./cmd/bench -prove -vectors                   # prove the KATs, not a generated input
//	go run ./cmd/bench -workload message                 # sponge over 32768 elements (~1 MiB)
//	go run ./cmd/bench -workload merkle                  # 2:1 tree over 16384 leaves
//	go run ./cmd/bench -workload merkle -scale -10       # 1/1024 the tree: a quick sanity run
//	go run ./cmd/bench -workload message -size 4,64,1024 # a sweep, cheapest first
//	go run ./cmd/bench -workload merkle -node sponge     # same mode for every construction
//	go run ./cmd/bench -instance griffin-bn254-t3 -out g.csv
//
// The workload sizes are large by default (16384 leaves is ~8.9M R1CS constraints
// for the heaviest construction, 13.5 GB peak): -scale shifts every size by 2^n, so
// -scale -10 is a sanity run that finishes in seconds and -scale 2 a bigger machine.
// Cost is affine in the call count, so a small sweep recovers the marginal cost
// exactly and the big sizes only confirm the extrapolation.
//
// gnark logs to stderr, as does this command's progress, so the CSV on stdout stays
// clean under `> results.csv`.
package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/registry"
)

func main() {
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, "bench:", err)
		os.Exit(1)
	}
}

func run() error {
	var (
		workloadFlag = flag.String("workload", "call", "circuit to measure: call, message (long sponge), or merkle (2:1 tree)")
		sizeFlag     = flag.String("size", "", "comma-separated workload sizes: elements (message) or leaves (merkle); default 32768 / 16384")
		scale        = flag.Int("scale", 0, "shift every workload size by 2^scale: -10 divides by 1024 (quick run), 2 quadruples it")
		nodeFlag     = flag.String("node", "auto", "merkle node: auto (arity-2 compression if declared, else sponge) or sponge")

		all          = flag.Bool("all", false, "measure every instance in the registry instead of the plan")
		instanceFlag = flag.String("instance", "", "comma-separated instances to measure instead of the plan")
		construction = flag.String("construction", "", "restrict to constructions (e.g. gmimc,anemoi)")
		field        = flag.String("field", "", "restrict to fields (e.g. bn254,bls12_381)")
		kindFlag     = flag.String("kind", "", "restrict to kinds (permutation,sponge,compression)")
		modeFlag     = flag.String("mode", "", "restrict to modes (permutation,sponge,jive-2,...)")

		backendFlag = flag.String("backend", "both", "constraint system(s): r1cs (Groth16), plonk, or both")
		prove       = flag.Bool("prove", false, "also run setup/prove/verify (slow and memory-hungry at workload sizes)")
		vectors     = flag.Bool("vectors", false, "with -prove: prove the reference vectors instead of a generated input")
		runs        = flag.Int("runs", 3, "prove/verify repetitions per measurement; reported as their mean and fastest")
		seed        = flag.String("seed", "bench", "seed for the generated witness input (reproducible)")

		failFast = flag.Bool("fail-fast", false, "stop at the first failure instead of reporting it and continuing")
		list     = flag.Bool("list", false, "print the resolved measurements and exit")
		out      = flag.String("out", "", "write the CSV here instead of stdout")
	)
	flag.Parse()

	wl, err := parseWorkload(*workloadFlag)
	if err != nil {
		return err
	}
	backends, err := parseBackends(*backendFlag)
	if err != nil {
		return err
	}
	kinds, err := parseKinds(*kindFlag)
	if err != nil {
		return err
	}
	sizes, err := sizesFor(wl, *sizeFlag, *scale)
	if err != nil {
		return err
	}
	spongeNodes, err := parseNode(*nodeFlag, wl)
	if err != nil {
		return err
	}
	if *vectors && wl != call {
		return fmt.Errorf("-vectors applies to -workload call only: a reference vector is one hash call")
	}

	// The work list: the plan, the whole registry, or what -instance names.
	insts := plan
	switch {
	case *instanceFlag != "":
		if insts, err = named(parseList(*instanceFlag)); err != nil {
			return err
		}
	case *all:
		insts = registry.Instances()
	}

	filter := harness.Filter{
		Constructions: parseList(*construction),
		Fields:        parseList(*field),
		Modes:         parseList(*modeFlag),
		Kinds:         kinds,
	}
	jobs, err := build(insts, wl, sizes, spongeNodes, filter)
	if err != nil {
		return err
	}
	if len(jobs) == 0 {
		return fmt.Errorf("nothing to measure (the filters exclude everything)")
	}

	if *list {
		for _, j := range jobs {
			fmt.Println(j.Case.Name)
		}
		fmt.Fprintf(os.Stderr, "bench: %d jobs x %d backend(s) = %d measurements\n",
			len(jobs), len(backends), len(jobs)*len(backends))
		return nil
	}

	return measure(jobs, backends, options{
		prove:    *prove,
		vectors:  *vectors,
		failFast: *failFast,
		out:      *out,
		proving:  harness.Options{Runs: *runs, Seed: *seed},
	})
}

// named resolves instance names typed on the command line. A plan does not come
// through here — its entries are typed references (harness.Pick) — so this is only
// for user input, where an unknown name is an error rather than a skip: a typo must
// not quietly shorten the table.
func named(names []string) ([]harness.Instance, error) {
	all := registry.Instances()
	out := make([]harness.Instance, 0, len(names))
	for _, name := range names {
		found := false
		for _, inst := range all {
			if inst.InstanceName() == name {
				out = append(out, inst)
				found = true
				break
			}
		}
		if !found {
			return nil, fmt.Errorf("-instance %q: unknown instance (see `go run ./cmd/bench -all -list`)", name)
		}
	}
	return out, nil
}

func parseBackends(s string) ([]harness.Backend, error) {
	if strings.TrimSpace(s) == "both" {
		return []harness.Backend{harness.R1CS, harness.PLONK}, nil
	}
	var out []harness.Backend
	for _, name := range parseList(s) {
		switch name {
		case "r1cs", "groth16":
			out = append(out, harness.R1CS)
		case "plonk":
			out = append(out, harness.PLONK)
		default:
			return nil, fmt.Errorf("unknown backend %q (want r1cs, plonk, or both)", name)
		}
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("no backend selected")
	}
	return out, nil
}

func parseKinds(s string) ([]harness.Kind, error) {
	var kinds []harness.Kind
	for _, name := range parseList(s) {
		k, err := harness.ParseKind(name)
		if err != nil {
			return nil, err
		}
		kinds = append(kinds, k)
	}
	return kinds, nil
}

func parseList(s string) []string {
	if strings.TrimSpace(s) == "" {
		return nil
	}
	return strings.Split(s, ",")
}
