package main

// The mechanics: turning instances into jobs for the chosen workload, and measuring
// them into CSV rows.

import (
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"math/bits"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/consensys/gnark/logger"
	"github.com/zkhash-sok/gnark-hashes/harness"
	"github.com/zkhash-sok/gnark-hashes/registry"
)

// workloadKind is which circuit to measure.
type workloadKind int

const (
	call    workloadKind = iota // one hash call: the target as the harness builds it
	message                     // a sponge absorbing n elements
	merkle                      // a 2:1 Merkle tree over n leaves
)

func (w workloadKind) String() string {
	switch w {
	case message:
		return "message"
	case merkle:
		return "merkle"
	default:
		return "call"
	}
}

func parseWorkload(s string) (workloadKind, error) {
	switch s {
	case "call":
		return call, nil
	case "message":
		return message, nil
	case "merkle":
		return merkle, nil
	default:
		return 0, fmt.Errorf("unknown -workload %q (want call, message, or merkle)", s)
	}
}

// job is one row of work: a target whose Case has been replaced by the
// workload's circuit, plus the size that workload was built at.
//
// Reusing the target is what keeps this command to one code path: the labels
// (construction, instance, field, kind, mode) and the label filter are the harness's
// for every workload, and only the circuit differs.
type job struct {
	harness.Target
	wl   workloadKind
	size int
}

// build resolves the instances into jobs, smallest size first, and applies the label
// filter. Sizes are the outer loop because rows are flushed as they are produced: a
// sweep that dies on the largest size still leaves every smaller one on disk.
func build(insts []harness.Instance, wl workloadKind, sizes []int, spongeNodes bool, filter harness.Filter) ([]job, error) {
	var out []job
	for _, n := range sizes {
		for _, inst := range insts {
			js, err := jobsFor(inst, wl, n, spongeNodes)
			if err != nil {
				return nil, err
			}
			for _, j := range js {
				if len(filter.Select([]harness.Target{j.Target})) == 1 {
					out = append(out, j)
				}
			}
		}
	}
	return out, nil
}

// jobsFor builds one instance's jobs: every mode it declares for a single call, and
// the one mode the workload uses otherwise.
func jobsFor(inst harness.Instance, wl workloadKind, n int, spongeNodes bool) ([]job, error) {
	targets := harness.Targets(inst)
	if wl == call {
		out := make([]job, len(targets))
		for i, tg := range targets {
			out[i] = job{Target: tg, wl: call, size: 1}
		}
		return out, nil
	}

	// A workload runs in exactly one mode, so find the target for that mode and swap
	// in the workload's circuit: same labels, different case.
	var mode string
	var c harness.Case
	var err error
	switch wl {
	case message:
		// The instance's own sponge (Anemoi's Hirose, else plain).
		if c, err = harness.Message(inst, n); err != nil {
			return nil, err
		}
		for _, tg := range targets {
			if tg.Kind == harness.Sponge {
				mode = tg.Mode
				break
			}
		}
	case merkle:
		newNode := harness.NewNode
		if spongeNodes {
			newNode = harness.NewSpongeNode
		}
		node, nerr := newNode(inst)
		if nerr != nil {
			return nil, nerr
		}
		if c, err = harness.Merkle(inst, node, n); err != nil {
			return nil, err
		}
		mode = node.Name()
	}
	for _, tg := range targets {
		if tg.Mode == mode {
			tg.Case = c
			return []job{{Target: tg, wl: wl, size: n}}, nil
		}
	}
	return nil, fmt.Errorf("%s: no %s target for the %s workload", inst.InstanceName(), mode, wl)
}

// sizesFor resolves -size and -scale into the sizes to measure, ascending so a sweep
// runs cheapest first. The defaults are the headline sizes: 32768 elements is 1 MiB
// of field elements, 16384 leaves is 16383 hash calls.
func sizesFor(wl workloadKind, sizeFlag string, scale int) ([]int, error) {
	if wl == call {
		if sizeFlag != "" {
			return nil, fmt.Errorf("-size applies to -workload message/merkle; a single call has no size")
		}
		if scale != 0 {
			return nil, fmt.Errorf("-scale applies to -workload message/merkle; a single call has no size")
		}
		return []int{1}, nil
	}

	sizes := []int{32768}
	if wl == merkle {
		sizes = []int{16384}
	}
	if sizeFlag != "" {
		sizes = nil
		for _, tok := range parseList(sizeFlag) {
			n, err := strconv.Atoi(strings.TrimSpace(tok))
			if err != nil || n < 1 {
				return nil, fmt.Errorf("invalid -size %q (want positive integers)", tok)
			}
			sizes = append(sizes, n)
		}
	}
	for i := 1; i < len(sizes); i++ { // insertion sort: the lists are tiny
		for j := i; j > 0 && sizes[j] < sizes[j-1]; j-- {
			sizes[j], sizes[j-1] = sizes[j-1], sizes[j]
		}
	}
	if scale == 0 {
		return sizes, nil
	}
	// Every size here is a power of two, so scaling keeps it one — which the Merkle
	// workload requires of its leaf count. A division that would truncate means the
	// flags disagree, and silently benchmarking a different size is worse than an
	// error.
	for i, n := range sizes {
		switch {
		case scale > 0:
			if scale >= 63 || n > math.MaxInt>>scale {
				return nil, fmt.Errorf("-scale %d overflows size %d", scale, n)
			}
			sizes[i] = n << scale
		default:
			d := -scale
			if d >= 63 || n>>d == 0 {
				return nil, fmt.Errorf("-scale %d shrinks size %d to nothing", scale, n)
			}
			if (n>>d)<<d != n {
				return nil, fmt.Errorf("-scale %d does not divide size %d evenly", scale, n)
			}
			sizes[i] = n >> d
		}
	}
	return sizes, nil
}

// parseNode reads -node. Forcing the sponge everywhere gives a same-mode column
// across all constructions, at the cost of not using the better primitive where one
// exists; it is meaningless for the other workloads.
func parseNode(s string, wl workloadKind) (spongeNodes bool, err error) {
	switch s {
	case "auto":
		return false, nil
	case "sponge":
		if wl != merkle {
			return false, fmt.Errorf("-node applies to -workload merkle only")
		}
		return true, nil
	default:
		return false, fmt.Errorf("unknown -node %q (want auto or sponge)", s)
	}
}

// options is what measure needs beyond the jobs themselves.
type options struct {
	prove    bool
	vectors  bool
	failFast bool
	out      string
	proving  harness.Options
}

// measure runs every job under every backend and writes one CSV row per pair.
func measure(jobs []job, backends []harness.Backend, o options) error {
	// gnark's compile/setup progress logs default to stdout; move them aside so the
	// CSV is the only thing there.
	logger.SetOutput(os.Stderr)

	sink := io.Writer(os.Stdout)
	if o.out != "" {
		f, err := os.Create(o.out)
		if err != nil {
			return err
		}
		defer f.Close()
		sink = f
	}
	w := csv.NewWriter(sink)
	defer w.Flush()
	if err := w.Write(header(o.prove)); err != nil {
		return err
	}

	total := len(jobs) * len(backends)
	done, failed := 0, 0
	started := time.Now()

	for _, j := range jobs {
		opts := o.proving
		if o.vectors {
			// One lookup per job, not per backend: the witness source is the same.
			v, ok := vectorFor(j.Target)
			if !ok {
				return fmt.Errorf("-vectors: no reference vector for %s", j.Case.Name)
			}
			opts.Vector = &v
		}
		for _, b := range backends {
			done++
			fmt.Fprintf(os.Stderr, "bench: [%d/%d] %s %s\n", done, total, b.Name, j.Case.Name)

			row, err := measureOne(j, b, o.prove, opts)
			if err != nil {
				if o.failFast {
					return err
				}
				// Keep going: one construction failing should not cost the whole
				// table. The row is skipped and reported at the end.
				failed++
				fmt.Fprintf(os.Stderr, "bench: FAILED %s %s: %v\n", b.Name, j.Case.Name, err)
				continue
			}
			if err := w.Write(row); err != nil {
				return err
			}
			// Flush per row so a long run's partial results survive an interrupt.
			w.Flush()
			if err := w.Error(); err != nil {
				return err
			}
		}
	}

	fmt.Fprintf(os.Stderr, "bench: %d/%d measured in %s", total-failed, total, time.Since(started).Round(time.Second))
	if failed > 0 {
		fmt.Fprintf(os.Stderr, " (%d FAILED)", failed)
	}
	fmt.Fprintln(os.Stderr)

	// PLONK's universal setup, reported once: it is one SRS per curve, shared by
	// every row over that curve, so the srs_ms column can only say how far the run
	// had grown it by that row. This is the whole of it. Empty without -prove, and
	// empty for a curve nothing measured.
	for _, u := range harness.DefaultSRS().Universal() {
		fmt.Fprintf(os.Stderr, "bench: universal KZG SRS, %s: %d powers of tau "+
			"(domain 2^%d) in %s — one per curve, not per row\n",
			u.Curve, u.Points, bits.TrailingZeros64(u.Domain),
			u.Generate.Round(time.Millisecond))
	}
	if failed > 0 {
		return fmt.Errorf("%d of %d measurements failed", failed, total)
	}
	return nil
}

// measureOne produces one row: the cost alone, or the whole pipeline with -prove.
func measureOne(j job, b harness.Backend, prove bool, opts harness.Options) ([]string, error) {
	labels := []string{
		j.Construction, j.Instance, j.FieldName, j.Kind.String(), j.Mode,
		strconv.Itoa(j.Width),
		j.wl.String(), strconv.Itoa(j.size), b.Name, b.System,
	}
	if !prove {
		start := time.Now()
		n, err := harness.Count(b, j.Case)
		if err != nil {
			return nil, err
		}
		return append(labels, strconv.Itoa(n), ms(time.Since(start))), nil
	}
	r, err := harness.Measure(b, j.Case, opts)
	if err != nil {
		return nil, err
	}
	return append(labels,
		strconv.Itoa(r.Constraints), ms(r.Compile),
		strconv.Itoa(r.PublicWires), strconv.Itoa(r.SecretWires),
		ms(r.SRS), ms(r.Setup), ms(r.Witness),
		ms(r.Prove), ms(r.ProveMin), ms(r.Verify), ms(r.VerifyMin), strconv.Itoa(r.Runs),
		strconv.FormatInt(r.ProofBytes, 10), strconv.FormatInt(r.PKBytes, 10),
		strconv.FormatInt(r.VKBytes, 10)), nil
}

// header names the CSV columns. width is the permutation's state size t (see
// harness.Target.Width — the number, never the reference's spelling of it);
// constraints is the backend's own unit (R1CS constraints or PLONK gates);
// durations are milliseconds; sizes are bytes of gnark's point-compressed encoding.
func header(prove bool) []string {
	h := []string{
		"construction", "instance", "field", "kind", "mode", "width",
		"workload", "size", "backend", "system",
		"constraints", "compile_ms",
	}
	if prove {
		h = append(h, "public_wires", "secret_wires",
			"srs_ms", "setup_ms", "witness_ms",
			"prove_ms", "prove_min_ms", "verify_ms", "verify_min_ms", "runs",
			"proof_bytes", "pk_bytes", "vk_bytes")
	}
	return h
}

// vectorFor finds a reference vector for a target — the first one it matches, by the
// same rule the vector tests use.
func vectorFor(tg harness.Target) (harness.Vector, bool) {
	for _, v := range registry.Vectors() {
		if tg.MatchesVector(v) {
			return v, true
		}
	}
	return harness.Vector{}, false
}

func ms(d time.Duration) string {
	return strconv.FormatFloat(float64(d)/float64(time.Millisecond), 'f', 3, 64)
}
