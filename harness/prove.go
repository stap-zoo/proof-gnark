// This file adds the *proving* pipeline to the harness. Everything else here
// stops at constraint satisfaction (see Solve, which runs no prover) or at a
// static count (Count); this is the part that runs a real setup, produces a real
// proof, and verifies it — the source of the SoK's timing and proof-size tables.
//
// It is written once, generically, for the same reason Count is: a construction
// contributes a [Case] and nothing else, so both proof systems apply to all eleven
// constructions and every mode without a line of per-construction code.
//
// Two pieces of machinery make that possible:
//
//   - The proof system rides on [Backend]. R1CS already meant "the Groth16
//     arithmetization" and PLONK "the SCS one", so each Backend simply carries
//     its prover, and the concrete gnark types (groth16.ProvingKey vs
//     plonk.ProvingKey, which share no interface) stay behind the closures in
//     setup. Callers see only [Prover], [Proof], and durations.
//   - A witness has to *satisfy* the circuit, and every circuit here asserts a
//     computed result against a witnessed one, so a generated input is useless
//     until something supplies the matching output. [Evaluate] does that by
//     running Case.Eval under gnark's native test engine — the same evaluator
//     the vector tests use — which means proving needs no out-of-circuit
//     reimplementation of any hash.
package harness

import (
	"fmt"
	"io"
	"math/big"
	"time"

	"github.com/consensys/gnark/backend/groth16"
	"github.com/consensys/gnark/backend/plonk"
	"github.com/consensys/gnark/backend/witness"
	"github.com/consensys/gnark/constraint"
	"github.com/consensys/gnark/frontend"
	"github.com/consensys/gnark/test"
	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Proof is an opaque handle on a produced proof. The harness never inspects a
// proof — it times it and measures its serialized size — so the only behaviour
// it needs is io.WriterTo, which every gnark proof provides. The concrete type
// is groth16.Proof or plonk.Proof, recovered inside the verifier closure that
// was created alongside it.
type Proof interface{ io.WriterTo }

// Options tunes the proving pipeline. The zero value is valid: one run, a fixed
// default seed, a generated input, and no on-disk SRS cache.
type Options struct {
	// Runs is how many times Prove and Verify are repeated per measurement;
	// the reported duration is the mean, alongside the fastest run. 0 means 1.
	Runs int

	// Seed makes the generated input reproducible: the same seed and Case
	// always give the same witness. Ignored when Vector is set.
	Seed string

	// Vector, if non-nil, is used as the witness instead of a generated input,
	// which makes the proof a proof of a *reference* vector — a stronger check
	// (it re-validates the KAT through the real prover) at the cost of being
	// limited to the inputs the reference publishes.
	Vector *Vector

	// SRS is the universal KZG setup PLONK draws from. nil means the
	// process-wide store ([DefaultSRS]), which is what keeps the zero value
	// valid and what lets every measurement in a run share one SRS per curve;
	// pass one to isolate a measurement from the rest of the process.
	SRS *SRS
}

// srs is the store to use: the caller's, or the process-wide one.
func (o Options) srs() *SRS {
	if o.SRS != nil {
		return o.SRS
	}
	return defaultSRS
}

func (o Options) runs() int {
	if o.Runs < 1 {
		return 1
	}
	return o.Runs
}

// keys is a completed setup, erased of its proof system: the two serializable
// keys plus closures that already hold the typed ones. This is what lets one
// Prover serve both Groth16 and PLONK.
type keys struct {
	pk, vk io.WriterTo

	// setup is the circuit-specific preprocessing: the proof system's own setup,
	// plus — for PLONK — the Lagrange-basis SRS derived for this circuit's
	// domain. Timed on this side of the interface rather than by the caller
	// because only this side knows which of the work belongs to the circuit and
	// which to the curve.
	setup time.Duration

	// universal is PLONK's universal KZG SRS as far as the run has grown it: one
	// object per curve (see srs.go), so it is *not* this circuit's cost — every
	// PLONK row over the curve carries the same number, and the largest of them
	// is the whole of it. Zero for Groth16, which has no universal setup: its
	// whole setup is per circuit.
	universal time.Duration

	prove  func(constraint.ConstraintSystem, witness.Witness) (Proof, error)
	verify func(Proof, witness.Witness) error
}

// proofSystem is the backend-specific half of the pipeline. Implementations are
// empty structs, which keeps Backend comparable.
type proofSystem interface {
	setup(ccs constraint.ConstraintSystem, o Options) (keys, error)
}

// groth16System is Groth16 over R1CS: a per-circuit trusted setup, then Prove
// and Verify.
type groth16System struct{}

func (groth16System) setup(ccs constraint.ConstraintSystem, _ Options) (keys, error) {
	start := time.Now()
	pk, vk, err := groth16.Setup(ccs)
	if err != nil {
		return keys{}, fmt.Errorf("groth16 setup: %w", err)
	}
	return keys{
		pk:    pk,
		vk:    vk,
		setup: time.Since(start),
		prove: func(ccs constraint.ConstraintSystem, full witness.Witness) (Proof, error) {
			return groth16.Prove(ccs, pk, full)
		},
		verify: func(p Proof, public witness.Witness) error {
			proof, ok := p.(groth16.Proof)
			if !ok {
				return fmt.Errorf("harness: %T is not a groth16 proof", p)
			}
			return groth16.Verify(proof, vk, public)
		},
	}, nil
}

// plonkSystem is PLONK over gnark's sparse constraint system: a universal KZG
// SRS (one per curve, from the store in srs.go — fine for benchmarking, a
// deployment needs a ceremony), then per-circuit preprocessing.
type plonkSystem struct{}

func (plonkSystem) setup(ccs constraint.ConstraintSystem, o Options) (keys, error) {
	// The canonical SRS is the curve's, not this circuit's, so what comes back
	// timed as *setup* is only the part that is this circuit's: the Lagrange
	// basis over its domain, which cannot be shared with a circuit of another
	// size. See srs.go for why that split is where it is.
	s, err := o.srs().For(ccs)
	if err != nil {
		return keys{}, err
	}
	start := time.Now()
	pk, vk, err := plonk.Setup(ccs, s.canonical, s.lagrange)
	if err != nil {
		return keys{}, fmt.Errorf("plonk setup: %w", err)
	}
	return keys{
		pk:        pk,
		vk:        vk,
		setup:     time.Since(start) + s.derive,
		universal: s.universal,
		prove: func(ccs constraint.ConstraintSystem, full witness.Witness) (Proof, error) {
			return plonk.Prove(ccs, pk, full)
		},
		verify: func(p Proof, public witness.Witness) error {
			proof, ok := p.(plonk.Proof)
			if !ok {
				return fmt.Errorf("harness: %T is not a plonk proof", p)
			}
			return plonk.Verify(proof, vk, public)
		},
	}, nil
}

// Witness is one assignment of a Case, in the two forms the pipeline needs: the
// full witness the prover consumes and the public part the verifier does. Input
// and Output are kept so a caller can report or re-check what was proved.
type Witness struct {
	Full   witness.Witness
	Public witness.Witness
	Input  []*big.Int
	Output []*big.Int
}

// newWitness builds the witness for an explicit input/output pair. out must be
// what the circuit computes from in, or proving fails when the solver hits the
// circuit's assertion — use [Evaluate] (or [GeneratedWitness]) to obtain it.
func newWitness(c Case, in, out []*big.Int) (Witness, error) {
	assignment := c.Build(variables(in), variables(out))
	full, err := frontend.NewWitness(assignment, c.Field)
	if err != nil {
		return Witness{}, fmt.Errorf("harness: witness for %s: %w", c.Name, err)
	}
	public, err := full.Public()
	if err != nil {
		return Witness{}, fmt.Errorf("harness: public witness for %s: %w", c.Name, err)
	}
	return Witness{Full: full, Public: public, Input: in, Output: out}, nil
}

// GeneratedWitness builds a witness for a freshly generated input of the Case's
// canonical size, deriving the matching output with [Evaluate]. This is the
// default benchmark witness: it needs no reference vector, so it works for any
// target, and it is deterministic in seed, so a measured row is reproducible.
func GeneratedWitness(c Case, seed string) (Witness, error) {
	in := GenerateInput(c.Field, c.InputSize, seed)
	out, err := Evaluate(c, in)
	if err != nil {
		return Witness{}, err
	}
	if len(out) != c.OutputSize {
		return Witness{}, fmt.Errorf("harness: %s evaluated to %d outputs, case declares %d",
			c.Name, len(out), c.OutputSize)
	}
	return newWitness(c, in, out)
}

// vectorWitness builds the witness for a reference vector, so the resulting
// proof is a proof of a known-answer test. The vector's own lengths win over the
// Case's declared sizes, exactly as in Solve — a sponge vector may absorb more
// than one rate-sized block.
func vectorWitness(c Case, v Vector) (Witness, error) {
	in, err := ParseElements(v.Input)
	if err != nil {
		return Witness{}, err
	}
	out, err := ParseElements(v.Output)
	if err != nil {
		return Witness{}, err
	}
	return newWitness(c, in, out)
}

// GenerateInput returns n deterministic pseudo-random field elements, drawn from
// the repository's shared XOF sampler under a seed that is namespaced so it can
// never collide with a parameter-derivation stream. Deterministic rather than
// random because a benchmark row should be reproducible; the values themselves
// do not affect cost (the circuit is fixed) — only that they are a valid input.
func GenerateInput(field *big.Int, n int, seed string) []*big.Int {
	return sampler.SHAKE256Mod([]byte("gnark-hashes/harness/prove:"+seed), field, n)
}

// Evaluate runs the Case's in-circuit computation natively and returns its
// output, which is what turns an arbitrary input into a satisfying witness.
//
// The evaluator is gnark's test engine: it executes Define over big.Int instead
// of emitting constraints, and it resolves hints by calling them directly — so
// the type-2/3/4 designs (Anemoi, Rescue-Prime, Arion, Griffin, Skyscraper,
// Polocolo), whose S-boxes are witnessed by a hint and checked in-circuit,
// evaluate here just as they do in the vector tests. No hash is reimplemented
// out of circuit.
func Evaluate(c Case, in []*big.Int) ([]*big.Int, error) {
	if c.Eval == nil {
		return nil, fmt.Errorf("harness: case %s has no Eval, so only a reference-vector witness is available", c.Name)
	}
	spec := &evalSpec{fn: c.Eval, field: c.Field}

	// The stub carries the capture target; the test engine shallow-copies the
	// struct (so spec, a pointer, is shared) and fills In from the assignment.
	stub := &evalCircuit{In: make([]frontend.Variable, len(in)), Spec: spec}
	assignment := &evalCircuit{In: variables(in), Spec: spec}
	if err := test.IsSolved(stub, assignment, c.Field); err != nil {
		return nil, fmt.Errorf("harness: evaluating %s: %w", c.Name, err)
	}
	if spec.err != nil {
		return nil, fmt.Errorf("harness: evaluating %s: %w", c.Name, spec.err)
	}
	return spec.out, nil
}

// evalSpec is the in-circuit function plus the captured result. It is held by
// *pointer* from the circuit struct on purpose: the test engine clones a circuit
// with reflect.DeepEqual as a self-check, and DeepEqual considers two non-nil
// func values unequal — a func stored directly in the struct would fail the
// clone. A shared pointer compares equal, and doubles as the channel the
// computed output comes back on.
type evalSpec struct {
	fn    func(frontend.API, []frontend.Variable) []frontend.Variable
	field *big.Int
	out   []*big.Int
	err   error
}

// evalCircuit exposes Case.Eval as a circuit so the test engine can run it. It
// asserts nothing: it exists only to observe what the computation produces.
type evalCircuit struct {
	In []frontend.Variable

	Spec *evalSpec `gnark:"-"`
}

// Define implements frontend.Circuit.
func (c *evalCircuit) Define(api frontend.API) error {
	got := c.Spec.fn(api, c.In)
	out := make([]*big.Int, len(got))
	for i, v := range got {
		b, err := nativeValue(v, c.Spec.field)
		if err != nil {
			c.Spec.err = fmt.Errorf("output %d: %w", i, err)
			return nil
		}
		out[i] = b
	}
	c.Spec.out = out
	return nil
}

// nativeValue reduces a value produced by the test engine to a field element.
// The engine represents variables as *big.Int, but a circuit may also pass a
// literal straight through to its output, so the constant forms gnark accepts in
// an assignment are handled too.
func nativeValue(v frontend.Variable, field *big.Int) (*big.Int, error) {
	var b *big.Int
	switch x := v.(type) {
	case *big.Int:
		b = new(big.Int).Set(x)
	case big.Int:
		b = new(big.Int).Set(&x)
	case int:
		b = big.NewInt(int64(x))
	case int64:
		b = big.NewInt(x)
	case uint64:
		b = new(big.Int).SetUint64(x)
	case string:
		var ok bool
		if b, ok = new(big.Int).SetString(x, 0); !ok {
			return nil, fmt.Errorf("harness: invalid field element %q", x)
		}
	default:
		return nil, fmt.Errorf("harness: cannot read a %T as a field element", v)
	}
	return b.Mod(b, field), nil
}

// Prover is one Case compiled under one Backend with its setup already done —
// everything needed to produce and check proofs for that target. Build it once
// (setup is by far the most expensive step) and Prove repeatedly.
type Prover struct {
	Backend Backend
	Case    Case
	CCS     constraint.ConstraintSystem

	// Compile and Setup are how long those one-time steps took. SRS is PLONK's
	// universal setup, which is per *curve* and so excluded from Setup — one
	// number for the whole run rather than for this circuit (see keys.universal).
	Compile, Setup, SRS time.Duration

	keys keys
}

// NewProver compiles the Case for the Backend's constraint system and runs the
// proof system's setup over it.
func NewProver(b Backend, c Case, o Options) (*Prover, error) {
	start := time.Now()
	ccs, err := compile(b, c)
	if err != nil {
		return nil, err
	}
	compileTime := time.Since(start)

	// The setup times itself: which part of it is the circuit's and which the
	// curve's is the proof system's business, not this function's.
	k, err := b.proofs.setup(ccs, o)
	if err != nil {
		return nil, fmt.Errorf("harness: setup %s (%s): %w", c.Name, b.Name, err)
	}

	return &Prover{
		Backend: b,
		Case:    c,
		CCS:     ccs,
		Compile: compileTime,
		Setup:   k.setup,
		SRS:     k.universal,
		keys:    k,
	}, nil
}

// Prove produces a proof for the given witness.
func (p *Prover) Prove(w Witness) (Proof, error) {
	proof, err := p.keys.prove(p.CCS, w.Full)
	if err != nil {
		return nil, fmt.Errorf("harness: prove %s (%s): %w", p.Case.Name, p.Backend.Name, err)
	}
	return proof, nil
}

// Verify checks a proof against the witness's public part.
func (p *Prover) Verify(proof Proof, w Witness) error {
	if err := p.keys.verify(proof, w.Public); err != nil {
		return fmt.Errorf("harness: verify %s (%s): %w", p.Case.Name, p.Backend.Name, err)
	}
	return nil
}

// ProvingKey and VerifyingKey expose the setup output for sizing (see [Size]).
func (p *Prover) ProvingKey() io.WriterTo   { return p.keys.pk }
func (p *Prover) VerifyingKey() io.WriterTo { return p.keys.vk }

// Size is the serialized byte size of a proof or key, in gnark's canonical
// (point-compressed) encoding — the number to report as "proof size".
func Size(w io.WriterTo) (int64, error) {
	if w == nil {
		return 0, nil
	}
	return w.WriteTo(io.Discard)
}

// Result is one measured row: what a single (Case, Backend) pair costs from
// compilation through verification. Durations for the repeated steps are means
// over Options.Runs, with the fastest run alongside; the one-time steps are
// measured once, which is all they happen.
type Result struct {
	Case    string
	Backend string // "r1cs" or "plonk" — the constraint system
	System  string // "groth16" or "plonk" — the proof system over it

	Constraints int // identical to what Count reports for this backend
	PublicWires int
	SecretWires int

	Compile time.Duration
	SRS     time.Duration // universal KZG SRS; PLONK only, and per curve — see below
	Setup   time.Duration // circuit-specific preprocessing (universal SRS excluded)
	Witness time.Duration // input generation + native evaluation + witness build

	Prove     time.Duration
	ProveMin  time.Duration
	Verify    time.Duration
	VerifyMin time.Duration
	Runs      int

	ProofBytes int64
	PKBytes    int64
	VKBytes    int64
}

// Measure runs the whole pipeline for one Case under one Backend and reports it.
//
// Note on SRS: it is one object per curve, not per row (srs.go). Result.SRS is that
// curve's SRS as far as the run has grown it, so it is the same number for every
// PLONK row over the curve and the largest of them is the whole cost; it is not
// added to Setup, which keeps Setup comparable across rows and across proof
// systems. Report the universal setup once per run — [SRS.Universal] is that
// number — rather than beside the constructions, which never differ in it.
func Measure(b Backend, c Case, o Options) (Result, error) {
	res := Result{Case: c.Name, Backend: b.Name, System: b.System, Runs: o.runs()}

	start := time.Now()
	w, err := witnessFor(c, o)
	if err != nil {
		return res, err
	}
	res.Witness = time.Since(start)

	p, err := NewProver(b, c, o)
	if err != nil {
		return res, err
	}
	res.Compile, res.Setup, res.SRS = p.Compile, p.Setup, p.SRS
	res.Constraints = p.CCS.GetNbConstraints()
	res.PublicWires = p.CCS.GetNbPublicVariables()
	res.SecretWires = p.CCS.GetNbSecretVariables()

	var proof Proof
	for i := 0; i < res.Runs; i++ {
		start = time.Now()
		if proof, err = p.Prove(w); err != nil {
			return res, err
		}
		res.Prove, res.ProveMin = accumulate(res.Prove, res.ProveMin, time.Since(start), i)

		start = time.Now()
		if err = p.Verify(proof, w); err != nil {
			return res, err
		}
		res.Verify, res.VerifyMin = accumulate(res.Verify, res.VerifyMin, time.Since(start), i)
	}
	res.Prove /= time.Duration(res.Runs)
	res.Verify /= time.Duration(res.Runs)

	if res.ProofBytes, err = Size(proof); err != nil {
		return res, fmt.Errorf("harness: sizing proof for %s: %w", c.Name, err)
	}
	if res.PKBytes, err = Size(p.ProvingKey()); err != nil {
		return res, fmt.Errorf("harness: sizing proving key for %s: %w", c.Name, err)
	}
	if res.VKBytes, err = Size(p.VerifyingKey()); err != nil {
		return res, fmt.Errorf("harness: sizing verifying key for %s: %w", c.Name, err)
	}
	return res, nil
}

// witnessFor picks the witness source: the reference vector if one was given,
// otherwise a generated input.
func witnessFor(c Case, o Options) (Witness, error) {
	if o.Vector != nil {
		return vectorWitness(c, *o.Vector)
	}
	return GeneratedWitness(c, o.Seed)
}

// accumulate adds d to a running total and tracks the minimum. The first sample
// seeds the minimum, since a zero min would never be beaten.
func accumulate(total, min, d time.Duration, i int) (time.Duration, time.Duration) {
	if i == 0 || d < min {
		min = d
	}
	return total + d, min
}

func variables(xs []*big.Int) []frontend.Variable {
	out := make([]frontend.Variable, len(xs))
	for i, x := range xs {
		out[i] = x
	}
	return out
}
