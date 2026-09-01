// The universal setup: one KZG SRS per *curve*, not one per circuit. Package doc:
// harness.go.
//
// PLONK's SRS is powers of a secret tau, so the SRS for the largest circuit
// *contains* every smaller one as a prefix — gnark asserts only
// `len(srs.Pk.G1) >= domain+3` and then takes exactly that prefix
// (`pk.Kzg.G1 = srs.Pk.G1[:vk.Size+3]`, backend/plonk/*/setup.go). One SRS per
// curve therefore serves every construction, width, mode and workload measured
// over that curve. That is what "universal" means here, and it is why the
// benchmark reports this cost once per run instead of as a column beside the
// constructions: per row it would only say which circuit happened to be the one
// that grew it.
//
// What is *not* universal is the second SRS plonk.Setup takes: the same tau
// re-expressed over the circuit's evaluation domain, required at exactly that
// domain's size (`len(srsLagrange.Pk.G1) != domain` is a hard error), so it can
// neither be truncated nor shared between two circuits of different size. It is
// derived data though — one inverse FFT and one fixed-base multiplication, no new
// randomness and no ceremony — so it is per-circuit preprocessing, and [SRS.For]
// times it into the setup rather than into the universal number. See
// [curveKZG.split] for how it is derived and what the alternative would cost.
//
// gnark ships test/unsafekzg for all of this, and the harness used to call it once
// per circuit. It memoizes on (curve, domain) and *regenerates* for the next domain
// rather than extending, so a run paid for as many SRSs as it had distinct domains
// and every row after the first of a size measured a warm memo of a few
// microseconds. Neither of those is the cost of anything. This file replaces it
// with one store per process, generated on first use and grown in place.
package harness

import (
	"crypto/rand"
	"fmt"
	"math/big"
	"sort"
	"sync"
	"time"

	"github.com/consensys/gnark-crypto/ecc"
	bls12381 "github.com/consensys/gnark-crypto/ecc/bls12-381"
	fr_bls12381 "github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	fft_bls12381 "github.com/consensys/gnark-crypto/ecc/bls12-381/fr/fft"
	kzg_bls12381 "github.com/consensys/gnark-crypto/ecc/bls12-381/kzg"
	"github.com/consensys/gnark-crypto/ecc/bn254"
	fr_bn254 "github.com/consensys/gnark-crypto/ecc/bn254/fr"
	fft_bn254 "github.com/consensys/gnark-crypto/ecc/bn254/fr/fft"
	kzg_bn254 "github.com/consensys/gnark-crypto/ecc/bn254/kzg"
	"github.com/consensys/gnark-crypto/kzg"
	"github.com/consensys/gnark/constraint"
)

// srsBlinding is how many powers of tau plonk.Setup needs beyond the circuit's
// domain: three, for the KZG opening of a blinded polynomial.
const srsBlinding = 3

// SRS is the universal KZG setup of one process: at most one canonical SRS per
// curve, generated when the first circuit over that curve asks for it and extended
// if a later circuit is larger. It is never regenerated, so a run builds each
// curve's powers of tau exactly once end to end.
//
// The zero value is ready to use, and use is safe from several goroutines.
type SRS struct {
	mu     sync.Mutex
	curves map[ecc.ID]*curveSRS
}

// defaultSRS is the store [Options] uses when a caller names none, which is what
// makes the zero Options valid and what lets every PLONK measurement in a process
// share a curve's SRS without threading one through. Process scope is the right
// scope: it is exactly what makes the SRS universal, and it is what gnark's own
// unsafekzg does (a package-level cache) minus the per-domain regeneration.
var defaultSRS = &SRS{}

// DefaultSRS returns that store, so a caller who ran measurements with the zero
// [Options] can still report what the universal setup cost — see [SRS.Universal],
// which is how cmd/bench reports it once per run.
func DefaultSRS() *SRS { return defaultSRS }

// curveSRS is one curve's canonical SRS, plus what it took to build.
//
// tau is the SRS's trapdoor, kept because extending an SRS means continuing its
// powers of tau. Keeping it is precisely what makes this setup *unsafe* — anyone
// holding tau can forge a proof — and that is the same bargain gnark's own
// test/unsafekzg makes: these SRSs exist to time a prover, never to secure one.
// A deployment needs a real ceremony; nothing here is written to disk.
type curveSRS struct {
	srs    kzg.SRS  // *kzg_bn254.SRS or *kzg_bls12381.SRS
	tau    *big.Int // the trapdoor, kept so the SRS can be extended
	points int      // powers of tau it holds

	// generate is the total time spent building it: the first generation plus
	// every extension, which together is what generating it at its final size
	// would have cost.
	generate time.Duration
}

// setupSRS is what [SRS.For] hands the PLONK backend: the two SRSs plonk.Setup
// takes, and the two durations that have to be reported differently.
type setupSRS struct {
	canonical, lagrange kzg.SRS

	// derive is the Lagrange-basis derivation: per circuit, so it belongs in the
	// setup measurement.
	derive time.Duration

	// universal is the curve's universal SRS as far as the run has grown it. It is
	// not this circuit's cost — every PLONK row over the curve reports the same
	// number, and the largest of them is the whole of it.
	universal time.Duration
}

// For returns what plonk.Setup needs for ccs: the canonical SRS cut to the
// circuit's domain, and the Lagrange basis over that domain.
//
// Nothing is compiled or sized ahead of time. The domain comes from the constraint
// system the caller has already compiled, the store grows only when this circuit is
// bigger than everything before it, and a curve that no measurement touches never
// gets an SRS at all — so commenting a field out of the work list removes its setup
// with it.
func (s *SRS) For(ccs constraint.ConstraintSystem) (setupSRS, error) {
	id, ops, err := kzgFor(ccs.Field())
	if err != nil {
		return setupSRS{}, err
	}
	// gnark's own rule for the circuit's evaluation domain: the placeholder
	// constraints for the public inputs count towards it (backend/plonk/*/setup.go).
	domain := ecc.NextPowerOfTwo(uint64(ccs.GetNbConstraints() + ccs.GetNbPublicVariables()))

	entry, err := s.canonical(id, ops, domain+srsBlinding)
	if err != nil {
		return setupSRS{}, err
	}
	start := time.Now()
	canonical, lagrange := ops.split(entry.srs, domain, entry.tau)
	return setupSRS{
		canonical: canonical,
		lagrange:  lagrange,
		derive:    time.Since(start),
		universal: entry.generate,
	}, nil
}

// canonical returns the curve's SRS with at least points powers of tau, generating
// or extending it first if it is short.
//
// The lock is held across generation, which serializes two goroutines that want the
// same curve. That is the point: the alternative is both of them building an SRS,
// which is the waste this store exists to remove. Measurements run sequentially
// anyway.
func (s *SRS) canonical(id ecc.ID, ops curveKZG, points uint64) (curveSRS, error) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.curves == nil {
		s.curves = make(map[ecc.ID]*curveSRS)
	}

	e, ok := s.curves[id]
	if !ok {
		tau, err := rand.Int(rand.Reader, id.ScalarField())
		if err != nil {
			return curveSRS{}, fmt.Errorf("harness: sampling an SRS trapdoor for %s: %w", id, err)
		}
		start := time.Now()
		srs, err := ops.generate(points, tau)
		if err != nil {
			return curveSRS{}, fmt.Errorf("harness: %s SRS of %d powers: %w", id, points, err)
		}
		e = &curveSRS{srs: srs, tau: tau, points: int(points), generate: time.Since(start)}
		s.curves[id] = e
		return *e, nil
	}

	if e.points < int(points) {
		start := time.Now()
		srs, err := ops.extend(e.srs, points, e.tau)
		if err != nil {
			return curveSRS{}, fmt.Errorf("harness: growing the %s SRS from %d to %d powers: %w",
				id, e.points, points, err)
		}
		e.srs, e.points = srs, int(points)
		e.generate += time.Since(start)
	}
	return *e, nil
}

// UniversalSetup is one curve's universal KZG SRS, as a run built it.
type UniversalSetup struct {
	Curve  ecc.ID
	Domain uint64 // the largest circuit domain it covers
	Points int    // powers of tau held: Domain + srsBlinding

	// Generate is the whole cost of the SRS: one number per curve for the run,
	// against the per-circuit preprocessing that lands in [Result.Setup].
	Generate time.Duration
}

// Universal reports the universal SRS of every curve this store has served, in
// curve order. This is the setup cost to report *once* for a run; a curve nothing
// measured is absent rather than zero.
func (s *SRS) Universal() []UniversalSetup {
	s.mu.Lock()
	defer s.mu.Unlock()
	out := make([]UniversalSetup, 0, len(s.curves))
	for id, e := range s.curves {
		out = append(out, UniversalSetup{
			Curve:    id,
			Domain:   uint64(e.points - srsBlinding),
			Points:   e.points,
			Generate: e.generate,
		})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Curve.String() < out[j].Curve.String() })
	return out
}

// ---------------------------------------------------------------------------
// The curve-specific half
// ---------------------------------------------------------------------------

// curveKZG is everything above that has to name a concrete curve type. The store's
// own logic — when to grow, what to hand plonk.Setup, what to time — is written
// once against this interface, and each implementation is the dozen lines that
// cannot be: gnark-crypto's KZG types are generated per curve and share no
// interface beyond opaque serialization, which is why gnark's own plonk.Setup is a
// type switch too.
type curveKZG interface {
	// generate returns a fresh SRS holding `points` powers of tau.
	generate(points uint64, tau *big.Int) (kzg.SRS, error)

	// extend returns srs with its powers of tau continued to `points`, reusing
	// the ones it already holds. points must exceed them.
	extend(srs kzg.SRS, points uint64, tau *big.Int) (kzg.SRS, error)

	// split cuts srs down to one circuit: the canonical prefix plonk.Setup
	// slices, and the Lagrange basis it demands at exactly this domain.
	//
	// The Lagrange basis is [L_0(tau)]G .. [L_{domain-1}(tau)]G, built from the
	// trapdoor: inverse-FFT the powers of tau as *field* elements, then one
	// fixed-base batch multiplication. That is what gnark's unsafekzg does, and it
	// costs about what generating the same number of canonical points costs, since
	// it is the same batch multiplication.
	//
	// gnark-crypto also exposes kzg.ToLagrangeG1, which needs no trapdoor and is
	// therefore the route a real deployment has to take. It is not used here: it
	// runs the FFT over curve *points* and then one full scalar multiplication per
	// point, which is ~10x slower (bn254, measured on one machine):
	//
	//	domain    from tau    ToLagrangeG1
	//	  2^12        36ms           158ms
	//	  2^16       428ms         5.653s
	//	  2^19      5.290s        51.488s
	//
	// Either way this is per-circuit work and it is charged to the setup, but note
	// what it is a function of: (curve, domain) and nothing else. Two constructions
	// whose gate counts round to the same power of two pay exactly the same for it,
	// so it is the part of the setup column that ranks circuit *sizes* rather than
	// constructions. That is why the cheap route is the one to use — at 2^19,
	// plonk.Setup itself is ~3s, so ToLagrangeG1 would leave that column ~94% a
	// power-of-two readout and ~6% anything a construction did.
	split(srs kzg.SRS, domain uint64, tau *big.Int) (canonical, lagrange kzg.SRS)
}

// curves is the field policy in executable form (README, "Fields"): a circuit is
// natively arithmetized only in its own curve's scalar field, so these two are the
// only fields anything in this repository compiles over, and a third would need an
// entry here rather than field emulation.
var curves = map[ecc.ID]curveKZG{
	ecc.BN254:     bn254KZG{},
	ecc.BLS12_381: bls12381KZG{},
}

// kzgFor identifies a constraint system's curve from its scalar field.
func kzgFor(field *big.Int) (ecc.ID, curveKZG, error) {
	for id, ops := range curves {
		if field.Cmp(id.ScalarField()) == 0 {
			return id, ops, nil
		}
	}
	return ecc.UNKNOWN, nil, fmt.Errorf("harness: no KZG setup for the field with modulus %s "+
		"(this repository implements the BN254 and BLS12-381 scalar fields only)", field)
}

type bn254KZG struct{}

func (bn254KZG) generate(points uint64, tau *big.Int) (kzg.SRS, error) {
	return kzg_bn254.NewSRS(points, tau)
}

func (bn254KZG) extend(srs kzg.SRS, points uint64, tau *big.Int) (kzg.SRS, error) {
	have := srs.(*kzg_bn254.SRS)
	n := uint64(len(have.Pk.G1))

	// The powers this SRS is missing, tau^n .. tau^(points-1), continued from
	// where the ones it has stop. This is the whole reason the trapdoor is kept:
	// powers of tau are incremental, so growing an SRS costs only the new points.
	var alpha fr_bn254.Element
	alpha.SetBigInt(tau)
	scalars := make([]fr_bn254.Element, points-n)
	scalars[0].Exp(alpha, new(big.Int).SetUint64(n))
	for i := 1; i < len(scalars); i++ {
		scalars[i].Mul(&scalars[i-1], &alpha)
	}
	_, _, g1, _ := bn254.Generators()
	tail := bn254.BatchScalarMultiplicationG1(&g1, scalars)

	// A fresh array rather than append: plonk.Setup keeps a sub-slice of this one
	// inside every proving key already produced from it, and a live prover must
	// not have its SRS move underneath it.
	grown := make([]bn254.G1Affine, points)
	copy(grown, have.Pk.G1)
	copy(grown[n:], tail)
	return &kzg_bn254.SRS{Vk: have.Vk, Pk: kzg_bn254.ProvingKey{G1: grown}}, nil
}

func (bn254KZG) split(srs kzg.SRS, domain uint64, tau *big.Int) (kzg.SRS, kzg.SRS) {
	full := srs.(*kzg_bn254.SRS)
	canonical := &kzg_bn254.SRS{
		Vk: full.Vk,
		Pk: kzg_bn254.ProvingKey{G1: full.Pk.G1[:domain+srsBlinding]},
	}

	// [L_0(tau),..,L_{n-1}(tau)] = FFT_inv(sum_j tau^j X^j), so the basis is an
	// inverse FFT of the powers of tau followed by one fixed-base multiplication.
	var alpha fr_bn254.Element
	alpha.SetBigInt(tau)
	powers := make([]fr_bn254.Element, domain)
	powers[0].SetOne()
	for i := 1; i < len(powers); i++ {
		powers[i].Mul(&powers[i-1], &alpha)
	}
	fft_bn254.NewDomain(domain).FFTInverse(powers, fft_bn254.DIF)
	fft_bn254.BitReverse(powers) //nolint:staticcheck // the replacement is generic over a type this call cannot name
	_, _, g1, _ := bn254.Generators()

	return canonical, &kzg_bn254.SRS{
		Vk: full.Vk,
		Pk: kzg_bn254.ProvingKey{G1: bn254.BatchScalarMultiplicationG1(&g1, powers)},
	}
}

type bls12381KZG struct{}

func (bls12381KZG) generate(points uint64, tau *big.Int) (kzg.SRS, error) {
	return kzg_bls12381.NewSRS(points, tau)
}

func (bls12381KZG) extend(srs kzg.SRS, points uint64, tau *big.Int) (kzg.SRS, error) {
	have := srs.(*kzg_bls12381.SRS)
	n := uint64(len(have.Pk.G1))

	var alpha fr_bls12381.Element
	alpha.SetBigInt(tau)
	scalars := make([]fr_bls12381.Element, points-n)
	scalars[0].Exp(alpha, new(big.Int).SetUint64(n))
	for i := 1; i < len(scalars); i++ {
		scalars[i].Mul(&scalars[i-1], &alpha)
	}
	_, _, g1, _ := bls12381.Generators()
	tail := bls12381.BatchScalarMultiplicationG1(&g1, scalars)

	grown := make([]bls12381.G1Affine, points)
	copy(grown, have.Pk.G1)
	copy(grown[n:], tail)
	return &kzg_bls12381.SRS{Vk: have.Vk, Pk: kzg_bls12381.ProvingKey{G1: grown}}, nil
}

func (bls12381KZG) split(srs kzg.SRS, domain uint64, tau *big.Int) (kzg.SRS, kzg.SRS) {
	full := srs.(*kzg_bls12381.SRS)
	canonical := &kzg_bls12381.SRS{
		Vk: full.Vk,
		Pk: kzg_bls12381.ProvingKey{G1: full.Pk.G1[:domain+srsBlinding]},
	}

	var alpha fr_bls12381.Element
	alpha.SetBigInt(tau)
	powers := make([]fr_bls12381.Element, domain)
	powers[0].SetOne()
	for i := 1; i < len(powers); i++ {
		powers[i].Mul(&powers[i-1], &alpha)
	}
	fft_bls12381.NewDomain(domain).FFTInverse(powers, fft_bls12381.DIF)
	fft_bls12381.BitReverse(powers) //nolint:staticcheck // as above
	_, _, g1, _ := bls12381.Generators()

	return canonical, &kzg_bls12381.SRS{
		Vk: full.Vk,
		Pk: kzg_bls12381.ProvingKey{G1: bls12381.BatchScalarMultiplicationG1(&g1, powers)},
	}
}
