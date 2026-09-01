// Package poseidon implements the Poseidon permutation (a Hades-strategy design:
// external rounds with a full-width power-map S-box, internal rounds with a
// single-branch power map, one shared MDS matrix) as a gnark circuit. See
// ref/hades for the reference (PoseidonParams / the Poseidon class).
package poseidon

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/algebra"
	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific instantiation of Poseidon —
// the single source of truth the permutation reads. Everything below the chosen
// numbers (field, t, alpha, R_ext, R_int, r/c/d, version) is derived by
// NewParameters via the Grain LFSR, matching hades/params.py PoseidonParams. A
// SNARK only ever evaluates the permutation forward, so everything here is
// forward data.
type Parameters struct {
	Field   *big.Int // field characteristic p
	T       int      // state size (branches)
	Alpha   int      // power-map S-box exponent
	RExt    int      // external (full) rounds
	RInt    int      // internal (partial) rounds
	RExtBeg int      // external rounds before the internal block (R_ext/2)
	RExtEnd int      // external rounds after
	R       int      // RExt + RInt
	U       int      // branches the internal S-box hits (1)
	Version string   // Grain seeding version: "isec" or "circom"

	Rate     int
	Capacity int
	Digest   int

	M              [][]*big.Int // single MDS matrix, used for both external and internal rounds
	RoundConstants [][]*big.Int // R rows, full width: added to every branch before each S-box

	// Internal is the internal-round block rewritten with sparse matrices — the
	// Poseidon paper's Appendix B factorisation, the optimisation circomlib ships
	// and the one thing this construction was measured without. Each internal
	// round's dense MDS product, t(t-1) PLONK gates, becomes 2t-2, and the
	// residual costs nothing because it folds into the MDS matrix of the external
	// round that feeds the block (algebra.PartialRounds explains the rewrite and
	// where the residual can go). It saves R_int*(t-1)(t-2) gates: 114 of the 639
	// at t=3, and quadratically more as t grows.
	//
	// nil at t=2, where the rewritten layer costs exactly what the dense 2x2
	// product costs, so there is nothing to win. Nothing changes in R1CS, and the
	// KAT vectors are unmoved: the rewrite is an identity.
	Internal *algebra.PartialRounds
}

// NewParameters builds and validates Poseidon parameters for the given field,
// deriving the round constants (Grain rejection sampling) and the Cauchy MDS
// matrix (Grain reduction sampling, off the same stream). version selects the
// Grain seed layout: "isec" (the reference implementation) or "circom".
func NewParameters(field *big.Int, t, alpha, rExt, rInt, rate, capacity, digest int, version string) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("poseidon: field modulus must be a positive integer")
	}
	if t < 2 {
		return nil, fmt.Errorf("poseidon: state size t must be >= 2, got %d", t)
	}
	if alpha < 3 {
		return nil, fmt.Errorf("poseidon: alpha must be >= 3, got %d", alpha)
	}
	if rExt < 2 || rExt%2 != 0 {
		return nil, fmt.Errorf("poseidon: R_ext must be even and >= 2, got %d", rExt)
	}
	if rInt < 1 {
		return nil, fmt.Errorf("poseidon: R_int must be >= 1, got %d", rInt)
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("poseidon: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("poseidon: digest size must be >= 1, got %d", digest)
	}
	if version != "isec" && version != "circom" {
		return nil, fmt.Errorf("poseidon: unknown version %q (want isec or circom)", version)
	}

	p := &Parameters{
		Field: new(big.Int).Set(field),
		T:     t, Alpha: alpha,
		RExt: rExt, RInt: rInt,
		RExtBeg: rExt / 2, RExtEnd: rExt - rExt/2,
		R: rExt + rInt, U: 1, Version: version,
		Rate: rate, Capacity: capacity, Digest: digest,
	}

	// Draw order must match the reference: full-width round constants first
	// (rejection sampling), then the MDS xs/ys (reduction sampling) off the same
	// Grain stream. A single MDS serves both external and internal rounds.
	g := sampler.NewPoseidonGrain(field, t, rExt, rInt, version)
	p.RoundConstants = g.Grid(p.R, t)
	g.SwitchToMod()
	p.M = sampledCauchyMDS(g, field, t)

	if err := p.initInternal(); err != nil {
		return nil, err
	}

	return p, nil
}

// initInternal rewrites the internal-round block with sparse matrices, keeping
// the result only where it beats applying the MDS matrix directly in every one
// of those rounds (see Parameters.Internal and algebra.PartialRounds).
//
// The block is the R_int rounds from RExtBeg on, each "S-box on branch 0, then
// M*x + the next round's constants"; the layer that feeds it is the last
// external round before it, which applies the same M and the block's first
// constants. R_ext >= 2 guarantees that round exists.
func (p *Parameters) initInternal() error {
	block := algebra.PartialBlock{
		M:         p.M,
		U:         p.U,
		Consts:    p.RoundConstants[p.RExtBeg+1 : p.RExtBeg+p.RInt+1],
		Pre:       p.M,
		PreConsts: p.RoundConstants[p.RExtBeg],
		Field:     p.Field,
	}
	f, err := algebra.SelectPartialRounds(block, algebra.DenseGates(p.M))
	if err != nil {
		return fmt.Errorf("poseidon: internal-round factorisation: %w", err)
	}
	p.Internal = f
	return nil
}

func (p *Parameters) isInternal(r int) bool {
	return p.RExtBeg <= r && r < p.RExtBeg+p.RInt
}

// ---------------------------------------------------------------------------
// Cauchy MDS (matches cauchy_mds_matrix with a sampler: the isec "sampled" strategy)
// ---------------------------------------------------------------------------

// sampledCauchyMDS draws xs, ys as a 2-row Grain grid and builds the Cauchy
// matrix M[i][j] = 1 / (xs[i] + ys[j]) over the field, resampling the whole
// batch until the 2t values are distinct and no denominator vanishes. Matches
// cauchy_mds_matrix(..., sampler=...) in utils/matrix.py.
func sampledCauchyMDS(g *sampler.Grain, field *big.Int, t int) [][]*big.Int {
	for {
		grid := g.Grid(2, t)
		xs, ys := grid[0], grid[1]
		if !allDistinct(xs, ys) {
			continue
		}
		if m, ok := cauchyMatrix(xs, ys, field); ok {
			return m
		}
	}
}

// cauchyMatrix returns M[i][j] = (xs[i] + ys[j])^{-1} mod p, or ok=false if any
// denominator is zero mod p (the Cauchy matrix is then undefined for these xs, ys).
func cauchyMatrix(xs, ys []*big.Int, field *big.Int) ([][]*big.Int, bool) {
	m := make([][]*big.Int, len(xs))
	for i, x := range xs {
		m[i] = make([]*big.Int, len(ys))
		for j, y := range ys {
			denom := new(big.Int).Add(x, y)
			denom.Mod(denom, field)
			if denom.Sign() == 0 {
				return nil, false
			}
			m[i][j] = new(big.Int).ModInverse(denom, field)
		}
	}
	return m, true
}

// allDistinct reports whether the concatenation xs || ys has no repeated value.
func allDistinct(xs, ys []*big.Int) bool {
	seen := make(map[string]struct{}, len(xs)+len(ys))
	for _, v := range append(append([]*big.Int{}, xs...), ys...) {
		k := v.String()
		if _, dup := seen[k]; dup {
			return false
		}
		seen[k] = struct{}{}
	}
	return true
}
