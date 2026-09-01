package gmimc

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific instantiation of S-GMiMC:
// the single source of truth every consumer (permutation, modes, tests) reads.
// This mirrors GMiMCParams in the Python reference (params.py).
//
// Only a handful of values are chosen per instance (Field, T, Rounds, and the
// sponge split Rate/Capacity/Digest — see instances.go); everything else is
// *derived* from them by NewParameters, matching the reference's _init_* helpers:
//
//   - Alpha: smallest S-box exponent coprime to p-1.
//   - M:     the cyclic-shift linear layer, determined by T alone. Stored as a
//     general t*t matrix of *big.Int (works for any construction, and is
//     free in-circuit since a 0/1 matrix folds into wire coefficients).
//   - RoundConstants: R elements drawn deterministically from SHAKE256, so they
//     reproduce from (p, t, R) exactly (see deriveRoundConstants).
//
// M_inv / permutation_inv from the reference are intentionally omitted: a SNARK
// only ever evaluates the permutation forward, so the inverse would be dead code.
type Parameters struct {
	Field          *big.Int     // field characteristic p
	T              int          // state size (branches)
	Rounds         int          // number of rounds R
	Alpha          int          // power-map S-box exponent
	Rate           int          // sponge rate r
	Capacity       int          // sponge capacity c
	Digest         int          // digest size d
	M              [][]*big.Int // t*t linear layer
	RoundConstants []*big.Int   // >= R affine round constants
}

// NewParameters builds and validates S-GMiMC parameters for the given field,
// deriving Alpha, M and RoundConstants.
func NewParameters(field *big.Int, t, rounds, rate, capacity, digest int) (*Parameters, error) {
	// Hard checks (the recommendation-only warnings from params.py are dropped:
	// they guard against toy fields and only matter when evaluating in Python).
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("gmimc: field modulus must be a positive integer")
	}
	if field.Bit(0) == 0 {
		return nil, fmt.Errorf("gmimc: field modulus must be odd, got an even value")
	}
	if t <= 1 {
		return nil, fmt.Errorf("gmimc: state size t must be > 1, got %d", t)
	}
	if rounds < 1 {
		return nil, fmt.Errorf("gmimc: rounds must be >= 1, got %d", rounds)
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("gmimc: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("gmimc: digest size must be >= 1, got %d", digest)
	}

	rc, err := deriveRoundConstants(field, t, rounds)
	if err != nil {
		return nil, fmt.Errorf("gmimc: deriving round constants: %w", err)
	}

	p := &Parameters{
		Field:          new(big.Int).Set(field),
		T:              t,
		Rounds:         rounds,
		Alpha:          defaultAlpha(field),
		Rate:           rate,
		Capacity:       capacity,
		Digest:         digest,
		M:              shiftMatrix(t),
		RoundConstants: rc,
	}

	// Post-override validation.
	if p.Alpha < 2 {
		return nil, fmt.Errorf("gmimc: alpha must be >= 2, got %d", p.Alpha)
	}
	if err := checkSquareMatrix(p.M, t); err != nil {
		return nil, fmt.Errorf("gmimc: matrix M: %w", err)
	}
	if len(p.RoundConstants) < p.Rounds {
		return nil, fmt.Errorf("gmimc: expected at least %d round constants, got %d", p.Rounds, len(p.RoundConstants))
	}
	return p, nil
}

// defaultAlpha returns the smallest exponent >= 2 coprime to p-1, matching
// _init_alpha in the reference.
func defaultAlpha(field *big.Int) int {
	pm1 := new(big.Int).Sub(field, big.NewInt(1))
	one := big.NewInt(1)
	g := new(big.Int)
	for a := int64(2); ; a++ {
		g.GCD(nil, nil, big.NewInt(a), pm1)
		if g.Cmp(one) == 0 {
			return int(a)
		}
	}
}

// shiftMatrix builds the t*t cyclic-shift permutation matrix (out[i]=state[i+1],
// out[t-1]=state[0]), matching _init_M in the reference.
func shiftMatrix(t int) [][]*big.Int {
	m := make([][]*big.Int, t)
	for i := range m {
		m[i] = make([]*big.Int, t)
		for j := range m[i] {
			m[i][j] = big.NewInt(0)
		}
	}
	for i := 0; i < t-1; i++ {
		m[i][i+1] = big.NewInt(1)
	}
	m[t-1][0] = big.NewInt(1)
	return m
}

// deriveRoundConstants reproduces _init_rcons: R field elements drawn from
// SHAKE256("GMiMC(p,t,R)aff") with the shared "mod" sampler. Verified byte-exact
// against the Python reference (see params_test.go and sampler/sampler_test.go).
func deriveRoundConstants(field *big.Int, t, rounds int) ([]*big.Int, error) {
	seed := fmt.Appendf(nil, "GMiMC(%s,%d,%d)aff", field.String(), t, rounds)
	return sampler.SHAKE256Mod(seed, field, rounds), nil
}

func checkSquareMatrix(m [][]*big.Int, t int) error {
	if len(m) != t {
		return fmt.Errorf("expected %d rows, got %d", t, len(m))
	}
	for i, row := range m {
		if len(row) != t {
			return fmt.Errorf("row %d has %d columns, expected %d", i, len(row), t)
		}
	}
	return nil
}
