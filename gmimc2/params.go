package gmimc2

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific instantiation of GMiMC2: the
// single source of truth every consumer (permutation, modes, tests) reads. It
// mirrors GMiMC2Params in gmimc2_ref.py, this construction's Python oracle.
//
// Only the numbers the specification fixes per instance are chosen (Field, T,
// Rounds, Alpha, and the sponge split Rate/Capacity/Digest — see instances.go);
// everything else is derived from them:
//
//   - M:    the cyclic-shift linear layer, determined by T alone. Stored as a
//     general t*t matrix of *big.Int (free in-circuit — a 0/1 matrix folds
//     into wire coefficients).
//   - MIO:  the input/output matrix, determined by T and Capacity (see ioMatrix).
//   - RoundConstants: R elements drawn deterministically from SHAKE256, so they
//     reproduce from (p, t, R, alpha) exactly (see deriveRoundConstants).
//
// Two things differ from the GMiMC parameters next door, and both are deliberate:
//
// Alpha is PASSED, not derived. GMiMC takes the smallest exponent coprime to p-1;
// GMiMC2 fixes alpha = 2^k, which is never coprime to p-1 over an odd field, so
// the S-box x -> x^(2^k) is not a bijection. That is sound because the round
// function is an expanding-round-function Feistel: the branch function is only
// ever evaluated forwards — a round is undone by subtracting F(x_0) from the
// branches it was added to — so F need not be invertible. alpha and the round
// count travel together, and both come from the instance.
//
// Rounds must be a multiple of T. Every round number the specification gives is
// the smallest multiple of t exceeding its security bound, and the improved t>4
// circuit relies on it (after R rounds the cyclic shift composes to the identity).
// Rejecting anything else keeps a typo from silently producing an instance the
// specification does not describe.
type Parameters struct {
	Field          *big.Int     // field characteristic p
	T              int          // state size (branches)
	Rounds         int          // number of rounds R
	Alpha          int          // power-map S-box exponent, a power of two
	Rate           int          // sponge rate r
	Capacity       int          // sponge capacity c
	Digest         int          // digest size d
	M              [][]*big.Int // t*t round linear layer (cyclic shift)
	MIO            [][]*big.Int // t*t input/output matrix, applied before and after the rounds
	RoundConstants []*big.Int   // >= R round constants
}

// NewParameters builds and validates GMiMC2 parameters for the given field,
// deriving M, MIO and RoundConstants.
func NewParameters(field *big.Int, t, rounds, alpha, rate, capacity, digest int) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("gmimc2: field modulus must be a positive integer")
	}
	if field.Bit(0) == 0 {
		return nil, fmt.Errorf("gmimc2: field modulus must be odd, got an even value")
	}
	if t <= 1 {
		return nil, fmt.Errorf("gmimc2: state size t must be > 1, got %d", t)
	}
	if rounds < 1 {
		return nil, fmt.Errorf("gmimc2: rounds must be >= 1, got %d", rounds)
	}
	if rounds%t != 0 {
		return nil, fmt.Errorf("gmimc2: rounds must be a multiple of t, got R=%d t=%d", rounds, t)
	}
	if alpha < 2 || alpha&(alpha-1) != 0 {
		return nil, fmt.Errorf("gmimc2: alpha must be 2^k with k >= 1, got %d", alpha)
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("gmimc2: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("gmimc2: digest size must be >= 1, got %d", digest)
	}

	mio, err := ioMatrix(t, capacity)
	if err != nil {
		return nil, fmt.Errorf("gmimc2: %w", err)
	}

	p := &Parameters{
		Field:          new(big.Int).Set(field),
		T:              t,
		Rounds:         rounds,
		Alpha:          alpha,
		Rate:           rate,
		Capacity:       capacity,
		Digest:         digest,
		M:              shiftMatrix(t),
		MIO:            mio,
		RoundConstants: deriveRoundConstants(field, t, rounds, alpha),
	}

	if err := checkSquareMatrix(p.M, t); err != nil {
		return nil, fmt.Errorf("gmimc2: matrix M: %w", err)
	}
	if err := checkSquareMatrix(p.MIO, t); err != nil {
		return nil, fmt.Errorf("gmimc2: matrix M_IO: %w", err)
	}
	if len(p.RoundConstants) < p.Rounds {
		return nil, fmt.Errorf("gmimc2: expected at least %d round constants, got %d", p.Rounds, len(p.RoundConstants))
	}
	return p, nil
}

// shiftMatrix builds the t*t cyclic-shift permutation matrix (out[i]=state[i+1],
// out[t-1]=state[0]) — the erf round's linear layer, unchanged from GMiMC.
func shiftMatrix(t int) [][]*big.Int {
	m := zeroMatrix(t)
	for i := 0; i < t-1; i++ {
		m[i][i+1] = big.NewInt(1)
	}
	m[t-1][0] = big.NewInt(1)
	return m
}

// ioMatrix returns M_IO, the matrix applied to the state before the first round
// and after the last: pi(x) = M_IO . R^(R) . ... . R^(1) (M_IO . x). It exists so
// that a DRL Groebner basis can be computed with a change of coordinates, and the
// specification picks it per capacity with minimal cost to performance:
//
//   - c = 1: the identity.
//   - t in {8,16}, c = r = t/2, compression mode: circ(1,0,...,0,2,0,...,0) with
//     the 2 in position t/2.
//   - t in {12,24}, c = t/3: circ(1,0,...,0,2,0,...,0,2,0,...,0) with the 2s in
//     positions t/3 and t/2.
//
// Only the first case is reachable here, and this returns an error for the others
// rather than building them: the two circulants belong to the specification's
// 32-bit and 64-bit round-number tables, whose instances need a Goldilocks or
// Mersenne31 field. Those are out of scope by the repo's field policy (see
// README, "Fields"), so a right-circulant builder here would be code no instance
// can reach and no vector can test. Adding one is the natural first step if that
// policy ever changes.
func ioMatrix(t, capacity int) ([][]*big.Int, error) {
	if capacity != 1 {
		return nil, fmt.Errorf("no M_IO defined for t=%d c=%d: only the c=1 (identity) case "+
			"is in scope; the specification's circulant M_IO belong to the 32/64-bit "+
			"instances, whose fields this repo does not implement over", t, capacity)
	}
	m := zeroMatrix(t)
	for i := 0; i < t; i++ {
		m[i][i] = big.NewInt(1)
	}
	return m, nil
}

// deriveRoundConstants draws R field elements from SHAKE256("GMiMC2(p,t,R,alpha)aff")
// with the shared "mod" sampler.
//
// The specification does not fix a constant derivation, so this follows the
// reference framework's GMiMC scheme ("GMiMC(p,t,R)aff", ref/gmimc/params.py) with
// alpha added to the seed. alpha is needed to identify an instance: the
// specification's 32-bit table gives one round number for all three alpha, so
// (p, t, R) alone would have three instances sharing a constant stream. The
// sampler itself is byte-exact against the reference's XOFFieldElementSampler —
// gmimc2_ref.py asserts that by reproducing GMiMC's own pinned constants.
func deriveRoundConstants(field *big.Int, t, rounds, alpha int) []*big.Int {
	seed := fmt.Appendf(nil, "GMiMC2(%s,%d,%d,%d)aff", field.String(), t, rounds, alpha)
	return sampler.SHAKE256Mod(seed, field, rounds)
}

func zeroMatrix(t int) [][]*big.Int {
	m := make([][]*big.Int, t)
	for i := range m {
		m[i] = make([]*big.Int, t)
		for j := range m[i] {
			m[i][j] = big.NewInt(0)
		}
	}
	return m
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
