// Package neptune implements the Neptune permutation (a Hades-strategy design:
// external rounds with a Lai-Massey pair-wise S-box, internal rounds with a
// partial power map) as a gnark circuit. See ref/hades for the reference.
package neptune

import (
	"crypto/sha3"
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/algebra"
)

// Parameters is the fully-expanded, field-specific instantiation of Neptune —
// the single source of truth the permutation reads. Everything below the chosen
// numbers (field, t, alpha, R_ext, R_int, r/c/d) is derived by NewParameters via
// the SHAKE128 sampler, matching hades/params.py NeptuneParams. A SNARK only ever
// evaluates the permutation forward, so everything here is forward data.
type Parameters struct {
	Field   *big.Int // field characteristic p
	T       int      // state size (branches), even
	Alpha   int      // internal-round power-map exponent
	RExt    int      // external (full) rounds
	RInt    int      // internal (partial) rounds
	RExtBeg int      // external rounds before the internal block (R_ext/2)
	RExtEnd int      // external rounds after
	R       int      // RExt + RInt
	U       int      // branches the internal S-box hits

	Rate     int
	Capacity int
	Digest   int

	MExt           [][]*big.Int // external-round matrix (t x t)
	MInt           [][]*big.Int // internal-round matrix J + diag(mu)
	RoundConstants [][]*big.Int // R+1 rows: row 0 is zero, rows 1..R-1 are ARK'd in the loop, row R whitens the output

	// SlpInt is the shared-sum addition program for M_int, verified against it by
	// algebra.SelectSLP and nil at t=2 where the dense product is already as
	// cheap. Forming s = sum_j x_j once gives every output as s + (mu_i - 1)*x_i,
	// 2t-1 additions against t(t-1). No program for M_ext: its even/odd split
	// leaves only t/2 nonzeros per row, which the dense product already realizes
	// at t in {2, 4}.
	//
	// It is the fallback, not what the registered instances run: where Internal
	// applies, the internal rounds go through that instead and this program is
	// only what they would have cost.
	SlpInt *algebra.SLP

	// Internal is the internal-round block rewritten with sparse matrices, the
	// factorisation algebra.PartialRounds derives: one dense row plus t-1 two-term
	// rows per internal round, 2t-2 gates against the 2t-1 of SlpInt and the
	// t(t-1) of the dense product. The residual folds into M_ext for the external
	// round that feeds the block, which is the one place it costs anything —
	// M_ext's rows have only t/2 nonzeros, so the product with the residual is
	// (t-1)t/2 gates dearer, 6 at t=4 against the 68 the rounds save.
	//
	// nil at t=2, where the rewritten layer costs exactly what the dense 2x2
	// product costs and only the residual would be left to pay. Nothing changes in
	// R1CS, and the KAT vectors are unmoved: the rewrite is an identity.
	Internal *algebra.PartialRounds

	// Lai-Massey external S-box constants (open Flystel, alpha=beta=1 in the spec).
	LMAlpha  *big.Int
	LMBeta   *big.Int
	LMGamma  *big.Int
	LMMatrix [][]*big.Int // 2x2
}

// NewParameters builds and validates Neptune parameters for the given field,
// deriving the matrices, round constants and Lai-Massey gamma from SHAKE128.
func NewParameters(field *big.Int, t, alpha, rExt, rInt, rate, capacity, digest int) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("neptune: field modulus must be a positive integer")
	}
	if t < 2 || t%2 != 0 {
		return nil, fmt.Errorf("neptune: state size t must be even and >= 2, got %d", t)
	}
	if alpha < 3 {
		return nil, fmt.Errorf("neptune: alpha must be >= 3, got %d", alpha)
	}
	if rExt < 2 || rExt%2 != 0 {
		return nil, fmt.Errorf("neptune: R_ext must be even and >= 2, got %d", rExt)
	}
	if rInt < 1 {
		return nil, fmt.Errorf("neptune: R_int must be >= 1, got %d", rInt)
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("neptune: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("neptune: digest size must be >= 1, got %d", digest)
	}

	p := &Parameters{
		Field: new(big.Int).Set(field),
		T:     t, Alpha: alpha,
		RExt: rExt, RInt: rInt,
		RExtBeg: rExt / 2, RExtEnd: rExt - rExt/2,
		R: rExt + rInt, U: 1,
		Rate: rate, Capacity: capacity, Digest: digest,
		LMAlpha:  big.NewInt(1),
		LMBeta:   big.NewInt(1),
		LMMatrix: [][]*big.Int{{big.NewInt(2), big.NewInt(1)}, {big.NewInt(1), big.NewInt(3)}},
	}

	// Draw order must match the reference: round constants, then M_ext, then the
	// internal-matrix diagonal, then the Lai-Massey gamma — off one shared stream.
	s := newSampler(field)
	p.RoundConstants = buildRoundConstants(s, p.R, t)
	p.MExt = buildMExt(s, t)
	p.MInt = buildMInt(s, t)
	p.LMGamma = s.nextNonzero()

	slp, ok := algebra.JPlusDiagSLP(p.MInt, p.Field)
	if !ok {
		return nil, fmt.Errorf("neptune: internal matrix is not of the form J+diag")
	}
	var err error
	if p.SlpInt, err = algebra.SelectSLP(slp, p.MInt, p.Field); err != nil {
		return nil, fmt.Errorf("neptune: internal matrix program: %w", err)
	}
	if err := p.initInternal(); err != nil {
		return nil, err
	}

	return p, nil
}

// initInternal rewrites the internal-round block with sparse matrices, keeping
// the result only where it beats applying M_int in every one of those rounds —
// as the shared-sum program where there is one, as the dense product where there
// is not (see Parameters.Internal and algebra.PartialRounds).
//
// The block is the R_int rounds from RExtBeg on, each "power map on branch 0,
// then M_int*x + the next round's constants"; the layer that feeds it is the last
// external round before it, which applies M_ext and the block's first constants.
// R_ext >= 2 guarantees that round exists.
func (p *Parameters) initInternal() error {
	perRound := algebra.DenseGates(p.MInt)
	if p.SlpInt != nil {
		perRound = p.SlpInt.Gates()
	}
	block := algebra.PartialBlock{
		M:         p.MInt,
		U:         p.U,
		Consts:    p.RoundConstants[p.RExtBeg+1 : p.RExtBeg+p.RInt+1],
		Pre:       p.MExt,
		PreConsts: p.RoundConstants[p.RExtBeg],
		Field:     p.Field,
	}
	f, err := algebra.SelectPartialRounds(block, perRound)
	if err != nil {
		return fmt.Errorf("neptune: internal-round factorisation: %w", err)
	}
	p.Internal = f
	return nil
}

// ---------------------------------------------------------------------------
// Derivations (consume the sampler in the reference's order)
// ---------------------------------------------------------------------------

// buildRoundConstants returns R+1 rows: a leading zero row (round 0 adds nothing
// before its S-box) followed by R sampled rows; the last is the output-whitening
// constant. Matches NeptuneParams._init_rcons.
func buildRoundConstants(s *sampler, r, t int) [][]*big.Int {
	rc := make([][]*big.Int, r+1)
	rc[0] = zeroRow(t)
	sampled := s.grid(r, t)
	for i := 0; i < r; i++ {
		rc[i+1] = sampled[i]
	}
	return rc
}

// buildMExt returns the even/odd "split" external matrix: M[2i][2j]=M'[i][j],
// M[2i+1][2j+1]=M”[i][j], else zero. M', M” are fixed circulants for t in
// {4,8} and otherwise sampled (M' then M”). Matches NeptuneParams._init_M_ext.
func buildMExt(s *sampler, t int) [][]*big.Int {
	half := t / 2
	var mp, mpp [][]*big.Int
	switch t {
	case 4:
		mp, mpp = circulant([]int64{2, 1}), circulant([]int64{1, 2})
	case 8:
		mp, mpp = circulant([]int64{3, 2, 1, 1}), circulant([]int64{1, 1, 2, 3})
	default:
		mp = s.grid(half, half)
		mpp = s.grid(half, half)
	}
	m := make([][]*big.Int, t)
	for i := range m {
		m[i] = zeroRow(t)
	}
	for i := 0; i < half; i++ {
		for j := 0; j < half; j++ {
			m[2*i][2*j] = mp[i][j]
			m[2*i+1][2*j+1] = mpp[i][j]
		}
	}
	return m
}

// buildMInt returns the internal matrix J + diag(mu-1), i.e. all-ones with the
// diagonal replaced by the sampled nonzero mu. Matches NeptuneParams._init_M_int
// (ones_plus_diag_matrix(mu-1)).
func buildMInt(s *sampler, t int) [][]*big.Int {
	mu := make([]*big.Int, t)
	for i := range mu {
		mu[i] = s.nextNonzero()
	}
	m := make([][]*big.Int, t)
	for i := 0; i < t; i++ {
		m[i] = make([]*big.Int, t)
		for j := 0; j < t; j++ {
			if i == j {
				m[i][j] = new(big.Int).Set(mu[i])
			} else {
				m[i][j] = big.NewInt(1)
			}
		}
	}
	return m
}

// circulant builds the right-circulant matrix from its first row:
// M[i][j] = row[(j-i) mod n]. Matches utils.matrix.circulant.
func circulant(row []int64) [][]*big.Int {
	n := len(row)
	m := make([][]*big.Int, n)
	for i := 0; i < n; i++ {
		m[i] = make([]*big.Int, n)
		for j := 0; j < n; j++ {
			m[i][j] = big.NewInt(row[((j-i)%n+n)%n])
		}
	}
	return m
}

func zeroRow(n int) []*big.Int {
	r := make([]*big.Int, n)
	for i := range r {
		r[i] = big.NewInt(0)
	}
	return r
}

// ---------------------------------------------------------------------------
// Deterministic field-element sampler: SHAKE128 with "bitmask" rejection
// sampling, reproducing utils.sampler.XOFFieldElementSampler for Neptune
// (seed = "Neptune" || p little-endian, padded to a multiple of 8 bytes).
// ---------------------------------------------------------------------------

type sampler struct {
	shake  *sha3.SHAKE
	field  *big.Int
	nBytes int
	mask   byte
	buf    []byte
}

func newSampler(field *big.Int) *sampler {
	bits := field.BitLen()
	// seed: "Neptune" || p as little-endian bytes, width rounded up to 64-bit words.
	nSeed := ((bits + 63) / 64) * 8
	be := make([]byte, nSeed)
	field.FillBytes(be)
	seed := append([]byte("Neptune"), reversed(be)...)

	h := sha3.NewSHAKE128()
	h.Write(seed)

	mask := byte(0xFF)
	if r := bits % 8; r != 0 {
		mask = byte((1 << r) - 1)
	}
	return &sampler{shake: h, field: field, nBytes: (bits + 7) / 8, mask: mask, buf: make([]byte, (bits+7)/8)}
}

// draw reads one raw candidate: nBytes of stream, little-endian, top byte masked
// to the field's bit length.
func (s *sampler) draw() *big.Int {
	if _, err := s.shake.Read(s.buf); err != nil {
		panic(fmt.Sprintf("neptune: shake read: %v", err)) // SHAKE never errors
	}
	b := make([]byte, s.nBytes)
	copy(b, s.buf)
	b[s.nBytes-1] &= s.mask
	return new(big.Int).SetBytes(reversed(b)) // little-endian -> big.Int
}

// next returns the next field element, rejecting candidates >= p.
func (s *sampler) next() *big.Int {
	for {
		if v := s.draw(); v.Cmp(s.field) < 0 {
			return v
		}
	}
}

// nextNonzero returns the next nonzero field element.
func (s *sampler) nextNonzero() *big.Int {
	for {
		if v := s.next(); v.Sign() != 0 {
			return v
		}
	}
}

// grid samples a rows x cols grid, row-major.
func (s *sampler) grid(rows, cols int) [][]*big.Int {
	out := make([][]*big.Int, rows)
	for i := range out {
		out[i] = make([]*big.Int, cols)
		for j := range out[i] {
			out[i][j] = s.next()
		}
	}
	return out
}

func reversed(b []byte) []byte {
	out := make([]byte, len(b))
	for i, x := range b {
		out[len(b)-1-i] = x
	}
	return out
}
