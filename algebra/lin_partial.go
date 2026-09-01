// This file is the parameter-time factorisation of a *partial-round block* — a run
// of rounds that share one matrix and whose S-box touches only the first few
// branches — plus the two-line applier that runs the result in circuit. It is the
// counterpart of lin_slp.go for a family of layers that cannot be cheapened one
// layer at a time: what is exploitable here is not the matrix but the block, so
// the object a construction stores is the whole rewritten run of rounds rather
// than a program for one layer.
//
// Everything here runs natively over *big.Int before the circuit exists, so none
// of it costs a constraint or a gate; it only decides how many the rounds will
// cost. Nothing changes in R1CS either — a constant linear layer folds into wire
// coefficients whatever its shape — so this moves the PLONK column alone.

package algebra

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark/frontend"
)

// PartialBlock is the block to factorise: len(Consts) rounds of "S-box on the
// first U branches, then y = M*x + Consts[k]", preceded by the layer
// y = Pre*x + PreConsts that feeds it. Every field is parameter data a
// construction has already derived.
//
// The S-box itself is absent, and does not need to be there: the rewrite holds
// for any map that touches only the first U branches — whatever it is, and even
// if it differs from round to round.
type PartialBlock struct {
	M         [][]*big.Int // the matrix every round in the block applies (t x t)
	U         int          // branches the S-box touches: the first U of t
	Consts    [][]*big.Int // one row per round: the constants that round's layer adds
	Pre       [][]*big.Int // the matrix of the layer that feeds the block (t x t)
	PreConsts []*big.Int   // and the constants that layer adds (nil for none)
	Field     *big.Int     // field characteristic p
}

// PartialRounds is the sparse factorisation of a [PartialBlock]: the same block
// of rounds, rewritten so each round's matrix is sparse. It is the Poseidon
// paper's Appendix B rewrite (circomlib's "sparse" matrices S with its
// pre-multiplied P), derived from the matrix here rather than transcribed, and it
// applies to any Hades-style design — Poseidon's internal rounds, Neptune's,
// Poseidon2's.
//
// # The rewrite
//
// The block computes y_k = M*S(y_{k-1}) + d_k for k = 1..R, with S the partial
// S-box. A matrix of the form N = diag(I_u, C) fixes every branch S touches, so
// N*S = S*N holds *exactly* — the two act on disjoint coordinates. Peeling one
// such factor out of each round, from the back,
//
//	N_R = I,   W_k = N_k*M,   N_{k-1} = diag(I_u, (W_k)_{u:,u:}),   G_k = W_k*N_{k-1}^{-1}
//
// rewrites the block over the changed coordinates v_k = N_k*y_k:
//
//	v_0 = N_0*y_0,   v_k = G_k*S(v_{k-1}) + N_k*d_k,   y_R = v_R
//
// by induction, since G_k*S(v_{k-1}) = G_k*N_{k-1}*S(y_{k-1}) = N_k*M*S(y_{k-1}).
// It is the same shape — one S-box and one affine layer per round — so nothing
// but the constants and the matrices change. The last coordinate change is the
// identity, so the block hands the rest of the permutation the true state and the
// rewrite stays inside the block.
//
// # What it costs
//
// G_k's bottom-right (t-u)x(t-u) block is the identity by construction, leaving u
// dense rows and t-u rows of u+1 terms:
//
//	u(t-1) + u(t-u) gates  =  2t-2 at u = 1
//
// against t(t-1) for a dense M, or 2t-1 for the shared-sum program of a J+diag M
// ([JPlusDiagSLP]). No program is needed to realize it: the zeros are in the
// matrix, and [MatVecMulConst] already skips them.
//
// The one thing the rewrite has to pay for is the residual N_0 at the block
// entry. It fixes the S-box's branches, so it cannot commute past the *full*
// S-box of the round before the block — but it can be multiplied into that
// round's matrix at parameter time, which is what [FactorPartialRounds] returns
// as Pre/PreConsts and what circomlib's P matrix is. That costs
// DenseGates(N_0*Pre) - DenseGates(Pre): nothing at all where the feeding layer
// is already dense (Poseidon, whose MDS matrix has no zero entry), and (t-1)t/2
// where it is the even/odd split (Neptune's M_ext, whose rows have t/2 nonzeros).
// The alternative placement — applying N_0 = diag(I_u, C) as its own layer, u
// free rows and t-u dense ones, (t-u)(t-u-1) gates — is not implemented because
// it never wins for the two constructions here: it ties Neptune at t=4 and loses
// everywhere else.
type PartialRounds struct {
	// Pre and PreConsts replace the feeding layer's own matrix and constants for
	// that one round: N_0*PartialBlock.Pre and N_0*PartialBlock.PreConsts.
	Pre       [][]*big.Int
	PreConsts []*big.Int

	// M[k] and Consts[k] are the layer of block round k, in round order:
	// G_{k+1} and N_{k+1}*d_{k+1} above.
	M      [][][]*big.Int
	Consts [][]*big.Int

	// n holds N_0..N_R, the coordinate change at each round boundary. Kept because
	// it is what [CheckPartialRounds] verifies the rest against.
	n [][][]*big.Int
}

// Rounds reports how many rounds the block has.
func (f *PartialRounds) Rounds() int { return len(f.M) }

// Gates is the modelled PLONK cost of the rewritten block: the feeding layer with
// the residual folded in, plus one layer per round. Compare it against
// DenseGates(Pre) + R * (what one round's layer costs today) — which is what
// [SelectPartialRounds] does.
func (f *PartialRounds) Gates() int {
	n := DenseGates(f.Pre)
	for _, m := range f.M {
		n += DenseGates(m)
	}
	return n
}

// Apply applies block round k's linear layer: the rewritten matrix, carrying the
// constants that follow it. Round constants stay free — after the rewrite every
// branch still has a gate to ride on (see MatVecMulConst).
func (f *PartialRounds) Apply(api frontend.API, k int, state []frontend.Variable) []frontend.Variable {
	return MatVecMulConst(api, f.M[k], state, f.Consts[k])
}

// ApplyPre applies the layer that feeds the block, in place of the matrix and
// constants the construction would otherwise apply in that round.
func (f *PartialRounds) ApplyPre(api frontend.API, state []frontend.Variable) []frontend.Variable {
	return MatVecMulConst(api, f.Pre, state, f.PreConsts)
}

// SelectPartialRounds is what a parameter constructor calls for a partial-round
// block: it factorises the block and keeps the result only when it is strictly
// cheaper in PLONK than what the construction does today — perRound gates for
// each of the block's rounds, plus the feeding layer as it stands.
//
// It returns nil with no error when the rewrite does not apply (a singular pivot
// block, so no factorisation exists) or does not pay: at t = 2 the rewritten
// layer costs 2t-2 = 2, exactly what the dense 2x2 product costs, so there is
// nothing to win and the residual can only lose. An error means the factorisation
// fails its own identities — always a bug, never a cost question — so no
// construction can ship a rewrite that was never checked against the block it
// replaces.
func SelectPartialRounds(b PartialBlock, perRound int) (*PartialRounds, error) {
	f, err := FactorPartialRounds(b)
	if err != nil || f == nil {
		return nil, err
	}
	if f.Gates() >= perRound*len(b.Consts)+DenseGates(b.Pre) {
		return nil, nil
	}
	return f, nil
}

// FactorPartialRounds derives the factorisation described on [PartialRounds] and
// verifies it with [CheckPartialRounds] before returning it. A nil result with a
// nil error means no factorisation exists: some (N_k*M)_{u:,u:} is singular, and
// the caller keeps applying M directly.
func FactorPartialRounds(b PartialBlock) (*PartialRounds, error) {
	t := len(b.M)
	p := b.Field
	if p == nil || p.Sign() <= 0 {
		return nil, fmt.Errorf("field modulus must be a positive integer")
	}
	if t < 2 {
		return nil, fmt.Errorf("matrix must be at least 2x2, got %d rows", t)
	}
	if b.U < 1 || b.U >= t {
		return nil, fmt.Errorf("S-box width u must be in [1, %d), got %d", t, b.U)
	}
	if err := checkSquare("matrix", b.M, t); err != nil {
		return nil, err
	}
	if err := checkSquare("feeding matrix", b.Pre, t); err != nil {
		return nil, err
	}
	r := len(b.Consts)
	if r == 0 {
		return nil, fmt.Errorf("block has no rounds")
	}
	for k, c := range b.Consts {
		if len(c) != t {
			return nil, fmt.Errorf("round %d has %d constants, want %d", k, len(c), t)
		}
	}
	if b.PreConsts != nil && len(b.PreConsts) != t {
		return nil, fmt.Errorf("feeding layer has %d constants, want %d", len(b.PreConsts), t)
	}

	f := &PartialRounds{
		M:      make([][][]*big.Int, r),
		Consts: make([][]*big.Int, r),
		n:      make([][][]*big.Int, r+1),
	}
	f.n[r] = identityMod(t)
	for k := r - 1; k >= 0; k-- {
		w := matMulMod(f.n[k+1], b.M, p)
		c := lowerBlock(w, b.U)
		cinv, ok := invertMod(c, p)
		if !ok {
			return nil, nil // singular pivot block: no factorisation
		}
		f.n[k] = blockDiagMod(b.U, c)
		f.M[k] = matMulMod(w, blockDiagMod(b.U, cinv), p)
		f.Consts[k] = matVecMod(f.n[k+1], b.Consts[k], p)
	}
	f.Pre = matMulMod(f.n[0], b.Pre, p)
	f.PreConsts = matVecMod(f.n[0], b.PreConsts, p)

	if err := CheckPartialRounds(f, b); err != nil {
		return nil, err
	}
	return f, nil
}

// CheckPartialRounds reports an error unless f rewrites exactly the block b. Four
// identities are enough, and together they are the derivation on [PartialRounds]
// read backwards:
//
//   - every coordinate change fixes the S-box's branches, N = diag(I_u, C), which
//     is exactly the condition for N*S = S*N to hold for *every* map S on those
//     branches — so the commutation the rewrite turns on is real;
//   - G_k*N_{k-1} = N_k*M, so a round in the changed coordinates computes what the
//     original round computed;
//   - the last change is the identity, so the block returns the true state;
//   - the feeding layer and every round's constants are their originals under the
//     change that applies there.
func CheckPartialRounds(f *PartialRounds, b PartialBlock) error {
	p := b.Field
	r := len(b.Consts)
	if len(f.M) != r || len(f.Consts) != r || len(f.n) != r+1 {
		return fmt.Errorf("factorisation has %d matrices and %d constant rows for %d rounds", len(f.M), len(f.Consts), r)
	}
	t := len(b.M)
	if !matEqualMod(f.n[r], identityMod(t), p) {
		return fmt.Errorf("the last coordinate change is not the identity")
	}
	for k := 0; k <= r; k++ {
		if err := fixesBranches(f.n[k], b.U); err != nil {
			return fmt.Errorf("coordinate change %d does not commute with the S-box: %w", k, err)
		}
	}
	for k := 0; k < r; k++ {
		if !matEqualMod(matMulMod(f.M[k], f.n[k], p), matMulMod(f.n[k+1], b.M, p), p) {
			return fmt.Errorf("round %d: G*N != N'*M", k)
		}
		if !vecEqualMod(f.Consts[k], matVecMod(f.n[k+1], b.Consts[k], p), p) {
			return fmt.Errorf("round %d: constants are not the originals under the coordinate change", k)
		}
	}
	if !matEqualMod(f.Pre, matMulMod(f.n[0], b.Pre, p), p) {
		return fmt.Errorf("the feeding layer's matrix does not carry the residual")
	}
	if !vecEqualMod(f.PreConsts, matVecMod(f.n[0], b.PreConsts, p), p) {
		return fmt.Errorf("the feeding layer's constants do not carry the residual")
	}
	return nil
}

// fixesBranches reports an error unless n is diag(I_u, C): the identity on the
// first u coordinates and blind to them everywhere else.
func fixesBranches(n [][]*big.Int, u int) error {
	for i, row := range n {
		for j, v := range row {
			if i >= u && j >= u {
				continue
			}
			want := 0
			if i == j {
				want = 1
			}
			if v.Cmp(big.NewInt(int64(want))) != 0 {
				return fmt.Errorf("entry [%d][%d] is %s, want %d", i, j, v, want)
			}
		}
	}
	return nil
}

// ---------------------------------------------------------------------------
// Native matrix arithmetic mod p. Small and exact: the largest inverse taken
// here is (t-1)x(t-1) and it happens once per round at parameter time.
// ---------------------------------------------------------------------------

func checkSquare(what string, m [][]*big.Int, t int) error {
	if len(m) != t {
		return fmt.Errorf("%s has %d rows, want %d", what, len(m), t)
	}
	for i, row := range m {
		if len(row) != t {
			return fmt.Errorf("%s row %d has width %d, want %d", what, i, len(row), t)
		}
	}
	return nil
}

func identityMod(t int) [][]*big.Int {
	m := make([][]*big.Int, t)
	for i := range m {
		m[i] = make([]*big.Int, t)
		for j := range m[i] {
			if i == j {
				m[i][j] = big.NewInt(1)
			} else {
				m[i][j] = big.NewInt(0)
			}
		}
	}
	return m
}

// lowerBlock returns the trailing (t-u)x(t-u) block m[u:][u:], copied.
func lowerBlock(m [][]*big.Int, u int) [][]*big.Int {
	out := make([][]*big.Int, len(m)-u)
	for i := range out {
		out[i] = make([]*big.Int, len(m)-u)
		for j := range out[i] {
			out[i][j] = new(big.Int).Set(m[u+i][u+j])
		}
	}
	return out
}

// blockDiagMod returns diag(I_u, c).
func blockDiagMod(u int, c [][]*big.Int) [][]*big.Int {
	m := identityMod(u + len(c))
	for i := range c {
		for j := range c[i] {
			m[u+i][u+j] = new(big.Int).Set(c[i][j])
		}
	}
	return m
}

func matMulMod(a, b [][]*big.Int, p *big.Int) [][]*big.Int {
	n, k := len(a), len(b)
	cols := len(b[0])
	out := make([][]*big.Int, n)
	for i := 0; i < n; i++ {
		out[i] = make([]*big.Int, cols)
		for j := 0; j < cols; j++ {
			acc := new(big.Int)
			for l := 0; l < k; l++ {
				acc.Add(acc, new(big.Int).Mul(a[i][l], b[l][j]))
			}
			out[i][j] = acc.Mod(acc, p)
		}
	}
	return out
}

// matVecMod returns m*v, and nil for a nil v (no constants stay no constants).
func matVecMod(m [][]*big.Int, v []*big.Int, p *big.Int) []*big.Int {
	if v == nil {
		return nil
	}
	out := make([]*big.Int, len(m))
	for i, row := range m {
		acc := new(big.Int)
		for j, c := range row {
			acc.Add(acc, new(big.Int).Mul(c, v[j]))
		}
		out[i] = acc.Mod(acc, p)
	}
	return out
}

// invertMod inverts m over the prime field p by Gauss-Jordan elimination on
// [m | I], reporting ok == false if m is singular.
func invertMod(m [][]*big.Int, p *big.Int) ([][]*big.Int, bool) {
	n := len(m)
	aug := make([][]*big.Int, n)
	for i := range aug {
		aug[i] = make([]*big.Int, 2*n)
		for j := 0; j < n; j++ {
			aug[i][j] = new(big.Int).Mod(m[i][j], p)
		}
		for j := 0; j < n; j++ {
			aug[i][n+j] = big.NewInt(0)
		}
		aug[i][n+i] = big.NewInt(1)
	}

	for col := 0; col < n; col++ {
		pivot := -1
		for row := col; row < n; row++ {
			if aug[row][col].Sign() != 0 {
				pivot = row
				break
			}
		}
		if pivot < 0 {
			return nil, false
		}
		aug[col], aug[pivot] = aug[pivot], aug[col]

		inv := new(big.Int).ModInverse(aug[col][col], p)
		if inv == nil {
			return nil, false // p not prime: no field to work in
		}
		for j := col; j < 2*n; j++ {
			aug[col][j].Mod(new(big.Int).Mul(aug[col][j], inv), p)
		}
		for row := 0; row < n; row++ {
			if row == col || aug[row][col].Sign() == 0 {
				continue
			}
			factor := new(big.Int).Set(aug[row][col])
			for j := col; j < 2*n; j++ {
				term := new(big.Int).Mul(factor, aug[col][j])
				aug[row][j].Mod(new(big.Int).Sub(aug[row][j], term), p)
			}
		}
	}

	out := make([][]*big.Int, n)
	for i := range out {
		out[i] = aug[i][n:]
	}
	return out, true
}

func matEqualMod(a, b [][]*big.Int, p *big.Int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !vecEqualMod(a[i], b[i], p) {
			return false
		}
	}
	return true
}

func vecEqualMod(a, b []*big.Int, p *big.Int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if new(big.Int).Mod(a[i], p).Cmp(new(big.Int).Mod(b[i], p)) != 0 {
			return false
		}
	}
	return true
}
