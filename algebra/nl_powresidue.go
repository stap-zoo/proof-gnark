package algebra

import (
	"fmt"
	"math/big"
	"sync"

	"github.com/consensys/gnark/constraint/solver"
	"github.com/consensys/gnark/frontend"
	"github.com/consensys/gnark/std/lookup/logderivlookup"
)

func init() {
	solver.RegisterHint(powerResidueHint)
}

// GrMode selects how the PowerResidueSbox obtains g^r in-circuit (the two ways
// differ only in cost, never in the result):
//
//   - GrLookup: a second lookup table G[j] = g^j (with G[m] = 0). One extra
//     lookup query per S-box, but the table adds ~m once-per-circuit
//     constraints under gnark's logderiv emulation — cheap for small m.
//   - GrBits: decompose r into log2(m) bits and recompose g^r as the product
//     prod_i (1 + b_i*(g^(2^i)-1)). ~2*log2(m) constraints per S-box, no extra
//     table — cheaper than GrLookup once m is large relative to the number of
//     S-box evaluations (few wide tables vs. many small products).
type GrMode int

const (
	GrLookup GrMode = iota
	GrBits
)

// PowerResidueSbox evaluates Polocolo's power-residue S-box (eprint 2025/926)
//
//	S(x) = x^{-1} * T[x^((p-1)/m)],   S(0) = 0,
//
// where m | p-1 is a power of two, g generates F_p^*, and sigma permutes
// {0,...,m-1}. Writing x = g^(qm+r), S(x) = g^(-qm + rm + sigma(r)).
//
// The "annihilator" power residue x^((p-1)/m) is NEVER evaluated in-circuit
// (a ~254-bit exponent): it exists only to select the residue class r, so the
// circuit witnesses r directly and verifies it with the paper's own Plonk
// constraint system (Section 3.3, Constraints (2)):
//
//	wm = w^m                    log2(m) squarings, w from a hint
//	wm * G_j == x               binds j to x's residue class: g^e is an m-th
//	                            power iff m | e, so for j < m the equation is
//	                            satisfiable iff j == r (Lemma 3)
//	winv = wm^{-1}              forces wm != 0 (the zero case cannot cheat)
//	y = winv * V_j              V_r = g^(rm + sigma(r)), winv = g^(-qm)
//
// with (G_j, V_j) from lookup tables keyed by the witnessed index j. Row m of
// both tables is 0, which makes S(0) = 0 hold with NO special-case constraints:
// at x = 0 the membership constraint forces G_j = 0 (since wm != 0), hence
// j = m, V_j = 0, y = 0; at x != 0 the zero row is unusable.
//
// Construct once per circuit (the tables are shared by every Apply call), like
// the Skyscraper Bar engine.
type PowerResidueSbox struct {
	api   frontend.API
	m     int
	l     int // log2(m)
	g     *big.Int
	mode  GrMode
	tblV  logderivlookup.Table
	tblG  logderivlookup.Table // GrLookup only
	gPow2 []*big.Int           // GrBits only: g^(2^i) - 1 for i = 0..l-1
}

// NewPowerResidueSbox builds the shared lookup tables for one circuit. m must
// be a power of two dividing p-1 (p the compile-time field modulus), g a
// generator of F_p^*, sigma a permutation of {0,...,m-1}.
func NewPowerResidueSbox(api frontend.API, m int, g *big.Int, sigma []int, mode GrMode) *PowerResidueSbox {
	if m < 2 || m&(m-1) != 0 {
		panic(fmt.Sprintf("algebra: PowerResidueSbox m must be a power of two >= 2, got %d", m))
	}
	if len(sigma) != m {
		panic(fmt.Sprintf("algebra: PowerResidueSbox sigma has length %d, want m=%d", len(sigma), m))
	}
	p := api.Compiler().Field()
	pMinus1 := new(big.Int).Sub(p, big.NewInt(1))
	if new(big.Int).Mod(pMinus1, big.NewInt(int64(m))).Sign() != 0 {
		panic(fmt.Sprintf("algebra: PowerResidueSbox m=%d does not divide p-1", m))
	}

	s := &PowerResidueSbox{api: api, m: m, g: new(big.Int).Set(g), mode: mode}
	for mm := m; mm > 1; mm >>= 1 {
		s.l++
	}

	// V[r] = g^(r*m + sigma(r)) for r < m; V[m] = 0 (the S(0) = 0 row).
	s.tblV = logderivlookup.New(api)
	e := new(big.Int)
	for r := 0; r < m; r++ {
		e.SetInt64(int64(r)*int64(m) + int64(sigma[r]))
		s.tblV.Insert(new(big.Int).Exp(g, e, p))
	}
	s.tblV.Insert(0)

	switch mode {
	case GrLookup:
		// G[r] = g^r for r < m; G[m] = 0.
		s.tblG = logderivlookup.New(api)
		for r := 0; r < m; r++ {
			s.tblG.Insert(new(big.Int).Exp(g, e.SetInt64(int64(r)), p))
		}
		s.tblG.Insert(0)
	case GrBits:
		// g^(2^i) - 1 constants for the bit-recomposition product.
		s.gPow2 = make([]*big.Int, s.l)
		acc := new(big.Int).Set(g)
		for i := 0; i < s.l; i++ {
			s.gPow2[i] = new(big.Int).Sub(acc, big.NewInt(1))
			acc.Mul(acc, acc).Mod(acc, p)
		}
	default:
		panic(fmt.Sprintf("algebra: unknown GrMode %d", mode))
	}
	return s
}

// Apply returns y = S(x). Cost per call: log2(m) squarings + 3 multiplications
// + 1 lookup into V, plus either 1 lookup into G (GrLookup) or ~log2(m)
// constraints for the bit recomposition (GrBits — one PLONK row per bit, see
// MulAffine, where the obvious spelling costs two).
func (s *PowerResidueSbox) Apply(x frontend.Variable) frontend.Variable {
	api := s.api

	// Hint: the residue-class index r (0 for x = 0), the zero flag z, and an
	// m-th root w of x*g^(-r) (1 for x = 0). All three are pinned below.
	outs, err := api.NewHint(powerResidueHint, 3, x, s.m, s.g)
	if err != nil {
		panic(fmt.Sprintf("algebra: PowerResidueSbox hint: %v", err))
	}
	r, z, w := outs[0], outs[1], outs[2]
	j := api.Add(r, api.Mul(z, s.m)) // table row: r for x != 0, m (the zero row) for x = 0

	// wm = w^m via log2(m) squarings.
	wm := w
	for i := 0; i < s.l; i++ {
		wm = api.Mul(wm, wm)
	}

	// G_j = g^r (0 on the zero row).
	var G frontend.Variable
	switch s.mode {
	case GrLookup:
		G = s.tblG.Lookup(j)[0]
	default: // GrBits
		api.AssertIsBoolean(z)
		bits := api.ToBinary(r, s.l) // also binds r to [0, m)
		G = frontend.Variable(1)
		for i, b := range bits {
			// G * (1 + b*(g^(2^i)-1)): a product plus a linear term in the same
			// wires, so ONE PLONK gate rather than the two that materialising the
			// bracket first costs (algebra.MulAffine). That halves this loop, which
			// is l = log2(m) steps long and runs once per S-box evaluation.
			G = MulAffine(api, G, api.Mul(b, s.gPow2[i]), bigOne)
		}
		G = MulAffine(api, G, api.Neg(z), bigOne) // G * (1 - z), likewise one gate
	}

	// Membership: wm * g^r == x pins j to x's residue class (and j = m iff x = 0,
	// because the next line forces wm != 0). The product is not materialised — it
	// shares a row with the equality (AssertProductIsEqual).
	AssertProductIsEqual(api, wm, G, x)
	winv := api.Inverse(wm) // wm = g^(qm) != 0 always; also what makes the zero row sound

	// y = g^(-qm) * g^(rm + sigma(r)) = S(x); the V lookup binds j to [0, m].
	return api.Mul(winv, s.tblV.Lookup(j)[0])
}

// powerResidueHint witnesses (r, z, w) for the S-box input x: inputs are
// [x, m, g]. For x = 0 it returns (0, 1, 1); otherwise (r, 0, w) with
// x = g^(qm+r), 0 <= r < m, and w ANY m-th root of x*g^(-r) (all m roots
// satisfy w^m = x*g^(-r); the circuit only ever uses w^m).
//
// This runs outside the constraint system, but it is the most expensive native
// computation of any construction here and it is not amortized by anything: it
// is called once per S-box evaluation — t*R times per permutation, ~15k times
// for a 1024-leaf Merkle tree — and every proof pays it twice, once to evaluate
// the circuit natively for the witness and again inside gnark's solver. Both
// halves are therefore table-driven; see [residueTables] for the precomputation
// and [residueTables.root] for the m-th root, which is the part that matters.
func powerResidueHint(field *big.Int, inputs []*big.Int, outputs []*big.Int) error {
	if len(inputs) != 3 || len(outputs) != 3 {
		return fmt.Errorf("algebra: powerResidueHint expects 3 inputs and 3 outputs, got %d and %d", len(inputs), len(outputs))
	}
	x, g := inputs[0], inputs[2]
	m := inputs[1].Int64()

	if x.Sign() == 0 {
		outputs[0].SetInt64(0)
		outputs[1].SetInt64(1)
		outputs[2].SetInt64(1)
		return nil
	}

	t, err := residueTablesFor(field, m, g)
	if err != nil {
		return err
	}

	// r: which of the m classes x^ann = g^(r*ann) lands in. One table lookup,
	// not the linear sweep over the m class representatives the definition
	// suggests.
	r, ok := t.class[string(new(big.Int).Exp(x, t.ann, field).Bytes())]
	if !ok {
		return fmt.Errorf("algebra: powerResidueHint: no residue class found (g not a generator, or m does not divide p-1)")
	}

	// w: an m-th root of u = x*g^(-r), which exists because u is an m-th power.
	u := new(big.Int).Mul(x, t.gPowInv[r])
	w := t.root(u.Mod(u, field))

	// Defensive self-check: w^m * g^r == x.
	chk := new(big.Int).Exp(w, big.NewInt(t.m), field)
	chk.Mul(chk, t.gPow[r]).Mod(chk, field)
	if chk.Cmp(new(big.Int).Mod(x, field)) != 0 {
		return fmt.Errorf("algebra: powerResidueHint: self-check failed")
	}

	outputs[0].SetInt64(r)
	outputs[1].SetInt64(0)
	outputs[2].Set(w)
	return nil
}

// residueTables is everything powerResidueHint can precompute for one (p, m, g),
// which is everything except the two exponentiations that depend on x.
//
// Writing p-1 = 2^s*q with q odd (s = 28 for BN254, 32 for BLS12-381) and
// m = 2^l, the m-th root of an m-th power u is taken by splitting the exponent
// with the CRT rather than by l iterated square roots — see [residueTables.root].
// The naive route costs l full Tonelli-Shanks runs, each O(s^2) multiplications
// in the 2-Sylow; the split costs one such sweep for all l levels together, and
// at Polocolo's m = 1024 that is where a factor ~12 of the hint sits.
type residueTables struct {
	p   *big.Int
	ann *big.Int // (p-1)/m, the "annihilator" exponent
	m   int64
	l   int // log2(m)
	s   int // 2-adicity of p-1

	class   map[string]int64 // g^(r*ann), big-endian bytes -> r
	gPow    []*big.Int       // g^r,    r = 0..m-1
	gPowInv []*big.Int       // g^(-r), r = 0..m-1

	// The m-th root's fixed parts. eOdd and eTwo are the CRT idempotents scaled
	// by 2^-l: raising to eOdd keeps the odd part of the exponent and divides it
	// by 2^l, raising to eTwo keeps the 2-part untouched.
	eOdd    *big.Int   // == 2^-l (mod q), == 0 (mod 2^s)
	eTwo    *big.Int   // == 1    (mod 2^s), == 0 (mod q)
	z       *big.Int   // g^q, a generator of the 2-Sylow (order 2^s)
	zNegPow []*big.Int // z^(-2^i), i = 0..s-1
	negOne  *big.Int   // z^(2^(s-1)) = p-1, the only element of order 2
}

// root returns some w with w^m == u, for u an m-th power (which x*g^(-r) always
// is). It never verifies that premise: the caller's self-check does.
//
// F_p^* is cyclic of order 2^s*q, so u = g^a and we want g^(a/2^l). Splitting
// that by the CRT into its part mod q and its part mod 2^s makes each half easy
// for a different reason:
//
//	odd part   raising to the 2^l-th power is invertible mod q, so u^eOdd is
//	           one exponentiation
//	2-part     it is not invertible there (squaring is 2-to-1 on a 2-group, which
//	           is exactly why no single exponent can extract this root), so
//	           u^eTwo is written as z^j by a Pohlig-Hellman discrete log — s-l
//	           steps of shrinking squaring chains — and z^(j/2^l) taken instead.
//	           j's low l bits are zero precisely because u is an m-th power.
//
// The product of the two halves has the right exponent modulo q and modulo 2^s,
// hence modulo p-1. Which of the m roots comes out is unspecified and does not
// matter: the circuit only ever squares w back up to w^m.
func (t *residueTables) root(u *big.Int) *big.Int {
	w := new(big.Int).Exp(u, t.eOdd, t.p)

	cur := new(big.Int).Exp(u, t.eTwo, t.p)
	j, e := new(big.Int), new(big.Int)
	for i := t.l; i < t.s; i++ {
		// cur has had bits < i cleared, so cur^(2^(s-1-i)) is +-1 and is -1
		// exactly when bit i of j is set.
		e.Set(cur)
		for k := i; k < t.s-1; k++ {
			e.Mul(e, e).Mod(e, t.p)
		}
		if e.Cmp(t.negOne) == 0 {
			j.SetBit(j, i, 1)
			cur.Mul(cur, t.zNegPow[i]).Mod(cur, t.p)
		}
	}
	w.Mul(w, new(big.Int).Exp(t.z, j.Rsh(j, uint(t.l)), t.p))
	return w.Mod(w, t.p)
}

// residueKey identifies one table set. A hint sees no gadget state — only the
// values passed to api.NewHint — so the field, m and g are all it can key on,
// and all three are what the tables depend on.
type residueKey struct {
	p, g string
	m    int64
}

var (
	residueMu    sync.RWMutex
	residueCache = map[residueKey]*residueTables{}
)

// residueTablesFor returns the cached tables for (p, m, g), building them on
// first use. The cache is process-wide because the hint is: gnark resolves hints
// by a global registry, and the same instance is measured under several
// backends. Two goroutines racing to build the same entry is harmless — the
// tables are a pure function of the key, so either copy will do.
func residueTablesFor(p *big.Int, m int64, g *big.Int) (*residueTables, error) {
	k := residueKey{p: p.String(), g: g.String(), m: m}

	residueMu.RLock()
	t, ok := residueCache[k]
	residueMu.RUnlock()
	if ok {
		return t, nil
	}

	t, err := newResidueTables(p, m, g)
	if err != nil {
		return nil, err
	}
	residueMu.Lock()
	residueCache[k] = t
	residueMu.Unlock()
	return t, nil
}

func newResidueTables(p *big.Int, m int64, g *big.Int) (*residueTables, error) {
	if m < 2 || m&(m-1) != 0 {
		return nil, fmt.Errorf("algebra: powerResidueHint: m must be a power of two >= 2, got %d", m)
	}
	pMinus1 := new(big.Int).Sub(p, big.NewInt(1))
	if new(big.Int).Mod(pMinus1, big.NewInt(m)).Sign() != 0 {
		return nil, fmt.Errorf("algebra: powerResidueHint: m=%d does not divide p-1", m)
	}

	t := &residueTables{p: new(big.Int).Set(p), m: m}
	t.ann = new(big.Int).Div(pMinus1, big.NewInt(m))
	for mm := m; mm > 1; mm >>= 1 {
		t.l++
	}
	q := new(big.Int).Set(pMinus1)
	for q.Bit(0) == 0 {
		q.Rsh(q, 1)
		t.s++
	}
	if t.l > t.s {
		return nil, fmt.Errorf("algebra: powerResidueHint: m=%d exceeds the 2-adicity 2^%d of p-1", m, t.s)
	}

	// The class representatives g^(r*ann) and the g^(+-r) tables, in one sweep.
	step := new(big.Int).Exp(g, t.ann, p)
	t.class = make(map[string]int64, m)
	t.gPow = make([]*big.Int, m)
	t.gPowInv = make([]*big.Int, m)
	cur, gr := big.NewInt(1), big.NewInt(1)
	for r := int64(0); r < m; r++ {
		t.class[string(cur.Bytes())] = r
		cur = new(big.Int).Mod(new(big.Int).Mul(cur, step), p)
		t.gPow[r] = new(big.Int).Set(gr)
		t.gPowInv[r] = new(big.Int).ModInverse(gr, p)
		gr = new(big.Int).Mod(new(big.Int).Mul(gr, g), p)
	}
	if len(t.class) != int(m) {
		return nil, fmt.Errorf("algebra: powerResidueHint: g^((p-1)/%d) has order < %d (g is not a generator)", m, m)
	}

	// eOdd = 2^s * ((2^-l * 2^-s) mod q): == 2^-l mod q and == 0 mod 2^s.
	// eTwo = q * (q^-1 mod 2^s):          == 1    mod 2^s and == 0 mod q.
	twoL := new(big.Int).Lsh(big.NewInt(1), uint(t.l))
	twoS := new(big.Int).Lsh(big.NewInt(1), uint(t.s))
	t.eOdd = new(big.Int).Mul(new(big.Int).ModInverse(twoL, q), new(big.Int).ModInverse(twoS, q))
	t.eOdd.Mod(t.eOdd, q).Mul(t.eOdd, twoS)
	t.eTwo = new(big.Int).Mul(q, new(big.Int).ModInverse(q, twoS))

	t.z = new(big.Int).Exp(g, q, p)
	t.zNegPow = make([]*big.Int, t.s)
	acc := new(big.Int).Set(t.z)
	two := big.NewInt(2)
	for i := 0; i < t.s; i++ {
		t.zNegPow[i] = new(big.Int).ModInverse(acc, p)
		t.negOne = new(big.Int).Set(acc)
		acc = new(big.Int).Exp(acc, two, p)
	}
	if t.negOne.Cmp(pMinus1) != 0 {
		return nil, fmt.Errorf("algebra: powerResidueHint: g^((p-1)/2^%d) is not -1 (g does not generate the 2-Sylow)", t.s)
	}
	return t, nil
}
