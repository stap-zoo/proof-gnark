// This file holds the parameter-time half of the straight-line program: given the
// matrix a construction's linear layer computes, derive an SLP for it. Everything
// here runs natively over *big.Int before the circuit exists, so none of it costs a
// constraint or a gate — it only decides how many the layer will cost. The SLP type
// itself, and what applying one costs, are in lin_slp.go.
//
// Two kinds live here. JPlusDiagSLP and SharedFormSLP read the program *out of* the
// matrix and report ok == false for a matrix of the wrong shape, so they cannot be
// applied to a layer they do not compute. DLM44SLP, BlockCirculantM4SLP and
// DLM3352SLP build a fixed published program for a named matrix family instead.
// Either way the caller passes the result through SelectSLP, which checks it
// entry-by-entry against the matrix — see lin_slp.go.

package algebra

import (
	"fmt"
	"math/big"
)

// JPlusDiagSLP derives the shared-sum program for a matrix of the form
// J + diag(d) — every off-diagonal entry equal to 1 — and reports ok == false
// for any other matrix, so it cannot be applied to a layer it does not compute.
// The program forms the running sum s = sum_j x_j once (t-1 gates) and reuses it
// for every output: out_i = s + (m[i][i] - 1)*x_i, free where m[i][i] == 1. Cost
// is t-1 + |{i : m[i][i] != 1}| gates, at most 2t-1, against t(t-1) dense — the
// dominant saving in Poseidon2's and Neptune's internal rounds, and it also
// covers circ(2,1,...,1) = J + I (Griffin's matrix, Poseidon2's t=3 external
// matrix).
//
// Deriving the program from the matrix rather than transcribing it is what makes
// this safe at any t: the diagonal it uses is read straight out of the layer it
// replaces. It is still checked with Matrix(p) by the caller.
func JPlusDiagSLP(m [][]*big.Int, p *big.Int) (SLP, bool) {
	t := len(m)
	if t < 2 || p == nil {
		return SLP{}, false
	}
	one := big.NewInt(1)
	for i, row := range m {
		if len(row) != t {
			return SLP{}, false
		}
		for j, v := range row {
			if i != j && new(big.Int).Mod(v, p).Cmp(one) != 0 {
				return SLP{}, false
			}
		}
	}

	b := &slpBuilder{t: t}
	sum := 0
	for j := 1; j < t; j++ {
		sum = b.add(one, sum, one, j)
	}
	out := make([]int, t)
	for i := 0; i < t; i++ {
		d := new(big.Int).Sub(m[i][i], one)
		d.Mod(d, p)
		if d.Sign() == 0 {
			out[i] = sum // row i is the bare sum
			continue
		}
		out[i] = b.add(one, sum, d, i)
	}
	return b.slp(out), true
}

// SharedFormSLP derives the shared-linear-form program for a matrix whose every
// row is a scalar multiple of ONE common linear form, off at most a single
// coordinate:
//
//	m[i] = c_i * v + e_i * delta_{k_i}
//
// The program computes V = <v, x> once and then reads every output off it:
// out_i = c_i*V + e_i*x_{k_i}. Cost is nnz(v) - 1 gates for V plus one per row
// that is not V itself, i.e. at most 2t-1 against t(t-1) dense.
//
// This is the general form of the shared sum: J + diag(d) is the case v = 1
// (every c_i = 1, e_i = d_i - 1), which [JPlusDiagSLP] derives directly because
// there v is known a priori. Searching for v is what covers a layer whose shared
// form is *not* the plain sum — Arion's circ(1,2,...,t) at t=3, where
// v = (6, 2, 3) up to scale gives
//
//	V  = 6*x0 + 2*x1 + 3*x2        (2 gates)
//	y0 = V - 5*x0, y1 = V/2 + x2/2, y2 = V/3 + 7*x1/3
//
// i.e. 5 gates against 6 dense, while the plain sum gives 8 (the residual
// circ(0,1,2) still has two nonzeros per row, so each output costs the same two
// gates as dense and the sum is pure overhead).
//
// Why a search over v is cheap and complete: the row a whose correction is e_a is
// c_a*v off coordinate k_a alone, so v is — up to a scale the c_i absorb — some
// row of m with at most one entry changed, and the changed entry is pinned by
// whichever other row does not correct that coordinate. Enumerating (row,
// coordinate) pairs therefore enumerates every candidate; each is validated in
// full, and the cheapest valid one wins. ok == false means no such v exists (the
// t=4 circulant is an example — with t-1 coordinates to agree on per row, four
// rows over-determine v).
//
// Like JPlusDiagSLP this derives the program *from* the matrix, and the caller
// still checks it with SelectSLP/CheckSLP, so it can never stand in for a layer
// it does not compute.
func SharedFormSLP(m [][]*big.Int, p *big.Int) (SLP, bool) {
	t := len(m)
	if t < 2 || p == nil || p.Sign() <= 0 {
		return SLP{}, false
	}
	for _, row := range m {
		if len(row) != t {
			return SLP{}, false
		}
	}
	mm := make([][]*big.Int, t)
	for i, row := range m {
		mm[i] = make([]*big.Int, t)
		for j, v := range row {
			mm[i][j] = new(big.Int).Mod(v, p)
		}
	}

	var best SLP
	found := false
	for _, v := range sharedFormCandidates(mm, p) {
		s, ok := sharedFormProgram(mm, v, p)
		if !ok {
			continue
		}
		if !found || s.Gates() < best.Gates() {
			best, found = s, true
		}
	}
	return best, found
}

// sharedFormCandidates enumerates every possible shared form v (up to scale):
// each row of m, and each row with one entry replaced by the value another row
// forces there.
func sharedFormCandidates(mm [][]*big.Int, p *big.Int) [][]*big.Int {
	t := len(mm)
	out := make([][]*big.Int, 0, t+t*t*t*t)
	for a := 0; a < t; a++ {
		out = append(out, mm[a]) // row a needs no correction: v is the row itself
		for k := 0; k < t; k++ { // row a's correction lands on coordinate k
			for b := 0; b < t; b++ {
				if b == a {
					continue
				}
				for kb := 0; kb < t; kb++ { // row b's correction lands on kb
					if kb == k {
						continue
					}
					if v, ok := pinSharedForm(mm, a, k, b, kb, p); ok {
						out = append(out, v)
					}
				}
			}
		}
	}
	return out
}

// pinSharedForm builds v = m[a] with coordinate k replaced by the value row b
// forces there: neither row corrects the coordinates outside {k, kb}, so c_b is
// read off those, and then v[k] = m[b][k]/c_b.
func pinSharedForm(mm [][]*big.Int, a, k, b, kb int, p *big.Int) ([]*big.Int, bool) {
	t := len(mm)
	var c *big.Int
	for j := 0; j < t; j++ {
		if j == k || j == kb {
			continue
		}
		if mm[a][j].Sign() == 0 {
			// v is zero here, so row b must be too — it has no correction left.
			if mm[b][j].Sign() != 0 {
				return nil, false
			}
			continue
		}
		r := divMod(mm[b][j], mm[a][j], p)
		if r == nil {
			return nil, false
		}
		if c == nil {
			c = r
		} else if c.Cmp(r) != 0 {
			return nil, false
		}
	}
	if c == nil || c.Sign() == 0 {
		return nil, false // nothing pins the scale (t = 2), or row b is degenerate
	}
	u := divMod(mm[b][k], c, p)
	if u == nil {
		return nil, false
	}
	v := make([]*big.Int, t)
	copy(v, mm[a])
	v[k] = u
	return v, true
}

// sharedFit is one row written as c*v + e*delta_k; k < 0 means no correction.
type sharedFit struct {
	c, e *big.Int
	k    int
}

// sharedFormProgram fits every row of mm to the form c*v + e*delta_k and emits
// the program, or reports ok == false if any row does not fit.
func sharedFormProgram(mm [][]*big.Int, v []*big.Int, p *big.Int) (SLP, bool) {
	t := len(mm)
	nz := make([]int, 0, t)
	for j, x := range v {
		if x.Sign() != 0 {
			nz = append(nz, j)
		}
	}
	// With fewer than two nonzeros V is a bare (scaled) input rather than a gate
	// output, which an SLP cannot name; such a matrix has no diffusion anyway.
	if len(nz) < 2 {
		return SLP{}, false
	}

	fits := make([]sharedFit, t)
	for i, row := range mm {
		f, ok := fitSharedForm(row, v, nz, p)
		if !ok {
			return SLP{}, false
		}
		fits[i] = f
	}

	one := big.NewInt(1)
	b := &slpBuilder{t: t}
	// V = <v, x>, as a chain over v's nonzero coordinates.
	cur := b.add(v[nz[0]], nz[0], v[nz[1]], nz[1])
	for _, j := range nz[2:] {
		cur = b.add(one, cur, v[j], j)
	}

	out := make([]int, t)
	for i, f := range fits {
		if f.k < 0 && f.c.Cmp(one) == 0 {
			out[i] = cur // row i is exactly V — free
			continue
		}
		k, e := f.k, f.e
		if k < 0 {
			// A pure rescaling c*V still needs a gate; the second term is zero.
			k, e = nz[0], new(big.Int)
		}
		out[i] = b.add(f.c, cur, e, k)
	}
	return b.slp(out), true
}

// fitSharedForm writes row as c*v + e*delta_k with at most one correction,
// preferring an exact multiple (k < 0). Any coordinate where row and c*v agree
// and v is nonzero pins c, so anchoring on each in turn finds c if it exists.
func fitSharedForm(row, v []*big.Int, nz []int, p *big.Int) (sharedFit, bool) {
	var best sharedFit
	found := false
	for _, anchor := range nz {
		c := divMod(row[anchor], v[anchor], p)
		if c == nil || c.Sign() == 0 {
			continue
		}
		mismatch := -1
		ok := true
		for j := range row {
			want := new(big.Int).Mul(c, v[j])
			want.Mod(want, p)
			if row[j].Cmp(want) == 0 {
				continue
			}
			if mismatch >= 0 {
				ok = false
				break
			}
			mismatch = j
		}
		if !ok {
			continue
		}
		f := sharedFit{c: c, k: mismatch, e: new(big.Int)}
		if mismatch >= 0 {
			f.e.Mul(c, v[mismatch])
			f.e.Sub(row[mismatch], f.e)
			f.e.Mod(f.e, p)
		}
		if !found || (best.k >= 0 && f.k < 0) {
			best, found = f, true
		}
	}
	return best, found
}

// divMod returns a/b mod p, or nil if b is not invertible.
func divMod(a, b, p *big.Int) *big.Int {
	inv := new(big.Int).ModInverse(b, p)
	if inv == nil {
		return nil
	}
	r := new(big.Int).Mul(a, inv)
	return r.Mod(r, p)
}

// DLM44SLP is the Duval-Leurent M^{8,4}_{4,4} program (DL18 Fig. 13, ref
// utils/matrix.py dl_m44_84_apply) computing scale * M4 * x in 8 additions,
// against 12 dense. M4 is the 4x4 diffusion block Poseidon2 and Griffin use for
// state sizes that are multiples of 4; at alpha = 2 (the only instantiation in
// use, MDS for every p > 2^31) it is [[5,7,1,3],[4,6,1,1],[1,3,5,7],[1,1,4,6]].
//
// The scale factor costs nothing: scale * M4 * x == M4 * (scale * x), and every
// step that reads an input already carries a coefficient, so the factor folds
// into those selectors. Poseidon2's t=4 external matrix is 2*M4.
func DLM44SLP(alpha, scale int64) SLP {
	b := &slpBuilder{t: 4}
	out := emitM4(b, [4]int{0, 1, 2, 3}, alpha, scale)
	return b.slp(out[:])
}

// BlockCirculantM4SLP is the program for the t x t block-circulant matrix
// circ(2,1,...,1) (x) M4 with t a multiple of 4 (ref utils/matrix.py
// m4_to_block_circulant_matrix) — Poseidon2's external matrix for t >= 4 and
// Griffin's for t a multiple of 4. Block (I,J) is 2*M4 on the block diagonal and
// M4 off it, so with x_I the I-th 4-element block the output block is
// M4 * (x_I + sum_J x_J): the block combination happens before M4, on t values,
// instead of after it.
//
// Cost by block count nb = t/4: nb = 1 is 8 gates (2*M4, see DLM44SLP); nb = 2
// takes 8 gates for the two combinations 2*x_0 + x_1 and x_0 + 2*x_1 plus 8 per
// block, i.e. 24 against 56 dense; nb >= 3 forms the block sum once (4*(nb-1)
// gates), adds x_I to it (4*nb) and runs M4 per block (8*nb), i.e. 16*nb - 4.
func BlockCirculantM4SLP(t int, alpha int64) SLP {
	if t <= 0 || t%4 != 0 {
		panic(fmt.Sprintf("algebra: BlockCirculantM4SLP needs t a positive multiple of 4, got %d", t))
	}
	nb := t / 4
	if nb == 1 {
		return DLM44SLP(alpha, 2)
	}

	one, two := big.NewInt(1), big.NewInt(2)
	b := &slpBuilder{t: t}
	// z[I] = 2*x_I + sum_{J != I} x_J, one 4-element block at a time.
	z := make([][4]int, nb)
	for k := 0; k < 4; k++ {
		if nb == 2 {
			z[0][k] = b.add(two, k, one, 4+k)
			z[1][k] = b.add(one, k, two, 4+k)
			continue
		}
		s := k
		for I := 1; I < nb; I++ {
			s = b.add(one, s, one, 4*I+k)
		}
		for I := 0; I < nb; I++ {
			z[I][k] = b.add(one, s, one, 4*I+k)
		}
	}
	out := make([]int, 0, t)
	for I := 0; I < nb; I++ {
		blk := emitM4(b, z[I], alpha, 1)
		out = append(out, blk[:]...)
	}
	return b.slp(out)
}

// DLM3352SLP is the Duval-Leurent M^{5,2}_{3,3} program (DL18 Fig. 6, ref
// utils/matrix.py dl_m33_52_apply) computing the 3x3 matrix
// [[1+a, 1, 1+a], [1, 1, a], [a, 1, 1]] in 5 additions against 6 dense —
// Anemoi's M_3 at a = g^i.
func DLM3352SLP(a, p *big.Int) SLP {
	one := big.NewInt(1)
	am := new(big.Int).Mod(a, p)
	b := &slpBuilder{t: 3}
	t0 := b.add(one, 0, am, 2)    // t  = x0 + a*x2
	u := b.add(one, 2, one, 1)    // u  = x2 + x1
	z2 := b.add(one, u, am, 0)    // z2 = u + a*x0
	y0 := b.add(one, t0, one, z2) // y0 = t + z2
	y1 := b.add(one, 1, one, t0)  // y1 = x1 + t
	return b.slp([]int{y0, y1, z2})
}

// emitM4 appends the 8-addition Duval-Leurent M4 program over the four values
// at indices `in`, folding `scale` into the input coefficients, and returns the
// indices of the four outputs.
func emitM4(b *slpBuilder, in [4]int, alpha, scale int64) [4]int {
	s := big.NewInt(scale)
	as := big.NewInt(alpha * scale)
	a2 := big.NewInt(alpha * alpha)
	one := big.NewInt(1)

	t0 := b.add(s, in[0], s, in[1]) // t0 = s*x0 + s*x1
	t1 := b.add(s, in[2], s, in[3]) // t1 = s*x2 + s*x3
	t2 := b.add(as, in[1], one, t1) // t2 = a*s*x1 + t1
	t3 := b.add(as, in[3], one, t0) // t3 = a*s*x3 + t0
	t4 := b.add(a2, t1, one, t3)    // t4 = a^2*t1 + t3
	t5 := b.add(a2, t0, one, t2)    // t5 = a^2*t0 + t2
	return [4]int{b.add(one, t3, one, t5), t5, b.add(one, t2, one, t4), t4}
}

// slpBuilder accumulates steps and hands back the index of each result, so
// generated programs can compose (feed one program's outputs into another)
// without hand-tracking offsets.
type slpBuilder struct {
	t     int
	steps []SLPStep
}

// add appends w = c1*values[a] + c2*values[b] and returns w's index.
func (b *slpBuilder) add(c1 *big.Int, a int, c2 *big.Int, bi int) int {
	b.steps = append(b.steps, SLPStep{C1: c1, A: a, C2: c2, B: bi})
	return b.t + len(b.steps) - 1
}

func (b *slpBuilder) slp(out []int) SLP { return SLP{T: b.t, Steps: b.steps, Out: out} }
