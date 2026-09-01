package algebra

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark/frontend"
)

// SLPStep is one straight-line-program gate w = C1*a + C2*b, where a and b
// reference earlier values by index: 0..T-1 are the input vector entries, and
// T+k is the output of Steps[k]. Coefficients are field constants, so in PLONK
// each step is exactly one fan-in-2 addition gate (constant multiplications
// fold into the gate's selectors) regardless of how large they are.
type SLPStep struct {
	C1 *big.Int
	A  int
	C2 *big.Int
	B  int
}

// SLP is a straight-line program computing a constant matrix-vector product
// M*x with a fixed sequence of two-term additions — the representation behind
// "low-addition" linear layers (Duval-Leurent ToSC 2018(2); Polocolo eprint
// 2025/926 Appendix A) and behind shared-sum forms such as J+diag. Applying it
// in-circuit realizes the published PLONK addition-gate counts, which a dense
// MatVecMul cannot: the SCS compiler charges one gate per nonzero past the
// first in every row and never discovers cross-row sharing (minimal SLP
// extraction is NP-hard). Under R1CS both forms are free (linear layers fold
// into wire coefficients), so an SLP only changes the PLONK column.
//
// Out lists which values (same indexing as SLPStep) form the output vector.
//
// Every SLP in this repo is validated against the matrix it replaces by
// comparing Matrix(p) entry-by-entry — see algebra/lin_slp_test.go and the
// per-construction parameter tests. A program is only ever used in place of a
// dense product once that comparison passes.
type SLP struct {
	T     int // input vector length
	Steps []SLPStep
	Out   []int
}

// NewSLP builds a program from a compact step table, each row {c1, a, c2, b}
// for w = c1*values[a] + c2*values[b]. This is the transcription form for
// programs copied out of a paper — their coefficients are always small, and a
// table of rows diffs line-for-line against the published listing. Use the
// struct form directly when a coefficient is a full field element.
func NewSLP(t int, steps [][4]int64, out []int) SLP {
	s := SLP{T: t, Steps: make([]SLPStep, len(steps)), Out: out}
	for i, r := range steps {
		s.Steps[i] = SLPStep{C1: big.NewInt(r[0]), A: int(r[1]), C2: big.NewInt(r[2]), B: int(r[3])}
	}
	return s
}

// Gates reports the number of PLONK addition gates Apply emits: one per step.
func (s SLP) Gates() int { return len(s.Steps) }

// Apply evaluates the program over circuit variables: len(in) must be T, and
// the result has len(Out) entries — one PLONK addition gate per step.
func (s SLP) Apply(api frontend.API, in []frontend.Variable) []frontend.Variable {
	return s.ApplyConst(api, in, nil)
}

// ApplyConst evaluates the program and adds consts[i] to output i (a nil slice, a
// short slice, or a nil/zero entry means "no constant there"). Constants cost
// nothing: each rides on the gate that produces its output, since a PLONK gate
// carries a constant selector — `Add(a, b, c)` is one gate where
// `Add(Add(a, b), c)` is two.
//
// Shifting a step that later steps also read would corrupt them, so the shift is
// compensated where it lands: a step reading a value shifted by d with coefficient
// C absorbs -C*d into its own constant selector, which is free because that gate
// exists anyway. So shifts never propagate past one level and sharing costs
// nothing — the shared sum of a J+diag program and the M4 program's t4/t5 (both
// outputs *and* inputs to later steps) fold like anything else.
//
// Two cases still cost a gate, exactly what a separate addition would: an output
// that is a bare input (no gate of ours produces it) and a second output pointing
// at the same step as an earlier one (one step cannot carry two shifts).
func (s SLP) ApplyConst(api frontend.API, in []frontend.Variable, consts []*big.Int) []frontend.Variable {
	if len(in) != s.T {
		panic(fmt.Sprintf("algebra: SLP expects %d inputs, got %d", s.T, len(in)))
	}

	// shift[v] is the constant the value at index v carries relative to what the
	// bare program would compute. Only steps claimed by an output are shifted.
	shift := make([]*big.Int, s.T+len(s.Steps))
	for i, idx := range s.Out {
		c := constAt(consts, i)
		if c == nil || c.Sign() == 0 || idx < s.T || shift[idx] != nil {
			continue
		}
		shift[idx] = c
	}

	vals := make([]frontend.Variable, s.T, s.T+len(s.Steps))
	copy(vals, in)
	for k, st := range s.Steps {
		// The gate's constant: the shift this step must carry, minus what its
		// inputs already carry into it.
		fold := new(big.Int)
		if want := shift[s.T+k]; want != nil {
			fold.Set(want)
		}
		fold.Sub(fold, mulShift(st.C1, shift[st.A]))
		fold.Sub(fold, mulShift(st.C2, shift[st.B]))

		a, b := scale(api, st.C1, vals[st.A]), scale(api, st.C2, vals[st.B])
		if fold.Sign() != 0 {
			vals = append(vals, api.Add(a, b, fold))
		} else {
			vals = append(vals, api.Add(a, b))
		}
	}

	out := make([]frontend.Variable, len(s.Out))
	for i, idx := range s.Out {
		out[i] = vals[idx]
		// What this output still owes: its constant less the shift it already carries.
		owed := new(big.Int)
		if c := constAt(consts, i); c != nil {
			owed.Set(c)
		}
		if got := shift[idx]; got != nil {
			owed.Sub(owed, got)
		}
		if owed.Sign() != 0 {
			out[i] = api.Add(out[i], owed)
		}
	}
	return out
}

// mulShift returns c*d, treating a nil shift as zero.
func mulShift(c, d *big.Int) *big.Int {
	if d == nil {
		return new(big.Int)
	}
	return new(big.Int).Mul(c, d)
}

func scale(api frontend.API, c *big.Int, v frontend.Variable) frontend.Variable {
	if c.Cmp(bigOne) == 0 {
		return v
	}
	return api.Mul(c, v)
}

// Matrix evaluates the program symbolically and returns the matrix M it
// computes, one coefficient row per output. Entries are reduced mod `mod`; pass
// nil to evaluate over the integers. This is how a program is asserted to equal
// the matrix it replaces — transcribed ones may contain a typo, and generated
// ones may encode the wrong shape, so nothing is applied in-circuit until the
// comparison passes.
func (s SLP) Matrix(mod *big.Int) [][]*big.Int {
	rows := make([][]*big.Int, 0, s.T+len(s.Steps))
	for i := 0; i < s.T; i++ {
		row := make([]*big.Int, s.T)
		for j := range row {
			row[j] = big.NewInt(0)
		}
		row[i] = big.NewInt(1)
		rows = append(rows, row)
	}
	for _, st := range s.Steps {
		row := make([]*big.Int, s.T)
		for j := 0; j < s.T; j++ {
			row[j] = new(big.Int).Mul(st.C1, rows[st.A][j])
			row[j].Add(row[j], new(big.Int).Mul(st.C2, rows[st.B][j]))
			if mod != nil {
				row[j].Mod(row[j], mod)
			}
		}
		rows = append(rows, row)
	}
	out := make([][]*big.Int, len(s.Out))
	for i, idx := range s.Out {
		out[i] = rows[idx]
	}
	return out
}

// PermuteInputs returns the program that reads input slot perm[i] wherever s
// reads input i, i.e. s'(x) = s(x[perm[0]], ..., x[perm[T-1]]). In matrix terms
// the columns are permuted: s'.Matrix()[o][perm[j]] == s.Matrix()[o][j]. Anemoi
// uses this for the y-lane, whose matrix M_y is M_x with every row rotated
// right — equivalently M_x applied to the left-rotated input.
func (s SLP) PermuteInputs(perm []int) SLP {
	if len(perm) != s.T {
		panic(fmt.Sprintf("algebra: PermuteInputs expects %d indices, got %d", s.T, len(perm)))
	}
	remap := func(i int) int {
		if i < s.T {
			return perm[i]
		}
		return i
	}
	out := SLP{T: s.T, Steps: make([]SLPStep, len(s.Steps)), Out: make([]int, len(s.Out))}
	for i, st := range s.Steps {
		out.Steps[i] = SLPStep{C1: st.C1, A: remap(st.A), C2: st.C2, B: remap(st.B)}
	}
	for i, idx := range s.Out {
		out.Out[i] = remap(idx)
	}
	return out
}

// DenseGates models the PLONK addition-gate cost of MatVecMul(m, ·): one gate
// per nonzero entry past the first in each row (the zero accumulator the first
// term lands on is folded away, and constant multiplications are free). It is
// the yardstick a candidate SLP has to beat.
//
// Treat it as a close upper bound, not an exact count: the SCS compiler
// deduplicates identical linear combinations, so a matrix with repeating rows
// measures cheaper than this (the t=8 block-circulant costs 53, not 56, and a
// J+diag whose diagonal is 1 in three places costs 5, not 12). Every selection
// made with this model is confirmed against harness.Count(harness.PLONK, ·) in
// the tests for that construction.
func DenseGates(m [][]*big.Int) int {
	n := 0
	for _, row := range m {
		nz := 0
		for _, v := range row {
			if v.Sign() != 0 {
				nz++
			}
		}
		if nz > 1 {
			n += nz - 1
		}
	}
	return n
}

// ApplyMatrix computes m*in, using the straight-line program when one is
// supplied (slp != nil) and the dense product otherwise. Constructions store a
// program only when it is strictly cheaper than dense in PLONK, so this is the
// single place a linear layer is applied and the cost choice stays in the
// parameters where it can be tested.
func ApplyMatrix(api frontend.API, m [][]*big.Int, slp *SLP, in []frontend.Variable) []frontend.Variable {
	return ApplyMatrixConst(api, m, slp, in, nil)
}

// ApplyMatrixConst computes m*in + consts through the program when one is supplied
// and the dense product otherwise, folding each constant into the gate that
// produces its output wherever that is free (see MatVecMulConst / SLP.ApplyConst).
// A construction calls this with the round constants that follow the layer, so ARK
// stops costing gates of its own.
func ApplyMatrixConst(api frontend.API, m [][]*big.Int, slp *SLP, in []frontend.Variable, consts []*big.Int) []frontend.Variable {
	if slp != nil {
		return slp.ApplyConst(api, in, consts)
	}
	return MatVecMulConst(api, m, in, consts)
}

// SelectSLP is what a parameter constructor calls for each linear layer: it
// verifies that slp computes exactly m over the field p, then returns the
// program if it is cheaper than the dense product in PLONK and nil if it is not
// (t = 2 layers, where the shared sum costs more than the two products it
// replaces). Pass the result to ApplyMatrix.
//
// An error means the program does not compute m — always a bug, never a cost
// question — so no construction can ship a program that was never checked
// against its matrix, whether or not a test covers it.
func SelectSLP(slp SLP, m [][]*big.Int, p *big.Int) (*SLP, error) {
	if err := CheckSLP(slp, m, p); err != nil {
		return nil, err
	}
	if slp.Gates() >= DenseGates(m) {
		return nil, nil
	}
	return &slp, nil
}

// CheckSLP reports an error unless slp computes exactly the matrix m over the
// field p.
func CheckSLP(slp SLP, m [][]*big.Int, p *big.Int) error {
	got := slp.Matrix(p)
	if len(got) != len(m) {
		return fmt.Errorf("program has %d outputs, matrix has %d rows", len(got), len(m))
	}
	for i := range m {
		if len(got[i]) != len(m[i]) {
			return fmt.Errorf("row %d: program width %d, matrix width %d", i, len(got[i]), len(m[i]))
		}
		for j := range m[i] {
			want := new(big.Int).Mod(m[i][j], p)
			if got[i][j].Cmp(want) != 0 {
				return fmt.Errorf("program computes M[%d][%d] = %s, matrix has %s", i, j, got[i][j], want)
			}
		}
	}
	return nil
}
