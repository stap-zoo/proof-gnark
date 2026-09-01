package algebra

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark/frontend"
)

var bigOne = big.NewInt(1)

// MatVecMul computes the matrix-vector product m·v for an arbitrary square
// matrix of constant field coefficients. It is the generic linear layer used by
// every construction; the specific matrix (cyclic shift, MDS, ...) comes from
// each hash's parameters.
//
// It is essentially free in-circuit: multiplying a variable by a constant is a
// scalar operation that folds into wire coefficients, so no R1CS constraints are
// produced. Zero coefficients are skipped and unit coefficients add the variable
// directly, keeping the linear combinations minimal for sparse matrices.
func MatVecMul(api frontend.API, m [][]*big.Int, v []frontend.Variable) []frontend.Variable {
	return MatVecMulConst(api, m, v, nil)
}

// MatVecMulConst computes m·v + consts, adding consts[i] to row i's result (a nil
// slice, a short slice, or a nil/zero entry means "no constant there"). Each
// constant rides along on the row's last addition gate instead of costing one of
// its own — a PLONK gate has a constant selector, so `Add(a, b, c)` is one gate
// where `Add(Add(a, b), c)` is two. A row with a single nonzero has no addition to
// ride on and pays a gate, exactly as a separate constant addition would; a zero
// row returns the constant itself.
//
// This is how round constants stop costing gates: a construction hands the round
// constants that follow a linear layer to the layer itself. It changes nothing
// mathematically — the KAT vectors are the proof — and nothing in R1CS.
func MatVecMulConst(api frontend.API, m [][]*big.Int, v []frontend.Variable, consts []*big.Int) []frontend.Variable {
	out := make([]frontend.Variable, len(m))
	for i, row := range m {
		if len(row) != len(v) {
			panic(fmt.Sprintf("algebra: matrix row %d has width %d, vector has length %d", i, len(row), len(v)))
		}
		terms := make([]frontend.Variable, 0, len(row))
		for j, coeff := range row {
			switch {
			case coeff.Sign() == 0:
				// skip zero term
			case coeff.Cmp(bigOne) == 0:
				terms = append(terms, v[j])
			default:
				terms = append(terms, api.Mul(coeff, v[j]))
			}
		}
		out[i] = sumWithConstant(api, terms, constAt(consts, i))
	}
	return out
}

// sumWithConstant returns sum(terms) + c: one addition gate per term past the
// first, with c folded into the last of them for free (a PLONK gate carries a
// constant selector). With fewer than two terms there is no gate to fold into, so
// a nonzero c costs one — the same as adding it separately.
func sumWithConstant(api frontend.API, terms []frontend.Variable, c *big.Int) frontend.Variable {
	pending := c != nil && c.Sign() != 0
	if len(terms) == 0 {
		if pending {
			return c
		}
		return 0
	}
	acc := terms[0]
	for i := 1; i < len(terms); i++ {
		if pending && i == len(terms)-1 {
			acc, pending = api.Add(acc, terms[i], c), false
			continue
		}
		acc = api.Add(acc, terms[i])
	}
	if pending {
		acc = api.Add(acc, c)
	}
	return acc
}

// constAt reads consts[i], treating a nil slice or a short one as "no constant".
func constAt(consts []*big.Int, i int) *big.Int {
	if i < len(consts) {
		return consts[i]
	}
	return nil
}
