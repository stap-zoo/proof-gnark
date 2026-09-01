package rescueprime

import (
	"fmt"
	"math/big"

	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific Rescue-Prime instance the
// permutation reads, mirroring RescuePrimeParams in the reference
// marvellous/params.py. The chosen numbers (field, t, alpha, g, R, r/c/d) are
// passed to NewParameters; the derived tables — the dense MDS matrix M and the
// 2*R round-constant rows — are rebuilt from them.
//
// No inverse data (M_inv, alpha_inv) is stored: the in-circuit permutation only
// evaluates forward, and its inverse S-box is handled by algebra.InvPow, whose
// hint derives the inverse exponent itself.
type Parameters struct {
	Field *big.Int // field characteristic p
	T     int      // state size
	R     int      // number of (double) rounds; the permutation runs 2*R half-rounds
	Alpha int      // S-box power-map exponent

	M              [][]*big.Int // t x t MDS matrix (Vandermonde-derived, transposed)
	RoundConstants [][]*big.Int // 2*R rows of t constants, one per half-round

	Rate     int
	Capacity int
	Digest   int
}

// NewParameters builds and validates Rescue-Prime parameters for the given
// field. g is a multiplicative generator (used to build the MDS matrix), R the
// round count, and r/c/d the sponge parameters — all taken explicitly per
// instance, as in the reference instances.py. The MDS matrix and round
// constants are derived.
func NewParameters(field *big.Int, t, alpha int, g *big.Int, R, rate, capacity, digest int) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("rescueprime: field modulus must be a positive integer")
	}
	if field.Bit(0) == 0 {
		return nil, fmt.Errorf("rescueprime: field modulus must be odd, got an even value")
	}
	if t < 1 {
		return nil, fmt.Errorf("rescueprime: t must be >= 1, got %d", t)
	}
	if alpha < 3 {
		return nil, fmt.Errorf("rescueprime: alpha must be >= 3, got %d", alpha)
	}
	if new(big.Int).GCD(nil, nil, big.NewInt(int64(alpha)), new(big.Int).Sub(field, big.NewInt(1))).Cmp(big.NewInt(1)) != 0 {
		return nil, fmt.Errorf("rescueprime: x^%d is not a permutation (gcd(alpha, p-1) != 1)", alpha)
	}
	if g == nil || g.Sign() <= 0 {
		return nil, fmt.Errorf("rescueprime: generator g must be a positive integer")
	}
	if rate+capacity != t {
		return nil, fmt.Errorf("rescueprime: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if R < 1 {
		return nil, fmt.Errorf("rescueprime: R must be >= 1, got %d", R)
	}
	if digest < 1 {
		return nil, fmt.Errorf("rescueprime: digest size must be >= 1, got %d", digest)
	}

	m, err := vandermondeMDSTranspose(field, t, new(big.Int).Mod(g, field))
	if err != nil {
		return nil, fmt.Errorf("rescueprime: MDS matrix: %w", err)
	}

	p := &Parameters{
		Field:          new(big.Int).Set(field),
		T:              t,
		R:              R,
		Alpha:          alpha,
		M:              m,
		RoundConstants: deriveRoundConstants(field, t, capacity, R),
		Rate:           rate,
		Capacity:       capacity,
		Digest:         digest,
	}
	return p, nil
}

// ---------------------------------------------------------------------------
// Round constants (ref params.py RescuePrimeParams._init_rcons): a 2*R x t grid
// drawn from SHAKE256 with "mod" sampling, seeded by the instance's
// "Rescue-XLIX(p,t,c,kappa)" label. kappa is fixed at 128 (the recommended
// level, and the only one the scalar-field instances use).
// ---------------------------------------------------------------------------

const kappa = 128

func deriveRoundConstants(field *big.Int, t, capacity, R int) [][]*big.Int {
	seed := fmt.Appendf(nil, "Rescue-XLIX(%s,%d,%d,%d)", field.String(), t, capacity, kappa)
	return sampler.Grid(seed, field, 2*R, t)
}

// ---------------------------------------------------------------------------
// MDS matrix (ref utils.matrix.vandermonde_mds_matrix with transpose=True).
// ---------------------------------------------------------------------------

// vandermondeMDSTranspose builds Rescue-Prime's t x t MDS matrix: the right half
// of the reduced row-echelon form of the t x 2t Vandermonde matrix
// V[i][j] = g^(i*j) (0 <= i < t, 0 <= j < 2t), transposed. g must be a primitive
// element, so the left t x t block is an invertible Vandermonde matrix and the
// echelon form is [I | A^{-1}B]; the returned matrix is (A^{-1}B)^T.
func vandermondeMDSTranspose(p *big.Int, t int, g *big.Int) ([][]*big.Int, error) {
	// V[i][j] = g^(i*j) mod p.
	v := make([][]*big.Int, t)
	for i := 0; i < t; i++ {
		v[i] = make([]*big.Int, 2*t)
		for j := 0; j < 2*t; j++ {
			v[i][j] = new(big.Int).Exp(g, big.NewInt(int64(i*j)), p)
		}
	}

	if err := rowReduceMod(v, p); err != nil {
		return nil, err
	}

	// Right half M0 = V[:, t:2t], then transpose: M[i][j] = M0[j][i].
	m := make([][]*big.Int, t)
	for i := 0; i < t; i++ {
		m[i] = make([]*big.Int, t)
		for j := 0; j < t; j++ {
			m[i][j] = new(big.Int).Set(v[j][t+i])
		}
	}
	return m, nil
}

// rowReduceMod reduces m to reduced row-echelon form in place over GF(p),
// mirroring Sage's Matrix.echelon_form. The matrix here has full row rank with
// an invertible leading square block, so every pivot column is a leading column
// (no column pivoting is needed).
func rowReduceMod(m [][]*big.Int, p *big.Int) error {
	rows := len(m)
	if rows == 0 {
		return nil
	}
	cols := len(m[0])
	pivotRow := 0
	for col := 0; col < cols && pivotRow < rows; col++ {
		// Find a nonzero pivot in this column at or below pivotRow.
		sel := -1
		for r := pivotRow; r < rows; r++ {
			if m[r][col].Sign() != 0 {
				sel = r
				break
			}
		}
		if sel == -1 {
			continue // free column
		}
		m[pivotRow], m[sel] = m[sel], m[pivotRow]

		// Normalize the pivot row so the pivot becomes 1.
		inv := new(big.Int).ModInverse(m[pivotRow][col], p)
		if inv == nil {
			return fmt.Errorf("pivot %s not invertible mod p", m[pivotRow][col])
		}
		for j := col; j < cols; j++ {
			m[pivotRow][j].Mod(new(big.Int).Mul(m[pivotRow][j], inv), p)
		}

		// Eliminate this column from every other row.
		for r := 0; r < rows; r++ {
			if r == pivotRow || m[r][col].Sign() == 0 {
				continue
			}
			factor := new(big.Int).Set(m[r][col])
			for j := col; j < cols; j++ {
				term := new(big.Int).Mul(factor, m[pivotRow][j])
				m[r][j].Mod(new(big.Int).Sub(m[r][j], term), p)
			}
		}
		pivotRow++
	}
	return nil
}
