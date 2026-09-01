package algebra

import (
	"math/big"

	"github.com/consensys/gnark/frontend"
)

// ClosedFlystel is Anemoi's open Flystel S-box H:(x,y)->(u,v) evaluated through
// its closed (verification) form, so the circuit only ever uses the cheap
// forward power map. See https://eprint.iacr.org/2022/840, Fig. 3.
//
// The open Flystel computes v = y - (x - Q_gamma(y))^{1/alpha}, an expensive
// inverse power. Using H(x,y)=(u,v) iff V(y,v)=(x,u), we recover the same (u,v)
// with a single InvPow (one hint + one degree-alpha check) plus the two quadratic
// Q maps:
//
//	Q_gamma(y) = beta*y^quad + gamma,   Q_delta(v) = beta*v^quad + delta
//	t = x - Q_gamma(y)
//	w = t^{1/alpha}                     (InvPow: hinted, checked by w^alpha == t)
//	v = y - w
//	u = t + Q_delta(v)
//
// quad is the Q-map exponent (2 in odd characteristic); beta, gamma, delta are
// field constants. Cost: log2(alpha) muls for the InvPow check + 2 squares.
//
// Neither Q-map is materialised: each one's additive constant rides the qC selector
// of the row that consumes it — the subtraction for gamma, the addition for delta —
// instead of costing a PLONK gate to add to beta*y^quad first. That is the same
// trick MatVecMulConst plays with round constants, and the reason it is available
// here is that a quadratic's constant term cannot ride the quadratic's OWN row (see
// algebra.MulAffine). Two gates per Flystel, and nothing changes in R1CS.
func ClosedFlystel(api frontend.API, x, y frontend.Variable, alpha, quad int, beta, gamma, delta *big.Int) (u, v frontend.Variable) {
	// t = x - Q_gamma(y) = x - beta*y^quad - gamma
	t := api.Sub(x, api.Mul(beta, Pow(api, y, quad)), gamma)
	w := InvPow(api, t, alpha)
	v = api.Sub(y, w)
	// u = t + Q_delta(v) = t + beta*v^quad + delta
	u = api.Add(t, api.Mul(beta, Pow(api, v, quad)), delta)
	return u, v
}
