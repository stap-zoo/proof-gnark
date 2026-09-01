package algebra

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark/constraint/solver"
	"github.com/consensys/gnark/frontend"
)

// Registering the inverse-power hint at package load makes it available to the
// real R1CS solver (the test engine calls the hint value directly, but proving
// resolves hints from this global registry).
func init() {
	solver.RegisterHint(invPowHint)
}

// invPowHint computes w = x^(1/alpha) natively, out of circuit. inputs = [x, alpha];
// the single output is x raised to alpha^{-1} mod (p-1), the field-specific inverse
// exponent. Everything the hint needs — including alpha^{-1} — is derived from the
// field modulus, so one hint serves every field and every alpha.
func invPowHint(field *big.Int, inputs []*big.Int, outputs []*big.Int) error {
	if len(inputs) != 2 || len(outputs) != 1 {
		return fmt.Errorf("algebra: invPowHint expects 2 inputs and 1 output, got %d and %d", len(inputs), len(outputs))
	}
	x, alpha := inputs[0], inputs[1]
	pMinus1 := new(big.Int).Sub(field, big.NewInt(1))
	alphaInv := new(big.Int).ModInverse(alpha, pMinus1)
	if alphaInv == nil {
		return fmt.Errorf("algebra: alpha=%s is not invertible mod p-1 (x^alpha is not a permutation)", alpha)
	}
	outputs[0].Exp(x, alphaInv, field)
	return nil
}

// InvPow returns the value w satisfying w^alpha == x, i.e. the inverse power map
// x^{1/alpha}. The high-degree inverse exponent is never evaluated in-circuit:
// a hint computes w natively (out of circuit), and the returned wire is pinned by
// the single low-degree constraint w^alpha == x. This is the standard "verify a
// high-degree inverse with a cheap forward map" trick, shared by the Anemoi
// Flystel and Rescue-Prime's inverse S-box. alpha must be the field's power-map
// exponent (gcd(alpha, p-1) == 1); otherwise the hint fails at solving time.
//
// The forward check does not materialise w^alpha: its last multiplication shares a
// PLONK row with the equality (AssertPowIsEqual), so the whole thing costs
// ceil(log2(alpha)) + popcount(alpha) - 1 gates rather than one more. That matters
// more than it looks — this is the only nonlinear work Rescue-Prime's backward
// half-rounds do, t times per round.
func InvPow(api frontend.API, x frontend.Variable, alpha int) frontend.Variable {
	if alpha < 1 {
		panic("algebra: InvPow alpha must be >= 1")
	}
	out, err := api.NewHint(invPowHint, 1, x, alpha)
	if err != nil {
		panic(fmt.Sprintf("algebra: InvPow hint: %v", err))
	}
	w := out[0]
	AssertPowIsEqual(api, w, alpha, x)
	return w
}
