package algebra

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark/constraint/solver"
	"github.com/consensys/gnark/frontend"
)

// Registering the decomposition hints at package load makes them available to the
// real solver (the test engine calls hints directly, but proving resolves them
// from this global registry).
func init() {
	solver.RegisterHint(decomposeBEHint)
}

// decomposeBEHint splits inputs[1] into len(outputs) big-endian words of
// inputs[0] bits each, most-significant word first. The solver passes the
// canonical representative of the value (in [0, field)), so the words are its
// unique base-2^wordBits digits. Purely a witness generator — every output is
// unconstrained until the decomposition below pins it.
func decomposeBEHint(_ *big.Int, inputs []*big.Int, outputs []*big.Int) error {
	if len(inputs) != 2 {
		return fmt.Errorf("algebra: decomposeBEHint expects 2 inputs, got %d", len(inputs))
	}
	wordBits := uint(inputs[0].Uint64())
	v := inputs[1]
	numWords := len(outputs)
	mask := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), wordBits), big.NewInt(1))
	for i := 0; i < numWords; i++ {
		shift := wordBits * uint(numWords-1-i)
		w := new(big.Int).Rsh(v, shift)
		outputs[i].And(w, mask)
	}
	return nil
}

// CanonicalDecomposeExternalWordRange is CanonicalDecompose with the per-word range
// checks omitted: it decomposes v, ties the words back to v (recomposition) and
// enforces canonicity, but does NOT itself constrain each word to wordBits bits.
//
// The CALLER MUST range-bind every returned word to [0, 2^wordBits) — e.g. by
// looking each one up in a 2^wordBits-entry table. This exists for the split-and-
// lookup S-box (Skyscraper's Bar), where each word is immediately fed through the
// S-box lookup table: that lookup already binds the word to the table's index
// domain [0, 2^wordBits), so a separate range check would be a redundant second
// lookup on the same wire. If the caller does NOT bind every word, BOTH results are
// unsound — an unbounded word breaks the recomposition (the prover picks any words
// summing to v mod p) and lets the canonicity limbs hi/lo exceed 2^halfBits, which
// assertLessThanModulus assumes they do not.
//
// rchk is still required: canonicity needs one halfBits-wide range check on a
// derived difference (see assertLessThanModulus), which is not a word of v and so is
// not covered by the caller's per-word lookups. Skyscraper passes gnark's
// std/rangecheck; skyscraper.bar's doc comment says why that beats routing the check
// back through the caller's own table.
func CanonicalDecomposeExternalWordRange(api frontend.API, rchk frontend.Rangechecker, v frontend.Variable, wordBits, numWords int) []frontend.Variable {
	return canonicalDecompose(api, rchk, v, wordBits, numWords, false)
}

func canonicalDecompose(api frontend.API, rchk frontend.Rangechecker, v frontend.Variable, wordBits, numWords int, checkWords bool) []frontend.Variable {
	if wordBits < 1 {
		panic("algebra: canonical decompose: wordBits must be >= 1")
	}
	if numWords < 2 || numWords%2 != 0 {
		panic("algebra: canonical decompose: numWords must be even and >= 2")
	}
	p := api.Compiler().Field()
	if wordBits*numWords < p.BitLen() {
		panic(fmt.Sprintf("algebra: canonical decompose: %d words of %d bits do not cover the %d-bit field",
			numWords, wordBits, p.BitLen()))
	}

	out, err := api.NewHint(decomposeBEHint, numWords, wordBits, v)
	if err != nil {
		panic(fmt.Sprintf("algebra: canonical decompose hint: %v", err))
	}
	words := out
	if checkWords {
		for _, w := range words {
			rchk.Check(w, wordBits)
		}
	}

	// The words are folded back to v *through* the two canonicity limbs rather than
	// by a recomposition of their own: hi and lo are needed by assertLessThanModulus
	// anyway, and v == hi*2^halfBits + lo says the same thing as
	// v == sum(w_i * radix^i) for a fraction of the PLONK rows — one two-term row
	// instead of a second full numWords-step accumulation. (R1CS is indifferent;
	// every recomposition is a free linear combination there.)
	//
	// That last row does NOT collapse into a single PLONK constraint the way the
	// per-word accumulation steps do: 2^halfBits is 2^128 here and gnark exposes gate
	// selectors as ints, so the coefficient cannot be a selector and the addition and
	// the equality stay two rows.
	half := numWords / 2
	halfBits := half * wordBits
	hi := RecomposeBE(api, words[:half], wordBits)
	lo := RecomposeBE(api, words[half:], wordBits)
	shift := new(big.Int).Lsh(big.NewInt(1), uint(halfBits))
	api.AssertIsEqual(api.Add(api.Mul(hi, shift), lo), v)

	assertLessThanModulus(api, rchk, hi, lo, p, halfBits)
	return words
}

// RecomposeBE folds big-endian words (most-significant first) back into a single
// value with the constant radix 2^wordBits — the inverse of the decomposition
// the decomposition produces. Multiplying by a constant folds into wire
// coefficients, so this adds no multiplication constraints.
func RecomposeBE(api frontend.API, words []frontend.Variable, wordBits int) frontend.Variable {
	radix := new(big.Int).Lsh(big.NewInt(1), uint(wordBits))
	var acc frontend.Variable = 0
	for _, w := range words {
		acc = api.Add(api.Mul(acc, radix), w)
	}
	return acc
}

// assertLessThanModulus asserts that the two-limb value (hi, lo) — with each limb
// < 2^halfBits, which the caller must already have enforced — is strictly less than
// the field modulus p = P_hi*2^halfBits + P_lo. It generalizes reilabs' hard-coded
// 128-bit BN254 split to any field whose low limb (p mod 2^halfBits) is non-zero,
// true for the curve scalar fields here.
//
// (hi, lo) < (P_hi, P_lo) is a lexicographic comparison: either hi < P_hi, or
// hi == P_hi and lo < P_lo. Which case applies is decided by eq = [hi == P_hi], and
// the two cases are the same shape — a difference that must be non-negative — so a
// select feeds ONE of them to ONE range check:
//
//	eq == 0:  P_hi - 1 - hi  in [0, 2^halfBits)  =>  hi <= P_hi - 1, and lo < 2^halfBits
//	                                                 gives (hi,lo) <= P_hi*2^halfBits - 1 < p
//	eq == 1:  P_lo - 1 - lo  in [0, 2^halfBits)  =>  hi == P_hi and lo <= P_lo - 1 < P_lo
//
// A failing case is caught because the difference wraps: if hi > P_hi then
// P_hi - 1 - hi is p minus something under 2^halfBits, which is nowhere near
// [0, 2^halfBits) as long as p >> 2^halfBits — it is ~2^(2*halfBits) here.
//
// The point of the select is that it halves the range-check bill. Long subtraction
// (a hinted borrow bit carrying between the limbs, then a checked difference *per
// limb*, which is what this used to do) checks 2*halfBits bits where the lexicographic
// form checks halfBits; a range check is the expensive part of the Bar's decomposition,
// so that is a fifth of the whole Skyscraper permutation — see skyscraper.bar.
//
// eq comes from api.IsZero, not a hint: IsZero is an iff, so hi == P_hi in the eq == 1
// branch needs no separate constraint pinning it, and eq needs no separate boolean
// constraint either.
func assertLessThanModulus(api frontend.API, rchk frontend.Rangechecker, hi, lo frontend.Variable, p *big.Int, halfBits int) {
	pow := new(big.Int).Lsh(big.NewInt(1), uint(halfBits))
	modHi := new(big.Int).Rsh(p, uint(halfBits))
	modLo := new(big.Int).Mod(p, pow)
	if modLo.Sign() == 0 {
		panic("algebra: canonical decompose: canonicity split requires a non-zero low limb of the modulus")
	}

	eq := api.IsZero(api.Sub(modHi, hi))
	hiCase := api.Sub(new(big.Int).Sub(modHi, big.NewInt(1)), hi)
	loCase := api.Sub(new(big.Int).Sub(modLo, big.NewInt(1)), lo)
	rchk.Check(api.Select(eq, loCase, hiCase), halfBits)
}
