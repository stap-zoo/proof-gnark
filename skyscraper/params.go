package skyscraper

import (
	"encoding/binary"
	"fmt"
	"math/big"
	"math/bits"

	"github.com/zkhash-sok/gnark-hashes/sampler"
)

// Parameters is the fully-expanded, field-specific Skyscraper instance the
// permutation reads, mirroring SkyscraperParams in the reference skyscraper/params.py.
// The chosen numbers (field, extension degree, R, r/c/d, wordBits) are passed to
// NewParameters; the derived data — the Montgomery squaring constant, the SHA-256
// round constants, the word S-box lookup table, and the Bar decomposition layout —
// are rebuilt from them.
//
// The extension degree n is the design's only width knob (state t = 2n); what it
// changes, and what it does not, is the package doc's "Width" section in
// skyscraper.go. Here it decides three things: the branch field itself (Beta), the
// shape of the round constants (R rows of N), and the length of the digit list a Bar
// rotates as one (N*NumWords, not NumWords).
type Parameters struct {
	Field *big.Int // BASE field characteristic p; a branch lives in GF(p^N)
	N     int      // extension degree n; a Feistel branch is one GF(p^n) element
	T     int      // Feistel state width in BASE-field elements, t = 2n

	// Beta is the extension modulus' constant term: GF(p^n) = F_p[X]/(X^n + Beta),
	// so X^n reduces to -Beta and that is all the Square layer needs of it. Every
	// modulus the reference instantiates is such a binomial (its fmod is [Beta, 0,
	// ..., 1]). Nil exactly when N == 1, where there is no extension.
	Beta *big.Int

	R         int          // number of Feistel rounds
	BarRounds map[int]bool // round indices using the Bar S-box (all others use Square)
	Rcons     [][]*big.Int // R rounds x N coordinates (the first and last rows are zero)
	MontRInv  *big.Int     // Montgomery constant (2^machineBits)^{-1}; the Square scale x -> x^2 * MontRInv

	// Bar layer decomposition layout: each of the N coordinates becomes NumWords
	// big-endian words of WordBits bits; the N*NumWords words are concatenated,
	// cyclically rotated left by Rot, and each is passed through SboxLUT.
	WordBits int
	NumWords int      // words per COORDINATE (the flattened list is N*NumWords long)
	Rot      int      // rotation of the FLATTENED word list, in words
	SboxLUT  []uint64 // word S-box table, SboxLUT[i] = byte S-box applied to each byte of i, 2^WordBits entries

	Rate     int
	Capacity int
	Digest   int
}

// NewParameters builds and validates Skyscraper parameters over the given base
// field. n is the extension degree of the branch field GF(p^n) and beta its modulus'
// constant term (X^n + beta; nil at n = 1), R the Feistel round count and r/c/d the
// sponge parameters — all taken explicitly per instance, as in the reference
// instances.py. The Montgomery constant, round constants, S-box table, and Bar layout
// are derived deterministically from (p, n, R).
//
// wordBits sets the Bar decomposition granularity: each COORDINATE of the branch is
// split into wordBits-bit words, each passed through a 2^wordBits-entry S-box table.
// It must be a whole number of bytes (the word S-box composes the byte S-box over each
// byte of a word) and is capped at 16 (a 2^16 table); 8 is the reference byte-word
// Bar, 16 the two-byte-word variant. Both produce byte-identical output, so they share
// the same known-answer vectors.
func NewParameters(field *big.Int, n int, beta *big.Int, R, rate, capacity, digest, wordBits int) (*Parameters, error) {
	if field == nil || field.Sign() <= 0 {
		return nil, fmt.Errorf("skyscraper: field modulus must be a positive integer")
	}
	if field.Bit(0) == 0 {
		return nil, fmt.Errorf("skyscraper: field modulus must be odd, got an even value")
	}
	if R < 2 {
		return nil, fmt.Errorf("skyscraper: R must be >= 2 (first and last rounds are constant-free), got %d", R)
	}
	if n < 1 {
		return nil, fmt.Errorf("skyscraper: extension degree n must be >= 1, got %d", n)
	}
	beta, err := checkExtension(field, n, beta)
	if err != nil {
		return nil, err
	}
	t := 2 * n // 2-branch Feistel over GF(p^n)
	if rate+capacity != t {
		return nil, fmt.Errorf("skyscraper: sponge invariant violated: r+c=%d != t=%d", rate+capacity, t)
	}
	if digest < 1 {
		return nil, fmt.Errorf("skyscraper: digest size must be >= 1, got %d", digest)
	}

	// Bar decomposition into wordBits-bit big-endian words, PER COORDINATE. The
	// reference splits each coordinate into m = ceil(bits/8) byte digits and rotates
	// the concatenated n*m digits by m/2; a wordBits-word version tiles the same
	// layout when wordBits is a whole number of bytes. The per-coordinate word count
	// must be even (for the half-word rotation), and the reference's (m/2)-byte
	// rotation must be a whole number of words so the byte permutation is identical —
	// i.e. Rot = numWords/2 words must equal (m/2) bytes in bits.
	if wordBits < 8 || wordBits%8 != 0 {
		return nil, fmt.Errorf("skyscraper: wordBits must be a positive multiple of 8 (byte-composite S-box), got %d", wordBits)
	}
	if wordBits > 16 {
		return nil, fmt.Errorf("skyscraper: wordBits=%d unsupported (2^%d-entry table is impractical; use 8 or 16)", wordBits, wordBits)
	}
	bitLen := field.BitLen()
	numWords := (bitLen + wordBits - 1) / wordBits
	if numWords%2 != 0 {
		return nil, fmt.Errorf("skyscraper: word length %d (wordBits=%d) must be even for the Bar half-rotation", numWords, wordBits)
	}
	byteLen := (bitLen + 7) / 8
	if byteLen%2 != 0 {
		return nil, fmt.Errorf("skyscraper: byte length %d must be even for the Bar half-rotation", byteLen)
	}
	if (numWords/2)*wordBits != (byteLen/2)*8 {
		return nil, fmt.Errorf("skyscraper: wordBits=%d does not tile the %d-byte Bar rotation evenly", wordBits, byteLen/2)
	}

	p := &Parameters{
		Field:     new(big.Int).Set(field),
		N:         n,
		T:         t,
		Beta:      beta,
		R:         R,
		BarRounds: map[int]bool{6: true, 7: true, 10: true, 11: true}, // reference default
		Rcons:     deriveRcons(field, R, n),
		MontRInv:  montgomeryRInv(field),
		WordBits:  wordBits,
		NumWords:  numWords,
		Rot:       numWords / 2,
		SboxLUT:   buildSboxLUT(wordBits),
		Rate:      rate,
		Capacity:  capacity,
		Digest:    digest,
	}
	return p, nil
}

// ---------------------------------------------------------------------------
// The branch field GF(p^n) = F_p[X]/(X^n + beta)
//
// The reference passes the modulus as a coefficient list (fmod = [beta, 0, ..., 1])
// and derives the coordinate polynomials of x -> x^2 by building the extension in
// Sage. Here the modulus is a BINOMIAL — which every reference instance's is — so
// the reduction is just X^n = -beta and the squaring is written out directly
// (skyscraper.go, engine.square); this function is the guard that beta really does
// define a field.
// ---------------------------------------------------------------------------

// checkExtension validates (n, beta) and returns the reduced beta to store.
//
// Irreducibility of a binomial is checked, not assumed: for n = 2, X^2 + beta is
// irreducible over F_p exactly when -beta is a quadratic non-residue, which Euler's
// criterion decides in one exponentiation. Degrees above 2 are rejected rather than
// waved through — the criterion is different for each degree (for prime n dividing
// p-1 it is that -beta is not an n-th power, and n coprime to p-1 is always
// reducible), and the Square layer would need its own reduction as well, so shipping
// an untested branch would be worse than not having one. The reference's N3
// instances are what that would be for.
func checkExtension(p *big.Int, n int, beta *big.Int) (*big.Int, error) {
	if n == 1 {
		if beta != nil {
			return nil, fmt.Errorf("skyscraper: n=1 is the prime field and takes no extension modulus, got beta=%v", beta)
		}
		return nil, nil
	}
	if n > 2 {
		return nil, fmt.Errorf("skyscraper: extension degree n=%d not implemented (n=1 and n=2 are); "+
			"a higher degree needs an irreducibility criterion of its own and a Square layer to match", n)
	}
	if beta == nil {
		return nil, fmt.Errorf("skyscraper: n=%d needs an extension modulus X^%d + beta, got beta=nil", n, n)
	}
	b := new(big.Int).Mod(beta, p)
	if b.Sign() == 0 {
		return nil, fmt.Errorf("skyscraper: beta must be non-zero mod p (X^%d + 0 is reducible)", n)
	}
	// X^2 + b is irreducible iff -b is a non-residue: (-b)^((p-1)/2) == -1.
	negB := new(big.Int).Sub(p, b)
	exp := new(big.Int).Rsh(new(big.Int).Sub(p, big.NewInt(1)), 1)
	if new(big.Int).Exp(negB, exp, p).Cmp(big.NewInt(1)) == 0 {
		return nil, fmt.Errorf("skyscraper: X^2 + %v is reducible over this field (-beta is a quadratic residue)", b)
	}
	return b, nil
}

// ---------------------------------------------------------------------------
// Montgomery squaring constant (ref params.py: mont_R = to_field(2^machineBit),
// mont_R_inv = mont_R^{-1}). machineBit is the field element's in-memory width
// rounded up to whole 64-bit limbs, so for the 254/255-bit curve scalar fields
// it is 256 and MontRInv = (2^256)^{-1} mod p. The Square layer scales x^2 by
// this constant, matching the reference implementation's Montgomery-space square.
// It is a BASE-field constant and does not depend on n: the reference scales every
// coordinate of the extension square by the same mont_R_inv.
// ---------------------------------------------------------------------------

func montgomeryRInv(field *big.Int) *big.Int {
	machineBit := ((field.BitLen() + 63) / 64) * 64
	montR := new(big.Int).Exp(big.NewInt(2), big.NewInt(int64(machineBit)), field)
	return new(big.Int).ModInverse(montR, field)
}

// ---------------------------------------------------------------------------
// Round constants (ref params.py _init_cons). One constant per (round, coordinate):
// constant j (0 <= j < (R-2)*n) is the SHA-256 digest of the 32-byte seed (j as
// big-endian uint32) || "Skyscraper" || zero-pad, read big-endian and reduced mod p.
// The stream is laid out round-major, n coordinates at a time, so the n = 1 chain is
// the prefix of every wider one. The first and last Feistel rounds add nothing (their
// constants are zero by construction). Matches the round constants of the Skyscraper
// reference (github.com/Skyscraper-Hash/skyscraper-sage).
// ---------------------------------------------------------------------------

func deriveRcons(field *big.Int, R, n int) [][]*big.Int {
	rcons := make([][]*big.Int, R)
	rcons[0] = zeroRow(n) // first round: no constant
	j := 0
	for i := 1; i < R-1; i++ {
		row := make([]*big.Int, n)
		for k := range row {
			row[k] = sampler.SHA256Mod(rconSeed(j), field)
			j++
		}
		rcons[i] = row
	}
	rcons[R-1] = zeroRow(n) // last round: no constant
	return rcons
}

func zeroRow(n int) []*big.Int {
	row := make([]*big.Int, n)
	for i := range row {
		row[i] = big.NewInt(0)
	}
	return row
}

func rconSeed(j int) []byte {
	seed := make([]byte, 32) // zero-initialized (the trailing pad)
	binary.BigEndian.PutUint32(seed[0:4], uint32(j))
	copy(seed[4:], []byte("Skyscraper"))
	return seed
}

// ---------------------------------------------------------------------------
// Byte S-box (ref utils/lut.py monolith_lut8). The reilabs reference expresses
// the same map compactly as sboxByte; the two are byte-for-byte identical (the
// reference's Daemen chi-landscape "001*" composed with a 1-bit rotation equals
// this bit-rotation formula), verified across all 256 inputs.
// ---------------------------------------------------------------------------

func sboxByte(b byte) byte {
	x := bits.RotateLeft8(^b, 1)
	y := bits.RotateLeft8(b, 2)
	z := bits.RotateLeft8(b, 3)
	return bits.RotateLeft8(b^(x&y&z), 1)
}

// buildSboxLUT builds the 2^wordBits-entry word S-box table. Each entry composes
// the byte S-box over the wordBits/8 bytes of the index, most-significant byte
// first: entry(i) = concat_b sboxByte(byte b of i). For wordBits=8 this is the
// plain byte table; for wordBits=16, entry(i) = sboxByte(i>>8)<<8 | sboxByte(i&0xff).
// Because the map is byte-wise, applying it to wordBits-words is byte-for-byte
// identical to applying the byte S-box to each byte — which is why the wider word
// size keeps the same known-answer vectors.
func buildSboxLUT(wordBits int) []uint64 {
	numBytes := wordBits / 8
	lut := make([]uint64, 1<<uint(wordBits))
	for i := range lut {
		var out uint64
		for b := 0; b < numBytes; b++ {
			shift := uint(8 * (numBytes - 1 - b)) // most-significant byte first
			by := byte((uint64(i) >> shift) & 0xff)
			out |= uint64(sboxByte(by)) << shift
		}
		lut[i] = out
	}
	return lut
}
