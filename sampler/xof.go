package sampler

import (
	"crypto/sha3"
	"math/big"
)

// XOFName selects the extendable-output function backing an XOFSampler.
type XOFName int

const (
	SHAKE128 XOFName = iota
	SHAKE256
)

// Sampling is the strategy mapping raw XOF bytes to a field element in [0, p),
// mirroring ref/utils/sampler.py's FieldElementSampler strategies.
type Sampling int

const (
	// Bitmask reads ceil(bitlen(p)/8) bytes, zeroes the bits above bitlen(p) in
	// the most-significant byte, and rejects (resamples) any value >= p. This is
	// the Rust field_element_from_shake strategy, used by Griffin with SHAKE128.
	Bitmask Sampling = iota
	// Mod reads ceil(bitlen(p)/8)+1 bytes and reduces mod p (no rejection). It is
	// the strategy behind the batch SHAKE256Mod / Grid helpers (Rescue-Prime, Arion).
	Mod
	// Bitshift reads ceil(bitlen(p)/8) bytes and right-SHIFTS the most-significant
	// byte by the excess bits (keeping its high bits, dropping its low bits) instead
	// of masking it, then rejects any value >= p. Matches the HashTape squeeze of
	// Polocolo's param_gen.sage (seq[0] >>= (-size % 8) on a big-endian read); used
	// by Polocolo with SHAKE128 and big-endian draws.
	Bitshift
)

// XOFSampler is a stateful, sequential field-element sampler over a seeded SHAKE
// stream, mirroring ref/utils/sampler.py's XOFFieldElementSampler. Draws are
// consumed one at a time via Next / NextNonzero.
//
// It complements the batch SHAKE256Mod / Grid helpers, which can pre-size their
// read because "mod" never rejects. A sequential reader is required when
// rejection breaks the stream alignment — e.g. Griffin interleaves its round
// constants and its rejection-sampled quadratic coefficients on a single stream,
// so the number of bytes consumed before the coefficients is not known up front.
//
// Reading n bytes at a time from Go's streaming SHAKE is byte-identical to the
// reference's grow-and-buffer digest(n) loop: an XOF's longer output extends its
// shorter one, so sequential reads land on the same stream positions.
type XOFSampler struct {
	h         *sha3.SHAKE
	p         *big.Int
	nBytes    int
	mask      byte
	shift     uint // right-shift of the most-significant byte (Bitshift only)
	reduce    bool // true for Mod (reduce), false for Bitmask/Bitshift (reject)
	bigEndian bool
}

// NewXOFSampler seeds a SHAKE stream and configures the per-draw width and
// reject/reduce rule for field p, with little-endian draws. It mirrors
// XOFFieldElementSampler(seed=seed, p=p, xof=xof, sampling=sampling).
func NewXOFSampler(seed []byte, p *big.Int, xof XOFName, sampling Sampling) *XOFSampler {
	return newXOFSampler(seed, p, xof, sampling, false)
}

// NewXOFSamplerBigEndian is NewXOFSampler with big-endian draws (each chunk's
// first byte is most significant), mirroring
// XOFFieldElementSampler(..., endianess="big") — Polocolo's round-constant
// derivation is (SHAKE128, Bitshift, big-endian).
func NewXOFSamplerBigEndian(seed []byte, p *big.Int, xof XOFName, sampling Sampling) *XOFSampler {
	return newXOFSampler(seed, p, xof, sampling, true)
}

func newXOFSampler(seed []byte, p *big.Int, xof XOFName, sampling Sampling, bigEndian bool) *XOFSampler {
	var h *sha3.SHAKE
	if xof == SHAKE256 {
		h = sha3.NewSHAKE256()
	} else {
		h = sha3.NewSHAKE128()
	}
	// SHAKE Write on an in-memory hash never errors; ignore per the hash.Hash contract.
	_, _ = h.Write(seed)

	bits := p.BitLen()
	s := &XOFSampler{h: h, p: new(big.Int).Set(p), mask: 0xFF, bigEndian: bigEndian}
	switch sampling {
	case Mod:
		s.nBytes = (bits+7)/8 + 1
		s.reduce = true
	case Bitshift:
		s.nBytes = (bits + 7) / 8
		s.shift = uint((8 - bits%8) % 8) // excess bits in the most significant byte, Python's (-bits) % 8
	default: // Bitmask
		s.nBytes = (bits + 7) / 8
		if m := bits % 8; m != 0 {
			s.mask = byte((1 << uint(m)) - 1)
		}
	}
	return s
}

// Next draws the next field element from the stream.
func (s *XOFSampler) Next() *big.Int {
	raw := make([]byte, s.nBytes)
	for {
		// SHAKE Read on an in-memory hash never errors and always fills raw.
		_, _ = s.h.Read(raw)
		var v *big.Int
		if s.bigEndian {
			raw[0] &= s.mask // big-endian: the first byte is most significant
			raw[0] >>= s.shift
			v = new(big.Int).SetBytes(raw)
		} else {
			raw[s.nBytes-1] &= s.mask // little-endian: the last byte is most significant
			raw[s.nBytes-1] >>= s.shift
			v = littleEndianInt(raw)
		}
		if s.reduce {
			return v.Mod(v, s.p)
		}
		if v.Cmp(s.p) < 0 { // rejection sampling
			return v
		}
	}
}

// NextNonzero draws the next field element, skipping zeros.
func (s *XOFSampler) NextNonzero() *big.Int {
	for {
		if v := s.Next(); v.Sign() != 0 {
			return v
		}
	}
}

// IsQuadraticNonResidue reports whether x is a quadratic non-residue mod p via
// Euler's criterion: x^((p-1)/2) == p-1. Both a residue (== 1) and x == 0 (== 0)
// return false, matching the reference legendre_symbol(x, p) == -1 test. It is
// the rejection predicate shared by the constructions that sample root-free
// quadratics (Arion's coeffs_g, Griffin's coeffs_G).
func IsQuadraticNonResidue(x, p *big.Int) bool {
	exp := new(big.Int).Rsh(new(big.Int).Sub(p, big.NewInt(1)), 1) // (p-1)/2
	ls := new(big.Int).Exp(x, exp, p)
	return ls.Cmp(new(big.Int).Sub(p, big.NewInt(1))) == 0
}
