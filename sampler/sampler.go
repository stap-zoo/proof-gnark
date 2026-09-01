// Package sampler is the deterministic field-element sampler shared by the
// parameter-derivation code of every construction, mirroring
// ref/utils/sampler.py (the XOFFieldElementSampler with SHAKE256 and "mod"
// sampling).
//
// It is native, out-of-circuit code (plain big.Int over a seeded XOF), not an
// in-circuit gadget — which is why it lives in its own package rather than in
// algebra/. Round constants (GMiMC, Rescue-Prime, Arion), Arion's linear-map
// coefficients, and Arion's rejection-sampled g-coefficients are all drawn from
// this one primitive, so it is implemented here once.
package sampler

import (
	"crypto/sha256"
	"crypto/sha3"
	"math/big"
)

// SHA256Mod draws one field element from a single SHA-256 digest of seed: the
// 32-byte digest is read big-endian and reduced mod field. This mirrors the
// reference XOFFieldElementSampler(xof="sha256", sampling="mod", n_bytes=32,
// endianess="big").next() — SHA-256 is a bounded 32-byte source, so "mod"
// sampling reads the whole digest and reduces once. Skyscraper seeds a fresh
// digest per round constant (index || label), so a single-element draw is the
// natural primitive.
func SHA256Mod(seed []byte, field *big.Int) *big.Int {
	sum := sha256.Sum256(seed)
	v := new(big.Int).SetBytes(sum[:]) // big-endian
	return v.Mod(v, field)
}

// SHAKE256Mod draws count field elements from SHAKE256(seed) using the "mod"
// strategy of the reference sampler: read ceil(bitlen(p)/8)+1 little-endian
// bytes per element and reduce mod p (no rejection). Because SHAKE is an XOF
// (its longer output extends the shorter one) and "mod" never rejects, reading
// the whole count*nBytes stream up front is byte-identical to the reference's
// grow-and-buffer loop that advances nBytes at a time.
func SHAKE256Mod(seed []byte, field *big.Int, count int) []*big.Int {
	nBytes := (field.BitLen()+7)/8 + 1 // "mod" sampling width

	stream := make([]byte, nBytes*count)
	h := sha3.NewSHAKE256()
	// SHAKE Write/Read on an in-memory hash never error; ignore per the stdlib
	// hash.Hash contract.
	_, _ = h.Write(seed)
	_, _ = h.Read(stream)

	out := make([]*big.Int, count)
	for k := 0; k < count; k++ {
		v := littleEndianInt(stream[k*nBytes : (k+1)*nBytes])
		out[k] = v.Mod(v, field)
	}
	return out
}

// Grid returns a rows x cols grid (row-major) of field elements, matching the
// reference sampler's grid(num_rows, num_cols): it is SHAKE256Mod over
// rows*cols elements reshaped, since grid draws column-by-column within each row
// off the same sequential stream.
func Grid(seed []byte, field *big.Int, rows, cols int) [][]*big.Int {
	flat := SHAKE256Mod(seed, field, rows*cols)
	out := make([][]*big.Int, rows)
	for r := 0; r < rows; r++ {
		out[r] = flat[r*cols : (r+1)*cols]
	}
	return out
}

// littleEndianInt interprets b as a little-endian unsigned integer (big.Int
// reads big-endian, so reverse first).
func littleEndianInt(b []byte) *big.Int {
	be := make([]byte, len(b))
	for i, x := range b {
		be[len(b)-1-i] = x
	}
	return new(big.Int).SetBytes(be)
}
