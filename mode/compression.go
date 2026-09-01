package mode

import (
	"fmt"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// DaviesMeyer is the Davies-Meyer compression: trunc(perm(x) + x) to the first
// digest elements, where x is the concatenation of the two input blocks. Defined
// when the permutation width is exactly 2*digest.
type DaviesMeyer struct {
	perm   permutation.Permutation
	digest int
}

var _ Compression = (*DaviesMeyer)(nil)

// NewDaviesMeyer builds a Davies-Meyer compression producing digest elements. It
// requires perm.Width() == 2*digest.
func NewDaviesMeyer(perm permutation.Permutation, digest int) (*DaviesMeyer, error) {
	if digest < 1 {
		return nil, fmt.Errorf("mode: davies-meyer digest must be >= 1, got %d", digest)
	}
	if perm.Width() != 2*digest {
		return nil, fmt.Errorf("mode: davies-meyer requires width == 2*digest, got width=%d digest=%d", perm.Width(), digest)
	}
	return &DaviesMeyer{perm: perm, digest: digest}, nil
}

func (d *DaviesMeyer) Name() string   { return "davies-meyer" }
func (d *DaviesMeyer) Arity() int     { return 2 }
func (d *DaviesMeyer) BlockSize() int { return d.digest }

// Compress maps two digest-sized blocks (message, chaining value) to one.
func (d *DaviesMeyer) Compress(api frontend.API, blocks [][]frontend.Variable) []frontend.Variable {
	if len(blocks) != 2 {
		panic(fmt.Sprintf("mode: davies-meyer expects 2 blocks, got %d", len(blocks)))
	}
	x := append(clone(blocks[0]), blocks[1]...) // x_m || x_c, length 2*digest = width
	y := d.perm.Permute(api, x)
	out := make([]frontend.Variable, d.digest)
	for i := range out {
		out[i] = api.Add(y[i], x[i]) // feed-forward, left-truncated to digest
	}
	return out
}

// Jive is Anemoi's Jive_b compression (eprint 2022/840, Sec. 3.2): viewing the
// state as b blocks of m = width/b elements,
//
//	Jive_b(x_1..x_b)_i = sum_j (x_j[i] + perm(x_1||..||x_b)_j[i])
//
// i.e. a b-to-1 compression producing m elements. Defined when width % b == 0.
type Jive struct {
	perm permutation.Permutation
	b    int
	m    int
}

var _ Compression = (*Jive)(nil)

// NewJive builds a Jive_b compression. It requires perm.Width() % b == 0.
func NewJive(perm permutation.Permutation, b int) (*Jive, error) {
	if b < 2 {
		return nil, fmt.Errorf("mode: jive requires b >= 2, got %d", b)
	}
	if perm.Width()%b != 0 {
		return nil, fmt.Errorf("mode: jive requires width divisible by b, got width=%d b=%d", perm.Width(), b)
	}
	return &Jive{perm: perm, b: b, m: perm.Width() / b}, nil
}

func (j *Jive) Name() string   { return fmt.Sprintf("jive-%d", j.b) }
func (j *Jive) Arity() int     { return j.b }
func (j *Jive) BlockSize() int { return j.m }

// Compress maps b blocks of m elements to one block of m.
func (j *Jive) Compress(api frontend.API, blocks [][]frontend.Variable) []frontend.Variable {
	if len(blocks) != j.b {
		panic(fmt.Sprintf("mode: jive expects %d blocks, got %d", j.b, len(blocks)))
	}
	state := make([]frontend.Variable, 0, j.perm.Width())
	for _, blk := range blocks {
		state = append(state, blk...) // x_1 || ... || x_b
	}
	out := j.perm.Permute(api, state)

	res := make([]frontend.Variable, j.m)
	for i := 0; i < j.m; i++ {
		var acc frontend.Variable = 0
		for jj := 0; jj < j.b; jj++ {
			acc = api.Add(acc, state[i+j.m*jj], out[i+j.m*jj])
		}
		res[i] = acc
	}
	return res
}
