package algebra

import "github.com/consensys/gnark/frontend"

// kvStore is gnark's per-compilation key-value store, declared structurally
// because gnark keeps the interface in an internal package. Every gnark builder
// implements it (and so does the test engine), which is how a gadget shares one
// instance of itself across a whole circuit — see std/rangecheck, which caches
// its checker exactly this way.
type kvStore interface {
	SetKeyValue(key, value any)
	GetKeyValue(key any) any
}

// Shared returns the per-compilation instance stored under key, calling build
// only the first time. It is how a gadget with *setup cost* — a lookup table, a
// range checker — is built once per circuit instead of once per use.
//
// This is what makes a lookup argument amortize. A logderivlookup table costs one
// constraint per entry, paid when the table is created, plus a small cost per
// query. Building it inside a round function would charge that table cost again
// for every hash call in the circuit, so a Merkle tree over 2^14 leaves would pay
// for 2^14 identical tables; with one shared table the entries are paid once and
// all 2^14 nodes' queries land in the same argument. Sharing is also what a real
// circuit would do, and it is sound for the same reason a single table is: the
// table is a fixed constant, and every query joins one lookup argument over it.
//
// key must be comparable and must identify the *configuration*, not just the
// gadget type: two instances with different tables (a different field, a
// different word size) must not collide. Use a small unexported struct type
// holding whatever distinguishes them.
//
// If the builder does not implement the store, build is called every time. That
// only costs gates — nothing becomes incorrect — so callers need no fallback.
func Shared[T any](api frontend.API, key any, build func() T) T {
	kv, ok := api.Compiler().(kvStore)
	if !ok {
		return build()
	}
	if got := kv.GetKeyValue(key); got != nil {
		if t, ok := got.(T); ok {
			return t
		}
	}
	t := build()
	kv.SetKeyValue(key, t)
	return t
}
