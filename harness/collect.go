package harness

import (
	"fmt"
	"math/big"
	"strings"

	"github.com/consensys/gnark/frontend"
	"github.com/zkhash-sok/gnark-hashes/mode"
	"github.com/zkhash-sok/gnark-hashes/permutation"
)

// Instance is the metadata every construction's instance exposes so the harness
// can label its targets. Every instance also has a permutation (so it always
// yields a permutation target); sponge/compression targets are added when the
// instance additionally implements mode.SpongeHash / mode.Compressor.
type Instance interface {
	Construction() string
	InstanceName() string
	FieldName() string
	Field() *big.Int
	Permutation() permutation.Permutation
}

// VectorSource is the optional interface for an instance whose reference vectors
// are labeled with a *different* instance name: an alternative arithmetization of
// a hash already registered, which computes the identical function and so
// satisfies the original instance's known-answer vectors.
//
// Skyscraper's 16-bit-word Bar is the case it exists for. It is byte-identical to
// the byte-word Bar, so it must satisfy every byte-word vector — but it needs an
// instance name of its own (…-n1-w16) to be a distinct row in the registry, in a
// CSV, and in a plan. Without this the two would have to share a name, which is
// exactly the ambiguity that made a plan entry unresolvable.
//
// An instance not implementing it has VectorInstance == InstanceName, which is
// every construction here bar that one.
type VectorSource interface {
	VectorInstance() string
}

// vectorInstance is the name a target's vectors are labeled with: the instance's
// own name unless it declares another via VectorSource.
func vectorInstance(inst Instance) string {
	if vs, ok := inst.(VectorSource); ok {
		return vs.VectorInstance()
	}
	return inst.InstanceName()
}

// Pick returns the named instance from a construction's own registry — the way a
// benchmark plan names one instance out of a construction's set:
//
//	Pick(skyscraper.Instances16, "skyscraper-bn254-n1-w16")
//
// The registry is a typed reference, so the *construction and variant* a plan
// entry means are a Go symbol (jump to it, and its parameters are right there)
// rather than something to be recovered from a name; only the choice within that
// registry is a string, and it is resolved once against the registry that must
// contain it.
//
// A name absent from the given registry panics rather than returning an error:
// every caller is package-level plan data, so a miss is a programming error that
// should surface at init, not a silently shorter table later.
func Pick[T Instance](insts []T, name string) Instance {
	for _, inst := range insts {
		if inst.InstanceName() == name {
			return inst
		}
	}
	have := make([]string, len(insts))
	for i, inst := range insts {
		have[i] = inst.InstanceName()
	}
	panic(fmt.Sprintf("harness: no instance %q in this registry (it has: %s)",
		name, strings.Join(have, ", ")))
}

// TargetsOf builds the targets of a construction's whole registry:
// TargetsOf(gmimc.Instances). It exists because a registry is a
// []<construction>.Instance and Targets takes []Instance (the interface), a
// conversion Go does not do implicitly — every construction used to carry the same
// five-line loop to do it.
func TargetsOf[T Instance](insts []T) []Target {
	return Targets(Lift(insts)...)
}

// Lift converts a construction's typed registry to []Instance. Needed wherever a
// registry is handed to something generic over instances (see TargetsOf, and
// compare's aggregation).
func Lift[T Instance](insts []T) []Instance {
	out := make([]Instance, len(insts))
	for i, inst := range insts {
		out[i] = inst
	}
	return out
}

// Targets builds the targets for the given instances: one permutation target
// each, plus one per sponge and per compression the instance declares.
func Targets(instances ...Instance) []Target {
	var out []Target
	for _, inst := range instances {
		out = append(out, permutationTarget(inst))
		if sh, ok := inst.(mode.SpongeHash); ok {
			for _, s := range sh.Sponges() {
				out = append(out, spongeTarget(inst, s))
			}
		}
		if c, ok := inst.(mode.Compressor); ok {
			for _, comp := range c.Compressions() {
				out = append(out, compressionTarget(inst, comp))
			}
		}
	}
	return out
}

func permutationTarget(inst Instance) Target {
	return Target{
		Construction:   inst.Construction(),
		Instance:       inst.InstanceName(),
		FieldName:      inst.FieldName(),
		Kind:           Permutation,
		Mode:           "permutation",
		Width:          inst.Permutation().Width(),
		Case:           PermutationCase(inst.InstanceName()+"/permutation", inst.Field(), inst.Permutation()),
		VectorInstance: vectorInstance(inst),
	}
}

// PermutationCase is the Case for a bare permutation circuit. It is what
// permutationTarget measures, exported so a construction can also cost a
// permutation it builds itself — e.g. the same instance with its linear-layer
// addition programs disabled, which is how the PLONK saving from each program is
// measured in the per-construction tests.
func PermutationCase(name string, field *big.Int, perm permutation.Permutation) Case {
	w := perm.Width()
	return Case{
		Name:       name,
		Field:      field,
		InputSize:  w,
		OutputSize: w,
		Build: func(in, out []frontend.Variable) frontend.Circuit {
			c := permutation.NewPermutationCircuit(perm)
			copy(c.Input, in)
			copy(c.Output, out)
			return c
		},
		Eval: perm.Permute,
	}
}

func spongeTarget(inst Instance, s mode.Sponge) Target {
	return Target{
		Construction:   inst.Construction(),
		Instance:       inst.InstanceName(),
		FieldName:      inst.FieldName(),
		Kind:           Sponge,
		Mode:           s.Name(),
		Width:          inst.Permutation().Width(),
		VectorInstance: vectorInstance(inst),
		Case: Case{
			Name:  inst.InstanceName() + "/" + s.Name(),
			Field: inst.Field(),
			// Canonical benchmark workload: one rate-sized absorption block.
			InputSize:  s.Rate(),
			OutputSize: s.DigestSize(),
			Build: func(in, out []frontend.Variable) frontend.Circuit {
				c := mode.NewSpongeCircuit(s, len(in))
				copy(c.Input, in)
				copy(c.Output, out)
				return c
			},
			Eval: s.Hash,
		},
	}
}

func compressionTarget(inst Instance, comp mode.Compression) Target {
	return Target{
		Construction:   inst.Construction(),
		Instance:       inst.InstanceName(),
		FieldName:      inst.FieldName(),
		Kind:           Compression,
		Mode:           comp.Name(),
		Width:          inst.Permutation().Width(),
		VectorInstance: vectorInstance(inst),
		Case: Case{
			Name:       inst.InstanceName() + "/" + comp.Name(),
			Field:      inst.Field(),
			InputSize:  comp.Arity() * comp.BlockSize(),
			OutputSize: comp.BlockSize(),
			Build: func(in, out []frontend.Variable) frontend.Circuit {
				c := mode.NewCompressionCircuit(comp)
				copy(c.Input, in)
				copy(c.Output, out)
				return c
			},
			// Compress takes the input blocks separately; the Case models input
			// as one flat slice (as CompressionCircuit.Define also does), so the
			// split is redone here.
			Eval: func(api frontend.API, in []frontend.Variable) []frontend.Variable {
				m := comp.BlockSize()
				blocks := make([][]frontend.Variable, comp.Arity())
				for j := range blocks {
					blocks[j] = in[j*m : (j+1)*m]
				}
				return comp.Compress(api, blocks)
			},
		},
	}
}
