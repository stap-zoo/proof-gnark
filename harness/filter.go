package harness

// Filter selects a subset of targets by their labels. Each field is a set of
// allowed values; an empty field means "no constraint on this dimension". So the
// zero Filter{} matches everything, Filter{Kinds: {Sponge}} keeps all sponges,
// and Filter{Fields: {"bn254"}, Modes: {"jive-2"}} keeps Jive on bn254.
type Filter struct {
	Constructions []string
	Fields        []string
	Modes         []string
	Kinds         []Kind
}

// Select returns the targets that match the filter, preserving order.
func (f Filter) Select(targets []Target) []Target {
	out := make([]Target, 0, len(targets))
	for _, t := range targets {
		if allowsString(f.Constructions, t.Construction) &&
			allowsString(f.Fields, t.FieldName) &&
			allowsString(f.Modes, t.Mode) &&
			allowsKind(f.Kinds, t.Kind) {
			out = append(out, t)
		}
	}
	return out
}

func allowsString(allowed []string, v string) bool {
	if len(allowed) == 0 {
		return true
	}
	for _, a := range allowed {
		if a == v {
			return true
		}
	}
	return false
}

func allowsKind(allowed []Kind, v Kind) bool {
	if len(allowed) == 0 {
		return true
	}
	for _, a := range allowed {
		if a == v {
			return true
		}
	}
	return false
}
