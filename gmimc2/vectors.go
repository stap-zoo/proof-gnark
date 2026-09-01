package gmimc2

import (
	"embed"

	"github.com/zkhash-sok/gnark-hashes/harness"
)

//go:embed testdata/*.json
var testdata embed.FS

// Vectors returns this construction's reference vectors (see harness.EmbeddedVectors).
// Unlike every other construction's, they do not come from the Python reference
// framework — ../ref has no GMiMC2 — but from gnark-hashes/gmimc2_ref.py, this
// repo's standalone oracle for it.
func Vectors() []harness.Vector { return harness.EmbeddedVectors("gmimc2", testdata) }
