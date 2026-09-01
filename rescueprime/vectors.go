package rescueprime

import (
	"embed"

	"github.com/zkhash-sok/gnark-hashes/harness"
)

//go:embed testdata/*.json
var testdata embed.FS

// Vectors returns this construction's reference vectors (see harness.EmbeddedVectors).
func Vectors() []harness.Vector { return harness.EmbeddedVectors("rescueprime", testdata) }
