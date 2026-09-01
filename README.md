## Benchmarked constructions

- GMiMC: `gmimc-bls12-t4`
- GMiMC2 with $\alpha=8$: `gmimc2-bls12-t4-a8`
- Neptune: `neptune-bls12-t4`
- Poseidon: `poseidon-bls12-t4`
- Poseidon2: `poseidon2-bls12-t4`
- Anemoi: `anemoi-bls12-381-scalar-t4`
- Rescue-Prime: `rescue-prime-bls12-t4`
- Arion: `arion-bls12-t4`
- Griffin: `griffin-bls12-t4`
- Skyscraper with 8-bit Bar words: `skyscraper-bls12-381-n2`
- Skyscraper with 16-bit Bar words: `skyscraper-bls12-381-n2-w16`
- Polocolo: `polocolo-bls12-381-scalar-t4`

All entries use the BLS12-381 scalar field and state width $t=4$. The benchmark plan measures the bare permutation, the construction's sponge mode, and, where defined, its compression mode.

## Prerequisites

- Go 1.26.2.
- Python 3.10 or newer; the table generator uses only the Python standard library.
- Network access for the initial Go module download.
- Sufficient memory and disk space for proving the full workloads. The largest default workloads are intentionally expensive.
- A LaTeX installation with the `booktabs` package to compile the generated tables.

## Commands

Run from the repository root:

```sh
go mod download

./tools/bench_table.py --workload permutation --tex-out permutation.tex
./tools/bench_table.py --workload many --tex-out many.tex
./tools/bench_table.py --workload merkle --tex-out merkle.tex
```
