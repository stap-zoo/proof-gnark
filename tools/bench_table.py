#!/usr/bin/env python3
"""Turn ./cmd/bench results into the paper's LaTeX table.

Runs the Go command with -prove and formats the CSV it prints. It takes the *same
flags as the Go command* -- anything it does not recognise is forwarded verbatim --
so -size, -scale, -node, -instance, -runs and the label filters all work as
documented there:

    tools/bench_table.py                                  # one permutation call each
    tools/bench_table.py --workload merkle -scale -10     # 16 leaves, a quick run
    tools/bench_table.py --workload many -size 4,64       # a sweep: one table per size
    tools/bench_table.py --workload merkle -node sponge   # same mode for every construction

Which instances are measured stays in cmd/bench/plan.go.

For a long run, measure and render as two steps instead: `./cmd/bench -prove -out
r.csv ...` writes the CSV row by row, flushing each, and `--csv-in r.csv` formats it.
Driving the command from here keeps the CSV in a pipe, so a run that ends non-zero --
which one failed measurement out of many is enough to do -- takes the whole table with
it. Two steps cost a re-read and leave every measured row on disk.

Three workloads, three different questions, each its own table:

    permutation   one bare permutation call, in no mode -- the unit the cost table
                  compares (the Go command's -workload call -mode permutation)
    many          a sponge absorbing n field elements
    merkle        a 2:1 Merkle tree over n leaves

The four reported quantities, and what each is read from:

    constraints         R1CS constraints / SCS gates, the primary SoK metric
    prove               prove_ms + witness_ms -- proving plus input generation,
                        native evaluation and witness construction, all paid once
                        per proof
    verify              verify_ms  -- same
    circuit setup       compile_ms + setup_ms: compiling the constraint system plus
                        the circuit-specific preprocessing -- Groth16's whole setup,
                        PLONK's preprocessing and its Lagrange-basis SRS, but not
                        the universal one. Read it knowing that the Lagrange half
                        is a function of the *padded* size alone (harness/srs.go),
                        so two constructions whose counts round to the same power
                        of two pay the same for it however far apart they are

PLONK's universal KZG SRS is deliberately *not* among them, because it is one object
per curve rather than one per row. It is powers of tau, so the SRS for the largest
circuit contains every smaller one as a prefix, and gnark takes exactly that prefix
(backend/plonk/*/setup.go asserts only len(srs.Pk.G1) >= domain+3 and then slices
it); harness/srs.go keeps one per curve and grows it. A column of it would only say
how far the run had grown it by each row, so it is reported once on stderr when the
run finishes -- see [srs_summary].

Verification time and proof size do not depend on the construction -- they are fixed
per (curve, proof system) up to a wire or two of public input -- so that column is a
check that nothing is anomalous rather than a ranking.

LaTeX requirements: booktabs.
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import shutil
import subprocess
import sys
import textwrap
from dataclasses import dataclass
from typing import Callable, Sequence

# --------------------------------------------------------------------------
# Naming
# --------------------------------------------------------------------------

# The paper's construction macros, keyed by the CSV's `construction` field
# (catalog.Instance.Construction()). A construction missing here is an error
# rather than a \texttt{} fallback: a table that silently spells one name
# differently from the rest of the paper is worse than one that fails to build.
MACROS = {
    "gmimc": r"\gmimc",
    # NOTE: \gmimctwo is the assumed analogue of \gmimc -- the paper's excerpts here
    # (gmimc2.tex, theoretical_performance.tex) only ever use \gmimchashtwo and
    # \gmimcERFtwo, so this is the one macro in this table not observed in the source.
    # If the preamble spells it differently, this line is the only place to change.
    "gmimc2": r"\gmimctwo",
    "poseidon": r"\poseidon",
    "poseidon2": r"\poseidontwo",
    "neptune": r"\neptune",
    "griffin": r"\griffin",
    "anemoi": r"\anemoi",
    "arion": r"\arion",
    "rescueprime": r"\rescueprime",
    "polocolo": r"\polocolo",
    "skyscraper": r"\skyscraper",
}

FIELDS = {
    "bn254": "BN254",
    "bls12_381": "BLS12-381",
    "bls12-381": "BLS12-381",
}

# (row-group heading, what the constraint column counts) per proof system.
SYSTEMS = {
    "groth16": ("R1CS", "constraints"),
    "plonk": ("PLONK", "gates"),
}

# The width token in an instance name: `-t3`, `-n1`. It is used only to find where
# the *qualifier* starts -- what follows it in
# `polocolo-bls12-381-scalar-t3-tight` or `skyscraper-bn254-n1-w16` -- because that
# is what goes beside the construction's name. Matching on the token rather than
# stripping a prefix is what makes `rescue-prime-bls12-t3` work, whose construction
# field ("rescueprime") is not a prefix of its name.
#
# The width itself is never read from here. It is the CSV's `width` column, which is
# the number the permutation actually has (harness.Target.Width): Skyscraper's own
# label for its state of two elements is "n1", so a column built out of these tokens
# would put it in a different row from every other t=2 instance.
VARIANT = re.compile(r"^[tn]\d+$")


def row_label(construction: str, instance: str) -> str:
    """The macro for a construction, plus the qualifier that tells two instances of
    the same width apart: Polocolo's tighter parameter set, Skyscraper's 16-bit-word
    arithmetization of the same hash."""
    try:
        macro = MACROS[construction]
    except KeyError:
        raise SystemExit(
            f"bench_table: no LaTeX macro for construction {construction!r}.\n"
            f"             add it to MACROS in {os.path.basename(__file__)} "
            f"(known: {', '.join(sorted(MACROS))})"
        )
    toks = instance.split("-")
    for i, tok in enumerate(toks):
        if VARIANT.match(tok):
            rest = toks[i + 1 :]
            return macro + (f" ({latex_escape('-'.join(rest))})" if rest else "")
    return macro


def latex_escape(s: str) -> str:
    for a, b in (("\\", r"\textbackslash{}"), ("_", r"\_"), ("&", r"\&"),
                 ("%", r"\%"), ("#", r"\#"), ("$", r"\$")):
        s = s.replace(a, b)
    return s


def field_label(name: str) -> str:
    return FIELDS.get(name, latex_escape(name))


def join_words(words: Sequence[str]) -> str:
    """One word, or two joined by "or", or a list whose last pair is joined by it."""
    if len(words) < 3:
        return " or ".join(words)
    return ", ".join(words[:-1]) + " or " + words[-1]


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Metric:
    header: str          # column head, unit appended by the renderer
    group: str           # spanning head above it, "" for none
    unit: str            # "" for a count, "ms" for a duration
    read: Callable[[dict, argparse.Namespace], float | int]


def _f(row: dict, col: str) -> float:
    return float(row[col])


def _prove(row: dict, args) -> float:
    return (_f(row, "prove_min_ms" if args.stat == "min" else "prove_ms")
            + _f(row, "witness_ms"))


def _verify(row: dict, args) -> float:
    return _f(row, "verify_min_ms" if args.stat == "min" else "verify_ms")


METRICS = (
    Metric("constraints", "", "",
           lambda row, args: int(row["constraints"])),
    Metric("prove", "", "ms", _prove),
    Metric("verify", "", "ms", _verify),
    Metric("circuit setup", "", "ms",
           lambda row, args: _f(row, "compile_ms") + _f(row, "setup_ms")),
)

# --sort: the order of the rows inside one proof system's block, all ascending.
# Alphabetical goes by construction and then instance, so a construction's widths
# stay together; the other three go by the number the reader is ranking on.
SORTS = {
    "alph": lambda row, args: (row["construction"], row["instance"]),
    "constraint": lambda row, args: int(row["constraints"]),
    "prover": _prove,
    "verifier": _verify,
}


# --------------------------------------------------------------------------
# Number formatting
# --------------------------------------------------------------------------


def scale_for(values: Sequence[float]) -> tuple[str, float]:
    """Pick one unit for a whole column, so its cells stay comparable by eye.

    Late thresholds on purpose: these columns span two orders of magnitude
    (griffin against skyscraper), and a unit chosen off the top of the column
    rounds the bottom of it away. Milliseconds hold until the largest value is
    ten seconds, by which point four digits of milliseconds is the wrong
    presentation anyway.
    """
    m = max(values, default=0.0)
    if m >= 600_000:
        return "min", 1 / 60_000
    if m >= 10_000:
        return "s", 1 / 1_000
    return "ms", 1.0


def decimals_for(values: Sequence[float]) -> int:
    """Decimals from the column's largest value: enough to keep two significant
    figures near the bottom, and no false precision at the top -- a third decimal
    of a millisecond is nanoseconds, which nothing here resolves."""
    m = max(values, default=0.0)
    if m >= 1_000:
        return 0
    if m >= 100:
        return 1
    if m >= 1:
        return 2
    return 3


def group_int(n: int) -> str:
    """1372 -> 1\\,372. A thin space works in and out of math mode."""
    s = str(n)
    out = []
    while len(s) > 3:
        out.append(s[-3:])
        s = s[:-3]
    out.append(s)
    return r"\,".join(reversed(out))


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------


@dataclass
class Column:
    """One rendered column: its heads, its alignment, and its formatted cells.

    unit gets a head row of its own rather than riding along in `head`: "prove" over
    "(ms)" is narrower than "prove (ms)", which is most of what shrinks the table,
    and it lines the column names up with each other.
    """
    head: str
    group: str
    align: str
    cells: list[str]
    unit: str = ""


def render_table(rows: list[dict], args, systems: list[str], caption: str,
                 label: str, common_field: str | None,
                 common_width: str | None) -> str:
    """One LaTeX table over `rows`, all of which share a workload size."""
    ordered: list[tuple[str, list[dict]]] = []
    key = SORTS[args.sort]
    for s in systems:
        block = [r for r in rows if r["system"] == s]
        block.sort(key=lambda r: key(r, args))
        ordered.append((s, block))
    flat = [r for _, block in ordered for r in block]

    # Label columns: construction, plus width and field only when the table mixes
    # them. A constant width or field belongs in the caption rather than every row.
    # What distinguishes two instances of the same width sits beside the name (see
    # [row_label]).
    labels = [
        Column("construction", "", "l",
               [row_label(r["construction"], r["instance"]) for r in flat]),
    ]
    if common_width is None:
        labels.append(Column("$t$", "", "c", [r["width"] for r in flat]))
    if common_field is None:
        labels.append(Column("field", "", "l", [field_label(r["field"]) for r in flat]))
    dedup(labels, flat, [len(block) for _, block in ordered])

    metrics: list[Column] = []
    for m in METRICS:
        vals = [m.read(r, args) for r in flat]
        if m.unit == "ms":
            unit, k = scale_for([float(v) for v in vals])
            d = decimals_for([float(v) * k for v in vals])
            cells = [f"{float(v) * k:.{d}f}" for v in vals]
            metrics.append(Column(m.header, m.group, "r", cells, unit=f"({unit})"))
        else:
            cells = [group_int(int(v)) for v in vals]
            metrics.append(Column(m.header, m.group, "r", cells))

    cols = labels + metrics
    tabular = _tabular(cols, _body(ordered, cols))
    if not args.float:
        return tabular

    out = io.StringIO()
    out.write("\\begin{table}[t]\n  \\centering\n")
    out.write(f"  \\caption{{{caption}}}\n  \\label{{{label}}}\n")
    out.write(indent(tabular, "  "))
    out.write("\\end{table}\n")
    return out.getvalue()


def dedup(labels: list[Column], rows: list[dict], block_sizes: list[int]) -> None:
    """Disambiguate rows that no label column tells apart, in place.

    Two different instances can reduce to the same construction label if one
    carries its distinguishing part somewhere [VARIANT] does not look. The other
    label columns often resolve it -- the same instance over two curves is two
    rows with different fields -- so the check is on the whole label tuple, and
    within one proof system's block, since the same instance appearing in both
    blocks is the table's own structure rather than a collision.
    """
    at = 0
    for n in block_sizes:
        seen: dict[tuple[str, ...], list[int]] = {}
        for i in range(at, at + n):
            seen.setdefault(tuple(c.cells[i] for c in labels), []).append(i)
        for key, idx in seen.items():
            if len(idx) > 1 and len({rows[i]["instance"] for i in idx}) > 1:
                print(f"bench_table: {len(idx)} rows share the label "
                      f"{' / '.join(key)!r}; appending instance names",
                      file=sys.stderr)
                for i in idx:
                    labels[0].cells[i] += \
                        f" \\texttt{{{latex_escape(rows[i]['instance'])}}}"
        at += n


def _tabular(cols: list[Column], body: str) -> str:
    spec = "".join(c.align for c in cols)
    out = io.StringIO()
    out.write(f"\\begin{{tabular}}{{{spec}}}\n  \\toprule\n")
    for line in _heads(cols):
        out.write(f"  {line}\n")
    out.write("  \\midrule\n")
    out.write(body)
    out.write("  \\bottomrule\n\\end{tabular}\n")
    return out.getvalue()


def _heads(cols: list[Column]) -> list[str]:
    """The head: a spanning row for the metric groups (when any), the column names,
    then the units on a row of their own (when any column has one)."""
    lines = []
    if any(c.group for c in cols):
        cells, rules, i = [], [], 0
        while i < len(cols):
            g = cols[i].group
            j = i
            while j < len(cols) and cols[j].group == g:
                j += 1
            if g:
                cells.append(f"\\multicolumn{{{j - i}}}{{c}}{{{g}}}")
                rules.append(f"\\cmidrule(lr){{{i + 1}-{j}}}")
            else:
                cells.extend([""] * (j - i))
            i = j
        lines.append(" & ".join(cells) + r" \\ " + " ".join(rules))
    lines.append(" & ".join(c.head for c in cols) + r" \\")
    if any(c.unit for c in cols):
        lines.append(" & ".join(c.unit for c in cols) + r" \\")
    return lines


def _body(ordered: list[tuple[str, list[dict]]], cols: list[Column]) -> str:
    """Rows, in per-proof-system blocks separated by their own heading.

    One block per proof system rather than side-by-side columns: a Groth16
    number and a PLONK number are not two measurements of one thing, and a row
    group makes that harder to misread than adjacent columns would.
    """
    out, at = io.StringIO(), 0
    for k, (system, block) in enumerate(ordered):
        if k:
            out.write("  \\midrule\n")
        if len(ordered) > 1:
            head, _ = SYSTEMS.get(system, (system, "constraints"))
            out.write(f"  \\multicolumn{{{len(cols)}}}{{l}}"
                      f"{{\\textit{{{head}}}}} \\\\\n")
            out.write("  \\midrule\n")
        for _ in block:
            out.write("  " + " & ".join(c.cells[at] for c in cols) + " \\\\\n")
            at += 1
    return out.getvalue()


def indent(s: str, pad: str) -> str:
    return "".join(pad + line if line.strip() else line for line in s.splitlines(True))


# --------------------------------------------------------------------------
# The universal setup, reported once
# --------------------------------------------------------------------------


def next_pow2(n: int) -> int:
    p = 1
    while p < n:
        p <<= 1
    return p


def srs_summary(rows: list[dict]) -> list[str]:
    """PLONK's universal KZG SRS: one number per curve, for stderr at the end.

    A KZG SRS is powers of tau, so the one for the largest circuit *contains* every
    smaller one as a prefix, and that is literally how gnark uses it
    (`pk.Kzg.G1 = srs.Pk.G1[:vk.Size+3]`). The harness keeps one per curve and grows
    it (harness/srs.go), so a row's srs_ms is that curve's SRS as far as the run had
    grown it by then -- the same number for every row until one of them grows it,
    and the largest is the whole of it. The domain comes back from the row that
    needed the most: gnark's own rule, nextPow2(constraints + public wires).

    Not included: the Lagrange-basis SRS, which plonk.Setup requires at *exactly*
    the circuit's domain and so cannot be shared between sizes. That one is derived
    data rather than a ceremony, so it is per-circuit preprocessing and is inside
    the setup column.
    """
    per: dict[str, tuple[int, float]] = {}
    for r in rows:
        if r["system"] != "plonk":
            continue
        domain = next_pow2(int(r["constraints"]) + int(r["public_wires"]))
        seen = per.get(r["field"], (0, 0.0))
        per[r["field"]] = (max(seen[0], domain), max(seen[1], _f(r, "srs_ms")))
    if not per:
        return []

    width = max(len(field_label(f)) for f in per)
    lines = ["universal KZG SRS -- PLONK only, one per curve, so not a table column:"]
    for field in sorted(per):
        domain, srs = per[field]
        lines.append(f"  {field_label(field):<{width}}  "
                     f"2^{domain.bit_length() - 1} domain  {srs:8.0f} ms")
    note = ("Generated once per curve, at the largest domain the run needed, and "
            "shared by every circuit over that curve -- plonk.Setup takes a prefix. "
            "Not a per-row cost, and not part of the setup columns. The Lagrange "
            "basis, which is per domain and cannot be shared, is in setup/per "
            "circuit where it belongs.")
    return lines + textwrap.wrap(note, width=76, initial_indent="  ",
                                 subsequent_indent="  ")


# --------------------------------------------------------------------------
# Running the Go command
# --------------------------------------------------------------------------


def module_root() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if not os.path.exists(os.path.join(root, "go.mod")):
        raise SystemExit(f"bench_table: no go.mod in {root} "
                         f"(these scripts live in <module>/tools/)")
    return root


def go_binary() -> str:
    go = shutil.which("go")
    if go:
        return go
    fallback = "/usr/local/go/bin/go"
    if os.path.exists(fallback):
        return fallback
    raise SystemExit("bench_table: no `go` on PATH (try /usr/local/go/bin)")


def has_flag(flags: Sequence[str], name: str) -> bool:
    """Whether `-name` or `--name` was forwarded, in either `-f v` or `-f=v` form."""
    return any(a in (f"-{name}", f"--{name}")
               or a.startswith((f"-{name}=", f"--{name}=")) for a in flags)


def run_bench(flags: list[str], args) -> str:
    """Run the Go command with -prove and return the CSV it printed.

    -prove is forced because every column but the constraint count needs it, and
    counting is the Go command's default. The command's stderr (the [k/n]
    progress) is left alone so a long run still reports itself.
    """
    flags = list(flags)
    if not has_flag(flags, "prove"):
        flags.append("-prove")
    for a in flags:
        if a in ("-prove=false", "--prove=false", "-prove=0"):
            raise SystemExit("bench_table: -prove=false leaves nothing to time")
    if has_flag(flags, "out"):
        raise SystemExit("bench_table: -out sends the CSV to a file instead of to this "
                         "script; run ./cmd/bench directly if that is what you want")

    cmd = [go_binary(), "run", "./cmd/bench"] + flags
    print("bench_table: " + " ".join(cmd), file=sys.stderr)
    if args.dry_run:
        return ""
    proc = subprocess.run(cmd, cwd=module_root(), stdout=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"bench_table: ./cmd/bench failed "
                         f"(exit {proc.returncode}); the CSV is incomplete")
    return proc.stdout


# --------------------------------------------------------------------------
# Front end
# --------------------------------------------------------------------------


# The three workloads: the flags each one runs the Go command with, and the
# \label{} stem it writes. `permutation` pins -mode as well as -workload, because
# the Go command's `call` workload measures every mode an instance declares and this
# table wants the bare permutation, which is the one target every construction has.
WORKLOADS = {
    "permutation": dict(flags=("-workload", "call", "-mode", "permutation"),
                        stem="cost"),
    "many": dict(flags=("-workload", "message"), stem="longdata"),
    "merkle": dict(flags=("-workload", "merkle"), stem="merkle"),
}


def caption_for(workload: str, common_field: str | None,
                common_width: str | None, kinds: list[str], row: dict) -> str:
    """The generated caption: what was measured, in which type of mode, and the
    field and width when the whole table shares them (mixed tables retain columns).

    The mode is named by its *type* and never by the individual modes: which mode a
    construction uses is the construction's own business -- Anemoi's Jive-2,
    Skyscraper's Davies-Meyer, Poseidon's plain sponge -- and a caption enumerating
    ten of them is noise. What a reader needs is the type, plus the fact that every
    construction used a mode its designers defined; the mode per row stays in the
    CSV.

    The type is the CSV's `kind`, printed verbatim, so these captions say
    "compression" exactly as harness.Kind, the `-kind` filter and the reference's own
    KAT labels do. Nothing here renames it: "truncation" would be wrong for Jive,
    which sums branches rather than dropping any of them, and a second vocabulary for
    the same three types is worth less than a shared one anyway.
    """
    over = f", over {field_label(common_field)}" if common_field else ""
    width = f", with $t={latex_escape(common_width)}$" if common_width else ""
    types = join_words([latex_escape(k) for k in kinds])
    a_types = join_words([f"a {latex_escape(k)}" for k in kinds])
    n = int(row["size"])
    if workload == "many":
        return (f"Cost of hashing {group_int(n)} field element"
                f"{'s' if n != 1 else ''} ({human_bytes(n * 32)}) as one absorption "
                f"inside a single circuit{over}{width}. Each construction uses the {types} "
                f"mode its designers define.")
    if workload == "merkle":
        return (f"Cost of proving a 2:1 Merkle root over {group_int(n)} leaves "
                f"({group_int(n - 1)} node {'hashes' if n != 2 else 'hash'}) inside "
                f"a single circuit{over}{width}. Each node is {a_types}, in the mode its "
                f"designers define.")
    return (f"Cost of a single {types} call inside a circuit{over}{width}.")


def parse_args(argv: list[str]) -> tuple[argparse.Namespace, list[str]]:
    p = argparse.ArgumentParser(
        prog="bench_table.py", allow_abbrev=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Run ./cmd/bench with -prove and write a LaTeX table of its "
                    "results.\n\n"
                    "Any flag not listed here is forwarded verbatim to ./cmd/bench, "
                    "so this script takes the same parameters it does (-size, "
                    "-scale, -node, -backend, -runs, -instance, -all, -seed, ...). "
                    "Use -- before them if one collides.",
        epilog="examples:\n"
               "  tools/bench_table.py --tex-out cost.tex\n"
               "  tools/bench_table.py --workload merkle -scale -10 -runs 10\n"
               "  tools/bench_table.py --workload many --sort prover\n")
    p.add_argument("--workload", choices=tuple(WORKLOADS), default="permutation",
                   help="which circuit to measure and caption: a bare permutation "
                        "call, a sponge over many elements, or a 2:1 Merkle tree "
                        "(default: permutation)")
    p.add_argument("--tex-out", metavar="PATH",
                   help="write the LaTeX here (default: stdout)")
    p.add_argument("--systems", default="groth16,plonk",
                   help="proof systems to report, in order (default: groth16,plonk)")
    p.add_argument("--only", metavar="N", dest="size",
                   help="render only these workload sizes (comma-separated); "
                        "default: one table per size in the CSV. Not to be confused "
                        "with -size, which is forwarded and decides what is measured")
    p.add_argument("--stat", choices=("mean", "min"), default="mean",
                   help="prove/verify: the mean over -runs, or the fastest run")
    p.add_argument("--sort", choices=tuple(SORTS), default="constraint",
                   help="row order within a proof system: alphabetically by "
                        "construction, or ascending by constraint count, prove time "
                        "or verify time (default: constraint)")
    p.add_argument("--no-float", dest="float", action="store_false",
                   help="emit the tabular only, without table/caption/label")
    p.add_argument("--csv-in", metavar="PATH", dest="csv_in",
                   help="render this CSV instead of running ./cmd/bench -- for a long "
                        "run measured separately with `./cmd/bench -prove -out PATH`, "
                        "whose rows are flushed as they are produced and so survive a "
                        "failure that would make the command exit non-zero and lose "
                        "the whole table here")
    p.add_argument("--dry-run", action="store_true",
                   help="print the benchmark command and exit")
    args, extras = p.parse_known_args(argv)
    args.systems = [s.strip() for s in args.systems.split(",") if s.strip()]
    for s in args.systems:
        if s not in SYSTEMS:
            raise SystemExit(f"bench_table: unknown proof system {s!r} "
                             f"(want {', '.join(SYSTEMS)})")
    if args.csv_in and extras:
        raise SystemExit(f"bench_table: --csv-in renders an existing CSV, so it "
                         f"measures nothing and {' '.join(extras)} would be silently "
                         f"ignored -- pass those to ./cmd/bench when producing the CSV")
    if args.workload == "permutation" and has_flag(extras, "mode"):
        raise SystemExit("bench_table: --workload permutation already pins -mode "
                         "permutation; a forwarded -mode would be overridden and the "
                         "caption would then describe the wrong measurement")
    return args, extras


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    args, extras = parse_args(argv)
    if args.csv_in:
        # --workload still matters: it picks the caption and the \label stem. Nothing
        # here checks that the CSV was measured with that workload's flags -- the
        # `workload` column says what was measured and render() reads the size from
        # it, but the caption is the caller's word.
        if args.dry_run:
            print(f"bench_table: would render {args.csv_in}", file=sys.stderr)
            return
        with open(args.csv_in) as f:
            raw = f.read()
    else:
        # The workload is this script's flag and the command's, so forward it -- last,
        # so it wins over anything the caller forwarded.
        raw = run_bench(extras + list(WORKLOADS[args.workload]["flags"]), args)
        if args.dry_run:
            return

    rows = list(csv.DictReader(io.StringIO(raw)))
    if not rows:
        raise SystemExit("bench_table: no rows in the benchmark output")
    missing = [c for c in ("public_wires", "srs_ms", "setup_ms", "witness_ms",
                           "prove_ms", "verify_ms") if c not in rows[0]]
    if missing:
        raise SystemExit(f"bench_table: the CSV has no {', '.join(missing)} -- it was "
                         f"produced without -prove, which is where the timings come "
                         f"from")
    if "width" not in rows[0]:
        raise SystemExit("bench_table: the CSV has no 'width' column, which is the "
                         "state size t -- it predates harness.Target.Width, so "
                         "re-run ./cmd/bench")

    # A sweep is one table per size: rows measured at different sizes do not belong
    # in the same column.
    tables = [render(args, group, size) for size, group in sizes(args, rows)]
    text = "\n".join(tables)
    if args.tex_out:
        with open(args.tex_out, "w") as f:
            f.write(text)
        print(f"bench_table: wrote {args.tex_out} "
              f"({len(tables)} table{'s' if len(tables) != 1 else ''})",
              file=sys.stderr)
    else:
        sys.stdout.write(text)

    for line in srs_summary(rows):
        print("bench_table: " + line, file=sys.stderr)


def sizes(args, rows: list[dict]) -> list[tuple[str, list[dict]]]:
    """Split the CSV into one table per workload size, smallest first."""
    if "size" not in rows[0]:
        raise SystemExit("bench_table: the CSV has no 'size' column "
                         "(is it from ./cmd/bench?)")
    wanted = None
    if args.size:
        wanted = {s.strip() for s in args.size.split(",") if s.strip()}
        unknown = wanted - {r["size"] for r in rows}
        if unknown:
            raise SystemExit(f"bench_table: no rows of size "
                             f"{', '.join(sorted(unknown))}")
    groups: dict[str, list[dict]] = {}
    for r in rows:
        if wanted is None or r["size"] in wanted:
            groups.setdefault(r["size"], []).append(r)
    return sorted(groups.items(), key=lambda kv: int(kv[0]))


def render(args, rows: list[dict], size: str) -> str:
    """Caption, label and field column for one size, then the table."""
    systems = [s for s in args.systems if any(r["system"] == s for r in rows)]
    if not systems:
        raise SystemExit(f"bench_table: no rows for proof system(s) "
                         f"{', '.join(args.systems)}")
    fields = {r["field"] for r in rows}
    common_field = rows[0]["field"] if len(fields) == 1 else None
    widths = {r["width"] for r in rows}
    common_width = rows[0]["width"] if len(widths) == 1 else None
    kinds = sorted({r["kind"] for r in rows})

    caption = caption_for(args.workload, common_field, common_width, kinds, rows[0])
    # With both proof systems the row-group headings say which is which; with one
    # there is no heading, so the caption has to carry it -- a column of times whose
    # proof system is nowhere stated is not a reportable number.
    if len(systems) == 1:
        head, counts = SYSTEMS[systems[0]]
        caption += f" {head}, so the constraint column counts {counts}."
    runs = rows[0].get("runs")
    if runs:
        caption += (f" Prove and verify are the "
                    f"{'fastest of' if args.stat == 'min' else 'mean over'}"
                    f" {runs} run{'s' if runs != '1' else ''}.")
    label = f"tab:{WORKLOADS[args.workload]['stem']}-{size}"
    return render_table(rows, args, systems, caption, label, common_field,
                        common_width)


def human_bytes(n: int) -> str:
    for unit, cut in (("GiB", 1 << 30), ("MiB", 1 << 20), ("KiB", 1 << 10)):
        if n >= cut:
            v = n / cut
            return f"{v:.0f}~{unit}" if v >= 10 or v == int(v) else f"{v:.1f}~{unit}"
    return f"{n}~B"


if __name__ == "__main__":
    main()
