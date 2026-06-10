# CodonForge

> Research preview: clean-room Rust CLI for protein-faithful mRNA coding-sequence optimization. Not suitable for production mRNA design or clinical use.

[![Rust](https://img.shields.io/badge/rust-1.96%2B-blue)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/raghuboi/codonforge/actions/workflows/ci.yml/badge.svg)](https://github.com/raghuboi/codonforge/actions/workflows/ci.yml)

CodonForge reads a protein sequence and generates a synonymous mRNA sequence. It supports:

- `greedy`: v0 highest-frequency human codon baseline
- `beam`: v1 deterministic beam search balancing CAI, GC%, U%, and adjacent-repeat penalty

It computes Codon Adaptation Index (CAI), reports simple composition metrics, and writes mRNA FASTA for downstream benchmarking.

CodonForge is intentionally small. It is useful as a deterministic baseline, benchmark integration point, and reproducible research tool, not as a complete mRNA design system.

## What this is not

CodonForge does not implement:

- LinearDesign algorithm parity
- O(N³) codon-constrained folding dynamic programming
- Turner energy secondary-structure scoring
- internal MFE prediction
- GPU acceleration
- UTR design
- clinical or therapeutic recommendations

Minimum free energy (MFE), structure-aware metrics, and downstream biological scores should be computed by external benchmark tools such as ViennaRNA-backed evaluators.

## Install

### From source

```bash
git clone https://github.com/raghuboi/codonforge.git
cd codonforge
cargo build --release
```

The binary will be at:

```bash
target/release/codonforge
```

### From crates.io

CodonForge is not published to crates.io yet. Until the first crates.io release is published, use the source checkout path above.

## Quickstart

From a source checkout:

```bash
echo "MNDTEAI" | cargo run --release --
```

Expected output shape:

```text
Input protein: MNDTEAI
mRNA sequence:  AUGAACGACACCGAGGCCAUC
mRNA structure: .....................
mRNA folding free energy: nan kcal/mol; mRNA CAI: 1.000
CodonForge strategy: greedy; GC%: 57.14; U%: 9.52; repeat penalty: 4.0
Runtime: 0.0000 seconds
```

Runtime varies by machine.

## Usage

### Protein sequence from stdin

```bash
echo "MNDTEAI" | codonforge
```

From a source checkout, use:

```bash
echo "MNDTEAI" | cargo run --release --
```

### Protein FASTA input

```bash
codonforge --input path/to/protein.fasta
```

### Write pure mRNA FASTA

```bash
codonforge \
  --input path/to/protein.fasta \
  --fasta-out path/to/output.fasta
```

### Beam search

```bash
codonforge \
  --input path/to/protein.fasta \
  --strategy beam \
  --beam-width 32 \
  --target-gc 55 \
  --target-u 20
```

Beam-search objective:

```text
score =
  + weight_cai * mean_log_cai_score
  - weight_gc * abs(GC% - target_gc)
  - weight_u * abs(U% - target_u)
  - weight_repeat * adjacent_repeat_penalty
```

The default weights are conservative and CPU-only. Increase GC/U weights when exploring composition realism; keep greedy as the baseline for comparisons.

### Custom codon usage table

```bash
codonforge \
  --input path/to/protein.fasta \
  --codonusage path/to/codon_table.csv
```

Codon usage tables must contain three columns:

```text
codon,amino_acid,frequency
```

Codons may use RNA `U` or DNA `T` notation. DNA notation is normalized to RNA.

## CLI flags

| Flag | Description |
| --- | --- |
| `-i, --input <FILE>` | Input protein FASTA file path. If omitted, CodonForge reads a plain protein sequence from stdin. |
| `-f, --fasta-out <FILE>` | Output pure mRNA FASTA file path. |
| `-c, --codonusage <FILE>` | Codon usage table CSV path. Default: `data/codon_usage_freq_table_human.csv`. |
| `--strategy <greedy|beam>` | Optimization strategy. Default: `greedy`. |
| `--beam-width <N>` | Beam width for `--strategy beam`. Default: `32`. |
| `--target-gc <PERCENT>` | Target GC percentage for beam search. Default: `55`. |
| `--target-u <PERCENT>` | Target U percentage for beam search. Default: `20`. |
| `--weight-cai <FLOAT>` | CAI objective weight for beam search. Default: `1.0`. |
| `--weight-gc <FLOAT>` | GC target penalty weight. Default: `0.02`. |
| `--weight-u <FLOAT>` | U target penalty weight. Default: `0.01`. |
| `--weight-repeat <FLOAT>` | Adjacent-repeat penalty weight. Default: `0.05`. |
| `-l, --lambda <FLOAT>` | Deprecated LinearDesign-style structure parameter; accepted for compatibility and ignored with a warning. |
| `-v, --verbose` | Print diagnostic information to stderr. |
| `-h, --help` | Print help. |

## Examples

### Toy example

```bash
cargo build --release
./examples/toy/run.sh
```

### RP mini-benchmark

```bash
cargo build --release
./examples/benchmark/run.sh
```

The mini-benchmark runs RHO, PRPF31, RPE65, and RPGR-ORF15 with both `greedy` and `beam` strategies. It is included for reproducible computational experiments only.

## Data provenance

The bundled human codon usage table is derived from the Kazusa Codon Usage Database Homo sapiens record:

- URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
- Organism: Homo sapiens `[gbpri]`
- Accessed: 2026-06-10

Frequencies in `data/codon_usage_freq_table_human.csv` are normalized from Kazusa raw counts within each amino-acid family. See `data/README.md` for conversion notes.

## Development

Run the full local gate:

```bash
cargo fmt --check
cargo test
cargo clippy --all-targets --all-features -- -D warnings
cargo doc --no-deps
cargo publish --dry-run
```

## Research use and citation

If you use CodonForge in research, cite the repository and include the exact git commit, CodonForge version, codon usage table source, and benchmark scripts used.

Citation metadata is provided in `CITATION.cff`.

## Limitations

- CAI depends on the selected codon usage table; compare tools only when using the same table.
- Beam search is a simple multi-objective heuristic, not a folding-aware optimizer.
- `mRNA folding free energy: nan` is a compatibility placeholder; CodonForge does not compute MFE.
- Composition metrics are useful for screening and benchmarking, not evidence of biological efficacy.

## License

MIT. See `LICENSE`.
