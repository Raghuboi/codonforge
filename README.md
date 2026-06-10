# CodonForge

> v0 research preview: greedy CAI optimizer only. Not suitable for production mRNA design or clinical use.

[![Rust](https://img.shields.io/badge/rust-1.96%2B-blue)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/raghuboi/codonforge/actions/workflows/ci.yml/badge.svg)](https://github.com/raghuboi/codonforge/actions/workflows/ci.yml)

CodonForge is a clean-room Rust CLI for protein-faithful mRNA coding-sequence optimization research. It reads a protein sequence, chooses synonymous human codons with a greedy highest-frequency strategy, computes Codon Adaptation Index (CAI), and writes mRNA FASTA for downstream benchmarking.

CodonForge v0 is intentionally small. It is useful as a deterministic baseline and benchmark integration point, not as a complete mRNA design system.

## What this is not

CodonForge v0 does not implement:

- LinearDesign algorithm parity
- O(N³) codon-constrained folding dynamic programming
- Turner energy secondary-structure scoring
- internal MFE prediction
- GPU acceleration
- UTR design
- clinical or therapeutic recommendations

Minimum free energy (MFE), GC%, repeat metrics, and other biological scores should be computed by external benchmark tools such as ViennaRNA-backed evaluators.

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

CodonForge is not published to crates.io yet. Until the first release is published, use the source checkout path above.

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
| `-l, --lambda <FLOAT>` | Structure weighting parameter reserved for v1; ignored in v0. |
| `-v, --verbose` | Print diagnostic information to stderr. |
| `-h, --help` | Print help. |
| `-V, --version` | Print version. |

## Validation

Run the local quality gates:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo doc --no-deps
```

Run a FASTA smoke test:

```bash
cargo run --release -- \
  --input examples/toy/input.fasta \
  --fasta-out /tmp/codonforge-toy.fasta
```

## Benchmark integration

CodonForge can be used by external benchmark pipelines that evaluate mRNA outputs with a consistent metric suite. A companion research workflow uses CodonForge on retinitis-pigmentosa-relevant proteins including RHO, PRPF31, RPE65, and RPGR-ORF15.

Expected integration pattern:

1. Build CodonForge with `cargo build --release`.
2. Run CodonForge on protein FASTA inputs with `--fasta-out`.
3. Evaluate the generated mRNA FASTA files with an external metric pipeline such as ViennaRNA-backed scripts.
4. Compare CAI, MFE, GC%, U%, repeat density, and protein-fidelity metrics against other methods.

This repository does not require the external research benchmark to run the core CLI.

## Limitations

- Greedy CAI optimizer only: v0 selects the highest-frequency codon for each amino acid independently.
- High GC content: pure greedy CAI optimization produced ~81-90% GC on the initial RP benchmark targets. v1 is expected to add GC/U/repeat-aware beam search.
- No internal MFE: stdout reports `nan` for folding free energy. Compute MFE externally.
- No LinearDesign parity: CodonForge does not implement LinearDesign's dynamic programming or energy model.
- No GPU in v0: optimization is CPU-only and deterministic.
- Computational research only: no clinical, diagnostic, or therapeutic claims.

## Roadmap

| Version | Scope |
| --- | --- |
| v0 | Greedy CAI optimizer, FASTA I/O, CAI, protein-fidelity validation |
| v1 | Beam search with CAI + GC/U/repeat constraints |
| v2 | GPU batch candidate scoring for large beams |
| v3 | Approximate or exact DP research if benchmarks justify it |

## Glossary

- CAI (Codon Adaptation Index): A measure of how closely codon usage matches a reference codon table. Values are in `(0, 1]`, with higher values indicating closer match to preferred codons.
- MFE (Minimum Free Energy): RNA secondary-structure stability score, usually computed by tools such as ViennaRNA. More negative values indicate more stable predicted structures.
- GC%: Percentage of guanine and cytosine nucleotides in a sequence.
- Synonymous codon: A different codon encoding the same amino acid.
- Clean-room implementation: Software written without copying another project's source code.

## Project files

- `CONTRIBUTING.md` — contribution workflow and required checks.
- `SECURITY.md` — security reporting policy.
- `CHANGELOG.md` — release history.
- `CITATION.cff` — machine-readable citation metadata.
- `docs/v1-beam-search-design.md` — planned GC-aware beam-search design.

## Citation

If you use CodonForge in research, cite this repository using `CITATION.cff`.

Foundational reference for CAI:

```bibtex
@article{sharp1987codon,
  title = {The codon adaptation index--a measure of directional synonymous codon usage bias},
  author = {Sharp, Paul M. and Li, Wen-Hsiung},
  journal = {Nucleic Acids Research},
  volume = {15},
  number = {5},
  pages = {1281--1295},
  year = {1987},
  doi = {10.1093/nar/15.5.1281}
}
```

CodonForge may be benchmarked alongside LinearDesign and related mRNA design tools, but it does not claim algorithmic equivalence to them.
