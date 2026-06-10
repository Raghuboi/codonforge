# CodonForge

A simple, clean-room Rust CLI codon optimizer for mRNA sequences.

## What is CodonForge?

CodonForge reads a protein sequence and generates a synonymous mRNA sequence using
greedy highest-frequency human codons. It computes the Codon Adaptation Index (CAI)
internally and outputs results in a format compatible with existing mRNA benchmarking
pipelines.

This is **not** a LinearDesign clone. It does not implement:
- O(N³) dynamic programming
- Turner energy model secondary structure prediction
- GPU-accelerated batch scoring
- Any code derived from LinearDesign source

## v0 Scope

- Greedy CAI optimization (highest-frequency codon per amino acid)
- CAI computation (Sharp & Li 1987 method)
- Protein sequence validation
- mRNA FASTA output for downstream benchmarking
- LinearDesign-compatible stdout format

## Install / Build

```bash
cargo build --release
```

The binary will be at `target/release/codonforge`.

## Usage

### Plain protein sequence from stdin

```bash
echo "MNDTEAI" | cargo run --release --
```

Output:
```
Input protein: MNDTEAI
mRNA sequence:  AUGAACGACACCGAGGCCAUC
mRNA structure: .....................
mRNA folding free energy: nan kcal/mol; mRNA CAI: 1.000
Runtime: 0.0000 seconds
```

### Protein FASTA input

```bash
cargo run --release -- \
  --input path/to/protein.fasta
```

### Output mRNA FASTA for benchmarking

```bash
cargo run --release -- \
  --input path/to/protein.fasta \
  --fasta-out path/to/output.fasta
```

### Custom codon usage table

```bash
cargo run --release -- \
  --input path/to/protein.fasta \
  --codonusage path/to/codon_table.csv
```

### CLI Flags

| Flag | Description |
|------|-------------|
| `-i, --input <FILE>` | Input protein FASTA file path |
| `-f, --fasta-out <FILE>` | Output mRNA FASTA file path |
| `-c, --codonusage <FILE>` | Codon usage table CSV path (default: `data/codon_usage_freq_table_human.csv`) |
| `-l, --lambda <FLOAT>` | Structure weighting parameter (parsed but ignored in v0) |
| `-v, --verbose` | Verbose output |

## Benchmark Integration

CodonForge integrates with the research benchmark pipeline:

```bash
# Run benchmark on all four RP-relevant targets
cd /home/raghuboi/Desktop/projects/research/experiments/gflownet-mrna-rp
../../.venv/bin/python scripts/run_codonforge_v0.py
```

This produces:
- `results/codonforge_v0/` - mRNA FASTA files for each target
- `results/codonforge_v0_summary.json` - Per-target metrics summary
- `results/benchmark_results_v7.csv` - Combined V5 + CodonForge results

## Limitations

- **Greedy CAI optimizer only**: Selects highest-frequency codon per amino acid. No consideration of secondary structure, GC content, or other constraints.
- **Not LinearDesign algorithm parity**: Does not implement Turner energy, DP, or GPU features.
- **No internal MFE**: Minimum free energy is reported as `nan` and must be computed externally (e.g., ViennaRNA).
- **No GPU in v0**: Single-threaded CPU-only optimization.
- **Computational research only**: No clinical or therapeutic claims. This tool is for computational benchmarking and prioritization only.

## Roadmap

| Version | Feature |
|---------|---------|
| v0 | Greedy CAI optimizer (current) |
| v1 | Beam search with CAI + GC constraints |
| v2 | GPU batch candidate scoring |
| v3 | Approximate/exact DP research |

## License

MIT

## Citation

When using CodonForge in research, please cite:

- Sharp PM, Li WH (1987). The codon adaptation index—a measure of directional synonymous codon usage bias. *Nucleic Acids Research*.
- LinearDesign / Ward et al. (2025) for mRNA design benchmarking context.

CodonForge is a clean-room implementation and does not claim algorithmic equivalence to LinearDesign.
