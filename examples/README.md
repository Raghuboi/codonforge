# CodonForge Examples

These examples verify CodonForge from a source checkout.

## Quick start

```bash
cargo build --release
cd examples/toy
./run.sh
```

## Examples

| Example | Description | Input | Output |
| --- | --- | --- | --- |
| [`toy/`](toy/) | Minimal seven-amino-acid protein (`MNDTEAI`) | Bundled FASTA | Human-readable stdout |
| [`rho/`](rho/) | RHO protein workflow for benchmark users | External FASTA | mRNA FASTA |
| [`benchmark/`](benchmark/) | Four RP-relevant proteins with greedy and beam strategies | Bundled protein FASTA | stdout + mRNA FASTA files |

## Benchmark

Run the 4-gene mini-benchmark:

```bash
cargo build --release
./examples/benchmark/run.sh
```

This runs CodonForge on RHO, PRPF31, RPE65, and RPGR-ORF15 using both `greedy` and `beam` strategies. The benchmark is a reproducibility smoke test for computational research, not a therapeutic design recommendation.
