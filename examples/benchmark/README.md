# CodonForge RP Mini-Benchmark

This directory contains a small public benchmark for exercising CodonForge on four retinitis-pigmentosa-relevant human proteins.

The benchmark is intended for reproducibility and CLI smoke testing. It is not a therapeutic design benchmark and does not make clinical claims.

## Run

From the repository root:

```bash
cargo build --release
./examples/benchmark/run.sh
```

The script runs both strategies:

- `greedy`: v0 highest-frequency codon baseline
- `beam`: v1 CAI + GC/U/repeat objective

Output FASTA files are written to `examples/benchmark/output/`, which is ignored by git.

## Targets

| File | Gene | Protein | Source |
| --- | --- | --- | --- |
| `targets/rho.fasta` | RHO | Rhodopsin | UniProt P08100 |
| `targets/prpf31.fasta` | PRPF31 | Pre-mRNA-processing factor 31 | UniProt Q9H3K4 |
| `targets/rpe65.fasta` | RPE65 | Retinoid isomerohydrolase RPE65 | UniProt Q9UBK9 |
| `targets/rpgr-orf15.fasta` | RPGR-ORF15 | Retinitis pigmentosa GTPase regulator isoform ORF15 | UniProt Q92834 / ORF15 isoform references |

Accessed: 2026-06-10.

## Metrics

CodonForge reports:

- CAI using the bundled human codon usage table
- GC percentage
- U percentage
- simple adjacent-repeat penalty
- runtime

MFE and downstream biological metrics should be computed by external benchmark tools such as ViennaRNA-backed evaluators.
