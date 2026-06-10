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

The RHO example intentionally does not bundle disease-target data. Pass your own protein FASTA file or use a benchmark dataset you are licensed to use.
