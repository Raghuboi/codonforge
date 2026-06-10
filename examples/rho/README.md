# RHO Protein Example

This example runs CodonForge on an externally supplied RHO protein FASTA file and writes mRNA FASTA output.

RHO is included here as a benchmark example because it is a retinitis-pigmentosa-relevant protein. This example is computational only and makes no clinical or therapeutic claims.

## Usage

```bash
cargo build --release
cd examples/rho
./run.sh /path/to/target-rho.fasta
```

The script writes:

```text
rho_output.fasta
```

You can also run CodonForge directly:

```bash
codonforge --input /path/to/target-rho.fasta --fasta-out rho_output.fasta --verbose
```

Expected v0 behavior:

- protein-faithful mRNA FASTA output
- high CAI from greedy highest-frequency codons
- high GC% is expected and documented as a v0 limitation
- MFE is not computed internally
