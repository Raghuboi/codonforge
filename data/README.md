# Codon Usage Data

## `codon_usage_freq_table_human.csv`

Human codon usage frequency table used by CodonForge v0 for greedy codon selection and CAI computation.

### Format

CSV with three columns:

```text
codon,amino_acid,frequency
```

Notes:

- Codons use RNA notation (`U` instead of `T`).
- Stop codons use `*` as the amino acid.
- The table contains 64 entries: 61 sense codons and 3 stop codons.
- Frequencies are normalized per amino-acid family in the bundled table.

### Provenance

The bundled file was copied from the local research benchmark fixture:

```text
experiments/gflownet-mrna-rp/methods/LinearDesign/codon_usage_freq_table_human.csv
```

That benchmark fixture originated from the LinearDesign-compatible research workflow used to evaluate CodonForge v0. The table is included to make the CLI runnable from a clean checkout and to preserve benchmark reproducibility.

The exact upstream biological database for these normalized frequencies is not independently verified in this repository. Before using CodonForge results in a publication, record the codon usage table source used in the experiment and consider replacing this bundled table with a directly cited database export such as Kazusa/Codon Usage Database or an NCBI-derived table.

### License note

The Rust code in this repository is MIT licensed. Codon usage data may be subject to the terms of its original database source. If you redistribute modified codon tables, include their source URL, access date, and license/terms in this directory.

### Implementation notes

- CodonForge strips a UTF-8 BOM if present.
- DNA-style tables using `T` are normalized to RNA `U` during parsing.
- Malformed non-codon rows are skipped; at least 61 codons must load successfully.
