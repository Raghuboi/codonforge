# Codon Usage Data

## `codon_usage_freq_table_human.csv`

Human codon usage frequency table used by CodonForge for greedy codon selection, beam-search scoring, and CAI computation.

## Format

CSV with three columns:

```text
codon,amino_acid,frequency
```

Notes:

- Codons use RNA notation (`U` instead of `T`).
- Stop codons use `*` as the amino acid.
- The table contains 64 entries: 61 sense codons and 3 stop codons.
- Frequencies are normalized per amino-acid family and sum to approximately 1.0 for each amino acid.

## Provenance

Source: Kazusa Codon Usage Database

- URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
- Organism: Homo sapiens `[gbpri]`
- Database record: 93,487 CDS entries, 40,662,582 codons
- Accessed: 2026-06-10
- Local raw capture: `data/kazusa-homo-sapiens-codon-usage.txt`

Conversion steps:

1. Parsed Kazusa fields: codon, amino acid, source fraction, per-thousand frequency, and raw count.
2. Kept RNA notation from Kazusa output.
3. Converted stop codons to amino acid `*`.
4. Recomputed `frequency` from raw counts normalized within each amino-acid family.
5. Wrote clean UTF-8 CSV without BOM.

The source fractions in Kazusa output are rounded to two decimals. CodonForge uses the raw counts to compute higher-precision within-family frequencies.

## License note

The Rust code in this repository is MIT licensed. Codon usage data comes from the publicly accessible Kazusa Codon Usage Database; cite the source URL and access date when using CodonForge results in publications.

## Implementation notes

- CodonForge strips a UTF-8 BOM if present.
- DNA-style tables using `T` are normalized to RNA `U` during parsing.
- Malformed non-codon rows are skipped; at least 61 codons must load successfully.
