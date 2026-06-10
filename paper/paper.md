---
title: 'CodonForge: A Reproducible Rust Framework for Multi-Objective mRNA Codon Optimization'
tags:
  - Rust
  - bioinformatics
  - codon optimization
  - mRNA
  - retinitis pigmentosa
  - beam search
authors:
  - name: Raghunath Prabhakar
    affiliation: "1"
affiliations:
  - name: Independent Researcher
    index: 1
date: 2026-06-10
bibliography: paper.bib
---

# Statement of need

mRNA codon optimization is a computational design step that balances expression proxies such as Codon Adaptation Index (CAI) with sequence-composition constraints such as GC content, uridine content, repeats, and downstream packaging or manufacturing considerations. Existing tools range from greedy codon optimizers to deep-learning models, but reproducible command-line baselines remain useful for comparing algorithmic trade-offs in disease-relevant benchmark sets.

CodonForge provides a clean-room Rust command-line tool for protein-faithful mRNA coding-sequence optimization research. It includes a greedy highest-frequency baseline and a deterministic beam-search strategy that balances CAI against GC%, U%, and adjacent-repeat penalties. The bundled mini-benchmark focuses on retinitis-pigmentosa-relevant proteins as computational examples only; CodonForge does not make clinical, diagnostic, or therapeutic claims.

# Summary

CodonForge reads a protein sequence and generates a synonymous mRNA sequence. The `greedy` strategy selects the highest-frequency codon for each amino acid from a cited human codon usage table. The `beam` strategy searches alternative synonymous codons with a configurable multi-objective score:

```text
score =
  + weight_cai * mean_log_cai_score
  - weight_gc * abs(GC% - target_gc)
  - weight_u * abs(U% - target_u)
  - weight_repeat * adjacent_repeat_penalty
```

Key features include:

- protein-fidelity validation by RNA translation
- FASTA input and output
- deterministic CPU-only optimization
- human codon usage data from the Kazusa Codon Usage Database
- RP mini-benchmark fixtures for RHO, PRPF31, RPE65, and RPGR-ORF15
- automated tests and CI for formatting, clippy, tests, docs, and package dry-run

# Example usage

```bash
# Greedy optimization
echo "MNDTEAI" | codonforge

# Beam search optimization
echo "MNDTEAI" | codonforge --strategy beam --target-gc 55 --target-u 20

# Mini-benchmark on RP-relevant targets
./examples/benchmark/run.sh
```

# Availability

CodonForge is available at <https://github.com/Raghuboi/codonforge> under the MIT license. Citation metadata is provided in `CITATION.cff`.

# Acknowledgements

This work was motivated by the author's personal connection to retinitis pigmentosa. All outputs are computational research artifacts and should not be interpreted as medical guidance.

# References
