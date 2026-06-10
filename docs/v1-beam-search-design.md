# CodonForge v1 Beam-Search Design

CodonForge v0 is a deterministic greedy CAI optimizer. It is useful as a baseline but produces high-GC sequences on RP benchmark targets because the highest-frequency human codons are often GC-rich.

v1 should add a deterministic CPU beam-search optimizer that trades a small amount of CAI for more realistic sequence composition.

## Goals

- Preserve exact protein fidelity.
- Keep greedy v0 behavior available.
- Add a beam strategy that reduces extreme GC%.
- Keep runtime practical for RHO, PRPF31, RPE65, and RPGR-ORF15.
- Avoid LinearDesign parity claims.
- Avoid O(N³) folding DP and Turner energy implementation.

## Proposed CLI

```bash
codonforge \
  --input protein.fasta \
  --strategy beam \
  --beam-width 32 \
  --target-gc 55 \
  --target-u 20 \
  --weight-cai 1.0 \
  --weight-gc 0.02 \
  --weight-u 0.01 \
  --weight-repeat 0.05 \
  --fasta-out output.fasta
```

Defaults:

```text
strategy = greedy
beam_width = 32
target_gc = 55.0
target_u = 20.0
weight_cai = 1.0
weight_gc = 0.02
weight_u = 0.01
weight_repeat = 0.05
```

## Scoring sketch

For each candidate partial mRNA sequence:

```text
score =
  + weight_cai * running_log_cai_score
  - weight_gc * abs(current_gc_percent - target_gc)
  - weight_u * abs(current_u_percent - target_u)
  - weight_repeat * homopolymer_repeat_penalty
```

Implementation detail: maintain cumulative counts rather than recomputing metrics from the full string at every beam expansion.

## Algorithm

```text
beam = [empty candidate]
for amino_acid in protein:
  next = []
  for candidate in beam:
    for codon in synonymous_codons(amino_acid):
      expanded = candidate + codon
      score expanded with cumulative metrics
      push expanded into next
  beam = top_k(next, beam_width)
return best candidate
```

## Verification target

Run v1 on the four RP benchmark targets and compare against v0:

- RHO
- PRPF31
- RPE65
- RPGR-ORF15

Success criteria for first v1 pass:

- 4/4 targets run.
- Protein translation exactly matches input for every target.
- GC% improves versus v0 for at least 3/4 targets.
- CAI remains above 0.65 for at least 3/4 targets.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test` pass.

## Out of scope for v1

- GPU kernels
- Turner energy model
- O(N³) folding DP
- UTR design
- Pareto-front output
- clinical interpretation
