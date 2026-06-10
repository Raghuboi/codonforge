# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Added an opt-in RL/GFlowNet research starting path: `docs/rl-research-start.md` plus runnable `examples/rl/` toy scaffold. No Rust API, dependencies, model weights, or default behavior changed.

## [0.2.0] - 2026-06-10

### Added
- v1 beam-search optimizer with configurable CAI, GC%, U%, and repeat-penalty objective.
- CLI flags for `--strategy`, `--beam-width`, `--target-gc`, `--target-u`, and objective weights.
- GC%, U%, and repeat-penalty metrics in stdout output.
- Public RP mini-benchmark in `examples/benchmark/` with four protein FASTA targets.
- Golden CLI tests and beam-search fidelity tests.
- JOSS paper scaffold and code of conduct.

### Changed
- Replaced the LinearDesign-derived bundled codon usage table with a Kazusa Codon Usage Database Homo sapiens export normalized from raw counts.
- Documented CAI table dependency and codon table provenance.
- `--lambda` now emits an explicit compatibility warning when supplied.

## [0.1.0] - 2026-06-10

### Added
- Greedy highest-frequency codon optimizer for protein-faithful human mRNA sequences.
- Codon Adaptation Index (CAI) computation using the Sharp and Li method.
- Protein sequence validation and RNA translation fidelity verification.
- Protein FASTA input and pure mRNA FASTA output.
- LinearDesign-style stdout for human-readable benchmarking logs.
- Human codon usage frequency table bundled under `data/`.
- Unit and integration tests for parser, optimizer, CAI, output formatting, and CLI behavior.
- MIT license.
