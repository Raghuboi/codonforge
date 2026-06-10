# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- GC-aware beam-search optimizer for v1.
- Optional GPU batch candidate scoring for later releases.

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
