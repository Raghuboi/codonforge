# Contributing

Thank you for considering a contribution to CodonForge.

CodonForge is research software. Contributions should prioritize correctness, reproducibility, and clear limitations over broad feature claims.

## Development setup

```bash
git clone https://github.com/raghuboi/codonforge.git
cd codonforge
cargo build
```

## Required checks

Run the full local gate before opening a pull request:

```bash
cargo fmt --check
cargo test
cargo clippy --all-targets --all-features -- -D warnings
cargo doc --no-deps
cargo publish --dry-run
```

Also smoke-test the CLI:

```bash
cargo run --release -- --help
cargo run --release -- --version
echo MNDTEAI | cargo run --release --
./examples/toy/run.sh
```

## Scientific and safety boundaries

Do not make clinical, diagnostic, or therapeutic claims. CodonForge is computational research software.

Do not copy source code from LinearDesign or other restricted projects. CodonForge must remain a clean-room implementation.

If you add data files, document their source, access date, license/terms, and any normalization steps.

## Commit style

Use conventional commits when practical:

```text
feat: add beam-search optimizer
fix: reject invalid amino acid input
docs: clarify codon table provenance
ci: add release workflow
```

## Pull request checklist

- [ ] Full local gate passes.
- [ ] CLI smoke tests pass.
- [ ] New behavior has tests.
- [ ] README/docs updated if user-visible behavior changed.
- [ ] No clinical or therapeutic claims added.
- [ ] Data provenance documented for any new data files.
