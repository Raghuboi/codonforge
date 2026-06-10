## Summary

<!-- What changed and why? -->

## Verification

Run before requesting review:

- [ ] `cargo fmt --check`
- [ ] `cargo test`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo doc --no-deps`
- [ ] `cargo publish --dry-run`
- [ ] `echo MNDTEAI | cargo run --release --`

## Research software checklist

- [ ] New behavior has tests.
- [ ] User-visible behavior is documented.
- [ ] New data files include source/provenance/license notes.
- [ ] No clinical, diagnostic, or therapeutic claims were added.
- [ ] No restricted-source code was copied into this repository.
