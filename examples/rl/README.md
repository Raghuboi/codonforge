# RL Toy Scaffold

This example shows how CodonForge can act as a tiny, deterministic substrate for external RL or GFlowNet experiments.

It is intentionally small:

- input protein: `MNDTEAI`
- no model weights
- no GPU dependency
- no therapeutic or biological claims
- Python is optional

Run from the repository root:

```bash
cargo build --release
./examples/rl/run.sh
```

The shell script first runs CodonForge's deterministic `greedy` and `beam` strategies. If `python3` is available, it also runs `python_loop.py`, a standard-library toy loop that samples synonymous codons and computes a shaped reward.

Real RL training should live in the research repo, not inside this Rust CLI.
