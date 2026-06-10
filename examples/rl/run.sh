#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CODONFORGE="${CODONFORGE:-$REPO_ROOT/target/release/codonforge}"
INPUT="$SCRIPT_DIR/input.fasta"
OUTPUT_DIR="$SCRIPT_DIR/output"

if [ ! -x "$CODONFORGE" ]; then
  echo "CodonForge binary not found at $CODONFORGE" >&2
  echo "Run: cargo build --release" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "== Greedy baseline =="
"$CODONFORGE" --input "$INPUT" --strategy greedy --fasta-out "$OUTPUT_DIR/toy_greedy.fasta"

echo
echo "== Beam baseline =="
"$CODONFORGE" --input "$INPUT" --strategy beam --fasta-out "$OUTPUT_DIR/toy_beam.fasta"

echo
if command -v python3 >/dev/null 2>&1; then
  echo "== Optional Python RL toy loop =="
  python3 "$SCRIPT_DIR/python_loop.py"     --input "$INPUT"     --codon-table "$REPO_ROOT/data/codon_usage_freq_table_human.csv"
else
  echo "python3 not found; skipping optional RL toy loop"
fi
