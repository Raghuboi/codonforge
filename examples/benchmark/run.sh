#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CODONFORGE="${CODONFORGE:-$REPO_ROOT/target/release/codonforge}"
TARGETS_DIR="$SCRIPT_DIR/targets"
OUTPUT_DIR="$SCRIPT_DIR/output"

if [ ! -x "$CODONFORGE" ]; then
  echo "CodonForge binary not found at $CODONFORGE" >&2
  echo "Run: cargo build --release" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "CodonForge Benchmark"
echo "=========================================="

for target_file in "$TARGETS_DIR"/*.fasta; do
  gene=$(basename "$target_file" .fasta)

  echo ""
  echo "--- $gene (greedy) ---"
  "$CODONFORGE" \
    --input "$target_file" \
    --fasta-out "$OUTPUT_DIR/${gene}_greedy.fasta" \
    --verbose

  echo ""
  echo "--- $gene (beam) ---"
  "$CODONFORGE" \
    --input "$target_file" \
    --fasta-out "$OUTPUT_DIR/${gene}_beam.fasta" \
    --strategy beam \
    --verbose
done

echo ""
echo "=========================================="
echo "Benchmark complete"
echo "Output FASTA files: $OUTPUT_DIR"
echo "=========================================="
