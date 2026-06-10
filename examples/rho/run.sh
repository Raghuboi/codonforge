#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <path-to-rho-protein.fasta>" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CODONFORGE="${CODONFORGE:-$REPO_ROOT/target/release/codonforge}"
INPUT="$1"
OUTPUT="$SCRIPT_DIR/rho_output.fasta"

if [ ! -x "$CODONFORGE" ]; then
  echo "CodonForge binary not found at $CODONFORGE" >&2
  echo "Run: cargo build --release" >&2
  exit 1
fi

"$CODONFORGE" --input "$INPUT" --fasta-out "$OUTPUT" --verbose

echo "Output written to: $OUTPUT"
echo "File size: $(wc -c < "$OUTPUT") bytes"
echo "Sequence length: $(tail -n +2 "$OUTPUT" | tr -d '\n' | wc -c) nucleotides"
