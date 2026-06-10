#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CODONFORGE="${CODONFORGE:-$REPO_ROOT/target/release/codonforge}"

if [ ! -x "$CODONFORGE" ]; then
  echo "CodonForge binary not found at $CODONFORGE" >&2
  echo "Run: cargo build --release" >&2
  exit 1
fi

"$CODONFORGE" --input "$SCRIPT_DIR/input.fasta"
