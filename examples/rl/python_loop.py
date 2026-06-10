#!/usr/bin/env python3
"""Toy RL/GFlowNet scaffold for CodonForge.

This is a standard-library demonstration only. It is not a training pipeline and
makes no biological or therapeutic claims.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path


def read_protein_fasta(path: Path) -> str:
    return "".join(
        line.strip().upper()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(">")
    )


def load_table(path: Path) -> dict[str, list[tuple[str, float]]]:
    by_aa: dict[str, list[tuple[str, float]]] = defaultdict(list)
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            aa = row["amino_acid"].strip().upper()
            if aa == "*":
                continue
            by_aa[aa].append((row["codon"].strip().upper(), float(row["frequency"])))
    return dict(by_aa)


def softmax(logits: list[float]) -> list[float]:
    max_logit = max(logits)
    exps = [math.exp(x - max_logit) for x in logits]
    total = sum(exps)
    return [x / total for x in exps]


def sample_codon(codons: list[tuple[str, float]], rng: random.Random) -> str:
    logits = [math.log(freq + 1e-12) for _, freq in codons]
    probs = softmax(logits)
    draw = rng.random()
    cumulative = 0.0
    for (codon, _), prob in zip(codons, probs):
        cumulative += prob
        if draw <= cumulative:
            return codon
    return codons[-1][0]


def metrics(seq: str) -> dict[str, float]:
    gc = sum(1 for base in seq if base in "GC")
    u = sum(1 for base in seq if base == "U")
    repeats = sum(1 for a, b in zip(seq, seq[1:]) if a == b)
    return {
        "gc_percent": 100.0 * gc / len(seq),
        "u_percent": 100.0 * u / len(seq),
        "repeat_penalty": float(repeats),
    }


def shaped_reward(seq: str, target_gc: float = 55.0, target_u: float = 20.0) -> float:
    m = metrics(seq)
    return (
        -0.02 * abs(m["gc_percent"] - target_gc)
        -0.01 * abs(m["u_percent"] - target_u)
        -0.05 * m["repeat_penalty"]
    )


def trajectory_balance_residual(reward: float, log_forward_prob: float) -> float:
    # Toy residual: real GFlowNets learn state/action flows. This just shows the
    # shape of a terminal consistency objective on one sampled trajectory.
    log_reward = reward if reward < 0 else math.log1p(reward)
    return abs(log_forward_prob - log_reward)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--codon-table", required=True, type=Path)
    parser.add_argument("--seed", default=7, type=int)
    args = parser.parse_args()

    protein = read_protein_fasta(args.input)
    table = load_table(args.codon_table)
    rng = random.Random(args.seed)

    sampled: list[str] = []
    log_forward_prob = 0.0
    for aa in protein:
        codons = table[aa]
        codon = sample_codon(codons, rng)
        sampled.append(codon)
        probs = softmax([math.log(freq + 1e-12) for _, freq in codons])
        idx = [c for c, _ in codons].index(codon)
        log_forward_prob += math.log(probs[idx] + 1e-12)

    mrna = "".join(sampled)
    reward = shaped_reward(mrna)
    m = metrics(mrna)

    print("research_preview: true")
    print("input_protein:", protein)
    print("sampled_mrna:", mrna)
    print(f"gc_percent: {m['gc_percent']:.2f}")
    print(f"u_percent: {m['u_percent']:.2f}")
    print(f"repeat_penalty: {m['repeat_penalty']:.1f}")
    print(f"toy_reward: {reward:.4f}")
    print(f"toy_trajectory_balance_residual: {trajectory_balance_residual(reward, log_forward_prob):.4f}")
    print("note: demonstration only; real RL training belongs in the research repo")


if __name__ == "__main__":
    main()
