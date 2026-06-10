# RL / GFlowNet Research Starting Path

CodonForge is a deterministic Rust baseline for protein-faithful synonymous mRNA sequence generation. This page explains how to use it as a small substrate for reinforcement-learning and GFlowNet experiments without turning the public CLI into a model-training framework.

## Goal

Use CodonForge to define a clean action space and baseline reward components for external RL/GFlowNet experiments on mRNA codon optimization.

## Non-goals

CodonForge does not claim therapeutic utility and does not include:

- clinical or patient recommendations
- in-repo model weights
- GPU training code
- structure-aware folding rewards
- UTR design
- wet-lab validation claims

Use this as computational research infrastructure only: benchmarking, prioritization, and hypothesis generation.

## Existing substrate

CodonForge already exposes the pieces an external RL loop needs:

- `src/codon.rs` defines amino-acid to synonymous-codon mappings and loads codon-usage tables.
- `src/optimize.rs` defines `Strategy`, `BeamConfig`, `OptimizationResult`, and the fidelity checks.
- `Strategy::Greedy` gives a deterministic highest-frequency baseline.
- `Strategy::Beam` gives a deterministic multi-objective baseline over CAI, GC%, U%, and adjacent-repeat penalty.
- `data/codon_usage_freq_table_human.csv` is the bundled Kazusa-derived human codon table.
- `examples/benchmark/` provides four RP/IRD-relevant protein targets: RHO, PRPF31, RPE65, and RPGR-ORF15.

The key constraint is translation fidelity: a generated mRNA coding sequence must translate back to the input protein sequence.

## RL formulation

A codon optimizer can be represented as a finite-horizon sequence decision problem:

```text
state_t  = input protein + partial mRNA up to amino acid t
action_t = choose one synonymous codon for protein[t]
terminal = full mRNA coding sequence
hard constraint = translate(mRNA) == input protein
reward = composition / expression / repeat / downstream proxy score
```

A simple shaped reward for starter experiments:

```text
reward(x) =
  + alpha * CAI(x)
  - beta  * abs(GC%(x) - target_gc)
  - gamma * abs(U%(x) - target_u)
  - delta * adjacent_repeat_penalty(x)
```

This matches the v1 beam-search objective closely enough to make beam search a deterministic baseline for RL experiments.

## Minimal experiment tracks

### 1. REINFORCE / policy-gradient toy loop

Use the CodonForge codon table as the action space. Start with a tiny policy over synonymous codons and update logit weights from terminal reward.

```text
for protein in targets:
  for episode in episodes:
    sample one codon per amino acid from policy
    compute CAI, GC%, U%, repeat penalty
    reward = shaped_reward(sequence)
    update logits with reward-normalized policy gradient
```

This is useful only as a smoke test. It is not a publishable method by itself.

### 2. PPO-style codon policy

Use the same state/action space but add a small actor/critic model outside CodonForge. Keep CodonForge as the deterministic evaluator/baseline and do not add heavy ML dependencies to the Rust CLI.

```text
actor(state_t) -> distribution over synonymous codons
critic(state_t) -> expected terminal reward
objective = PPO clipped policy loss + value loss + entropy bonus
```

### 3. GFlowNet trajectory balance

Treat protein-to-mRNA generation as a directed acyclic graph where each complete synonymous mRNA has reward `R(x)`.

```text
forward flow: partial_mRNA_t -> partial_mRNA_{t+1}
terminal reward: R(full_mRNA)
training loss: trajectory-balance consistency
sample objective: diverse high-reward mRNA sequences, not only the single best sequence
```

This is the most relevant path if the research goal is a Pareto front of diverse candidate sequences.

## RP/IRD-specific reward ideas

For the current mini-benchmark, useful terms include:

- RHO: composition realism plus translation fidelity for a compact photoreceptor protein target.
- PRPF31: splice-aware penalties are a potential extension, but require external splice predictors.
- RPE65: AAV cassette feasibility and sequence composition tradeoffs.
- RPGR-ORF15: purine-rich / low-cytosine constraints should be evaluated carefully because ORF15 has unusual sequence biology.

Do not encode disease or clinical efficacy as a reward unless backed by direct experimental or clinical evidence.

## Suggested workflow

1. Run deterministic baselines:

```bash
cargo build --release
./examples/benchmark/run.sh
```

2. Run the optional toy RL scaffold:

```bash
./examples/rl/run.sh
```

3. Move real RL training into the research repo, not into CodonForge:

```text
/home/raghuboi/Desktop/projects/research/experiments/gflownet-mrna-rp/
```

4. Cite CodonForge version, git commit, codon table source, benchmark inputs, and evaluator version in every report.

## Safety language

Allowed:

- computational benchmark
- in-silico candidate prioritization
- proxy metric
- sequence-design tradeoff
- reproducibility baseline
- hypothesis generation

Avoid without direct evidence:

- treatment works
- cure
- efficacy
- patient recommendation
- clinical candidate
- therapeutic claim
