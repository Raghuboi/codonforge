use std::collections::HashMap;

use crate::codon::{aa_to_codons, codon_frequency, max_freq_for_aa, translate_rna};

/// Optimization strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strategy {
    Greedy,
    Beam,
}

impl Strategy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Greedy => "greedy",
            Self::Beam => "beam",
        }
    }
}

/// Configuration for beam-search optimization.
#[derive(Debug, Clone, Copy)]
pub struct BeamConfig {
    pub beam_width: usize,
    pub target_gc: f64,
    pub target_u: f64,
    pub weight_cai: f64,
    pub weight_gc: f64,
    pub weight_u: f64,
    pub weight_repeat: f64,
}

impl Default for BeamConfig {
    fn default() -> Self {
        Self {
            beam_width: 32,
            target_gc: 55.0,
            target_u: 20.0,
            weight_cai: 1.0,
            weight_gc: 0.02,
            weight_u: 0.01,
            weight_repeat: 0.05,
        }
    }
}

#[derive(Debug, Clone)]
struct BeamCandidate {
    rna: String,
    gc_count: usize,
    u_count: usize,
    log_cai_sum: f64,
    repeat_penalty: f64,
    score: f64,
}

impl BeamCandidate {
    fn empty(capacity: usize) -> Self {
        Self {
            rna: String::with_capacity(capacity),
            gc_count: 0,
            u_count: 0,
            log_cai_sum: 0.0,
            repeat_penalty: 0.0,
            score: 0.0,
        }
    }
}

/// Greedy codon optimizer: for each amino acid, select the synonymous codon
/// with the highest frequency in the reference codon usage table.
///
/// This is a simple, fast optimization strategy that maximizes CAI locally
/// but does not consider secondary structure or other constraints.
pub fn greedy_optimize(
    protein: &str,
    codon_table: &HashMap<String, (char, f64)>,
) -> Result<String, anyhow::Error> {
    let mut rna = String::with_capacity(protein.len() * 3);

    for aa in protein.chars() {
        let Some(codons) = aa_to_codons(aa) else {
            return Err(anyhow::anyhow!("Invalid amino acid: {}", aa));
        };

        let mut best_codon = codons[0];
        let mut best_freq = codon_frequency(codon_table, codons[0]);

        for &codon in &codons[1..] {
            let freq = codon_frequency(codon_table, codon);
            if freq > best_freq {
                best_codon = codon;
                best_freq = freq;
            }
        }

        rna.push_str(best_codon);
    }

    Ok(rna)
}

/// Beam-search optimizer that balances CAI against sequence-composition targets.
///
/// The objective is intentionally simple and deterministic: preserve protein
/// fidelity, keep CAI reasonable, and avoid the extreme GC-rich/zero-U output
/// produced by pure greedy CAI optimization.
pub fn beam_optimize(
    protein: &str,
    codon_table: &HashMap<String, (char, f64)>,
    config: &BeamConfig,
) -> Result<String, anyhow::Error> {
    if config.beam_width == 0 {
        return Err(anyhow::anyhow!("beam_width must be greater than 0"));
    }

    let capacity = protein.len() * 3;
    let mut beam = vec![BeamCandidate::empty(capacity)];

    for aa in protein.chars() {
        let Some(codons) = aa_to_codons(aa) else {
            return Err(anyhow::anyhow!("Invalid amino acid: {}", aa));
        };

        let mut next = Vec::with_capacity(beam.len() * codons.len());
        for candidate in &beam {
            let prev_codon = candidate.rna.get(candidate.rna.len().saturating_sub(3)..);
            for &codon in &codons {
                next.push(expand_candidate(
                    candidate,
                    codon,
                    prev_codon,
                    aa,
                    codon_table,
                    config,
                ));
            }
        }

        next.sort_by(|a, b| b.score.total_cmp(&a.score).then_with(|| a.rna.cmp(&b.rna)));
        next.truncate(config.beam_width);
        beam = next;
    }

    let best = beam
        .into_iter()
        .max_by(|a, b| a.score.total_cmp(&b.score).then_with(|| b.rna.cmp(&a.rna)))
        .ok_or_else(|| anyhow::anyhow!("Beam search produced no candidates"))?;

    Ok(best.rna)
}

fn expand_candidate(
    candidate: &BeamCandidate,
    codon: &str,
    prev_codon: Option<&str>,
    aa: char,
    codon_table: &HashMap<String, (char, f64)>,
    config: &BeamConfig,
) -> BeamCandidate {
    let max_freq = max_freq_for_aa(codon_table, aa);
    let freq = codon_frequency(codon_table, codon);
    let cai_weight = if max_freq > 0.0 && freq > 0.0 {
        (freq / max_freq).ln()
    } else {
        f64::NEG_INFINITY
    };

    let mut rna = candidate.rna.clone();
    rna.push_str(codon);

    let mut expanded = BeamCandidate {
        rna,
        gc_count: candidate.gc_count + gc_in_codon(codon),
        u_count: candidate.u_count + u_in_codon(codon),
        log_cai_sum: candidate.log_cai_sum + cai_weight,
        repeat_penalty: candidate.repeat_penalty + repeat_penalty_increment(codon, prev_codon),
        score: 0.0,
    };
    expanded.score = score_candidate(&expanded, config);
    expanded
}

fn score_candidate(candidate: &BeamCandidate, config: &BeamConfig) -> f64 {
    let len = candidate.rna.len();
    if len == 0 {
        return 0.0;
    }

    let gc_pct = 100.0 * candidate.gc_count as f64 / len as f64;
    let u_pct = 100.0 * candidate.u_count as f64 / len as f64;

    let codon_count = (len / 3).max(1) as f64;
    let mean_log_cai = candidate.log_cai_sum / codon_count;
    let mean_repeat_penalty = candidate.repeat_penalty / codon_count;

    config.weight_cai * mean_log_cai
        - config.weight_gc * (gc_pct - config.target_gc).abs()
        - config.weight_u * (u_pct - config.target_u).abs()
        - config.weight_repeat * mean_repeat_penalty
}

fn gc_in_codon(codon: &str) -> usize {
    codon.chars().filter(|&c| c == 'G' || c == 'C').count()
}

fn u_in_codon(codon: &str) -> usize {
    codon.chars().filter(|&c| c == 'U').count()
}

fn repeat_penalty_increment(codon: &str, prev_codon: Option<&str>) -> f64 {
    let chars: Vec<char> = codon.chars().collect();
    let mut penalty = 0.0;

    for i in 1..chars.len() {
        if chars[i] == chars[i - 1] {
            penalty += 1.0;
        }
    }

    if let Some(prev) = prev_codon {
        if let (Some(last), Some(first)) = (prev.chars().last(), chars.first()) {
            if last == *first {
                penalty += 1.0;
            }
        }
    }

    penalty
}

pub fn nucleotide_percent(rna: &str, bases: &[char]) -> f64 {
    if rna.is_empty() {
        return 0.0;
    }
    let count = rna.chars().filter(|c| bases.contains(c)).count();
    100.0 * count as f64 / rna.len() as f64
}

pub fn repeat_penalty(rna: &str) -> f64 {
    let mut penalty = 0.0;
    let mut prev = None;
    for c in rna.chars() {
        if prev == Some(c) {
            penalty += 1.0;
        }
        prev = Some(c);
    }
    penalty
}

/// Verify that the optimized mRNA translates back to the original protein.
pub fn verify_fidelity(protein: &str, rna: &str) -> Result<(), anyhow::Error> {
    let translated = translate_rna(rna)?;

    if translated != protein {
        return Err(anyhow::anyhow!(
            "Fidelity check failed: expected '{}' but got '{}'",
            protein,
            translated
        ));
    }

    Ok(())
}

/// Optimization result.
#[derive(Debug)]
pub struct OptimizationResult {
    pub protein: String,
    pub rna: String,
    pub cai: f64,
    pub gc_percent: f64,
    pub u_percent: f64,
    pub repeat_penalty: f64,
    pub strategy: Strategy,
    pub runtime_secs: f64,
}

/// Run the full optimization pipeline with a selected strategy.
pub fn optimize_with_strategy(
    protein: &str,
    codon_table: &HashMap<String, (char, f64)>,
    strategy: Strategy,
    beam_config: &BeamConfig,
) -> Result<OptimizationResult, anyhow::Error> {
    let start = std::time::Instant::now();

    let rna = match strategy {
        Strategy::Greedy => greedy_optimize(protein, codon_table)?,
        Strategy::Beam => beam_optimize(protein, codon_table, beam_config)?,
    };

    verify_fidelity(protein, &rna)?;
    let cai = crate::cai::compute_cai(&rna, codon_table)?;
    let runtime = start.elapsed().as_secs_f64();

    Ok(OptimizationResult {
        protein: protein.to_string(),
        gc_percent: nucleotide_percent(&rna, &['G', 'C']),
        u_percent: nucleotide_percent(&rna, &['U']),
        repeat_penalty: repeat_penalty(&rna),
        rna,
        cai,
        strategy,
        runtime_secs: runtime,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codon::load_codon_usage;

    fn table() -> HashMap<String, (char, f64)> {
        load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap()
    }

    #[test]
    fn test_greedy_optimize_mndteai() {
        let codon_table = table();
        let rna = greedy_optimize("MNDTEAI", &codon_table).unwrap();

        assert_eq!(rna.len(), 21, "7 amino acids * 3 = 21 nucleotides");
        assert_eq!(translate_rna(&rna).unwrap(), "MNDTEAI");
        assert!(rna.chars().all(|c| "AUGC".contains(c)));
    }

    #[test]
    fn test_beam_optimize_mndteai() {
        let codon_table = table();
        let config = BeamConfig {
            beam_width: 8,
            target_gc: 50.0,
            target_u: 25.0,
            weight_cai: 1.0,
            weight_gc: 0.05,
            weight_u: 0.05,
            weight_repeat: 0.05,
        };
        let rna = beam_optimize("MNDTEAI", &codon_table, &config).unwrap();

        assert_eq!(rna.len(), 21);
        assert_eq!(translate_rna(&rna).unwrap(), "MNDTEAI");
        assert!(rna.chars().all(|c| "AUGC".contains(c)));
    }

    #[test]
    fn test_beam_rejects_zero_width() {
        let codon_table = table();
        let config = BeamConfig {
            beam_width: 0,
            ..BeamConfig::default()
        };
        let err = beam_optimize("MNDTEAI", &codon_table, &config).unwrap_err();
        assert!(err.to_string().contains("beam_width"));
    }

    #[test]
    fn test_greedy_optimize_longer() {
        let codon_table = table();
        let protein = "MNDTEAIGVHKR";
        let rna = greedy_optimize(protein, &codon_table).unwrap();
        assert_eq!(rna.len(), protein.len() * 3);
        assert_eq!(translate_rna(&rna).unwrap(), protein);
    }

    #[test]
    fn test_greedy_optimize_invalid_aa() {
        let codon_table = table();
        let err = greedy_optimize("MXN", &codon_table).unwrap_err();
        assert!(err.to_string().contains("Invalid amino acid"));
    }

    #[test]
    fn test_verify_fidelity_pass() {
        verify_fidelity("MNDTEAI", "AUGAACGACACCGAGGCCAUU").unwrap();
    }

    #[test]
    fn test_verify_fidelity_fail() {
        let err = verify_fidelity("MNDTEAI", "AUGAUGC").unwrap_err();
        assert!(
            err.to_string().contains("Fidelity check failed")
                || err.to_string().contains("not a multiple of 3")
        );
    }

    #[test]
    fn test_optimize_pipeline() {
        let codon_table = table();
        let result = optimize_with_strategy(
            "MNDTEAI",
            &codon_table,
            Strategy::Greedy,
            &BeamConfig::default(),
        )
        .unwrap();

        assert_eq!(result.protein, "MNDTEAI");
        assert_eq!(result.strategy, Strategy::Greedy);
        assert!(result.rna.len().is_multiple_of(3));
        assert!(result.cai > 0.0);
        assert!(result.cai <= 1.0);
        assert!((0.0..=100.0).contains(&result.gc_percent));
        assert!((0.0..=100.0).contains(&result.u_percent));
        assert!(result.runtime_secs >= 0.0);
    }

    #[test]
    fn test_beam_pipeline() {
        let codon_table = table();
        let result = optimize_with_strategy(
            "ACDEFGHIKLMNPQRSTVWY",
            &codon_table,
            Strategy::Beam,
            &BeamConfig::default(),
        )
        .unwrap();

        assert_eq!(translate_rna(&result.rna).unwrap(), result.protein);
        assert_eq!(result.strategy, Strategy::Beam);
        assert!(result.cai > 0.0);
        assert!(result.cai <= 1.0);
    }

    #[test]
    fn test_all_proteins_translate_back() {
        let codon_table = table();
        let test_proteins = [
            "M",
            "A",
            "R",
            "N",
            "D",
            "C",
            "E",
            "Q",
            "G",
            "H",
            "I",
            "L",
            "K",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
            "ACDEFGHIKLMNPQRSTVWY",
            "MNDTEAI",
            "MNDTEAIGVHKR",
        ];

        for protein in test_proteins {
            for strategy in [Strategy::Greedy, Strategy::Beam] {
                let result =
                    optimize_with_strategy(protein, &codon_table, strategy, &BeamConfig::default())
                        .unwrap_or_else(|e| {
                            panic!("optimize failed for {protein} with {strategy:?}: {e}")
                        });
                assert_eq!(translate_rna(&result.rna).unwrap(), protein);
                assert_eq!(result.rna.len(), protein.len() * 3);
                assert!(result.rna.chars().all(|c| "AUGC".contains(c)));
            }
        }
    }

    #[test]
    fn test_beam_reduces_gc_on_gc_rich_protein() {
        let codon_table = table();
        let protein = "MSLADELLADLEEAAEEEEGGSYGEEEEEPAIEDVQEETQLDLSGDSVKT";
        let greedy = optimize_with_strategy(
            protein,
            &codon_table,
            Strategy::Greedy,
            &BeamConfig::default(),
        )
        .unwrap();
        let beam = optimize_with_strategy(
            protein,
            &codon_table,
            Strategy::Beam,
            &BeamConfig::default(),
        )
        .unwrap();

        assert_eq!(translate_rna(&beam.rna).unwrap(), protein);
        assert!(
            beam.gc_percent < greedy.gc_percent,
            "beam GC% {} should be lower than greedy GC% {}",
            beam.gc_percent,
            greedy.gc_percent
        );
        assert!(
            beam.u_percent > greedy.u_percent,
            "beam U% {} should be higher than greedy U% {}",
            beam.u_percent,
            greedy.u_percent
        );
        assert!(
            beam.cai > 0.80,
            "beam CAI should stay useful, got {}",
            beam.cai
        );
    }

    #[test]
    fn test_nucleotide_percent() {
        assert_eq!(nucleotide_percent("", &['G', 'C']), 0.0);
        assert_eq!(nucleotide_percent("AUGC", &['G', 'C']), 50.0);
        assert_eq!(nucleotide_percent("UUUU", &['U']), 100.0);
    }
}
