use std::collections::HashMap;

use crate::codon::{codon_frequency, max_freq_for_aa};

/// Compute the Codon Adaptation Index (CAI) for an mRNA sequence.
///
/// CAI measures how similar the codon usage is to the usage in a reference set
/// of highly expressed genes. Values range from 0 (non-optimal) to 1 (optimal).
///
/// Algorithm (Sharp & Li 1987):
///   For each codon i in the sequence:
///     w_i = freq(codon_i) / max_freq(amino_acid(codon_i))
///   CAI = exp(mean(log(w_i)))
///
/// If any codon has zero frequency or max_freq is 0, that codon contributes w=0
/// and the overall CAI becomes 0.
pub fn compute_cai(
    rna: &str,
    codon_table: &HashMap<String, (char, f64)>,
) -> Result<f64, anyhow::Error> {
    if !rna.len().is_multiple_of(3) {
        return Err(anyhow::anyhow!(
            "RNA length {} is not a multiple of 3",
            rna.len()
        ));
    }

    let codon_count = rna.len() / 3;
    if codon_count == 0 {
        return Err(anyhow::anyhow!("Empty RNA sequence"));
    }

    let mut log_sum = 0.0_f64;
    let mut valid_count = 0usize;

    for chunk in rna.as_bytes().chunks(3) {
        let codon = std::str::from_utf8(chunk).unwrap();
        let (aa, _) = codon_table
            .get(codon)
            .ok_or_else(|| anyhow::anyhow!("Unknown codon: {}", codon))?;

        // Skip stop codons in CAI calculation
        if *aa == '*' {
            continue;
        }

        let max_f = max_freq_for_aa(codon_table, *aa);
        if max_f <= 0.0 {
            return Ok(0.0);
        }

        let freq = codon_frequency(codon_table, codon);
        let w = freq / max_f;

        // Handle zero frequency
        if w <= 0.0 {
            return Ok(0.0);
        }

        log_sum += w.ln();
        valid_count += 1;
    }

    if valid_count == 0 {
        return Err(anyhow::anyhow!("No valid codons found for CAI calculation"));
    }

    let cai = (log_sum / valid_count as f64).exp();

    // CAI should be in [0, 1]
    if !(0.0..=1.0).contains(&cai) {
        return Err(anyhow::anyhow!("CAI out of expected range [0, 1]: {}", cai));
    }

    Ok(cai)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codon::load_codon_usage;

    #[test]
    fn test_cai_mndteai() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

        // MNDTEAI optimized sequence (greedy highest-frequency codons)
        // M=AUG, N=AAC, D=GAC, T=ACC, E=GAG, A=GCC, I=AUU
        let rna = "AUGAACGACACCGAGGCCAUU";
        let cai = compute_cai(rna, &codon_table).unwrap();

        assert!(cai > 0.0, "CAI should be positive");
        assert!(cai <= 1.0, "CAI should be <= 1");
        // Greedy optimization should give CAI close to 1.0
        assert!(
            cai >= 0.9,
            "Greedy optimized sequence should have CAI >= 0.9, got {}",
            cai
        );
    }

    #[test]
    fn test_cai_invalid_sequence() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

        assert!(compute_cai("AUGA", &codon_table).is_err());
        assert!(compute_cai("", &codon_table).is_err());
    }

    #[test]
    fn test_cai_single_codon() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

        let cai = compute_cai("AUG", &codon_table).unwrap();
        assert!(cai > 0.0);
        assert!(cai <= 1.0);
    }
}
