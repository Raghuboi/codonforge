use std::collections::HashMap;

use crate::codon::{aa_to_codons, codon_frequency, translate_rna};

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

        // Find the codon with highest frequency
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

/// Verify that the optimized mRNA translates back to the original protein
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

/// Optimization result
#[derive(Debug)]
pub struct OptimizationResult {
    pub protein: String,
    pub rna: String,
    pub cai: f64,
    pub runtime_secs: f64,
}

/// Run the full optimization pipeline
pub fn optimize(
    protein: &str,
    codon_table: &HashMap<String, (char, f64)>,
) -> Result<OptimizationResult, anyhow::Error> {
    let start = std::time::Instant::now();

    // Optimize
    let rna = greedy_optimize(protein, codon_table)?;

    // Verify fidelity
    verify_fidelity(protein, &rna)?;

    // Compute CAI
    let cai = crate::cai::compute_cai(&rna, codon_table)?;

    let runtime = start.elapsed().as_secs_f64();

    Ok(OptimizationResult {
        protein: protein.to_string(),
        rna,
        cai,
        runtime_secs: runtime,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codon::load_codon_usage;

    #[test]
    fn test_greedy_optimize_mndteai() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

        let rna = greedy_optimize("MNDTEAI", &codon_table).unwrap();

        // Verify length
        assert_eq!(rna.len(), 21, "7 amino acids * 3 = 21 nucleotides");

        // Verify it translates back
        let protein = translate_rna(&rna).unwrap();
        assert_eq!(protein, "MNDTEAI");

        // Verify all codons are valid RNA
        assert!(rna.chars().all(|c| "AUGC".contains(c)));
    }

    #[test]
    fn test_greedy_optimize_longer() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

        let protein = "MNDTEAIGVHKR";
        let rna = greedy_optimize(protein, &codon_table).unwrap();
        assert_eq!(rna.len(), protein.len() * 3);

        let translated = translate_rna(&rna).unwrap();
        assert_eq!(translated, protein);
    }

    #[test]
    fn test_greedy_optimize_invalid_aa() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

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
        // Error could be from translate_rna (length not multiple of 3) or fidelity check
        assert!(
            err.to_string().contains("Fidelity check failed")
                || err.to_string().contains("not a multiple of 3")
        );
    }

    #[test]
    fn test_optimize_pipeline() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

        let result = optimize("MNDTEAI", &codon_table).unwrap();

        assert_eq!(result.protein, "MNDTEAI");
        assert!(result.rna.len() % 3 == 0);
        assert!(result.cai > 0.0);
        assert!(result.cai <= 1.0);
        assert!(result.runtime_secs >= 0.0);
    }

    #[test]
    fn test_optimize_all_amino_acids() {
        let codon_table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();

        // Test all 20 standard amino acids
        let protein = "ACDEFGHIKLMNPQRSTVWY";
        let result = optimize(protein, &codon_table).unwrap();

        let translated = translate_rna(&result.rna).unwrap();
        assert_eq!(translated, protein);
        assert!(result.cai > 0.0);
    }
}
