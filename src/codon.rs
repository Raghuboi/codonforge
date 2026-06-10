use std::collections::HashMap;

use anyhow::Context;

/// Standard genetic code: RNA codon (U) -> amino acid
pub fn codon_to_aa(codon: &str) -> Option<char> {
    match codon {
        // Stop codons
        "UAA" | "UAG" | "UGA" => Some('*'),
        // Alanine
        "GCU" | "GCC" | "GCA" | "GCG" => Some('A'),
        // Arginine
        "CGU" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => Some('R'),
        // Asparagine
        "AAU" | "AAC" => Some('N'),
        // Aspartic acid
        "GAU" | "GAC" => Some('D'),
        // Cysteine
        "UGU" | "UGC" => Some('C'),
        // Glutamic acid
        "GAA" | "GAG" => Some('E'),
        // Glutamine
        "CAA" | "CAG" => Some('Q'),
        // Glycine
        "GGU" | "GGC" | "GGA" | "GGG" => Some('G'),
        // Histidine
        "CAU" | "CAC" => Some('H'),
        // Isoleucine
        "AUU" | "AUC" | "AUA" => Some('I'),
        // Leucine
        "UUA" | "UUG" | "CUU" | "CUC" | "CUA" | "CUG" => Some('L'),
        // Lysine
        "AAA" | "AAG" => Some('K'),
        // Methionine (start)
        "AUG" => Some('M'),
        // Phenylalanine
        "UUU" | "UUC" => Some('F'),
        // Proline
        "CCU" | "CCC" | "CCA" | "CCG" => Some('P'),
        // Serine
        "UCU" | "UCC" | "UCA" | "UCG" | "AGU" | "AGC" => Some('S'),
        // Threonine
        "ACU" | "ACC" | "ACA" | "ACG" => Some('T'),
        // Tryptophan
        "UGG" => Some('W'),
        // Tyrosine
        "UAU" | "UAC" => Some('Y'),
        // Valine
        "GUU" | "GUC" | "GUA" | "GUG" => Some('V'),
        _ => None,
    }
}

/// All synonymous codons for each amino acid (RNA, U notation)
pub fn aa_to_codons(aa: char) -> Option<Vec<&'static str>> {
    match aa {
        'A' => Some(vec!["GCU", "GCC", "GCA", "GCG"]),
        'R' => Some(vec!["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]),
        'N' => Some(vec!["AAU", "AAC"]),
        'D' => Some(vec!["GAU", "GAC"]),
        'C' => Some(vec!["UGU", "UGC"]),
        'E' => Some(vec!["GAA", "GAG"]),
        'Q' => Some(vec!["CAA", "CAG"]),
        'G' => Some(vec!["GGU", "GGC", "GGA", "GGG"]),
        'H' => Some(vec!["CAU", "CAC"]),
        'I' => Some(vec!["AUU", "AUC", "AUA"]),
        'L' => Some(vec!["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"]),
        'K' => Some(vec!["AAA", "AAG"]),
        'M' => Some(vec!["AUG"]),
        'F' => Some(vec!["UUU", "UUC"]),
        'P' => Some(vec!["CCU", "CCC", "CCA", "CCG"]),
        'S' => Some(vec!["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]),
        'T' => Some(vec!["ACU", "ACC", "ACA", "ACG"]),
        'W' => Some(vec!["UGG"]),
        'Y' => Some(vec!["UAU", "UAC"]),
        'V' => Some(vec!["GUU", "GUC", "GUA", "GUG"]),
        '*' => Some(vec!["UAA", "UAG", "UGA"]),
        _ => None,
    }
}

/// Validate that a character is a valid amino acid (excluding stop)
pub fn is_valid_aa(c: char) -> bool {
    matches!(
        c,
        'A' | 'C'
            | 'D'
            | 'E'
            | 'F'
            | 'G'
            | 'H'
            | 'I'
            | 'K'
            | 'L'
            | 'M'
            | 'N'
            | 'P'
            | 'Q'
            | 'R'
            | 'S'
            | 'T'
            | 'V'
            | 'W'
            | 'Y'
    )
}

/// Translate an RNA sequence (string of codons) back to protein
pub fn translate_rna(rna: &str) -> Result<String, anyhow::Error> {
    if !rna.len().is_multiple_of(3) {
        return Err(anyhow::anyhow!(
            "RNA length {} is not a multiple of 3",
            rna.len()
        ));
    }

    let mut protein = String::with_capacity(rna.len() / 3);
    for chunk in rna.as_bytes().chunks(3) {
        let codon = std::str::from_utf8(chunk).unwrap();
        match codon_to_aa(codon) {
            Some('*') => break, // Stop codon
            Some(aa) => protein.push(aa),
            None => {
                return Err(anyhow::anyhow!("Invalid codon: {}", codon));
            }
        }
    }
    Ok(protein)
}

/// Load the bundled human codon usage table.
pub fn load_default_codon_usage() -> Result<HashMap<String, (char, f64)>, anyhow::Error> {
    parse_codon_usage_content(
        include_str!("../data/codon_usage_freq_table_human.csv"),
        "bundled data/codon_usage_freq_table_human.csv",
    )
}

/// Load codon usage table from CSV path.
/// CSV format: codon, aa, frequency (RNA U notation, may have BOM and CRLF line endings).
/// Returns HashMap of codon -> (amino_acid, frequency).
pub fn load_codon_usage(path: &str) -> Result<HashMap<String, (char, f64)>, anyhow::Error> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read codon table: {path}"))?;
    parse_codon_usage_content(&content, path)
}

fn parse_codon_usage_content(
    content: &str,
    source: &str,
) -> Result<HashMap<String, (char, f64)>, anyhow::Error> {
    // Strip UTF-8 BOM if present.
    let content = content.strip_prefix('\u{feff}').unwrap_or(content);

    // Strip CR from CRLF line endings.
    let content = content.replace('\r', "");

    let mut table: HashMap<String, (char, f64)> = HashMap::new();

    for (line_num, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let fields: Vec<&str> = trimmed.split(',').collect();
        if fields.len() < 3 {
            eprintln!(
                "Warning: skipping malformed line {} in {}: '{}'",
                line_num + 1,
                source,
                trimmed
            );
            continue;
        }

        let codon = fields[0].trim().to_uppercase();
        let aa_str = fields[1].trim();
        let freq_str = fields[2].trim();

        // Normalize T -> U for DNA-style tables.
        let codon = codon.replace('T', "U");

        // Skip header-like rows or malformed entries.
        if codon.len() != 3 || !codon.chars().all(|c| "AUGC".contains(c)) {
            continue;
        }

        let freq: f64 = freq_str.parse().with_context(|| {
            format!(
                "Failed to parse frequency '{}' on line {} in {}",
                freq_str,
                line_num + 1,
                source
            )
        })?;

        let aa = if aa_str == "*" {
            '*'
        } else {
            aa_str.chars().next().unwrap_or('*')
        };

        table.insert(codon, (aa, freq));
    }

    if table.len() < 61 {
        return Err(anyhow::anyhow!(
            "Codon table from {} only loaded {} codons, expected at least 61",
            source,
            table.len()
        ));
    }

    Ok(table)
}

/// Get max frequency for an amino acid from the codon usage table
pub fn max_freq_for_aa(table: &HashMap<String, (char, f64)>, aa: char) -> f64 {
    let Some(codons) = aa_to_codons(aa) else {
        return 0.0;
    };

    codons
        .iter()
        .filter_map(|&c| table.get(c).map(|(_, f)| *f))
        .fold(0.0_f64, f64::max)
}

/// Frequency of a specific codon from the table
pub fn codon_frequency(table: &HashMap<String, (char, f64)>, codon: &str) -> f64 {
    table.get(codon).map(|(_, f)| *f).unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_codon_to_aa() {
        assert_eq!(codon_to_aa("AUG"), Some('M'));
        assert_eq!(codon_to_aa("UAA"), Some('*'));
        assert_eq!(codon_to_aa("GGU"), Some('G'));
        assert_eq!(codon_to_aa("XXX"), None);
    }

    #[test]
    fn test_aa_to_codons() {
        let codons = aa_to_codons('A').unwrap();
        assert_eq!(codons.len(), 4);
        assert!(codons.contains(&"GCU"));

        let codons = aa_to_codons('M').unwrap();
        assert_eq!(codons.len(), 1);
        assert_eq!(codons[0], "AUG");

        assert!(aa_to_codons('X').is_none());
    }

    #[test]
    fn test_is_valid_aa() {
        assert!(is_valid_aa('M'));
        assert!(is_valid_aa('A'));
        assert!(!is_valid_aa('*'));
        assert!(!is_valid_aa('X'));
        assert!(!is_valid_aa('B'));
    }

    #[test]
    fn test_translate_rna() {
        let rna = "AUGAACGAUACGGAGGCGAUC";
        let protein = translate_rna(rna).unwrap();
        assert_eq!(protein, "MNDTEAI");
    }

    #[test]
    fn test_translate_rna_invalid_length() {
        assert!(translate_rna("AUGA").is_err());
    }

    #[test]
    fn test_translate_rna_invalid_codon() {
        assert!(translate_rna("XXX").is_err());
    }

    #[test]
    fn test_translate_rna_stop_codon() {
        let rna = "AUGUAA";
        let protein = translate_rna(rna).unwrap();
        assert_eq!(protein, "M");
    }

    #[test]
    fn test_load_codon_usage() {
        let table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();
        assert!(table.len() >= 61);
        // Check a known codon
        assert!(table.contains_key("AUG"));
    }

    #[test]
    fn test_max_freq_for_aa() {
        let table = load_codon_usage("data/codon_usage_freq_table_human.csv").unwrap();
        let max_a = max_freq_for_aa(&table, 'A');
        assert!(max_a > 0.0);
        assert!(max_a <= 1.0);
    }
}
