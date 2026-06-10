use anyhow::Context;
use std::fs;

use crate::codon::is_valid_aa;

/// Parsed FASTA entry
#[derive(Debug, Clone)]
pub struct FastaEntry {
    pub header: String,
    pub sequence: String,
}

/// Parse a protein FASTA file
pub fn parse_protein_fasta(path: &str) -> Result<Vec<FastaEntry>, anyhow::Error> {
    let content =
        fs::read_to_string(path).with_context(|| format!("Failed to read file: {}", path))?;
    parse_protein_fasta_str(&content)
}

/// Parse protein FASTA from a string
pub fn parse_protein_fasta_str(content: &str) -> Result<Vec<FastaEntry>, anyhow::Error> {
    let mut entries: Vec<FastaEntry> = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq: String = String::new();

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if let Some(stripped) = trimmed.strip_prefix('>') {
            // Save previous entry if exists
            if let Some(header) = current_header.take() {
                if !current_seq.is_empty() {
                    entries.push(FastaEntry {
                        header,
                        sequence: current_seq.clone(),
                    });
                    current_seq.clear();
                }
            }
            // New header
            current_header = Some(stripped.trim().to_string());
        } else {
            // Sequence line - strip whitespace and uppercase
            for c in trimmed.chars() {
                if !c.is_whitespace() {
                    current_seq.push(c.to_ascii_uppercase());
                }
            }
        }
    }

    // Don't forget the last entry
    if let Some(header) = current_header {
        if !current_seq.is_empty() {
            entries.push(FastaEntry {
                header,
                sequence: current_seq,
            });
        }
    }

    if entries.is_empty() {
        return Err(anyhow::anyhow!("No valid FASTA entries found"));
    }

    Ok(entries)
}

/// Parse a plain protein sequence from stdin
pub fn read_protein_from_stdin() -> Result<String, anyhow::Error> {
    let mut input = String::new();
    std::io::stdin()
        .read_line(&mut input)
        .context("Failed to read from stdin")?;

    let seq: String = input
        .chars()
        .filter(|c| !c.is_whitespace())
        .map(|c| c.to_ascii_uppercase())
        .collect();

    if seq.is_empty() {
        return Err(anyhow::anyhow!("Empty input from stdin"));
    }

    Ok(seq)
}

/// Validate that a protein sequence contains only valid amino acids
pub fn validate_protein_sequence(seq: &str) -> Result<(), anyhow::Error> {
    for (i, c) in seq.chars().enumerate() {
        if !is_valid_aa(c) {
            return Err(anyhow::anyhow!(
                "Invalid amino acid '{}' at position {}",
                c,
                i
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_fasta_single_entry() {
        let content = ">test_seq\nMNDTEAI";
        let entries = parse_protein_fasta_str(content).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].header, "test_seq");
        assert_eq!(entries[0].sequence, "MNDTEAI");
    }

    #[test]
    fn test_parse_fasta_multi_entry() {
        let content = ">seq1\nMNDTEAI\n>seq2\nABCDEFG";
        let entries = parse_protein_fasta_str(content).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence, "MNDTEAI");
        assert_eq!(entries[1].sequence, "ABCDEFG");
    }

    #[test]
    fn test_parse_fasta_with_whitespace() {
        let content = ">seq1\nMN DT\nEA I";
        let entries = parse_protein_fasta_str(content).unwrap();
        assert_eq!(entries[0].sequence, "MNDTEAI");
    }

    #[test]
    fn test_parse_fasta_lowercase() {
        let content = ">seq1\nmndteai";
        let entries = parse_protein_fasta_str(content).unwrap();
        assert_eq!(entries[0].sequence, "MNDTEAI");
    }

    #[test]
    fn test_parse_fasta_empty() {
        let entries = parse_protein_fasta_str("").unwrap_err();
        assert!(entries.to_string().contains("No valid FASTA entries"));
    }

    #[test]
    fn test_validate_protein_valid() {
        validate_protein_sequence("MNDTEAI").unwrap();
    }

    #[test]
    fn test_validate_protein_invalid() {
        let err = validate_protein_sequence("MNXTEAI").unwrap_err();
        assert!(err.to_string().contains("Invalid amino acid"));
    }
}
