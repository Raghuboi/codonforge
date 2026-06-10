use crate::optimize::OptimizationResult;

/// Format output in LinearDesign-compatible stdout style
pub fn format_stdout(result: &OptimizationResult) -> String {
    let structure = ".".repeat(result.rna.len());

    format!(
        "Input protein: {}\nmRNA sequence:  {}\nmRNA structure: {}\nmRNA folding free energy: nan kcal/mol; mRNA CAI: {:.3}\nRuntime: {:.4} seconds",
        result.protein,
        result.rna,
        structure,
        result.cai,
        result.runtime_secs
    )
}

/// Format mRNA sequence as FASTA with line wrapping
pub fn format_fasta(name: &str, rna: &str, line_width: usize) -> String {
    let mut fasta = format!(">{}\n", name);

    for (count, c) in rna.chars().enumerate() {
        if count > 0 && count % line_width == 0 {
            fasta.push('\n');
        }
        fasta.push(c);
    }

    fasta
}

/// Write FASTA output to a file
pub fn write_fasta(path: &str, name: &str, rna: &str) -> Result<(), anyhow::Error> {
    let fasta = format_fasta(name, rna, 80);
    std::fs::write(path, &fasta)
        .map_err(|e| anyhow::anyhow!("Failed to write FASTA to {}: {}", path, e))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_stdout() {
        let result = OptimizationResult {
            protein: "MNDTEAI".to_string(),
            rna: "AUGAACGACACCGAGGCCAUU".to_string(),
            cai: 0.95,
            runtime_secs: 0.001,
        };

        let output = format_stdout(&result);
        assert!(output.contains("Input protein: MNDTEAI"));
        assert!(output.contains("mRNA sequence:"));
        assert!(output.contains("AUGAACGACACCGAGGCCAUU"));
        assert!(output.contains("mRNA structure:"));
        // Structure should have 21 dots for 21 nucleotides
        let structure_line = output
            .lines()
            .find(|l| l.starts_with("mRNA structure:"))
            .unwrap();
        let dots: String = structure_line.chars().filter(|&c| c == '.').collect();
        assert_eq!(
            dots.len(),
            21,
            "Structure should have 21 dots for 21 nucleotides"
        );
        assert!(output.contains("mRNA folding free energy: nan kcal/mol"));
        assert!(output.contains("mRNA CAI: 0.950"));
        assert!(output.contains("Runtime:"));
    }

    #[test]
    fn test_format_fasta() {
        let rna = "AUGAACGACACCGAGGCCAUU";
        let fasta = format_fasta("test", rna, 80);

        assert!(fasta.starts_with(">test\n"));
        assert!(fasta.contains("AUGAACGACACCGAGGCCAUU"));
    }

    #[test]
    fn test_format_fasta_wrapping() {
        let rna = "A".repeat(100);
        let fasta = format_fasta("test", &rna, 10);

        // Check that lines are wrapped at 10 characters
        let lines: Vec<&str> = fasta.lines().collect();
        assert!(lines.len() > 2); // header + at least 2 sequence lines
        for line in &lines[1..] {
            assert!(line.len() <= 10, "Line too long: {}", line.len());
        }
    }

    #[test]
    fn test_write_fasta() {
        use std::fs;
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("test_codonforge.fasta");

        write_fasta(&path.to_string_lossy(), "test", "AUGAACGAC").unwrap();

        let content = fs::read_to_string(&path).unwrap();
        assert!(content.starts_with(">test\n"));
        assert!(content.contains("AUGAACGAC"));

        fs::remove_file(path).ok();
    }
}
