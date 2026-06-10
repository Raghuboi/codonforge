mod cai;
mod codon;
mod fasta;
mod optimize;
mod output;

use anyhow::Context;
use clap::Parser;
use std::collections::HashMap;

/// CodonForge v0 - A simple codon optimizer for mRNA sequences
///
/// Reads a protein sequence and generates an optimized mRNA sequence
/// using greedy highest-frequency human codons.
#[derive(Parser, Debug)]
#[command(name = "codonforge")]
#[command(version = "0.1.0")]
#[command(
    about = "Greedy codon optimizer for mRNA sequences",
    long_about = "CodonForge reads a protein sequence and generates a synonymous mRNA sequence \
using greedy highest-frequency human codons. It computes the Codon Adaptation Index \
(CAI) and can write pure FASTA output for downstream benchmarking.\n\n\
This is computational research software, not a clinical or therapeutic design tool."
)]
struct Cli {
    /// Input protein FASTA file path
    #[arg(short, long)]
    input: Option<String>,

    /// Output mRNA FASTA file path
    #[arg(short, long)]
    fasta_out: Option<String>,

    /// Codon usage table CSV path
    #[arg(short = 'c', long = "codonusage")]
    codonusage: Option<String>,

    /// Structure weighting parameter (reserved for v1; ignored in v0)
    #[arg(short = 'l', long = "lambda")]
    lambda: Option<f64>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn run() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    // Load codon usage table.
    let (codon_table, codon_source): (HashMap<String, (char, f64)>, String) =
        if let Some(codon_path) = cli.codonusage {
            (
                codon::load_codon_usage(&codon_path).context("Failed to load codon usage table")?,
                codon_path,
            )
        } else {
            (
                codon::load_default_codon_usage()
                    .context("Failed to load bundled codon usage table")?,
                "bundled data/codon_usage_freq_table_human.csv".to_string(),
            )
        };

    if cli.verbose {
        eprintln!("Loaded {} codons from {}", codon_table.len(), codon_source);
    }

    // Read protein input
    let (protein, gene_name) = if let Some(input_path) = &cli.input {
        let entries =
            fasta::parse_protein_fasta(input_path).context("Failed to parse input FASTA")?;

        if entries.len() > 1 && cli.verbose {
            eprintln!("Warning: Multiple FASTA entries found, using first");
        }

        let entry = &entries[0];
        let name = if entry.header.is_empty() {
            // Extract name from filename
            std::path::Path::new(input_path)
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "unknown".to_string())
        } else {
            entry.header.clone()
        };

        (entry.sequence.clone(), name)
    } else {
        // Read from stdin
        let seq = fasta::read_protein_from_stdin().context("Failed to read from stdin")?;
        (seq, "stdin".to_string())
    };

    // Validate protein sequence
    fasta::validate_protein_sequence(&protein).context("Invalid protein sequence")?;

    if cli.verbose {
        eprintln!(
            "Protein: {} ({} amino acids)",
            &protein[..std::cmp::min(50, protein.len())],
            protein.len()
        );
    }

    // Optimize
    let result = optimize::optimize(&protein, &codon_table).context("Optimization failed")?;

    // Print stdout output
    println!("{}", output::format_stdout(&result));

    // Write FASTA output if requested
    if let Some(fasta_path) = &cli.fasta_out {
        output::write_fasta(fasta_path, &gene_name, &result.rna)
            .context("Failed to write FASTA output")?;

        if cli.verbose {
            eprintln!("Wrote FASTA to {}", fasta_path);
        }
    }

    Ok(())
}

fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}
