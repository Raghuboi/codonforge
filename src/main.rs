mod cai;
mod codon;
mod fasta;
mod optimize;
mod output;

use anyhow::Context;
use clap::{Parser, ValueEnum};
use std::collections::HashMap;

/// CodonForge - a codon optimizer for mRNA sequence research.
///
/// Reads a protein sequence and generates a synonymous mRNA sequence using
/// either greedy highest-frequency codons or v1 beam search with composition
/// constraints.
#[derive(Parser, Debug)]
#[command(name = "codonforge")]
#[command(version)]
#[command(
    about = "Codon optimizer for mRNA sequence research",
    long_about = "CodonForge reads a protein sequence and generates a synonymous mRNA sequence. \
It supports the v0 greedy highest-frequency baseline and a v1 beam-search strategy \
that balances CAI with GC%, U%, and repeat penalties. It computes the Codon \
Adaptation Index (CAI) and can write pure FASTA output for downstream benchmarking.\n\n\
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

    /// Optimization strategy
    #[arg(long, value_enum, default_value_t = CliStrategy::Greedy)]
    strategy: CliStrategy,

    /// Beam width for --strategy beam
    #[arg(long, default_value_t = 32)]
    beam_width: usize,

    /// Target GC percentage for --strategy beam
    #[arg(long, default_value_t = 55.0)]
    target_gc: f64,

    /// Target U percentage for --strategy beam
    #[arg(long, default_value_t = 20.0)]
    target_u: f64,

    /// CAI weight for --strategy beam
    #[arg(long, default_value_t = 1.0)]
    weight_cai: f64,

    /// GC target penalty weight for --strategy beam
    #[arg(long, default_value_t = 0.02)]
    weight_gc: f64,

    /// U target penalty weight for --strategy beam
    #[arg(long, default_value_t = 0.01)]
    weight_u: f64,

    /// Homopolymer/repeat penalty weight for --strategy beam
    #[arg(long, default_value_t = 0.05)]
    weight_repeat: f64,

    /// Deprecated structure weighting parameter from LinearDesign-style workflows.
    /// Accepted for compatibility and ignored by CodonForge.
    #[arg(short = 'l', long = "lambda")]
    lambda: Option<f64>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum CliStrategy {
    Greedy,
    Beam,
}

impl From<CliStrategy> for optimize::Strategy {
    fn from(value: CliStrategy) -> Self {
        match value {
            CliStrategy::Greedy => Self::Greedy,
            CliStrategy::Beam => Self::Beam,
        }
    }
}

fn run() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    if cli.lambda.is_some() {
        eprintln!("Warning: --lambda is accepted for compatibility but ignored by CodonForge");
    }

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

    let (protein, gene_name) = if let Some(input_path) = &cli.input {
        let entries =
            fasta::parse_protein_fasta(input_path).context("Failed to parse input FASTA")?;

        if entries.len() > 1 && cli.verbose {
            eprintln!("Warning: Multiple FASTA entries found, using first");
        }

        let entry = &entries[0];
        let name = if entry.header.is_empty() {
            std::path::Path::new(input_path)
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "unknown".to_string())
        } else {
            entry.header.clone()
        };

        (entry.sequence.clone(), name)
    } else {
        let seq = fasta::read_protein_from_stdin().context("Failed to read from stdin")?;
        (seq, "stdin".to_string())
    };

    fasta::validate_protein_sequence(&protein).context("Invalid protein sequence")?;

    if cli.verbose {
        eprintln!(
            "Protein: {} ({} amino acids)",
            &protein[..std::cmp::min(50, protein.len())],
            protein.len()
        );
    }

    let beam_config = optimize::BeamConfig {
        beam_width: cli.beam_width,
        target_gc: cli.target_gc,
        target_u: cli.target_u,
        weight_cai: cli.weight_cai,
        weight_gc: cli.weight_gc,
        weight_u: cli.weight_u,
        weight_repeat: cli.weight_repeat,
    };

    let result =
        optimize::optimize_with_strategy(&protein, &codon_table, cli.strategy.into(), &beam_config)
            .context("Optimization failed")?;

    println!("{}", output::format_stdout(&result));

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
