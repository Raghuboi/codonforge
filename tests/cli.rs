use std::io::Write;
use std::process::Command;

fn codonforge_cmd() -> Command {
    // Use cargo run for testing
    let mut cmd = Command::new("cargo");
    cmd.args(["run", "--release", "--bin", "codonforge", "--"]);
    cmd
}

#[test]
fn test_cli_help() {
    let output = codonforge_cmd()
        .arg("--help")
        .output()
        .expect("Failed to execute codonforge");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("codonforge"));
    assert!(stdout.contains("--input"));
    assert!(stdout.contains("--fasta-out"));
}

#[test]
fn test_cli_stdin_input() {
    let mut child = codonforge_cmd()
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("Failed to spawn codonforge");

    // Write to stdin
    {
        let stdin = child.stdin.as_mut().expect("Failed to get stdin");
        stdin
            .write_all(b"MNDTEAI\n")
            .expect("Failed to write to stdin");
    }

    let output = child.wait_with_output().expect("Failed to get output");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Input protein: MNDTEAI"));
    assert!(stdout.contains("mRNA sequence:"));
    assert!(stdout.contains("mRNA CAI:"));
}

#[test]
fn test_cli_fasta_input() {
    // Create a temporary FASTA file with unique name
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join(format!("test_input_{}.fasta", std::process::id()));

    std::fs::write(&input_path, ">test_protein\nMNDTEAI\n").expect("Failed to write temp file");

    let output = codonforge_cmd()
        .arg("--input")
        .arg(&input_path)
        .output()
        .expect("Failed to execute codonforge");

    let success = output.status.success();
    let stderr = String::from_utf8_lossy(&output.stderr);
    let stdout = String::from_utf8_lossy(&output.stdout);

    std::fs::remove_file(&input_path).ok();

    if !success {
        eprintln!("STDERR: {}", stderr);
        eprintln!("STDOUT: {}", stdout);
    }

    assert!(success, "codonforge should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Input protein: MNDTEAI"));
}

#[test]
fn test_cli_fasta_output() {
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("test_input.fasta");
    let output_path = temp_dir.join("test_output.fasta");

    std::fs::write(&input_path, ">test\nMNDTEAI\n").expect("Failed to write input");

    let output = codonforge_cmd()
        .arg("--input")
        .arg(&input_path)
        .arg("--fasta-out")
        .arg(&output_path)
        .output()
        .expect("Failed to execute codonforge");

    std::fs::remove_file(&input_path).ok();

    assert!(output.status.success());

    // Check that output file exists and contains FASTA
    let content = std::fs::read_to_string(&output_path).expect("Failed to read output file");
    assert!(content.starts_with(">"));
    assert!(content.contains("AUG")); // Start codon

    std::fs::remove_file(&output_path).ok();
}

#[test]
fn test_cli_invalid_amino_acid() {
    let mut child = codonforge_cmd()
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("Failed to spawn codonforge");

    // Write invalid amino acid
    {
        let stdin = child.stdin.as_mut().expect("Failed to get stdin");
        stdin.write_all(b"MXN\n").expect("Failed to write to stdin");
    }

    let output = child.wait_with_output().expect("Failed to get output");
    assert!(!output.status.success());

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Invalid amino acid") || stderr.contains("Error:"));
}

#[test]
fn test_cli_real_fasta() {
    // Test with actual research target
    let target_path =
        "/home/raghuboi/Desktop/projects/research/experiments/gflownet-mrna-rp/target-rho.fasta";

    // Skip if file doesn't exist (e.g., in CI)
    if !std::path::Path::new(target_path).exists() {
        return;
    }

    let output = codonforge_cmd()
        .arg("--input")
        .arg(target_path)
        .output()
        .expect("Failed to execute codonforge");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Input protein:"));
    assert!(stdout.contains("mRNA CAI:"));
}
