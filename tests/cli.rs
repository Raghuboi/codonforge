use std::io::Write;
use std::process::Command;

fn codonforge_cmd() -> Command {
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
    assert!(stdout.contains("--strategy"));
    assert!(stdout.contains("beam"));
}

#[test]
fn test_cli_stdin_input() {
    let mut child = codonforge_cmd()
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("Failed to spawn codonforge");

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

    let content = std::fs::read_to_string(&output_path).expect("Failed to read output file");
    assert!(content.starts_with(">"));
    assert!(content.contains("AUG"));

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
fn test_cli_toy_fasta() {
    let output = codonforge_cmd()
        .arg("--input")
        .arg("examples/toy/input.fasta")
        .output()
        .expect("Failed to execute codonforge");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Input protein: MNDTEAI"));
    assert!(stdout.contains("mRNA CAI:"));
}

#[test]
fn test_golden_toy_output() {
    let output = codonforge_cmd()
        .arg("--input")
        .arg("examples/toy/input.fasta")
        .output()
        .expect("Failed to execute codonforge");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert!(
        stdout.contains("AUGAACGACACCGAGGCCAUC"),
        "Greedy optimizer should produce deterministic output for MNDTEAI"
    );
    assert!(
        stdout.contains("mRNA CAI: 1.000"),
        "Greedy optimizer should achieve CAI=1.000 with bundled table"
    );
    assert!(stdout.contains("CodonForge strategy: greedy"));
    assert!(stdout.contains("Runtime:"));
}

#[test]
fn test_cli_beam_strategy() {
    let mut child = codonforge_cmd()
        .args([
            "--strategy",
            "beam",
            "--target-gc",
            "50",
            "--target-u",
            "25",
        ])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("Failed to spawn codonforge");

    {
        let stdin = child.stdin.as_mut().expect("Failed to get stdin");
        stdin
            .write_all(b"MNDTEAI\n")
            .expect("Failed to write to stdin");
    }

    let output = child.wait_with_output().expect("Failed to get output");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("CodonForge strategy: beam"));
    assert!(stdout.contains("GC%:"));
    assert!(stdout.contains("U%:"));
}

#[test]
fn test_cli_lambda_warning() {
    let mut child = codonforge_cmd()
        .args(["--lambda", "0.5"])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("Failed to spawn codonforge");

    {
        let stdin = child.stdin.as_mut().expect("Failed to get stdin");
        stdin
            .write_all(b"MNDTEAI\n")
            .expect("Failed to write to stdin");
    }

    let output = child.wait_with_output().expect("Failed to get output");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--lambda"));
    assert!(stderr.contains("ignored"));
}
