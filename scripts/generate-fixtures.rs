// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Writes every synthetic BAM fixture to a directory.
//!
//! The test-suite builds fixtures in memory, so this binary exists for the times
//! a human needs one on disk: reproducing a bug with `samtools`, checking a hex
//! dump, or seeding a benchmark.
//!
//! ```console
//! cargo run --bin bamsplit-generate-fixtures -- tests/data/generated
//!
//! # One large, realistically compressible BAM for benchmarking.
//! cargo run --release --bin bamsplit-generate-fixtures -- \
//!     --benchmark 8 250000 /tmp/bench.bam
//! ```

fn main() -> std::process::ExitCode {
    let arguments: Vec<String> = std::env::args().skip(1).collect();

    // `--benchmark <references> <records-per-reference> <path>` writes one large
    // realistically compressible BAM instead of the catalogue.
    if arguments.first().map(String::as_str) == Some("--benchmark") {
        let references: usize = arguments.get(1).and_then(|v| v.parse().ok()).unwrap_or(8);
        let per_reference: usize = arguments
            .get(2)
            .and_then(|v| v.parse().ok())
            .unwrap_or(250_000);
        let path = std::path::PathBuf::from(
            arguments
                .get(3)
                .cloned()
                .unwrap_or_else(|| "bench.bam".to_string()),
        );
        let fixture = bamsplit_fixtures::Fixture::benchmark(references, per_reference);
        return match fixture.write(&path) {
            Ok(()) => {
                let size = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
                eprintln!(
                    "wrote {} records across {references} references to {} ({size} bytes)",
                    references * per_reference,
                    path.display()
                );
                std::process::ExitCode::SUCCESS
            }
            Err(error) => {
                eprintln!("cannot write {}: {error}", path.display());
                std::process::ExitCode::FAILURE
            }
        };
    }

    let directory = arguments
        .first()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::path::PathBuf::from("tests/data/generated"));

    match bamsplit_fixtures::write_catalogue(&directory) {
        Ok(paths) => {
            for path in &paths {
                println!("{}", path.display());
            }
            eprintln!("wrote {} fixtures to {}", paths.len(), directory.display());
            std::process::ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("cannot write fixtures to {}: {error}", directory.display());
            std::process::ExitCode::FAILURE
        }
    }
}
