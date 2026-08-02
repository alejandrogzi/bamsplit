// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Using `bamsplit` as a library rather than a command.
//!
//! ```console
//! cargo run --example library_split -- sample.bam out
//! ```
//!
//! The library returns the same validated manifest the CLI writes, so a Rust
//! pipeline can branch on per-output counts without shelling out or parsing
//! JSON.

use bamsplit_core::{ChromOptions, Classify as _, RunOptions};

fn main() -> std::process::ExitCode {
    let mut arguments = std::env::args_os().skip(1);
    let (Some(input), Some(out_dir)) = (arguments.next(), arguments.next()) else {
        eprintln!("usage: library_split <input.bam> <out-dir>");
        return std::process::ExitCode::from(2);
    };

    let options = ChromOptions {
        // Materialize every reference, even the empty ones, so a downstream
        // scatter has a predictable number of tasks.
        emit_empty: true,
        ..ChromOptions::new(std::path::PathBuf::from(out_dir))
    };
    let run = RunOptions {
        threads: std::thread::available_parallelism()
            .map(std::num::NonZero::get)
            .unwrap_or(1),
        manifest: bamsplit_core::ManifestFormat::Both,
        command_line: "library_split".to_string(),
        ..RunOptions::default()
    };

    match bamsplit_core::split_by_chromosome(std::path::PathBuf::from(input), &options, &run) {
        Ok(report) => {
            let manifest = &report.manifest;
            println!(
                "{} records -> {} outputs via the {} engine in {:.2}s",
                manifest.input_records,
                report.created_outputs(),
                manifest.selected_engine,
                manifest.elapsed_seconds
            );
            println!("{:<24} {:>10}  {:<12} digest", "key", "records", "index");
            for output in &manifest.outputs {
                println!(
                    "{:<24} {:>10}  {:<12} {}",
                    output.logical_key,
                    output.stats.record_count,
                    output.index_type.as_deref().unwrap_or("-"),
                    output.stats.raw_record_digest
                );
            }
            std::process::ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("bamsplit: {}", bamsplit_core::error::render_chain(&error));
            std::process::ExitCode::from(u8::try_from(error.exit_code().code()).unwrap_or(1))
        }
    }
}
