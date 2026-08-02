// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The command-line surface: exit codes, stream discipline, and help text.

use std::process::Command;

use bamsplit_fixtures::Fixture;

use crate::support::*;

fn bamsplit() -> Command {
    Command::new(binary())
}

#[test]
fn help_and_version_succeed() {
    for arguments in [vec!["--help"], vec!["--version"], vec!["chrom", "--help"]] {
        let output = bamsplit()
            .args(&arguments)
            .output()
            .unwrap_or_else(|error| panic!("{arguments:?}: {error}"));
        assert!(output.status.success(), "{arguments:?} failed");
        assert!(!output.stdout.is_empty(), "{arguments:?} printed nothing");
    }
}

#[test]
fn help_documents_every_subcommand_and_exit_code() {
    let output = bamsplit().arg("--help").output().expect("runs");
    let text = String::from_utf8_lossy(&output.stdout);
    for expected in ["chrom", "shard", "tag", "region", "inspect", "EXIT CODES"] {
        assert!(
            text.contains(expected),
            "`--help` omits {expected}:\n{text}"
        );
    }
}

#[test]
fn an_unknown_flag_exits_with_the_argument_code() {
    let status = bamsplit()
        .args(["chrom", "in.bam", "--nonesuch"])
        .status()
        .expect("runs");
    assert_eq!(status.code(), Some(2));
}

#[test]
fn a_successful_split_exits_zero_and_keeps_stdout_clean() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    let output = bamsplit()
        .args(["chrom"])
        .arg(&workspace.input)
        .arg("--out-dir")
        .arg(&out)
        .args(["--compression-level", "1"])
        .output()
        .expect("runs");

    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stdout.is_empty(),
        "a split must not write to stdout: {:?}",
        String::from_utf8_lossy(&output.stdout)
    );
    assert!(!output.stderr.is_empty(), "logs belong on stderr");
    assert_lossless(&fixture, &out);
}

#[test]
fn an_existing_output_directory_exits_with_the_conflict_code() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    for expected in [0, 4] {
        let status = bamsplit()
            .args(["chrom"])
            .arg(&workspace.input)
            .arg("--out-dir")
            .arg(&out)
            .args(["--compression-level", "1", "--quiet"])
            .status()
            .expect("runs");
        assert_eq!(status.code(), Some(expected));
    }
}

#[test]
fn malformed_input_exits_with_the_input_code() {
    let fixture = Fixture::invalid_reference_id();
    let workspace = Workspace::with(&fixture);
    let status = bamsplit()
        .args(["chrom"])
        .arg(&workspace.input)
        .arg("--out-dir")
        .arg(workspace.path("out"))
        .args(["--quiet"])
        .status()
        .expect("runs");
    assert_eq!(status.code(), Some(3));
}

#[test]
fn inspect_writes_its_report_to_stdout() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);

    let output = bamsplit()
        .args(["inspect"])
        .arg(&workspace.input)
        .args(["--full", "--quiet"])
        .output()
        .expect("runs");

    assert_eq!(output.status.code(), Some(0));
    let text = String::from_utf8_lossy(&output.stdout);
    assert!(text.contains("sort order       coordinate"), "{text}");
    assert!(text.contains("records          12"), "{text}");
    assert!(text.contains("unplaced-unmapped 2"), "{text}");
}

#[test]
fn inspect_json_is_parseable() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);

    let output = bamsplit()
        .args(["inspect"])
        .arg(&workspace.input)
        .args(["--full", "--json", "--quiet"])
        .output()
        .expect("runs");

    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("valid JSON on stdout");
    assert_eq!(value["reference_count"], 3);
    assert_eq!(value["scan"]["total_records"], 12);
}

#[test]
fn region_splits_by_generated_windows() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    let output = bamsplit()
        .args(["region"])
        .arg(&workspace.input)
        .args(["--window-size", "5000"])
        .arg("--out-dir")
        .arg(&out)
        .args(["--compression-level", "1", "--emit-empty"])
        .output()
        .expect("runs");

    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    // chr1 is 10 000 bases, so it tiles into two 5 000-base windows. The second
    // holds no reads, which is why `--emit-empty` is needed to materialize it.
    assert!(out.join("chr1%3A1-5000.bam").exists());
    assert!(out.join("chr1%3A5001-10000.bam").exists());
    assert!(out.join("chrX%3A1-5000.bam").exists());
}

#[test]
fn an_unknown_region_option_value_is_an_argument_error() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);

    for (option, value) in [
        ("--feature", "nonesuch"),
        ("--assignment", "nonesuch"),
        ("--type", "7"),
        ("--format", "vcf"),
        ("--alignment-geometry", "nonesuch"),
        ("--missing-feature", "nonesuch"),
    ] {
        let status = bamsplit()
            .args(["region"])
            .arg(&workspace.input)
            .args(["--window-size", "5000"])
            .arg("--out-dir")
            .arg(workspace.path("out"))
            .args([option, value, "--quiet"])
            .status()
            .expect("runs");
        assert_eq!(status.code(), Some(2), "{option} {value}");
    }
}

#[test]
fn region_requires_exactly_one_source() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let status = bamsplit()
        .args(["region"])
        .arg(&workspace.input)
        .arg("--out-dir")
        .arg(workspace.path("out"))
        .arg("--quiet")
        .status()
        .expect("runs");
    assert_eq!(
        status.code(),
        Some(2),
        "neither --regions nor --window-size"
    );
}

#[test]
fn the_filename_template_is_applied() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    let status = bamsplit()
        .args(["chrom"])
        .arg(&workspace.input)
        .arg("--out-dir")
        .arg(&out)
        .args([
            "--filename-template",
            "sample1_{key}",
            "--compression-level",
            "1",
            "--quiet",
        ])
        .status()
        .expect("runs");

    assert_eq!(status.code(), Some(0));
    assert!(out.join("sample1_chr1.bam").exists());
    assert!(out.join("sample1_unmapped.bam").exists());
}

#[test]
fn an_invalid_filename_template_is_an_argument_error() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let status = bamsplit()
        .args(["chrom"])
        .arg(&workspace.input)
        .arg("--out-dir")
        .arg(workspace.path("out"))
        .args(["--filename-template", "../{key}", "--quiet"])
        .status()
        .expect("runs");
    assert_eq!(status.code(), Some(2));
}

#[test]
fn no_pg_suppresses_the_program_record() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    bamsplit()
        .args(["chrom"])
        .arg(&workspace.input)
        .arg("--out-dir")
        .arg(&out)
        .args(["--no-pg", "--compression-level", "1", "--quiet"])
        .status()
        .expect("runs");

    let (header, _) = read_bam(&out.join("chr1.bam"));
    let text =
        String::from_utf8_lossy(&header.to_text().expect("the header re-serializes")).into_owned();
    assert!(
        !text.contains("ID:bamsplit"),
        "`--no-pg` must not add an @PG record:\n{text}"
    );
}

#[test]
fn the_program_record_is_added_and_chained_by_default() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    bamsplit()
        .args(["chrom"])
        .arg(&workspace.input)
        .arg("--out-dir")
        .arg(&out)
        .args(["--compression-level", "1", "--quiet"])
        .status()
        .expect("runs");

    let (header, _) = read_bam(&out.join("chr1.bam"));
    let text =
        String::from_utf8_lossy(&header.to_text().expect("the header re-serializes")).into_owned();
    assert!(text.contains("@PG\tID:bamsplit\tPN:bamsplit"), "{text}");
    assert!(text.contains("CL:"), "{text}");
    // The new record must chain onto the last existing one rather than
    // replacing it.
    assert!(text.contains("PP:"), "{text}");
}
