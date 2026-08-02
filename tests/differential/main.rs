// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Differential tests against `samtools`, the external correctness oracle.
//!
//! `bamsplit` and `samtools` share no code, so agreement between them is real
//! evidence rather than a tautology. Each test compares:
//!
//! * record counts (`samtools view -c`);
//! * record-body digests, normalized so the expected `@PG` difference does not
//!   matter;
//! * structural validity (`samtools quickcheck`);
//! * index usability (`samtools idxstats`, and a coordinate query).
//!
//! # Skipping
//!
//! `samtools` is not required to build or test `bamsplit`, so without it these
//! tests skip loudly rather than failing: a red suite on a laptop that simply
//! lacks an optional tool trains people to ignore red suites.
//!
//! `--features differential` inverts that. It is the CI switch, and it says "I
//! have installed samtools and I expect these to run" — so a missing binary
//! becomes a failure instead of a skip. Otherwise a broken CI image would turn
//! the whole suite into a no-op that still reports green.

use std::path::Path;
use std::process::Command;

use bamsplit_core::{ChromOptions, RunOptions};
use bamsplit_fixtures::Fixture;

/// Whether a usable `samtools` is on `PATH`.
fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .output()
        .is_ok_and(|output| output.status.success())
}

/// Skips the calling test unless `samtools` is present.
///
/// Under `--features differential` a missing `samtools` is a failure, not a
/// skip; see the module documentation.
macro_rules! require_samtools {
    () => {
        if !samtools_available() {
            assert!(
                !cfg!(feature = "differential"),
                "`--features differential` was requested but samtools is not on PATH"
            );
            eprintln!(
                "SKIP: samtools is not on PATH; run `cargo test --features differential` on a \
                 machine with samtools to exercise the differential suite"
            );
            return;
        }
    };
}

fn samtools(arguments: &[&str]) -> String {
    let output = Command::new("samtools")
        .args(arguments)
        .output()
        .unwrap_or_else(|error| panic!("samtools {arguments:?}: {error}"));
    assert!(
        output.status.success(),
        "samtools {arguments:?} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout).into_owned()
}

/// The record bodies `samtools` sees, as normalized SAM text.
///
/// Only the eleven mandatory fields are compared. Auxiliary data is excluded
/// because `samtools view` renders it in its own order, which would make the
/// comparison test the renderer rather than the split; the in-process integration
/// suite compares full record bytes instead.
fn samtools_records(bam: &Path, region: Option<&str>) -> Vec<String> {
    let path = bam.to_string_lossy().into_owned();
    let mut arguments = vec!["view", path.as_str()];
    if let Some(region) = region {
        arguments.push(region);
    }
    samtools(&arguments)
        .lines()
        .map(|line| line.split('\t').take(11).collect::<Vec<_>>().join("\t"))
        .collect()
}

fn split_fixture(fixture: &Fixture) -> (tempfile::TempDir, std::path::PathBuf, std::path::PathBuf) {
    let directory = tempfile::tempdir().expect("temp dir");
    let input = directory.path().join("input.bam");
    fixture.write(&input).expect("fixture written");
    let out = directory.path().join("out");
    bamsplit_core::split_by_chromosome(
        &input,
        &ChromOptions::new(&out),
        &RunOptions {
            compression_level: 1,
            command_line: "bamsplit differential".to_string(),
            ..RunOptions::default()
        },
    )
    .expect("split succeeds");
    (directory, input, out)
}

#[test]
fn samtools_accepts_every_output() {
    require_samtools!();
    let (_directory, _input, out) = split_fixture(&Fixture::coordinate_sorted());
    for stem in ["chr1", "chr2", "chrX", "unmapped"] {
        let bam = out.join(format!("{stem}.bam"));
        samtools(&["quickcheck", "-v", &bam.to_string_lossy()]);
    }
}

#[test]
fn per_chromosome_counts_match_samtools() {
    require_samtools!();
    let (_directory, input, out) = split_fixture(&Fixture::coordinate_sorted());
    // The input needs an index for a region query.
    samtools(&["index", &input.to_string_lossy()]);

    for stem in ["chr1", "chr2", "chrX"] {
        let expected = samtools(&["view", "-c", &input.to_string_lossy(), stem])
            .trim()
            .to_string();
        let actual = samtools(&[
            "view",
            "-c",
            &out.join(format!("{stem}.bam")).to_string_lossy(),
        ])
        .trim()
        .to_string();
        assert_eq!(actual, expected, "{stem} count differs");
    }
}

#[test]
fn per_chromosome_record_digests_match_samtools() {
    require_samtools!();
    let (_directory, input, out) = split_fixture(&Fixture::coordinate_sorted());
    samtools(&["index", &input.to_string_lossy()]);

    for stem in ["chr1", "chr2", "chrX"] {
        let expected = samtools_records(&input, Some(stem));
        let actual = samtools_records(&out.join(format!("{stem}.bam")), None);
        assert_eq!(actual, expected, "{stem} records differ");
    }
}

#[test]
fn generated_indexes_are_usable_by_samtools() {
    require_samtools!();
    let (_directory, _input, out) = split_fixture(&Fixture::coordinate_sorted());

    let bam = out.join("chr1.bam");
    // `idxstats` reads the index rather than the records, so it only succeeds if
    // the index `bamsplit` wrote is structurally valid.
    let stats = samtools(&["idxstats", &bam.to_string_lossy()]);
    let chr1 = stats
        .lines()
        .find(|line| line.starts_with("chr1\t"))
        .expect("chr1 row");
    let columns: Vec<&str> = chr1.split('\t').collect();
    assert_eq!(columns[1], "10000", "the full dictionary is retained");
    assert_eq!(columns[2], "4", "mapped records on chr1");
    assert_eq!(columns[3], "1", "the placed-unmapped record");

    // A coordinate query must return exactly the records in that window.
    let queried = samtools(&["view", "-c", &bam.to_string_lossy(), "chr1:1-250"]);
    assert_eq!(queried.trim(), "2");
}

#[test]
fn a_csi_index_is_usable_by_samtools() {
    require_samtools!();
    let (_directory, _input, out) = split_fixture(&Fixture::beyond_bai_limit());
    let bam = out.join("chrHuge.bam");
    assert!(out.join("chrHuge.bam.csi").exists());
    samtools(&["quickcheck", "-v", &bam.to_string_lossy()]);
    let stats = samtools(&["idxstats", &bam.to_string_lossy()]);
    assert!(stats.contains("chrHuge\t600000000\t1"), "{stats}");
}

#[test]
fn every_engine_agrees_with_samtools() {
    require_samtools!();
    use bamsplit_core::engine::EngineKind;

    let fixture = Fixture::coordinate_sorted();
    let directory = tempfile::tempdir().expect("temp dir");
    let input = directory.path().join("input.bam");
    fixture.write(&input).expect("fixture written");
    samtools(&["index", &input.to_string_lossy()]);

    for (name, engine) in [
        ("stream", EngineKind::Stream),
        ("spool", EngineKind::Spool),
        ("indexed", EngineKind::Indexed),
    ] {
        let out = directory.path().join(name);
        bamsplit_core::split_by_chromosome(
            &input,
            &ChromOptions::new(&out),
            &RunOptions {
                engine,
                threads: 4,
                compression_level: 1,
                command_line: format!("bamsplit differential {name}"),
                ..RunOptions::default()
            },
        )
        .unwrap_or_else(|error| panic!("{name}: {error}"));

        for stem in ["chr1", "chr2", "chrX"] {
            let expected = samtools_records(&input, Some(stem));
            let actual = samtools_records(&out.join(format!("{stem}.bam")), None);
            assert_eq!(actual, expected, "{name}/{stem} differs from samtools");
        }
        samtools(&["quickcheck", "-v", &out.join("chr1.bam").to_string_lossy()]);
    }
}

#[test]
fn the_long_cigar_record_survives_a_round_trip() {
    require_samtools!();
    let (_directory, input, out) = split_fixture(&Fixture::long_cigar());
    let expected = samtools_records(&input, None);
    let actual = samtools_records(&out.join("chr1.bam"), None);
    assert_eq!(actual, expected);
    samtools(&["quickcheck", "-v", &out.join("chr1.bam").to_string_lossy()]);
}

#[test]
fn every_auxiliary_type_survives_a_round_trip() {
    require_samtools!();
    let (_directory, input, out) = split_fixture(&Fixture::every_tag_type());
    // Here the auxiliary fields are the point, so compare the whole line.
    let expected = samtools(&["view", &input.to_string_lossy()]);
    let actual = samtools(&["view", &out.join("chr1.bam").to_string_lossy()]);
    assert_eq!(actual, expected);
}

#[test]
fn hostile_reference_names_round_trip_through_samtools() {
    require_samtools!();
    let (_directory, input, out) = split_fixture(&Fixture::hostile_reference_names());
    let mut total = 0usize;
    for entry in std::fs::read_dir(&out).expect("listable").flatten() {
        let path = entry.path();
        if path.extension().and_then(|extension| extension.to_str()) != Some("bam") {
            continue;
        }
        samtools(&["quickcheck", "-v", &path.to_string_lossy()]);
        total += samtools(&["view", "-c", &path.to_string_lossy()])
            .trim()
            .parse::<usize>()
            .expect("a count");
    }
    let expected: usize = samtools(&["view", "-c", &input.to_string_lossy()])
        .trim()
        .parse()
        .expect("a count");
    assert_eq!(total, expected);
}

#[test]
fn a_header_only_output_is_accepted_by_samtools() {
    require_samtools!();
    let fixture = Fixture::header_only();
    let directory = tempfile::tempdir().expect("temp dir");
    let input = directory.path().join("input.bam");
    fixture.write(&input).expect("fixture written");
    let out = directory.path().join("out");
    bamsplit_core::split_by_chromosome(
        &input,
        &ChromOptions {
            emit_empty: true,
            ..ChromOptions::new(&out)
        },
        &RunOptions {
            compression_level: 1,
            ..RunOptions::default()
        },
    )
    .expect("split succeeds");

    for stem in ["chr1", "chr2"] {
        let bam = out.join(format!("{stem}.bam"));
        samtools(&["quickcheck", "-v", &bam.to_string_lossy()]);
        assert_eq!(
            samtools(&["view", "-c", &bam.to_string_lossy()]).trim(),
            "0"
        );
    }
}
