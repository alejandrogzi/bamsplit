// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Malformed input, hostile names, and failure cleanup.
//!
//! Every case here must produce a typed error and an output directory that is
//! either complete or empty — never a directory of plausible-looking partial
//! BAMs.

use bamsplit_core::{ChromOptions, Classify as _, ExitCode, RunOptions};
use bamsplit_fixtures::{Damage, Fixture};

use crate::support::*;

fn expect_failure(fixture: &Fixture, expected: ExitCode) {
    let workspace = Workspace::with(fixture);
    let out = workspace.path("out");
    let error = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &run_options(),
    )
    .expect_err("must fail");
    assert_eq!(error.exit_code(), expected, "{error}");

    // Whatever went wrong, no output may survive.
    let remaining: Vec<_> = std::fs::read_dir(&out)
        .map(|entries| {
            entries
                .flatten()
                .map(|entry| entry.file_name().to_string_lossy().into_owned())
                .collect()
        })
        .unwrap_or_default();
    assert!(remaining.is_empty(), "left behind {remaining:?}");
}

#[test]
fn a_truncated_record_is_rejected() {
    expect_failure(
        &Fixture::coordinate_sorted().with_damage(Damage::TruncatedRecord),
        ExitCode::InvalidInput,
    );
}

#[test]
fn a_negative_block_size_is_rejected() {
    expect_failure(
        &Fixture::coordinate_sorted().with_damage(Damage::NegativeBlockSize),
        ExitCode::InvalidInput,
    );
}

#[test]
fn an_invalid_reference_id_is_rejected() {
    expect_failure(&Fixture::invalid_reference_id(), ExitCode::InvalidInput);
}

#[test]
fn an_invalid_mate_reference_id_is_rejected() {
    expect_failure(
        &Fixture::invalid_mate_reference_id(),
        ExitCode::InvalidInput,
    );
}

#[test]
fn a_missing_bgzf_eof_marker_does_not_lose_records() {
    // A missing marker is a warning sign, not corruption: every record is still
    // there, so `bamsplit` reads them all rather than refusing the file.
    let fixture = Fixture::coordinate_sorted().with_damage(Damage::MissingEof);
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");
    assert_lossless(&fixture, &out);

    // Every *output* does carry the marker.
    for stem in ["chr1", "chr2", "chrX", "unmapped"] {
        let bytes = std::fs::read(out.join(format!("{stem}.bam"))).expect("readable");
        let marker = bamsplit_core::output::BGZF_EOF;
        assert_eq!(&bytes[bytes.len() - marker.len()..], &marker[..], "{stem}");
    }
}

#[test]
fn a_malformed_auxiliary_section_only_fails_when_it_is_read() {
    // `chrom` routing never looks at auxiliary data, so a broken tag section is
    // irrelevant to it and the record is transferred untouched.
    let fixture = Fixture::malformed_tag();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("chrom");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("chrom does not read tags");
    assert_lossless(&fixture, &out);

    // Routing *on* a tag does read it, and then it is an input error.
    let error = bamsplit_core::split_by_tag(
        &workspace.input,
        &bamsplit_core::TagOptions::by_tag(workspace.path("tag"), "XX"),
        &run_options(),
    )
    .expect_err("must fail");
    assert_eq!(error.exit_code(), ExitCode::InvalidInput, "{error}");
}

#[test]
fn hostile_reference_names_produce_safe_files_inside_the_directory() {
    let fixture = Fixture::hostile_reference_names();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");

    assert_lossless(&fixture, &out);
    let canonical_out = out.canonicalize().expect("exists");
    for entry in std::fs::read_dir(&out).expect("listable").flatten() {
        let path = entry.path().canonicalize().expect("exists");
        assert!(path.starts_with(&canonical_out), "{path:?} escaped");
    }

    let manifest = read_manifest(&out);
    let encoded: Vec<&str> = manifest
        .outputs
        .iter()
        .map(|output| output.encoded_key.as_str())
        .collect();
    assert!(encoded.contains(&"chr1%2Falternate"), "{encoded:?}");
    assert!(encoded.contains(&"%2E%2E"), "{encoded:?}");
    assert!(encoded.contains(&"%43ON"), "{encoded:?}");

    // The manifest keeps the original bytes, so the mapping is reversible.
    for output in &manifest.outputs {
        let decoded = bamsplit_core::output::filename::decode(&output.encoded_key);
        assert_eq!(
            decoded.as_deref(),
            Some(output.logical_key.as_bytes()),
            "{}",
            output.encoded_key
        );
    }
}

#[test]
fn an_existing_output_is_refused_then_replaced_with_force() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("first run");

    let error = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &run_options(),
    )
    .expect_err("must refuse");
    assert_eq!(error.exit_code(), ExitCode::OutputConflict, "{error}");

    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &RunOptions {
            force: true,
            ..run_options()
        },
    )
    .expect("force replaces");
    assert_lossless(&fixture, &out);
}

#[test]
fn a_bam_with_no_references_is_handled() {
    let fixture = Fixture::no_references();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");
    let outputs = collect_outputs(&out);
    assert_eq!(outputs.len(), 1);
    assert_eq!(outputs["unmapped"].len(), 1);
}

#[test]
fn a_header_only_bam_produces_no_outputs_and_a_valid_manifest() {
    let fixture = Fixture::header_only();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");
    assert!(collect_outputs(&out).is_empty());
    read_manifest(&out).validate().expect("conservation holds");
}

#[test]
fn a_nonexistent_input_is_an_argument_error() {
    let directory = tempfile::tempdir().expect("temp dir");
    let error = bamsplit_core::split_by_chromosome(
        directory.path().join("absent.bam"),
        &ChromOptions::new(directory.path().join("out")),
        &run_options(),
    )
    .expect_err("must fail");
    assert_eq!(error.exit_code(), ExitCode::InvalidArguments, "{error}");
}

#[test]
fn a_non_bam_input_is_an_input_error() {
    let directory = tempfile::tempdir().expect("temp dir");
    let path = directory.path().join("not.bam");
    std::fs::write(&path, b"this is not a BAM at all").expect("seeded");
    let error = bamsplit_core::split_by_chromosome(
        &path,
        &ChromOptions::new(directory.path().join("out")),
        &run_options(),
    )
    .expect_err("must fail");
    assert_eq!(error.exit_code(), ExitCode::InvalidInput, "{error}");
}

#[test]
fn an_unknown_requested_reference_is_an_argument_error() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let error = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions {
            include: vec!["chr99".to_string()],
            ..ChromOptions::new(workspace.path("out"))
        },
        &run_options(),
    )
    .expect_err("must fail");
    assert_eq!(error.exit_code(), ExitCode::InvalidArguments, "{error}");
    assert!(
        error.to_string().contains("--ignore-missing-references"),
        "{error}"
    );
}

#[test]
fn an_invalid_compression_level_is_an_argument_error() {
    let fixture = Fixture::header_only();
    let workspace = Workspace::with(&fixture);
    let error = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(workspace.path("out")),
        &RunOptions {
            compression_level: 99,
            ..run_options()
        },
    )
    .expect_err("must fail");
    assert_eq!(error.exit_code(), ExitCode::InvalidArguments, "{error}");
}

#[test]
fn no_temporary_files_survive_a_successful_run() {
    let fixture = Fixture::unsorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");

    let leftovers: Vec<String> = std::fs::read_dir(&out)
        .expect("listable")
        .flatten()
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .filter(|name| name.contains(".part."))
        .collect();
    assert!(leftovers.is_empty(), "{leftovers:?}");
}
