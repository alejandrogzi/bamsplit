// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Losslessness and the manifest's conservation equations.

use bamsplit_core::bam::header::ReadGroupField;
use bamsplit_core::routing::chrom::UnplacedPolicy;
use bamsplit_core::routing::shard::ShardKeySource;
use bamsplit_core::{ChromOptions, ManifestFormat, ShardOptions, TagOptions};
use bamsplit_fixtures::Fixture;

use crate::support::*;

#[test]
fn a_chromosome_split_conserves_every_record_body() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    let report = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &run_options(),
    )
    .expect("split succeeds");

    assert_lossless(&fixture, &out);
    let manifest = read_manifest(&out);
    assert_eq!(manifest.input_records, fixture.records.len() as u64);
    assert_eq!(manifest.unique_emitted_records, manifest.input_records);
    assert_eq!(manifest.dropped_records, 0);
    assert_eq!(manifest.duplicate_emissions, 0);
    assert_eq!(report.created_outputs(), 4, "chr1, chr2, chrX, unmapped");
    manifest.validate().expect("conservation holds");
}

#[test]
fn placed_unmapped_records_stay_with_their_reference_by_default() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");

    let outputs = collect_outputs(&out);
    // `halfmapped` is UNMAPPED but carries chr1:400, so it belongs to chr1.
    assert_eq!(outputs["chr1"].len(), 5);
    assert_eq!(
        outputs["unmapped"].len(),
        2,
        "only the truly unplaced records"
    );
}

#[test]
fn placed_unmapped_records_can_be_diverted_to_the_unmapped_output() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions {
            placed_unmapped: bamsplit_core::routing::chrom::PlacedUnmapped::Unmapped,
            ..ChromOptions::new(&out)
        },
        &run_options(),
    )
    .expect("split succeeds");

    let outputs = collect_outputs(&out);
    assert_eq!(outputs["chr1"].len(), 4);
    assert_eq!(outputs["unmapped"].len(), 3);
    assert_lossless(&fixture, &out);
}

#[test]
fn dropping_unplaced_records_is_accounted_for() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions {
            unplaced: UnplacedPolicy::Drop,
            ..ChromOptions::new(&out)
        },
        &run_options(),
    )
    .expect("split succeeds");

    let manifest = read_manifest(&out);
    assert_eq!(manifest.dropped_records, 2);
    assert_eq!(
        manifest.unique_emitted_records + manifest.dropped_records,
        manifest.input_records
    );
    manifest.validate().expect("conservation holds");
    assert!(!out.join("unmapped.bam").exists());
}

#[test]
fn excluded_references_are_dropped_and_counted() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions {
            exclude: vec!["chrX".to_string()],
            ..ChromOptions::new(&out)
        },
        &run_options(),
    )
    .expect("split succeeds");

    assert!(!out.join("chrX.bam").exists());
    let manifest = read_manifest(&out);
    assert_eq!(manifest.dropped_records, 2, "the two chrX records");
    manifest.validate().expect("conservation holds");
}

#[test]
fn a_shard_split_conserves_records_and_keeps_templates_together() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    bamsplit_core::split_by_shard(
        &workspace.input,
        &ShardOptions::new(&out, 4),
        &run_options(),
    )
    .expect("split succeeds");

    assert_lossless(&fixture, &out);

    // Every record sharing a QNAME must be in exactly one shard.
    let mut shard_of: std::collections::HashMap<Vec<u8>, String> = std::collections::HashMap::new();
    for (stem, bodies) in collect_outputs(&out) {
        for body in bodies {
            let record = bamsplit_core::RawRecord::new(&body).expect("valid");
            let name = record.qname().expect("valid").to_vec();
            let previous = shard_of.insert(name.clone(), stem.clone());
            if let Some(previous) = previous {
                assert_eq!(
                    previous,
                    stem,
                    "{:?} was split across shards",
                    String::from_utf8_lossy(&name)
                );
            }
        }
    }
    read_manifest(&out).validate().expect("conservation holds");
}

#[test]
fn shard_assignment_is_reproducible() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);

    let run = |name: &str| {
        let out = workspace.path(name);
        bamsplit_core::split_by_shard(
            &workspace.input,
            &ShardOptions::new(&out, 8),
            &run_options(),
        )
        .expect("split succeeds");
        collect_outputs(&out)
    };
    assert_eq!(run("first"), run("second"));
}

#[test]
fn a_tag_split_by_read_group_field_conserves_records() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    bamsplit_core::split_by_tag(
        &workspace.input,
        &TagOptions::by_field(&out, ReadGroupField::Sample),
        &run_options(),
    )
    .expect("split succeeds");

    assert_lossless(&fixture, &out);
    let outputs = collect_outputs(&out);
    assert!(outputs.contains_key("sample0"));
    assert!(outputs.contains_key("sample1"));
    // Records with no RG, and the one naming a group the header does not
    // declare, land in the fallback output rather than vanishing.
    assert!(outputs.contains_key("no-tag"));
    read_manifest(&out).validate().expect("conservation holds");
}

#[test]
fn the_cardinality_guard_fires_and_leaves_nothing_behind() {
    let fixture = Fixture::high_tag_cardinality(50);
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    let error = bamsplit_core::split_by_tag(
        &workspace.input,
        &TagOptions {
            max_outputs: 5,
            ..TagOptions::by_tag(&out, "CB")
        },
        &run_options(),
    )
    .expect_err("must refuse");

    use bamsplit_core::Classify as _;
    assert_eq!(error.exit_code(), bamsplit_core::ExitCode::OutputConflict);
    assert!(
        error.to_string().contains("bamsplit shard --key tag:CB"),
        "{error}"
    );
    let remaining: Vec<_> = std::fs::read_dir(&out)
        .map(|entries| entries.flatten().map(|entry| entry.path()).collect())
        .unwrap_or_default();
    assert!(remaining.is_empty(), "{remaining:?}");
}

#[test]
fn the_limit_can_be_lifted_explicitly() {
    let fixture = Fixture::high_tag_cardinality(20);
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    bamsplit_core::split_by_tag(
        &workspace.input,
        &TagOptions {
            max_outputs: 1,
            allow_high_cardinality: true,
            ..TagOptions::by_tag(&out, "CB")
        },
        &run_options(),
    )
    .expect("split succeeds");

    assert_eq!(collect_outputs(&out).len(), 20);
    assert_lossless(&fixture, &out);
}

#[test]
fn a_tsv_manifest_is_directly_usable_by_a_workflow_manager() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &bamsplit_core::RunOptions {
            manifest: ManifestFormat::Both,
            ..run_options()
        },
    )
    .expect("split succeeds");

    let tsv = std::fs::read_to_string(out.join("bamsplit.manifest.tsv")).expect("written");
    let mut lines = tsv.lines();
    let header: Vec<&str> = lines.next().expect("header row").split('\t').collect();
    assert_eq!(
        &header[..5],
        ["logical_key", "resolved_key", "encoded_key", "bam", "index"]
    );
    let rows: Vec<&str> = lines.collect();
    assert_eq!(rows.len(), 4);
    for row in rows {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(columns.len(), header.len(), "{row}");
        assert!(out.join(columns[3]).exists(), "{}", columns[3]);
        assert!(out.join(columns[4]).exists(), "{}", columns[4]);
    }
}

#[test]
fn emit_empty_materializes_every_declared_reference() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);

    let without = workspace.path("without");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions {
            include: vec!["chr1".to_string(), "chr2".to_string(), "chrX".to_string()],
            ..ChromOptions::new(&without)
        },
        &run_options(),
    )
    .expect("split succeeds");

    let with = workspace.path("with");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions {
            emit_empty: true,
            include: vec!["chr1".to_string(), "chr2".to_string(), "chrX".to_string()],
            ..ChromOptions::new(&with)
        },
        &run_options(),
    )
    .expect("split succeeds");

    // Every reference has records in this fixture, so both runs create the same
    // set; the difference shows up on a fixture with an unused reference.
    let fixture = Fixture::header_only();
    let empty_workspace = Workspace::with(&fixture);
    let out = empty_workspace.path("out");
    bamsplit_core::split_by_chromosome(
        &empty_workspace.input,
        &ChromOptions {
            emit_empty: true,
            ..ChromOptions::new(&out)
        },
        &run_options(),
    )
    .expect("split succeeds");

    for stem in ["chr1", "chr2", "unmapped"] {
        let path = out.join(format!("{stem}.bam"));
        assert!(path.exists(), "{path:?}");
        let (header, bodies) = read_bam(&path);
        assert!(bodies.is_empty());
        assert_eq!(header.reference_count(), 2, "the full dictionary is kept");
    }
    let _ = (without, with);
}

#[test]
fn skipped_outputs_are_listed_in_the_manifest() {
    let fixture = Fixture::header_only();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");

    let manifest = read_manifest(&out);
    assert_eq!(manifest.input_records, 0);
    let skipped: Vec<&str> = manifest
        .outputs
        .iter()
        .filter(|output| output.skipped)
        .map(|output| output.logical_key.as_str())
        .collect();
    assert_eq!(skipped, ["chr1", "chr2", "unmapped"]);
    assert_eq!(
        std::fs::read_dir(&out).expect("listable").count(),
        1,
        "manifest only"
    );
}

#[test]
fn shard_key_sources_all_conserve_records() {
    for key in [
        ShardKeySource::QName,
        ShardKeySource::Record,
        ShardKeySource::Tag(*b"RG"),
        ShardKeySource::ReadGroupField(ReadGroupField::Sample),
    ] {
        let fixture = Fixture::coordinate_sorted();
        let workspace = Workspace::with(&fixture);
        let out = workspace.path("out");
        bamsplit_core::split_by_shard(
            &workspace.input,
            &ShardOptions {
                key: key.clone(),
                ..ShardOptions::new(&out, 4)
            },
            &run_options(),
        )
        .unwrap_or_else(|error| panic!("{}: {error}", key.as_string()));
        assert_lossless(&fixture, &out);
    }
}
