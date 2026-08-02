// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The three engines must produce identical results.
//!
//! This is the property that makes engine selection an optimization rather than
//! a behavioural choice, so it is asserted directly on the bytes.

use bamsplit_core::engine::{EngineKind, IoBackend};
use bamsplit_core::index::IndexMode;
use bamsplit_core::{ChromOptions, RunOptions};
use bamsplit_fixtures::Fixture;

use crate::support::*;

fn split_with(
    workspace: &Workspace,
    name: &str,
    engine: EngineKind,
    run: RunOptions,
) -> std::path::PathBuf {
    let out = workspace.path(name);
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &RunOptions { engine, ..run },
    )
    .unwrap_or_else(|error| panic!("{engine} engine: {error}"));
    out
}

#[test]
fn stream_and_spool_agree_on_a_sorted_input() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let stream = split_with(&workspace, "stream", EngineKind::Stream, run_options());
    let spool = split_with(&workspace, "spool", EngineKind::Spool, run_options());
    assert_eq!(collect_outputs(&stream), collect_outputs(&spool));
    assert_lossless(&fixture, &stream);
}

#[test]
fn every_engine_agrees_including_the_indexed_one() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);

    // The indexed engine needs an index; build one from the fixture itself.
    let index_path = {
        let mut path = workspace.input.clone().into_os_string();
        path.push(".bai");
        std::path::PathBuf::from(path)
    };
    build_index(&workspace.input, &index_path);

    let stream = split_with(&workspace, "stream", EngineKind::Stream, run_options());
    let indexed = split_with(
        &workspace,
        "indexed",
        EngineKind::Indexed,
        RunOptions {
            threads: 4,
            ..run_options()
        },
    );
    let spool = split_with(&workspace, "spool", EngineKind::Spool, run_options());

    let baseline = collect_outputs(&stream);
    assert_eq!(baseline, collect_outputs(&indexed), "indexed differs");
    assert_eq!(baseline, collect_outputs(&spool), "spool differs");

    // Digests must agree too, which is a stronger claim than equal bodies:
    // it also pins the order records were written in.
    let digest = |directory: &std::path::Path| {
        read_manifest(directory)
            .outputs
            .into_iter()
            .map(|output| (output.logical_key, output.stats.raw_record_digest))
            .collect::<std::collections::BTreeMap<_, _>>()
    };
    assert_eq!(digest(&stream), digest(&indexed));
    assert_eq!(digest(&stream), digest(&spool));
}

/// Builds a BAI for a fixture by streaming it through the index builder.
fn build_index(bam: &std::path::Path, destination: &std::path::Path) {
    use bamsplit_core::bam::RawRecordReader;
    use bamsplit_core::bam::header::BamHeader;
    use bamsplit_core::bam::raw_record::VirtualPositionSource;
    use bamsplit_core::index::{AlignmentContext, IndexBuilder, IndexKind};

    let file = std::fs::File::open(bam).expect("readable");
    let mut reader = noodles_bgzf::io::Reader::new(file);
    let header = BamHeader::read_from(&mut reader).expect("valid header");
    let mut builder = IndexBuilder::new(IndexKind::Bai, header.reference_count());
    let mut records = RawRecordReader::new(reader);

    let mut start = VirtualPositionSource::virtual_position(records.get_ref());
    while let Some(record) = records.read_record().expect("valid record") {
        let context = match (
            record.reference_sequence_id().expect("valid"),
            record.alignment_start().expect("valid"),
            record.alignment_end().expect("valid"),
        ) {
            (Some(reference_id), Some(begin), Some(end)) => Some(AlignmentContext {
                reference_id: usize::try_from(reference_id).expect("in range"),
                start: u64::try_from(begin).expect("in range") + 1,
                end: u64::try_from(end).expect("in range"),
                mapped: !record.is_unmapped().expect("valid"),
            }),
            _ => None,
        };
        let end = VirtualPositionSource::virtual_position(records.get_ref());
        builder.add(context, start, end).expect("indexable");
        start = end;
    }
    builder
        .finish(destination)
        .expect("written")
        .expect("built");
}

#[test]
fn an_interleaved_input_falls_back_from_stream_to_spool() {
    let fixture = Fixture::interleaved();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");

    // `--engine stream` is honoured, the violation is detected, and the run
    // completes via the spool engine rather than producing partial files.
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &RunOptions {
            engine: EngineKind::Stream,
            ..run_options()
        },
    )
    .expect("split succeeds after the fallback");

    assert_lossless(&fixture, &out);
    let manifest = read_manifest(&out);
    assert_eq!(manifest.selected_engine, "spool");
    assert!(
        manifest
            .notes
            .iter()
            .any(|note| note.contains("ordering violation")),
        "{:?}",
        manifest.notes
    );
}

#[test]
fn a_reversed_but_still_grouped_input_streams_fine() {
    // Reversing a coordinate-sorted stream leaves every routing key contiguous,
    // so the streaming engine can handle it even though the header no longer
    // claims a sort order.
    let fixture = Fixture::unsorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &RunOptions {
            engine: EngineKind::Stream,
            ..run_options()
        },
    )
    .expect("split succeeds");
    assert_eq!(read_manifest(&out).selected_engine, "stream");
    assert_lossless(&fixture, &out);
}

#[test]
fn auto_selects_spool_for_an_unsorted_input() {
    let fixture = Fixture::unsorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");
    let manifest = read_manifest(&out);
    assert_eq!(manifest.selected_engine, "spool");
    assert!(manifest.engine_selection_reason.contains("SO:unsorted"));
    assert_lossless(&fixture, &out);
}

#[test]
fn a_query_name_sorted_input_is_handled() {
    let fixture = Fixture::query_name_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");
    assert_lossless(&fixture, &out);
}

#[test]
fn requesting_the_indexed_engine_without_an_index_is_an_error() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let error = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(workspace.path("out")),
        &RunOptions {
            engine: EngineKind::Indexed,
            ..run_options()
        },
    )
    .expect_err("must refuse");

    use bamsplit_core::Classify as _;
    assert_eq!(error.exit_code(), bamsplit_core::ExitCode::InvalidArguments);
    assert!(error.to_string().contains("samtools index"), "{error}");
}

#[test]
fn a_thread_budget_does_not_change_the_output() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let single = split_with(&workspace, "one", EngineKind::Stream, run_options());
    let many = split_with(
        &workspace,
        "many",
        EngineKind::Stream,
        RunOptions {
            threads: 8,
            ..run_options()
        },
    );
    // Byte-identical, not merely equivalent: compression is deterministic.
    for stem in ["chr1", "chr2", "chrX", "unmapped"] {
        let name = format!("{stem}.bam");
        assert_eq!(
            std::fs::read(single.join(&name)).expect("readable"),
            std::fs::read(many.join(&name)).expect("readable"),
            "{name} differs between thread counts"
        );
    }
}

#[test]
fn parallel_and_serial_spool_finalization_agree() {
    // The spool engine fans finalization out across the task budget. Whichever
    // order the pool happens to schedule in, the outputs and the manifest must
    // be identical — that is what makes the parallelism an optimization rather
    // than a behavioural change.
    let fixture = Fixture::interleaved();
    let workspace = Workspace::with(&fixture);

    let serial = split_with(
        &workspace,
        "serial",
        EngineKind::Spool,
        RunOptions {
            threads: 1,
            ..run_options()
        },
    );
    let parallel = split_with(
        &workspace,
        "parallel",
        EngineKind::Spool,
        RunOptions {
            threads: 8,
            ..run_options()
        },
    );

    assert_eq!(collect_outputs(&serial), collect_outputs(&parallel));
    for stem in ["chr1", "chr2", "chrX"] {
        let name = format!("{stem}.bam");
        assert_eq!(
            std::fs::read(serial.join(&name)).expect("readable"),
            std::fs::read(parallel.join(&name)).expect("readable"),
            "{name} differs between finalization strategies"
        );
    }

    let keys = |directory: &std::path::Path| {
        read_manifest(directory)
            .outputs
            .into_iter()
            .map(|output| (output.logical_key, output.stats.raw_record_digest))
            .collect::<Vec<_>>()
    };
    assert_eq!(
        keys(&serial),
        keys(&parallel),
        "manifest order must be stable"
    );
    assert_lossless(&fixture, &parallel);
}

#[test]
fn parallel_spool_finalization_rolls_back_on_failure() {
    let fixture = Fixture::interleaved();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    // Block one of the outputs so its commit fails while the others succeed.
    std::fs::create_dir_all(&out).expect("created");
    std::fs::write(out.join("chr2.bam"), b"in the way").expect("seeded");

    let error = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &RunOptions {
            engine: EngineKind::Spool,
            threads: 8,
            ..run_options()
        },
    )
    .expect_err("must fail");

    use bamsplit_core::Classify as _;
    assert_eq!(error.exit_code(), bamsplit_core::ExitCode::OutputConflict);

    // Only the pre-existing file may survive; nothing the run committed.
    let remaining: Vec<String> = std::fs::read_dir(&out)
        .expect("listable")
        .flatten()
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(remaining, ["chr2.bam"], "{remaining:?}");
    assert_eq!(
        std::fs::read(out.join("chr2.bam")).expect("readable"),
        b"in the way"
    );
}

#[test]
fn a_descriptor_budget_of_one_still_completes() {
    let fixture = Fixture::many_contigs(50);
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &RunOptions {
            max_open_files: 1,
            ..run_options()
        },
    )
    .expect("split succeeds");
    assert_eq!(collect_outputs(&out).len(), 50);
    assert_lossless(&fixture, &out);
}

#[test]
fn buffered_and_mapped_input_agree() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);

    let mut outputs = Vec::new();
    for (name, io) in [("buffered", IoBackend::Buffered), ("mmap", IoBackend::Mmap)] {
        let out = workspace.path(name);
        bamsplit_core::split_by_chromosome(
            &workspace.input,
            &ChromOptions::new(&out),
            &RunOptions {
                io,
                ..run_options()
            },
        )
        .unwrap_or_else(|error| panic!("{name}: {error}"));
        outputs.push(collect_outputs(&out));
    }
    assert_eq!(outputs[0], outputs[1]);
}

#[test]
fn a_reference_beyond_the_bai_limit_gets_a_csi() {
    let fixture = Fixture::beyond_bai_limit();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(&workspace.input, &ChromOptions::new(&out), &run_options())
        .expect("split succeeds");

    assert!(out.join("chrHuge.bam.csi").exists(), "CSI should be chosen");
    assert!(!out.join("chrHuge.bam.bai").exists());
    let manifest = read_manifest(&out);
    assert!(manifest.index_selection_reason.contains("beyond BAI"));
}

#[test]
fn forcing_bai_beyond_its_limit_is_rejected() {
    let fixture = Fixture::beyond_bai_limit();
    let workspace = Workspace::with(&fixture);
    let error = bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(workspace.path("out")),
        &RunOptions {
            index: IndexMode::Bai,
            ..run_options()
        },
    )
    .expect_err("must refuse");
    assert!(error.to_string().contains("--index csi"), "{error}");
}

#[test]
fn index_none_writes_no_index() {
    let fixture = Fixture::coordinate_sorted();
    let workspace = Workspace::with(&fixture);
    let out = workspace.path("out");
    bamsplit_core::split_by_chromosome(
        &workspace.input,
        &ChromOptions::new(&out),
        &RunOptions {
            index: IndexMode::None,
            ..run_options()
        },
    )
    .expect("split succeeds");

    let indexes: Vec<_> = std::fs::read_dir(&out)
        .expect("listable")
        .flatten()
        .filter(|entry| {
            let name = entry.file_name();
            let name = name.to_string_lossy();
            name.ends_with(".bai") || name.ends_with(".csi")
        })
        .collect();
    assert!(indexes.is_empty());
}
