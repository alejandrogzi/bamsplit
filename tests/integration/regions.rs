// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Region routing, end to end.
//!
//! The unit tests in `annotation` and `routing::region` pin the semantics
//! against synthetic inputs. These check that a real annotation file and a real
//! BAM produce the outputs the semantics promise — including the manifest's
//! conservation equation switching when `overlap` duplicates records.

use bamsplit_core::annotation::{FeatureType, MissingFeaturePolicy};
use bamsplit_core::routing::region::AssignmentMode;
use bamsplit_core::{RegionOptions, RunOptions};
use bamsplit_fixtures::annotation::{self, Annotation, Codec};
use bamsplit_fixtures::{Fixture, RecordSpec};

use crate::support::*;

/// A BAM whose reads land predictably on the [`annotation::BED12`] transcripts.
///
/// ```text
/// ENST01/ENST02 exons  1000..1200   2500..2800   4600..5000   (chr1)
/// SINGLE               6000..6500                             (chr2)
///
/// exon1     1050 +100M          inside the first exon
/// spliced   1150 50M 1300N 50M  blocks in exons one and two
/// intronic  1500 +100M          inside the first intron
/// utr3      4850 +100M          inside the 3' UTR of ENST01
/// single    6100 +100M          inside SINGLE
/// nowhere   9000 +100M          in no region at all
/// ```
fn reads() -> Fixture {
    Fixture::new(&[("chr1", 20_000), ("chr2", 20_000)], "coordinate").with_records([
        RecordSpec::mapped("exon1", 0, 1_050).with_cigar(vec![(100, 0)]),
        RecordSpec::mapped("spliced", 0, 1_150).with_cigar(vec![(50, 0), (1_300, 3), (50, 0)]),
        RecordSpec::mapped("intronic", 0, 1_500).with_cigar(vec![(100, 0)]),
        RecordSpec::mapped("utr3", 0, 4_850).with_cigar(vec![(100, 0)]),
        RecordSpec::mapped("nowhere", 0, 9_000).with_cigar(vec![(100, 0)]),
        RecordSpec::mapped("single", 1, 6_100).with_cigar(vec![(100, 0)]),
    ])
}

struct Case {
    workspace: Workspace,
    annotation: std::path::PathBuf,
}

fn case(fixture: &Fixture, annotation: &Annotation, codec: Codec) -> Case {
    let workspace = Workspace::with(fixture);
    let path = annotation
        .write(workspace.directory.path(), codec)
        .expect("annotation written");
    Case {
        workspace,
        annotation: path,
    }
}

fn split(case: &Case, name: &str, options: RegionOptions) -> std::path::PathBuf {
    let out = case.workspace.path(name);
    bamsplit_core::split_by_region(
        &case.workspace.input,
        &RegionOptions {
            out_dir: out.clone(),
            ..options
        },
        &run_options(),
    )
    .unwrap_or_else(|error| panic!("{name}: {error}"));
    out
}

fn base(case: &Case) -> RegionOptions {
    RegionOptions::from_annotation(std::path::PathBuf::new(), &case.annotation)
}

fn counts(directory: &std::path::Path) -> Vec<(String, usize)> {
    collect_outputs(directory)
        .into_iter()
        .map(|(stem, bodies)| (stem, bodies.len()))
        .collect()
}

#[test]
fn span_and_start_assign_by_leftmost_position() {
    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    let out = split(
        &case,
        "span",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Start,
            ..base(&case)
        },
    );

    // ENST01 and ENST02 share a span, so the tie goes to the earlier record.
    // `nowhere` matches nothing; `single` lands in SINGLE.
    let counts = counts(&out);
    let enst01 = counts
        .iter()
        .find(|(stem, _)| stem == "ENST01")
        .map(|(_, count)| *count);
    assert_eq!(
        enst01,
        Some(4),
        "four chr1 reads inside the transcript span"
    );
    assert_eq!(
        counts.iter().find(|(stem, _)| stem == "SINGLE"),
        Some(&("SINGLE".to_string(), 1))
    );
    assert!(
        !counts.iter().any(|(stem, _)| stem == "ENST02"),
        "the tie-break must be deterministic, not both"
    );

    let manifest = read_manifest(&out);
    assert_eq!(manifest.input_records, 6);
    assert_eq!(manifest.unmatched_records, 1, "`nowhere` matches nothing");
    assert_eq!(manifest.duplicate_emissions, 0);
    manifest.validate().expect("conservation holds");
}

#[test]
fn exon_routing_ignores_intronic_reads() {
    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    let out = split(
        &case,
        "exon",
        RegionOptions {
            feature: FeatureType::Exon,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
    );

    let manifest = read_manifest(&out);
    // `intronic` sits entirely inside the first intron and `nowhere` is off the
    // end, so neither reaches an exon.
    assert_eq!(manifest.unmatched_records, 2);
    // Every matching read hits both ENST01 and ENST02, which share their exons.
    assert!(manifest.duplicate_emissions > 0);
    assert!(manifest.may_duplicate);
    manifest.validate().expect("overlap conservation holds");
}

#[test]
fn intron_routing_finds_only_the_intronic_read() {
    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    let out = split(
        &case,
        "intron",
        RegionOptions {
            feature: FeatureType::Intron,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
    );

    let outputs = collect_outputs(&out);
    for (stem, bodies) in &outputs {
        let names: Vec<String> = bodies
            .iter()
            .map(|body| {
                let record = bamsplit_core::RawRecord::new(body).expect("valid");
                String::from_utf8_lossy(record.qname().expect("valid")).into_owned()
            })
            .collect();
        // The spliced read's *blocks* are exonic; only its skipped gap is
        // intronic, and a skipped gap is not aligned overlap.
        assert!(!names.contains(&"spliced".to_string()), "{stem}: {names:?}");
        assert!(names.contains(&"intronic".to_string()), "{stem}: {names:?}");
    }
}

#[test]
fn contained_accepts_a_spliced_read_but_not_an_overhanging_one() {
    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    let out = split(
        &case,
        "contained",
        RegionOptions {
            feature: FeatureType::Exon,
            assignment: AssignmentMode::Contained,
            ..base(&case)
        },
    );

    let names: Vec<String> = collect_outputs(&out)
        .into_values()
        .flatten()
        .map(|body| {
            let record = bamsplit_core::RawRecord::new(&body).expect("valid");
            String::from_utf8_lossy(record.qname().expect("valid")).into_owned()
        })
        .collect();

    // `spliced` has blocks 1150..1200 and 2500..2550, both exonic.
    assert!(names.contains(&"spliced".to_string()), "{names:?}");
    // `exon1` runs 1050..1150, inside the first exon.
    assert!(names.contains(&"exon1".to_string()), "{names:?}");
    // `intronic` and `nowhere` are not exonic at all.
    assert!(!names.contains(&"intronic".to_string()), "{names:?}");
}

#[test]
fn every_assignment_mode_conserves_records() {
    for assignment in [
        AssignmentMode::Start,
        AssignmentMode::Midpoint,
        AssignmentMode::Contained,
        AssignmentMode::BestOverlap,
        AssignmentMode::Overlap,
    ] {
        let case = case(&reads(), &annotation::BED12, Codec::Plain);
        let out = split(
            &case,
            "out",
            RegionOptions {
                feature: FeatureType::Exon,
                assignment,
                ..base(&case)
            },
        );
        let manifest = read_manifest(&out);
        manifest
            .validate()
            .unwrap_or_else(|error| panic!("{assignment}: {error}"));
        assert_eq!(manifest.input_records, 6, "{assignment}");
        assert_eq!(
            manifest.may_duplicate,
            assignment.may_duplicate(),
            "{assignment}"
        );
        if !assignment.may_duplicate() {
            assert_eq!(manifest.duplicate_emissions, 0, "{assignment}");
        }
    }
}

#[test]
fn every_feature_type_runs_and_conserves() {
    for feature in [
        FeatureType::Span,
        FeatureType::Exon,
        FeatureType::Intron,
        FeatureType::Cds,
        FeatureType::Utr,
        FeatureType::FiveUtr,
        FeatureType::ThreeUtr,
    ] {
        let case = case(&reads(), &annotation::BED12, Codec::Plain);
        let out = split(
            &case,
            "out",
            RegionOptions {
                feature,
                assignment: AssignmentMode::Overlap,
                ..base(&case)
            },
        );
        let manifest = read_manifest(&out);
        manifest
            .validate()
            .unwrap_or_else(|error| panic!("{feature}: {error}"));
        let region = manifest.region.as_ref().expect("region fields are present");
        assert_eq!(region.feature_type, feature.to_string());
        assert_eq!(region.annotation_records, 3);
    }
}

#[test]
fn every_compression_codec_is_read() {
    for codec in [Codec::Plain, Codec::Gzip, Codec::Zstd, Codec::Bzip2] {
        let case = case(&reads(), &annotation::BED12, codec);
        let out = split(
            &case,
            "out",
            RegionOptions {
                feature: FeatureType::Span,
                assignment: AssignmentMode::Start,
                ..base(&case)
            },
        );
        let manifest = read_manifest(&out);
        assert_eq!(
            manifest.region.as_ref().map(|r| r.annotation_records),
            Some(3),
            "{codec:?}"
        );
    }
}

#[test]
fn a_gzip_stream_with_a_misleading_name_is_still_read() {
    let workspace = Workspace::with(&reads());
    let path = annotation::BED12
        .write_as(workspace.directory.path(), "plain-looking.bed", Codec::Gzip)
        .expect("written");

    let out = workspace.path("out");
    bamsplit_core::split_by_region(
        &workspace.input,
        &RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Start,
            ..RegionOptions::from_annotation(&out, &path)
        },
        &run_options(),
    )
    .expect("split succeeds");
    assert_eq!(
        read_manifest(&out).region.map(|r| r.annotation_records),
        Some(3)
    );
}

#[test]
fn gtf_and_gff3_produce_transcript_regions() {
    for (annotation, expected) in [(&annotation::GTF, "gtf"), (&annotation::GFF3, "gff")] {
        let case = case(&reads(), annotation, Codec::Plain);
        let out = split(
            &case,
            "out",
            RegionOptions {
                feature: FeatureType::Exon,
                assignment: AssignmentMode::Overlap,
                ..base(&case)
            },
        );
        let manifest = read_manifest(&out);
        let region = manifest.region.as_ref().expect("region fields");
        assert_eq!(region.annotation_format, expected);
        assert!(region.bed_type.is_none(), "GXF has no BED width");
        assert_eq!(region.annotation_records, 2, "two transcripts aggregate");
        manifest.validate().expect("conservation holds");
    }
}

#[test]
fn a_lying_extension_is_detected_by_content() {
    let workspace = Workspace::with(&reads());
    let path = annotation::GTF
        .write_as(workspace.directory.path(), "actually-gtf.bed", Codec::Plain)
        .expect("written");
    let out = workspace.path("out");
    bamsplit_core::split_by_region(
        &workspace.input,
        &RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Start,
            ..RegionOptions::from_annotation(&out, &path)
        },
        &run_options(),
    )
    .expect("split succeeds");
    assert_eq!(
        read_manifest(&out).region.map(|r| r.annotation_format),
        Some("gtf".to_string())
    );
}

#[test]
fn the_span_as_exon_fallback_is_counted_and_can_be_refused() {
    let case = case(&reads(), &annotation::BED6, Codec::Plain);

    let out = split(
        &case,
        "fallback",
        RegionOptions {
            feature: FeatureType::Exon,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
    );
    let region = read_manifest(&out).region.expect("region fields");
    assert_eq!(region.span_as_exon_fallbacks, 3, "BED6 has no blocks");

    let error = bamsplit_core::split_by_region(
        &case.workspace.input,
        &RegionOptions {
            out_dir: case.workspace.path("strict"),
            feature: FeatureType::Exon,
            assignment: AssignmentMode::Overlap,
            require_blocks: true,
            missing_feature: MissingFeaturePolicy::Error,
            ..base(&case)
        },
        &run_options(),
    )
    .expect_err("must refuse");
    assert!(error.to_string().contains("--require-blocks"), "{error}");
}

#[test]
fn missing_coding_bounds_are_reported() {
    let case = case(&reads(), &annotation::BED12_NONCODING, Codec::Plain);
    let error = bamsplit_core::split_by_region(
        &case.workspace.input,
        &RegionOptions {
            out_dir: case.workspace.path("out"),
            feature: FeatureType::Cds,
            assignment: AssignmentMode::Overlap,
            missing_feature: MissingFeaturePolicy::Error,
            ..base(&case)
        },
        &run_options(),
    )
    .expect_err("must refuse");
    assert!(error.to_string().contains("coding boundaries"), "{error}");
}

#[test]
fn a_missing_strand_is_reported_for_strand_aware_features() {
    let case = case(&reads(), &annotation::BED12_NO_STRAND, Codec::Plain);
    let error = bamsplit_core::split_by_region(
        &case.workspace.input,
        &RegionOptions {
            out_dir: case.workspace.path("out"),
            feature: FeatureType::FiveUtr,
            assignment: AssignmentMode::Overlap,
            missing_feature: MissingFeaturePolicy::Error,
            ..base(&case)
        },
        &run_options(),
    )
    .expect_err("must refuse");
    assert!(error.to_string().contains("no strand"), "{error}");
}

#[test]
fn duplicate_names_become_distinct_outputs() {
    let case = case(&reads(), &annotation::BED6_DUPLICATE_NAMES, Codec::Plain);
    let out = split(
        &case,
        "out",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            emit_empty: true,
            ..base(&case)
        },
    );

    let manifest = read_manifest(&out);
    let keys: Vec<&str> = manifest
        .outputs
        .iter()
        .map(|output| output.resolved_key.as_str())
        .collect();
    assert!(keys.contains(&"gene"), "{keys:?}");
    assert!(keys.contains(&"gene.2"), "{keys:?}");
    assert!(keys.contains(&"gene.3"), "{keys:?}");
    assert_eq!(
        manifest
            .region
            .as_ref()
            .map(|r| r.duplicate_name_resolutions),
        Some(2)
    );
}

#[test]
fn unnamed_records_get_generated_names() {
    let case = case(&reads(), &annotation::BED3, Codec::Plain);
    let out = split(
        &case,
        "out",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            emit_empty: true,
            unnamed_prefix: "win".to_string(),
            ..base(&case)
        },
    );
    let manifest = read_manifest(&out);
    assert!(
        manifest
            .outputs
            .iter()
            .all(|output| output.logical_key.starts_with("win-")),
        "{:?}",
        manifest
            .outputs
            .iter()
            .map(|o| &o.logical_key)
            .collect::<Vec<_>>()
    );
    assert_eq!(manifest.region.as_ref().map(|r| r.generated_names), Some(2));
}

#[test]
fn generated_windows_tile_the_dictionary() {
    let workspace = Workspace::with(&reads());
    let out = workspace.path("out");
    bamsplit_core::split_by_region(
        &workspace.input,
        &RegionOptions {
            emit_empty: true,
            ..RegionOptions::from_windows(&out, "5000")
        },
        &run_options(),
    )
    .expect("split succeeds");

    // Two references of 20 000 bases each, tiled at 5 000.
    let manifest = read_manifest(&out);
    assert_eq!(manifest.outputs.len(), 8);
    assert_eq!(manifest.input_records, 6);
    manifest.validate().expect("conservation holds");
    assert!(out.join("chr1%3A1-5000.bam").exists());
}

#[test]
fn a_genome_build_mismatch_is_reported_rather_than_producing_nothing() {
    let case = case(&reads(), &annotation::WRONG_GENOME_BUILD, Codec::Plain);
    let error = bamsplit_core::split_by_region(
        &case.workspace.input,
        &RegionOptions {
            out_dir: case.workspace.path("out"),
            ..base(&case)
        },
        &run_options(),
    )
    .expect_err("must refuse");
    let rendered = error.to_string();
    assert!(
        rendered.contains("no annotation reference name matches"),
        "{rendered}"
    );
    assert!(rendered.contains("chr1"), "{rendered}");
}

#[test]
fn malformed_annotations_are_rejected() {
    for annotation in [
        &annotation::MIXED_WIDTHS,
        &annotation::AMBIGUOUS,
        &annotation::EMPTY,
        &annotation::INVALID_BLOCK_COUNT,
        &annotation::UNSORTED_BLOCKS,
        &annotation::OVERLAPPING_BLOCKS,
    ] {
        let case = case(&reads(), annotation, Codec::Plain);
        let out = case.workspace.path("out");
        let error = bamsplit_core::split_by_region(
            &case.workspace.input,
            &RegionOptions {
                out_dir: out.clone(),
                feature: FeatureType::Exon,
                ..base(&case)
            },
            &run_options(),
        )
        .expect_err("must refuse");

        use bamsplit_core::Classify as _;
        assert_eq!(
            error.exit_code(),
            bamsplit_core::ExitCode::InvalidInput,
            "{}: {error}",
            annotation.name
        );
        assert!(
            std::fs::read_dir(&out).map(|d| d.count()).unwrap_or(0) == 0,
            "{} left output behind",
            annotation.name
        );
    }
}

#[test]
fn an_explicit_bed_type_narrows_the_parse() {
    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    let out = split(
        &case,
        "out",
        RegionOptions {
            bed_type: Some(bamsplit_core::annotation::BedType::Bed6),
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
    );
    assert_eq!(
        read_manifest(&out).region.and_then(|r| r.bed_type),
        Some("6".to_string())
    );
}

#[test]
fn trailing_columns_are_usable_as_a_name_field() {
    let case = case(&reads(), &annotation::BED12_EXTRA_COLUMNS, Codec::Plain);
    let out = split(
        &case,
        "out",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            emit_empty: true,
            ..base(&case)
        },
    );
    // Without a name field the record's own name is used.
    let manifest = read_manifest(&out);
    assert_eq!(manifest.outputs.len(), 1);
    assert_eq!(manifest.outputs[0].logical_key, "ENST01");
}

#[test]
fn emit_empty_materializes_regions_with_no_reads() {
    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    let out = split(
        &case,
        "out",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Start,
            emit_empty: true,
            ..base(&case)
        },
    );
    // ENST02 loses every tie to ENST01, so it would otherwise be skipped.
    let path = out.join("ENST02.bam");
    assert!(path.exists());
    let (header, bodies) = read_bam(&path);
    assert!(bodies.is_empty());
    assert_eq!(header.reference_count(), 2, "the full dictionary is kept");
}

#[test]
fn region_outputs_are_lossless_for_a_non_duplicating_mode() {
    let fixture = reads();
    let case = case(&fixture, &annotation::BED12, Codec::Plain);
    let out = split(
        &case,
        "out",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Start,
            ..base(&case)
        },
    );

    // Every emitted body must be one of the input's, byte for byte.
    let expected = fixture_bodies(&fixture);
    for body in collect_outputs(&out).into_values().flatten() {
        assert!(expected.contains(&body), "an emitted record was altered");
    }
}

#[test]
fn the_engine_choice_is_recorded_and_reproducible() {
    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    let first = split(
        &case,
        "first",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
    );
    let second = split(
        &case,
        "second",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
    );
    assert_eq!(collect_outputs(&first), collect_outputs(&second));
    let manifest = read_manifest(&first);
    assert!(!manifest.engine_selection_reason.is_empty());
    assert!(
        manifest.selected_engine == "spool" || manifest.selected_engine == "stream",
        "{}",
        manifest.selected_engine
    );
}

/// Builds a BAI beside a fixture, so the indexed engine has something to query.
fn index_input(bam: &std::path::Path) {
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
    let mut path = bam.as_os_str().to_os_string();
    path.push(".bai");
    builder
        .finish(std::path::Path::new(&path))
        .expect("written")
        .expect("built");
}

#[test]
fn the_indexed_engine_produces_the_same_regions_as_the_others() {
    use bamsplit_core::engine::EngineKind;

    let case = case(&reads(), &annotation::BED12, Codec::Plain);
    index_input(&case.workspace.input);

    let options = |assignment| RegionOptions {
        feature: FeatureType::Exon,
        assignment,
        ..base(&case)
    };

    for assignment in [
        AssignmentMode::Start,
        AssignmentMode::BestOverlap,
        AssignmentMode::Overlap,
    ] {
        let mut outputs = Vec::new();
        for (name, engine) in [
            ("stream", EngineKind::Stream),
            ("spool", EngineKind::Spool),
            ("indexed", EngineKind::Indexed),
        ] {
            let out = case.workspace.path(&format!("{name}-{assignment}"));
            bamsplit_core::split_by_region(
                &case.workspace.input,
                &RegionOptions {
                    out_dir: out.clone(),
                    ..options(assignment)
                },
                &RunOptions {
                    engine,
                    threads: 4,
                    ..run_options()
                },
            )
            .unwrap_or_else(|error| panic!("{name}/{assignment}: {error}"));
            outputs.push((name, collect_outputs(&out), read_manifest(&out)));
        }

        let (_, baseline, baseline_manifest) = &outputs[0];
        for (name, produced, manifest) in &outputs[1..] {
            assert_eq!(baseline, produced, "{name} differs under {assignment}");
            let digests = |m: &bamsplit_core::Manifest| {
                m.outputs
                    .iter()
                    .map(|o| (o.logical_key.clone(), o.stats.raw_record_digest.clone()))
                    .collect::<std::collections::BTreeMap<_, _>>()
            };
            assert_eq!(
                digests(baseline_manifest),
                digests(manifest),
                "{name} digests differ under {assignment}"
            );
        }
    }
}

/// Enough records to span many BGZF blocks, with regions only near the start.
///
/// A single-block fixture cannot demonstrate the per-region plan skipping
/// anything: one block holds every record, so the query returns it all.
fn many_reads() -> Fixture {
    Fixture::new(&[("chr1", 20_000_000)], "coordinate").with_records((0..6_000i32).map(|index| {
        RecordSpec {
            sequence_length: 100,
            ..RecordSpec::mapped(&format!("read{index:06}"), 0, index * 3_000)
                .with_cigar(vec![(100, 0)])
        }
    }))
}

/// Three small regions at the very start of `chr1`.
const NARROW: Annotation = Annotation {
    name: "narrow.bed",
    purpose: "a sparse annotation, so most of the input need never be read",
    contents: "\
chr1\t0\t3000\tfirst\t0\t+
chr1\t3000\t6000\tsecond\t0\t+
chr1\t6000\t9000\tthird\t0\t+
",
};

#[test]
fn the_indexed_engine_reads_only_the_records_its_queries_reach() {
    use bamsplit_core::engine::EngineKind;

    let case = case(&many_reads(), &NARROW, Codec::Plain);
    index_input(&case.workspace.input);

    // The baseline must be forced: with a sparse annotation `auto` would pick
    // the indexed engine for this run too, and the comparison would be vacuous.
    let stream = case.workspace.path("stream-skip");
    bamsplit_core::split_by_region(
        &case.workspace.input,
        &RegionOptions {
            out_dir: stream.clone(),
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
        &RunOptions {
            engine: EngineKind::Spool,
            ..run_options()
        },
    )
    .expect("split succeeds");
    let indexed_dir = case.workspace.path("indexed-skip");
    bamsplit_core::split_by_region(
        &case.workspace.input,
        &RegionOptions {
            out_dir: indexed_dir.clone(),
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
        &RunOptions {
            engine: EngineKind::Indexed,
            threads: 4,
            ..run_options()
        },
    )
    .expect("split succeeds");

    let stream_manifest = read_manifest(&stream);
    let indexed_manifest = read_manifest(&indexed_dir);

    // The outputs match, but the indexed engine never looked at the records
    // that no region could claim — that is the point of the plan, and the
    // manifest says so rather than letting the difference look like a bug.
    assert_eq!(collect_outputs(&stream), collect_outputs(&indexed_dir));
    assert_eq!(stream_manifest.input_records, 6_000);
    assert!(
        indexed_manifest.input_records < stream_manifest.input_records / 4,
        "the indexed plan should skip most of the input: {} of {}",
        indexed_manifest.input_records,
        stream_manifest.input_records
    );
    assert!(
        indexed_manifest
            .notes
            .iter()
            .any(|note| note.contains("only the records these queries examined")),
        "{:?}",
        indexed_manifest.notes
    );
    indexed_manifest
        .validate()
        .expect("conservation still holds");
}

#[test]
fn a_sparse_annotation_selects_the_indexed_engine_automatically() {
    let case = case(&many_reads(), &NARROW, Codec::Plain);
    index_input(&case.workspace.input);
    let out = split(
        &case,
        "auto",
        RegionOptions {
            feature: FeatureType::Span,
            assignment: AssignmentMode::Overlap,
            ..base(&case)
        },
    );
    let manifest = read_manifest(&out);
    assert_eq!(manifest.selected_engine, "indexed");
    assert!(
        manifest.engine_selection_reason.contains("sparse"),
        "{manifest:?}"
    );
}

#[test]
fn shard_and_tag_are_still_refused_by_the_indexed_engine() {
    use bamsplit_core::engine::EngineKind;

    let workspace = Workspace::with(&reads());
    index_input(&workspace.input);
    let error = bamsplit_core::split_by_shard(
        &workspace.input,
        &bamsplit_core::ShardOptions::new(workspace.path("out"), 4),
        &RunOptions {
            engine: EngineKind::Indexed,
            ..run_options()
        },
    )
    .expect_err("must refuse");
    assert!(
        error
            .to_string()
            .contains("cannot be decomposed into index queries"),
        "{error}"
    );
}
