// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `bamsplit region`.

use bamsplit_core::annotation::{AnnotationFormat, BedType, FeatureType, MissingFeaturePolicy};
use bamsplit_core::bam::AlignmentGeometry;
use bamsplit_core::error::ConfigError;
use bamsplit_core::routing::region::AssignmentMode;
use bamsplit_core::{RegionOptions, Result};

use crate::cli::{GlobalOptions, RegionArgs};

/// Parses one of the free-form string options into its enum.
///
/// These stay strings rather than `clap` enums because several accept values
/// the derive cannot express: `--type` mixes the keyword `auto` with integers,
/// and `--format auto` shadows a real format name.
fn parse<T>(
    value: &str,
    option: &'static str,
    constraint: &'static str,
    parser: impl Fn(&str) -> Option<T>,
) -> std::result::Result<T, ConfigError> {
    parser(value).ok_or_else(|| ConfigError::OutOfRange {
        option,
        constraint,
        value: value.to_string(),
    })
}

/// Runs a per-region split.
///
/// # Errors
///
/// Propagates every [`bamsplit_core::Error`], including annotation detection
/// and parse failures.
pub fn run(globals: &GlobalOptions, args: RegionArgs) -> Result<()> {
    let format = match args.format.as_str() {
        "auto" => None,
        other => Some(parse(
            other,
            "--format",
            "one of auto, bed, gtf, gff",
            AnnotationFormat::parse,
        )?),
    };
    let bed_type = match args.bed_type.as_str() {
        "auto" => None,
        other => Some(parse(
            other,
            "--type",
            "one of auto, 3, 4, 5, 6, 8, 9, 12",
            BedType::parse,
        )?),
    };
    let feature = parse(
        &args.feature,
        "--feature",
        "one of span, exon, intron, cds, utr, five-utr, three-utr",
        FeatureType::parse,
    )?;
    let assignment = parse(
        &args.assignment,
        "--assignment",
        "one of start, midpoint, contained, best-overlap, overlap",
        AssignmentMode::parse,
    )?;
    let missing_feature = parse(
        &args.missing_feature,
        "--missing-feature",
        "one of skip, empty, error",
        MissingFeaturePolicy::parse,
    )?;
    let geometry = parse(
        &args.alignment_geometry,
        "--alignment-geometry",
        "one of blocks, span",
        |value| match value {
            "blocks" => Some(AlignmentGeometry::Blocks),
            "span" => Some(AlignmentGeometry::Span),
            _ => None,
        },
    )?;

    let run = globals.to_run_options(crate::cli::command_line(), args.index, None);
    let options = RegionOptions {
        out_dir: args.out_dir,
        regions: args.regions,
        window_size: args.window_size,
        format,
        bed_type,
        feature,
        assignment,
        geometry,
        name_field: args.name_field,
        unnamed_prefix: args.unnamed_prefix,
        missing_feature,
        require_blocks: args.require_blocks,
        emit_empty: args.emit_empty,
    };

    let report = bamsplit_core::split_by_region(&args.input, &options, &run)?;
    super::report_summary(&report);
    if let Some(region) = &report.manifest.region {
        tracing::info!(
            "{} annotation record(s) -> {} region(s), {} segment(s), {} ambiguous assignment(s)",
            region.annotation_records,
            report.manifest.outputs.len(),
            region.derived_segment_count,
            region.ambiguous_assignments
        );
        if region.missing_feature_records > 0 {
            tracing::warn!(
                "{} record(s) lacked the information `--feature {}` needs",
                region.missing_feature_records,
                region.feature_type
            );
        }
    }
    Ok(())
}
