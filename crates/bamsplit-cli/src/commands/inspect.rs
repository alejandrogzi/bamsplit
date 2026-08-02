// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `bamsplit inspect`.
//!
//! Writes its report to **stdout**, so it can be piped, while logs stay on
//! stderr. `--json` emits the same data machine-readably.

use std::io::Write as _;

use bamsplit_core::planning::RoutingMode;
use bamsplit_core::{InspectOptions, InspectReport, Result};

use crate::cli::{ByArg, GlobalOptions, InspectArgs};

/// Inspects a BAM.
///
/// # Errors
///
/// Propagates every [`bamsplit_core::Error`], and any failure writing to stdout.
pub fn run(globals: &GlobalOptions, args: InspectArgs) -> Result<()> {
    let feature =
        bamsplit_core::annotation::FeatureType::parse(&args.feature).ok_or_else(|| {
            bamsplit_core::error::ConfigError::OutOfRange {
                option: "--feature",
                constraint: "one of span, exon, intron, cds, utr, five-utr, three-utr",
                value: args.feature.clone(),
            }
        })?;
    let run = globals.to_run_options(crate::cli::command_line(), crate::cli::IndexArg::None, None);
    let options = InspectOptions {
        full: args.full,
        tag: args.tag,
        by: args.by.map(|by| match by {
            ByArg::Chrom => RoutingMode::Chrom,
            ByArg::Shard => RoutingMode::Shard,
            ByArg::Tag => RoutingMode::Tag,
            ByArg::Region => RoutingMode::Region {
                dense: false,
                regions: 0,
            },
        }),
        regions: args.regions,
        feature,
    };
    let report = bamsplit_core::inspect(&args.input, &options, &run)?;

    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    if args.json {
        let json = serde_json::to_string_pretty(&report)
            .unwrap_or_else(|error| format!("{{\"error\": \"{error}\"}}"));
        writeln!(out, "{json}")?;
    } else {
        render(&mut out, &report)?;
    }
    Ok(())
}

fn render<W: std::io::Write>(out: &mut W, report: &InspectReport) -> std::io::Result<()> {
    writeln!(out, "input            {}", report.input_path)?;
    if let Some(size) = report.input_size {
        writeln!(out, "size             {size} bytes")?;
    }
    writeln!(out, "valid BAM        {}", yes_no(report.valid_bam))?;
    writeln!(out, "header size      {} bytes", report.header_bytes)?;
    writeln!(out, "sort order       {}", report.sort_order)?;
    writeln!(out, "references       {}", report.reference_count)?;
    writeln!(out, "read groups      {}", report.read_groups.len())?;
    writeln!(out, "seekable         {}", yes_no(report.seekable))?;
    writeln!(out, "BGZF EOF marker  {}", yes_no(report.bgzf_eof_present))?;
    writeln!(
        out,
        "index            {}",
        report.index_type.as_deref().unwrap_or("none")
    )?;
    if let Some(consistent) = report.index_consistent {
        writeln!(out, "index consistent {}", yes_no(consistent))?;
    }
    writeln!(out, "mmap eligible    {}", yes_no(report.mmap_eligible))?;
    writeln!(out, "predicted outputs {}", report.predicted_output_count)?;
    writeln!(out, "recommended engine {}", report.recommended_engine)?;
    writeln!(out, "  reason         {}", report.engine_reason)?;

    if !report.references_exceeding_bai.is_empty() {
        writeln!(
            out,
            "references beyond BAI limits: {} (CSI will be used)",
            report.references_exceeding_bai.join(", ")
        )?;
    }
    if !report.unsafe_output_names.is_empty() {
        writeln!(out, "reference names needing filename encoding:")?;
        for (logical, encoded) in &report.unsafe_output_names {
            writeln!(out, "  {logical}  ->  {encoded}")?;
        }
    }

    writeln!(out)?;
    writeln!(out, "{:<24} {:>12}", "reference", "length")?;
    for (name, length) in &report.references {
        writeln!(out, "{name:<24} {length:>12}")?;
    }

    if let Some(annotation) = &report.annotation {
        writeln!(out)?;
        writeln!(out, "annotation       {}", annotation.path)?;
        writeln!(out, "  format         {}", annotation.format)?;
        if let Some(bed_type) = &annotation.bed_type {
            writeln!(out, "  BED type       {bed_type}")?;
        }
        writeln!(out, "  compression    {}", annotation.compression)?;
        writeln!(out, "  detected via   {}", annotation.detection_reason)?;
        if annotation.trailing_columns > 0 {
            writeln!(out, "  extra columns  {}", annotation.trailing_columns)?;
        }
        writeln!(out, "  records        {}", annotation.annotation_records)?;
        writeln!(out, "  feature        {}", annotation.feature)?;
        writeln!(out, "  usable regions {}", annotation.usable_regions)?;
        writeln!(out, "  segments       {}", annotation.derived_segments)?;
        writeln!(out, "  covered bases  {}", annotation.covered_bases)?;
        if annotation.missing_feature_records > 0 {
            writeln!(
                out,
                "  records missing `--feature {}`: {}",
                annotation.feature, annotation.missing_feature_records
            )?;
            if annotation.missing_blocks > 0 {
                writeln!(out, "    no block structure {}", annotation.missing_blocks)?;
            }
            if annotation.missing_coding_bounds > 0 {
                writeln!(
                    out,
                    "    no coding bounds   {}",
                    annotation.missing_coding_bounds
                )?;
            }
            if annotation.missing_strand > 0 {
                writeln!(out, "    no strand          {}", annotation.missing_strand)?;
            }
        }
        if annotation.span_as_exon_fallbacks > 0 {
            writeln!(
                out,
                "  span-as-exon fallbacks {}",
                annotation.span_as_exon_fallbacks
            )?;
        }
        if annotation.duplicate_names > 0 {
            writeln!(out, "  duplicate names {}", annotation.duplicate_names)?;
        }
        if annotation.generated_names > 0 {
            writeln!(out, "  generated names {}", annotation.generated_names)?;
        }
        if annotation.empty_regions > 0 {
            writeln!(out, "  empty regions   {}", annotation.empty_regions)?;
        }
        writeln!(
            out,
            "  overlapping region pairs {}",
            annotation.overlapping_regions
        )?;
        if !annotation.references_absent_from_bam.is_empty() {
            writeln!(
                out,
                "  annotation references absent from the BAM: {}",
                annotation.references_absent_from_bam.join(", ")
            )?;
        }
        if !annotation.references_absent_from_annotation.is_empty() {
            writeln!(
                out,
                "  BAM references absent from the annotation: {}",
                annotation.references_absent_from_annotation.join(", ")
            )?;
        }
    }

    if let Some(scan) = &report.scan {
        writeln!(out)?;
        writeln!(out, "records          {}", scan.total_records)?;
        writeln!(out, "  primary        {}", scan.primary_records)?;
        writeln!(out, "  secondary      {}", scan.secondary_records)?;
        writeln!(out, "  supplementary  {}", scan.supplementary_records)?;
        writeln!(out, "  mapped         {}", scan.mapped_records)?;
        writeln!(out, "  placed-unmapped   {}", scan.placed_unmapped_records)?;
        writeln!(
            out,
            "  unplaced-unmapped {}",
            scan.unplaced_unmapped_records
        )?;
        writeln!(out, "  duplicate      {}", scan.duplicate_records)?;
        writeln!(out, "  QC fail        {}", scan.qc_fail_records)?;
        writeln!(
            out,
            "coordinate order violations {}",
            scan.coordinate_order_violations
        )?;
        writeln!(out, "indexable output {}", yes_no(scan.estimated_indexable))?;
        writeln!(
            out,
            "estimated temporary storage {} bytes",
            scan.estimated_temporary_bytes
        )?;
        if let Some(cardinality) = scan.tag_cardinality {
            writeln!(out, "tag cardinality  {cardinality}")?;
        }
        if scan.malformed_tags > 0 {
            writeln!(out, "malformed auxiliary sections {}", scan.malformed_tags)?;
        }
        writeln!(out)?;
        writeln!(out, "{:<24} {:>12}", "reference", "records")?;
        for (name, count) in &scan.per_reference_counts {
            writeln!(out, "{name:<24} {count:>12}")?;
        }
    }
    Ok(())
}

const fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "no" }
}
