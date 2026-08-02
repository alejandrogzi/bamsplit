// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Turning annotation records into logical regions.
//!
//! # One record, one output
//!
//! Each annotation record becomes exactly **one** output key. `--feature`
//! chooses which of that record's genomic segments participate in overlap
//! testing; it never splits a transcript into one BAM per exon. A transcript's
//! exons are the *segments* of a single region.
//!
//! ```text
//! BED12  chr1  1000  5000  ENST01  0  +  1200  4800  0,0,0  3  200,300,400  0,1500,3600
//!
//! --feature span     ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓   1000..5000
//! --feature exon     ▓▓        ▓▓▓          ▓▓▓▓            three exons
//! --feature intron     ▓▓▓▓▓▓▓▓   ▓▓▓▓▓▓▓▓▓▓                two introns
//! --feature cds       ▓        ▓▓▓          ▓▓▓             exonic ∩ 1200..4800
//! --feature utr      ▓                          ▓           exonic ∖ 1200..4800
//! ```
//!
//! All of them produce one region keyed `ENST01`.
//!
//! # What `genepred` does, and what this module has to add
//!
//! `GenePred::exons()` falls back to the whole span when a record has no block
//! structure, which is convenient but hides the distinction `--require-blocks`
//! needs. So block structure is checked here directly, and the fallback is
//! counted rather than silently applied.
//!
//! Likewise `coding_exons()` and the UTR accessors return an empty vector both
//! for "no coding bounds" and for "coding bounds that intersect nothing".
//! Those are different, and only the first should trigger
//! `--missing-feature`.

use std::collections::HashMap;

use genepred::GenePred;

use crate::annotation::detect::Detection;
use crate::error::AnnotationError;
use crate::interval::{Interval, merge_in_place, total_len};

/// Which genomic segments of an annotation record participate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FeatureType {
    /// The whole `start..end` span. Always available.
    #[default]
    Span,
    /// The record's exon blocks.
    Exon,
    /// The gaps between consecutive exons.
    Intron,
    /// The coding portions of exons — not the whole thick span.
    Cds,
    /// Every untranslated *exonic* segment.
    Utr,
    /// The strand-aware 5′ exonic UTR.
    FiveUtr,
    /// The strand-aware 3′ exonic UTR.
    ThreeUtr,
}

impl FeatureType {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Span => "span",
            Self::Exon => "exon",
            Self::Intron => "intron",
            Self::Cds => "cds",
            Self::Utr => "utr",
            Self::FiveUtr => "five-utr",
            Self::ThreeUtr => "three-utr",
        }
    }

    /// Parses an option value.
    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "span" => Self::Span,
            "exon" => Self::Exon,
            "intron" => Self::Intron,
            "cds" => Self::Cds,
            "utr" => Self::Utr,
            "five-utr" | "5utr" => Self::FiveUtr,
            "three-utr" | "3utr" => Self::ThreeUtr,
            _ => return None,
        })
    }

    /// Whether the feature needs block structure to be meaningful.
    #[must_use]
    pub const fn needs_blocks(self) -> bool {
        !matches!(self, Self::Span)
    }

    /// Whether the feature needs `thickStart`/`thickEnd`.
    #[must_use]
    pub const fn needs_coding_bounds(self) -> bool {
        matches!(self, Self::Cds | Self::Utr | Self::FiveUtr | Self::ThreeUtr)
    }

    /// Whether the feature needs a strand.
    #[must_use]
    pub const fn needs_strand(self) -> bool {
        matches!(self, Self::FiveUtr | Self::ThreeUtr)
    }
}

impl std::fmt::Display for FeatureType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// What to do with a record that lacks the information the feature needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MissingFeaturePolicy {
    /// Drop the record, and count it. The default.
    #[default]
    Skip,
    /// Keep it as a region with no segments.
    Empty,
    /// Fail the run.
    Error,
}

impl MissingFeaturePolicy {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Skip => "skip",
            Self::Empty => "empty",
            Self::Error => "error",
        }
    }

    /// Parses an option value.
    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "skip" => Self::Skip,
            "empty" => Self::Empty,
            "error" => Self::Error,
            _ => return None,
        })
    }
}

impl std::fmt::Display for MissingFeaturePolicy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A transcript's orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    /// `+`
    Forward,
    /// `-`
    Reverse,
    /// `.` or `?`
    Unknown,
}

impl From<genepred::Strand> for Strand {
    fn from(value: genepred::Strand) -> Self {
        match value {
            genepred::Strand::Forward => Self::Forward,
            genepred::Strand::Reverse => Self::Reverse,
            genepred::Strand::Unknown => Self::Unknown,
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Forward => "+",
            Self::Reverse => "-",
            Self::Unknown => ".",
        })
    }
}

/// How regions are derived and named.
#[derive(Debug, Clone)]
pub struct FeatureOptions {
    /// Which segments participate.
    pub feature: FeatureType,
    /// What to do about records missing the needed information.
    pub missing: MissingFeaturePolicy,
    /// Turn the BED span-as-single-exon fallback into an error.
    pub require_blocks: bool,
    /// An annotation attribute to prefer as the logical name.
    pub name_field: Option<Vec<u8>>,
    /// The prefix for generated names.
    pub unnamed_prefix: String,
}

impl Default for FeatureOptions {
    fn default() -> Self {
        Self {
            feature: FeatureType::Span,
            missing: MissingFeaturePolicy::Skip,
            require_blocks: false,
            name_field: None,
            unnamed_prefix: "region".to_string(),
        }
    }
}

/// The identity of one annotated region.
///
/// Ordering is by annotation input order, which is what every documented
/// tie-break falls back to. The name is included so the ordering is total even
/// if two records somehow shared an ordinal.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RegionKey {
    ordinal: u64,
    name: Vec<u8>,
}

impl RegionKey {
    /// Creates a key.
    #[must_use]
    pub fn new(ordinal: u64, name: Vec<u8>) -> Self {
        Self { ordinal, name }
    }

    /// The record's 0-based position in the annotation.
    #[must_use]
    pub const fn ordinal(&self) -> u64 {
        self.ordinal
    }

    /// The resolved, unique logical name.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }
}

impl crate::routing::RoutingKey for RegionKey {
    fn logical(&self) -> &[u8] {
        &self.name
    }

    fn manifest_fields(&self) -> Vec<(&'static str, serde_json::Value)> {
        vec![("annotation_ordinal", serde_json::Value::from(self.ordinal))]
    }
}

/// What is known about a region beyond its coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionMetadata {
    /// The name before duplicate resolution.
    pub original_name: Vec<u8>,
    /// Whether the name had to be invented.
    pub generated_name: bool,
    /// `0` for the first use of a name, `2` for `<name>.2`, and so on.
    pub duplicate_suffix: u32,
    /// Whether the whole span stood in for absent block structure.
    pub span_as_exon_fallback: bool,
    /// The transcript's orientation, when the input carries one.
    pub strand: Option<Strand>,
    /// How many segments the region ended up with.
    pub derived_segment_count: usize,
}

/// One annotation record, ready for routing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LogicalRegion {
    /// The output identity.
    pub key: RegionKey,
    /// The reference this region sits on.
    pub chrom: Vec<u8>,
    /// The segments that participate in overlap testing, sorted and merged.
    pub segments: Vec<Interval>,
    /// The smallest interval covering every segment.
    pub envelope: Interval,
    /// The record's 0-based position in the annotation.
    pub source_ordinal: u64,
    /// Everything else worth reporting.
    pub metadata: RegionMetadata,
}

impl LogicalRegion {
    /// The number of reference bases the segments cover.
    #[must_use]
    pub fn covered_bases(&self) -> i64 {
        total_len(&self.segments)
    }

    /// Whether the region has no segments at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }
}

/// What extraction observed, for the manifest.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExtractionStats {
    /// Annotation records read.
    pub annotation_records: u64,
    /// Segments across every region.
    pub derived_segments: u64,
    /// Records lacking the information the feature needed.
    pub missing_feature: u64,
    /// Records dropped because of that.
    pub skipped_records: u64,
    /// Regions that ended up with no segments.
    pub empty_regions: u64,
    /// Names that had to be invented.
    pub generated_names: u64,
    /// Duplicate names that were disambiguated.
    pub duplicate_names: u64,
    /// BED records whose whole span stood in for missing blocks.
    pub span_as_exon_fallbacks: u64,
}

/// Why a record could not produce its requested feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Missing {
    Blocks,
    Coding,
    Strand,
}

/// Turns annotation records into logical regions.
///
/// # Errors
///
/// Returns [`AnnotationError`] for invalid coordinates, invalid block
/// structure, or — under [`MissingFeaturePolicy::Error`] — a record that lacks
/// the information the feature needs.
pub fn extract_regions(
    records: &[GenePred],
    detection: &Detection,
    options: &FeatureOptions,
) -> Result<(Vec<LogicalRegion>, ExtractionStats), AnnotationError> {
    let mut stats = ExtractionStats {
        annotation_records: records.len() as u64,
        ..ExtractionStats::default()
    };
    let mut regions = Vec::with_capacity(records.len());
    let mut seen_names: HashMap<Vec<u8>, u32> = HashMap::with_capacity(records.len());

    for (ordinal, record) in records.iter().enumerate() {
        let ordinal = ordinal as u64;
        validate_coordinates(record, ordinal)?;

        let (original_name, generated) = resolve_name(record, ordinal, options);
        let outcome = derive_segments(record, detection, options, ordinal, &original_name)?;

        let (raw_segments, fallback) = match outcome {
            Ok(value) => value,
            Err(missing) => {
                stats.missing_feature += 1;
                match options.missing {
                    MissingFeaturePolicy::Skip => {
                        stats.skipped_records += 1;
                        continue;
                    }
                    MissingFeaturePolicy::Empty => (Vec::new(), false),
                    MissingFeaturePolicy::Error => {
                        return Err(missing_error(missing, &original_name, ordinal, options));
                    }
                }
            }
        };

        let mut segments: Vec<Interval> = raw_segments
            .into_iter()
            .map(|(start, end)| Interval::new(start as i64, end as i64))
            .collect();
        // Adjacent or overlapping segments of the *same* region are merged: two
        // abutting exons cover a contiguous stretch, and leaving them separate
        // would make overlap arithmetic double-count the boundary.
        merge_in_place(&mut segments);

        // `merge_in_place` leaves the segments sorted, so the envelope is just
        // the first start and the last end.
        let envelope = match (segments.first(), segments.last()) {
            (Some(first), Some(last)) => Interval::new(first.start, last.end),
            _ => Interval::empty_at(record.start as i64),
        };

        let (name, duplicate_suffix) = disambiguate(&mut seen_names, &original_name);
        if duplicate_suffix > 0 {
            stats.duplicate_names += 1;
        }
        if generated {
            stats.generated_names += 1;
        }
        if fallback {
            stats.span_as_exon_fallbacks += 1;
        }
        if segments.is_empty() {
            stats.empty_regions += 1;
        }
        stats.derived_segments += segments.len() as u64;

        regions.push(LogicalRegion {
            key: RegionKey::new(ordinal, name),
            chrom: record.chrom.clone(),
            metadata: RegionMetadata {
                original_name,
                generated_name: generated,
                duplicate_suffix,
                span_as_exon_fallback: fallback,
                strand: record.strand.map(Into::into),
                derived_segment_count: segments.len(),
            },
            envelope,
            segments,
            source_ordinal: ordinal,
        });
    }

    Ok((regions, stats))
}

/// The segments a record contributes, or why it cannot.
///
/// The outer `Result` is a hard error — malformed input. The inner `Result` is
/// the soft case that `--missing-feature` governs.
#[allow(clippy::type_complexity)]
fn derive_segments(
    record: &GenePred,
    detection: &Detection,
    options: &FeatureOptions,
    ordinal: u64,
    name: &[u8],
) -> Result<Result<(Vec<(u64, u64)>, bool), Missing>, AnnotationError> {
    let blocks = block_structure(record, ordinal, name)?;

    // `span` needs nothing beyond the record's own coordinates.
    if options.feature == FeatureType::Span {
        return Ok(Ok((vec![(record.start, record.end)], false)));
    }

    // Every other feature is exon-derived. A record with no block structure gets
    // the whole span as one exon — but only for BED, and only when the caller
    // has not asked for the stricter reading. An aggregated GTF/GFF record with
    // no exons means the aggregation produced nothing usable, which is an input
    // problem rather than a tolerable gap.
    let fallback = !blocks;
    if fallback {
        if options.require_blocks || detection.format.is_gxf() {
            return Ok(Err(Missing::Blocks));
        }
        if options.feature == FeatureType::Intron {
            // One exon has no introns; that is an empty region, not a failure.
            return Ok(Ok((Vec::new(), true)));
        }
    }

    if options.feature.needs_coding_bounds() && !has_coding_bounds(record) {
        return Ok(Err(Missing::Coding));
    }
    if options.feature.needs_strand()
        && !matches!(
            record.strand,
            Some(genepred::Strand::Forward | genepred::Strand::Reverse)
        )
    {
        return Ok(Err(Missing::Strand));
    }

    let segments = match options.feature {
        // Handled above.
        FeatureType::Span => unreachable!("span returns early"),
        FeatureType::Exon => record.exons(),
        FeatureType::Intron => record.introns(),
        FeatureType::Cds => record.coding_exons(),
        FeatureType::Utr => record.utr_exons(),
        FeatureType::FiveUtr => record.five_prime_utr(),
        FeatureType::ThreeUtr => record.three_prime_utr(),
    };
    Ok(Ok((segments, fallback)))
}

/// Whether a record carries usable block structure, validating it if so.
///
/// # Errors
///
/// Returns [`AnnotationError::InvalidBlocks`] for blocks that are unsorted,
/// overlapping, or zero-length — none of which can be interpreted, and all of
/// which would silently corrupt intron and CDS derivation.
fn block_structure(record: &GenePred, ordinal: u64, name: &[u8]) -> Result<bool, AnnotationError> {
    let (Some(count), Some(starts), Some(ends)) =
        (record.block_count, &record.block_starts, &record.block_ends)
    else {
        return Ok(false);
    };
    if count == 0 || starts.is_empty() || ends.is_empty() {
        return Ok(false);
    }

    let invalid = |reason: String| AnnotationError::InvalidBlocks {
        name: String::from_utf8_lossy(name).into_owned(),
        ordinal,
        reason,
    };

    if starts.len() != ends.len() {
        return Err(invalid(format!(
            "{} block start(s) but {} block end(s)",
            starts.len(),
            ends.len()
        )));
    }
    if starts.len() != count as usize {
        return Err(invalid(format!(
            "blockCount is {count} but {} block(s) were parsed",
            starts.len()
        )));
    }

    let mut previous_end = None;
    for (index, (start, end)) in starts.iter().zip(ends).enumerate() {
        if start >= end {
            return Err(invalid(format!(
                "block {index} is {start}..{end}, which is empty or reversed"
            )));
        }
        if let Some(previous) = previous_end
            && *start < previous
        {
            return Err(invalid(format!(
                "block {index} starts at {start}, before the previous block ended at {previous}"
            )));
        }
        previous_end = Some(*end);
    }
    Ok(true)
}

fn has_coding_bounds(record: &GenePred) -> bool {
    matches!((record.thick_start, record.thick_end), (Some(start), Some(end)) if start < end)
}

fn validate_coordinates(record: &GenePred, ordinal: u64) -> Result<(), AnnotationError> {
    if record.end <= record.start {
        return Err(AnnotationError::InvalidCoordinates {
            name: String::from_utf8_lossy(record.name.as_deref().unwrap_or(b"")).into_owned(),
            chrom: String::from_utf8_lossy(&record.chrom).into_owned(),
            ordinal,
            start: record.start,
            end: record.end,
            reason: "end must be greater than start",
        });
    }
    if record.chrom.is_empty() {
        return Err(AnnotationError::InvalidCoordinates {
            name: String::from_utf8_lossy(record.name.as_deref().unwrap_or(b"")).into_owned(),
            chrom: String::new(),
            ordinal,
            start: record.start,
            end: record.end,
            reason: "the reference name is empty",
        });
    }
    Ok(())
}

fn missing_error(
    missing: Missing,
    name: &[u8],
    ordinal: u64,
    options: &FeatureOptions,
) -> AnnotationError {
    let name = String::from_utf8_lossy(name).into_owned();
    let feature = options.feature.as_str();
    match missing {
        Missing::Blocks => AnnotationError::MissingBlocks {
            name,
            ordinal,
            feature,
            hint: if options.require_blocks {
                "; `--require-blocks` turned the span-as-exon fallback into this error"
            } else {
                ""
            },
        },
        Missing::Coding => AnnotationError::MissingCoding {
            name,
            ordinal,
            feature,
        },
        Missing::Strand => AnnotationError::MissingStrand {
            name,
            ordinal,
            feature,
        },
    }
}

/// The logical name of a record, and whether it had to be invented.
fn resolve_name(record: &GenePred, ordinal: u64, options: &FeatureOptions) -> (Vec<u8>, bool) {
    if let Some(field) = &options.name_field
        && let Some(value) = record
            .get_extra(field)
            .and_then(genepred::ExtraValue::first)
        && !value.is_empty()
    {
        return (value.to_vec(), false);
    }
    if let Some(name) = record.name.as_deref()
        && !name.is_empty()
        && name != b"."
    {
        return (name.to_vec(), false);
    }
    let generated = format!(
        "{}-{}-{}-{}-{}",
        options.unnamed_prefix,
        String::from_utf8_lossy(&record.chrom),
        record.start,
        record.end,
        ordinal
    );
    (generated.into_bytes(), true)
}

/// Makes a name unique: `name`, `name.2`, `name.3`, …
///
/// Two annotation records that genuinely share a name are common — a gene
/// present on a primary contig and an alt — and overwriting one with the other
/// would silently lose reads.
fn disambiguate(seen: &mut HashMap<Vec<u8>, u32>, name: &[u8]) -> (Vec<u8>, u32) {
    let count = seen.entry(name.to_vec()).or_insert(0);
    *count += 1;
    if *count == 1 {
        return (name.to_vec(), 0);
    }
    let mut unique = name.to_vec();
    unique.extend_from_slice(format!(".{count}").as_bytes());
    (unique, *count)
}

/// Parses a `--window-size` value.
///
/// Accepts a plain integer or one with a `k`, `M`, or `G` suffix, in **decimal**
/// SI units — `50M` is 50 000 000, not 52 428 800. Genomic coordinates are
/// counted in bases, not bytes, so powers of ten are what a user means.
///
/// # Errors
///
/// Returns [`AnnotationError::InvalidWindowSize`] for an unparseable value, an
/// unknown suffix, zero, or an overflowing product.
pub fn parse_window_size(text: &str) -> Result<u64, AnnotationError> {
    let invalid = |reason: &'static str| AnnotationError::InvalidWindowSize {
        value: text.to_string(),
        reason,
    };

    let trimmed = text.trim();
    if trimmed.is_empty() {
        return Err(invalid("the value is empty"));
    }
    let (digits, multiplier) = match trimmed.as_bytes()[trimmed.len() - 1] {
        b'k' | b'K' => (&trimmed[..trimmed.len() - 1], 1_000u64),
        b'm' | b'M' => (&trimmed[..trimmed.len() - 1], 1_000_000),
        b'g' | b'G' => (&trimmed[..trimmed.len() - 1], 1_000_000_000),
        byte if byte.is_ascii_digit() => (trimmed, 1),
        _ => return Err(invalid("the suffix must be one of k, M, or G")),
    };

    let value: u64 = digits
        .parse()
        .map_err(|_| invalid("expected a non-negative integer, optionally suffixed k, M, or G"))?;
    let size = value
        .checked_mul(multiplier)
        .ok_or_else(|| invalid("the value overflows a 64-bit integer"))?;
    if size == 0 {
        return Err(invalid("a window must be at least one base"));
    }
    Ok(size)
}

/// Generates fixed windows across a reference dictionary.
///
/// Windows are half-open and non-overlapping; the last window on each reference
/// is truncated to the reference length. Ordering is reference order, then
/// position, so ordinals — and therefore every tie-break — are deterministic.
#[must_use]
pub fn generate_windows(references: &[(Vec<u8>, u32)], window_size: u64) -> Vec<LogicalRegion> {
    let window_size = window_size.max(1);
    let mut regions = Vec::new();
    let mut ordinal = 0u64;

    for (chrom, length) in references {
        let length = u64::from(*length);
        let mut start = 0u64;
        while start < length {
            let end = start.saturating_add(window_size).min(length);
            // Named the way a user would type the region into `samtools view`:
            // 1-based, inclusive on both ends.
            let name =
                format!("{}:{}-{}", String::from_utf8_lossy(chrom), start + 1, end).into_bytes();
            let interval = Interval::new(start as i64, end as i64);
            regions.push(LogicalRegion {
                key: RegionKey::new(ordinal, name.clone()),
                chrom: chrom.clone(),
                segments: vec![interval],
                envelope: interval,
                source_ordinal: ordinal,
                metadata: RegionMetadata {
                    original_name: name,
                    generated_name: true,
                    duplicate_suffix: 0,
                    span_as_exon_fallback: false,
                    strand: None,
                    derived_segment_count: 1,
                },
            });
            ordinal += 1;
            start = end;
        }
    }
    regions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::detect::{AnnotationFormat, BedType, Compression};

    fn detection(format: AnnotationFormat) -> Detection {
        Detection {
            format,
            bed_type: (format == AnnotationFormat::Bed).then_some(BedType::Bed12),
            compression: Compression::None,
            trailing_columns: 0,
            reason: "test".to_string(),
        }
    }

    /// A three-exon transcript on `chr1`, coding from 1200 to 4800.
    ///
    /// ```text
    /// exons   1000..1200   2500..2800   4600..5000
    /// thick        1200..4800
    /// ```
    fn transcript(strand: Option<genepred::Strand>) -> GenePred {
        let mut record =
            GenePred::from_coords(b"chr1".to_vec(), 1000, 5000, genepred::Extras::new());
        record.set_name(Some(b"ENST01".to_vec()));
        record.set_strand(strand);
        record.set_thick_start(Some(1200));
        record.set_thick_end(Some(4800));
        record.set_block_count(Some(3));
        record.set_block_starts(Some(vec![1000, 2500, 4600]));
        record.set_block_ends(Some(vec![1200, 2800, 5000]));
        record
    }

    fn blockless(name: &[u8]) -> GenePred {
        let mut record = GenePred::from_coords(b"chr1".to_vec(), 100, 500, genepred::Extras::new());
        record.set_name(Some(name.to_vec()));
        record
    }

    fn extract(
        records: &[GenePred],
        format: AnnotationFormat,
        options: &FeatureOptions,
    ) -> Result<(Vec<LogicalRegion>, ExtractionStats), AnnotationError> {
        extract_regions(records, &detection(format), options)
    }

    fn segments(regions: &[LogicalRegion]) -> Vec<(i64, i64)> {
        regions
            .iter()
            .flat_map(|region| region.segments.iter().map(|s| (s.start, s.end)))
            .collect()
    }

    fn options(feature: FeatureType) -> FeatureOptions {
        FeatureOptions {
            feature,
            ..FeatureOptions::default()
        }
    }

    #[test]
    fn span_uses_the_whole_record() {
        let (regions, stats) = extract(
            &[transcript(Some(genepred::Strand::Forward))],
            AnnotationFormat::Bed,
            &options(FeatureType::Span),
        )
        .expect("extracted");
        assert_eq!(segments(&regions), [(1000, 5000)]);
        assert_eq!(regions[0].envelope, Interval::new(1000, 5000));
        assert_eq!(stats.annotation_records, 1);
        assert_eq!(stats.derived_segments, 1);
    }

    #[test]
    fn exon_uses_the_blocks() {
        let (regions, _) = extract(
            &[transcript(Some(genepred::Strand::Forward))],
            AnnotationFormat::Bed,
            &options(FeatureType::Exon),
        )
        .expect("extracted");
        assert_eq!(
            segments(&regions),
            [(1000, 1200), (2500, 2800), (4600, 5000)]
        );
        assert_eq!(regions[0].envelope, Interval::new(1000, 5000));
    }

    #[test]
    fn intron_uses_the_gaps() {
        let (regions, _) = extract(
            &[transcript(Some(genepred::Strand::Forward))],
            AnnotationFormat::Bed,
            &options(FeatureType::Intron),
        )
        .expect("extracted");
        assert_eq!(segments(&regions), [(1200, 2500), (2800, 4600)]);
    }

    #[test]
    fn cds_is_exonic_not_the_whole_thick_span() {
        let (regions, _) = extract(
            &[transcript(Some(genepred::Strand::Forward))],
            AnnotationFormat::Bed,
            &options(FeatureType::Cds),
        )
        .expect("extracted");
        // thick 1200..4800 intersected with the exons; the first exon ends
        // exactly at 1200 so it contributes nothing.
        assert_eq!(segments(&regions), [(2500, 2800), (4600, 4800)]);
        // The introns inside the thick span are deliberately absent.
        assert!(!segments(&regions).contains(&(1200, 4800)));
    }

    #[test]
    fn utr_is_exonic_and_excludes_introns() {
        let (regions, _) = extract(
            &[transcript(Some(genepred::Strand::Forward))],
            AnnotationFormat::Bed,
            &options(FeatureType::Utr),
        )
        .expect("extracted");
        assert_eq!(segments(&regions), [(1000, 1200), (4800, 5000)]);
    }

    #[test]
    fn five_and_three_prime_utrs_are_strand_aware() {
        let forward = extract(
            &[transcript(Some(genepred::Strand::Forward))],
            AnnotationFormat::Bed,
            &options(FeatureType::FiveUtr),
        )
        .expect("extracted")
        .0;
        assert_eq!(segments(&forward), [(1000, 1200)], "5' is upstream on +");

        let forward_three = extract(
            &[transcript(Some(genepred::Strand::Forward))],
            AnnotationFormat::Bed,
            &options(FeatureType::ThreeUtr),
        )
        .expect("extracted")
        .0;
        assert_eq!(
            segments(&forward_three),
            [(4800, 5000)],
            "3' is downstream on +"
        );

        let reverse = extract(
            &[transcript(Some(genepred::Strand::Reverse))],
            AnnotationFormat::Bed,
            &options(FeatureType::FiveUtr),
        )
        .expect("extracted")
        .0;
        assert_eq!(segments(&reverse), [(4800, 5000)], "5' is downstream on -");

        let reverse_three = extract(
            &[transcript(Some(genepred::Strand::Reverse))],
            AnnotationFormat::Bed,
            &options(FeatureType::ThreeUtr),
        )
        .expect("extracted")
        .0;
        assert_eq!(
            segments(&reverse_three),
            [(1000, 1200)],
            "3' is upstream on -"
        );
    }

    #[test]
    fn a_missing_strand_follows_the_policy() {
        for strand in [None, Some(genepred::Strand::Unknown)] {
            let record = transcript(strand);

            let (skipped, stats) = extract(
                std::slice::from_ref(&record),
                AnnotationFormat::Bed,
                &options(FeatureType::FiveUtr),
            )
            .expect("extracted");
            assert!(skipped.is_empty());
            assert_eq!(stats.missing_feature, 1);
            assert_eq!(stats.skipped_records, 1);

            let (kept, _) = extract(
                std::slice::from_ref(&record),
                AnnotationFormat::Bed,
                &FeatureOptions {
                    missing: MissingFeaturePolicy::Empty,
                    ..options(FeatureType::FiveUtr)
                },
            )
            .expect("extracted");
            assert_eq!(kept.len(), 1);
            assert!(kept[0].is_empty());

            let error = extract(
                std::slice::from_ref(&record),
                AnnotationFormat::Bed,
                &FeatureOptions {
                    missing: MissingFeaturePolicy::Error,
                    ..options(FeatureType::FiveUtr)
                },
            )
            .expect_err("must reject");
            assert!(
                matches!(error, AnnotationError::MissingStrand { .. }),
                "{error}"
            );
        }
    }

    #[test]
    fn missing_coding_bounds_follow_the_policy() {
        let mut record = transcript(Some(genepred::Strand::Forward));
        record.set_thick_start(None);
        record.set_thick_end(None);

        for feature in [
            FeatureType::Cds,
            FeatureType::Utr,
            FeatureType::FiveUtr,
            FeatureType::ThreeUtr,
        ] {
            let (regions, stats) = extract(
                std::slice::from_ref(&record),
                AnnotationFormat::Bed,
                &options(feature),
            )
            .expect("extracted");
            assert!(regions.is_empty(), "{feature}");
            assert_eq!(stats.missing_feature, 1, "{feature}");

            let error = extract(
                std::slice::from_ref(&record),
                AnnotationFormat::Bed,
                &FeatureOptions {
                    missing: MissingFeaturePolicy::Error,
                    ..options(feature)
                },
            )
            .expect_err("must reject");
            assert!(
                matches!(error, AnnotationError::MissingCoding { .. }),
                "{error}"
            );
        }
    }

    #[test]
    fn a_blockless_bed_record_falls_back_to_its_span() {
        let (regions, stats) = extract(
            &[blockless(b"g1")],
            AnnotationFormat::Bed,
            &options(FeatureType::Exon),
        )
        .expect("extracted");
        assert_eq!(segments(&regions), [(100, 500)]);
        assert!(regions[0].metadata.span_as_exon_fallback);
        assert_eq!(stats.span_as_exon_fallbacks, 1);
    }

    #[test]
    fn require_blocks_turns_the_fallback_into_an_error() {
        let error = extract(
            &[blockless(b"g1")],
            AnnotationFormat::Bed,
            &FeatureOptions {
                require_blocks: true,
                missing: MissingFeaturePolicy::Error,
                ..options(FeatureType::Exon)
            },
        )
        .expect_err("must reject");
        assert!(
            matches!(error, AnnotationError::MissingBlocks { .. }),
            "{error}"
        );
        assert!(error.to_string().contains("--require-blocks"), "{error}");
    }

    #[test]
    fn an_aggregated_gxf_record_with_no_exons_is_always_missing() {
        let (regions, stats) = extract(
            &[blockless(b"t1")],
            AnnotationFormat::Gtf,
            &options(FeatureType::Exon),
        )
        .expect("extracted");
        assert!(
            regions.is_empty(),
            "the GTF aggregation produced nothing usable"
        );
        assert_eq!(stats.missing_feature, 1);
    }

    #[test]
    fn a_single_exon_transcript_has_no_introns() {
        let mut record = transcript(Some(genepred::Strand::Forward));
        record.set_block_count(Some(1));
        record.set_block_starts(Some(vec![1000]));
        record.set_block_ends(Some(vec![5000]));

        let (regions, stats) = extract(
            &[record],
            AnnotationFormat::Bed,
            &options(FeatureType::Intron),
        )
        .expect("extracted");
        assert_eq!(regions.len(), 1);
        assert!(regions[0].is_empty());
        assert_eq!(stats.empty_regions, 1);
    }

    #[test]
    fn adjacent_segments_are_merged() {
        let mut record = transcript(Some(genepred::Strand::Forward));
        // Two exons that touch: 1000..1200 and 1200..1400.
        record.set_block_count(Some(2));
        record.set_block_starts(Some(vec![1000, 1200]));
        record.set_block_ends(Some(vec![1200, 1400]));
        record.set_end(1400);

        let (regions, _) = extract(
            &[record],
            AnnotationFormat::Bed,
            &options(FeatureType::Exon),
        )
        .expect("extracted");
        assert_eq!(segments(&regions), [(1000, 1400)], "abutting exons merge");
    }

    #[test]
    fn invalid_block_structure_is_rejected() {
        /// `(label, block starts, block ends, declared blockCount)`.
        type Case = (&'static str, Vec<u64>, Vec<u64>, Option<u32>);

        let cases: [Case; 3] = [
            ("overlapping", vec![1000, 1100], vec![1200, 1300], Some(2)),
            ("reversed", vec![1000, 900], vec![1200, 950], Some(2)),
            (
                "count mismatch",
                vec![1000, 2000],
                vec![1200, 2200],
                Some(5),
            ),
        ];
        for (label, starts, ends, count) in cases {
            let mut record = transcript(Some(genepred::Strand::Forward));
            record.set_block_count(count);
            record.set_block_starts(Some(starts));
            record.set_block_ends(Some(ends));
            let error = extract(
                &[record],
                AnnotationFormat::Bed,
                &options(FeatureType::Exon),
            )
            .expect_err("must reject");
            assert!(
                matches!(error, AnnotationError::InvalidBlocks { .. }),
                "{label}: {error}"
            );
        }
    }

    #[test]
    fn a_zero_length_block_is_rejected() {
        let mut record = transcript(Some(genepred::Strand::Forward));
        record.set_block_count(Some(2));
        record.set_block_starts(Some(vec![1000, 2000]));
        record.set_block_ends(Some(vec![1000, 2200]));
        let error = extract(
            &[record],
            AnnotationFormat::Bed,
            &options(FeatureType::Exon),
        )
        .expect_err("must reject");
        assert!(
            matches!(error, AnnotationError::InvalidBlocks { .. }),
            "{error}"
        );
    }

    #[test]
    fn reversed_record_coordinates_are_rejected() {
        let mut record = blockless(b"g1");
        record.set_start(500);
        record.set_end(100);
        let error = extract(
            &[record],
            AnnotationFormat::Bed,
            &options(FeatureType::Span),
        )
        .expect_err("must reject");
        assert!(
            matches!(error, AnnotationError::InvalidCoordinates { .. }),
            "{error}"
        );
    }

    #[test]
    fn duplicate_names_are_disambiguated_deterministically() {
        let records = vec![blockless(b"gene"), blockless(b"gene"), blockless(b"gene")];
        let (regions, stats) =
            extract(&records, AnnotationFormat::Bed, &options(FeatureType::Span))
                .expect("extracted");
        let names: Vec<String> = regions
            .iter()
            .map(|region| String::from_utf8_lossy(region.key.name()).into_owned())
            .collect();
        assert_eq!(names, ["gene", "gene.2", "gene.3"]);
        assert_eq!(stats.duplicate_names, 2);
        // The original is preserved for the manifest.
        assert_eq!(regions[1].metadata.original_name, b"gene");
        assert_eq!(regions[1].metadata.duplicate_suffix, 2);
    }

    #[test]
    fn unnamed_records_get_a_deterministic_generated_name() {
        let mut record = blockless(b"");
        record.set_name(None);
        let (regions, stats) = extract(
            std::slice::from_ref(&record),
            AnnotationFormat::Bed,
            &options(FeatureType::Span),
        )
        .expect("extracted");
        assert_eq!(
            String::from_utf8_lossy(regions[0].key.name()),
            "region-chr1-100-500-0"
        );
        assert!(regions[0].metadata.generated_name);
        assert_eq!(stats.generated_names, 1);

        let (custom, _) = extract(
            std::slice::from_ref(&record),
            AnnotationFormat::Bed,
            &FeatureOptions {
                unnamed_prefix: "win".to_string(),
                ..options(FeatureType::Span)
            },
        )
        .expect("extracted");
        assert!(String::from_utf8_lossy(custom[0].key.name()).starts_with("win-"));
    }

    #[test]
    fn a_dot_name_counts_as_unnamed() {
        let record = blockless(b".");
        let (regions, _) = extract(
            &[record],
            AnnotationFormat::Bed,
            &options(FeatureType::Span),
        )
        .expect("extracted");
        assert!(regions[0].metadata.generated_name);
    }

    #[test]
    fn a_name_field_takes_precedence() {
        let mut record = blockless(b"fallback");
        record.add_extra(b"gene_name".to_vec(), b"BRCA1".to_vec());
        let (regions, _) = extract(
            &[record],
            AnnotationFormat::Bed,
            &FeatureOptions {
                name_field: Some(b"gene_name".to_vec()),
                ..options(FeatureType::Span)
            },
        )
        .expect("extracted");
        assert_eq!(regions[0].key.name(), b"BRCA1");
        assert!(!regions[0].metadata.generated_name);
    }

    #[test]
    fn an_absent_name_field_falls_back_to_the_record_name() {
        let record = blockless(b"fallback");
        let (regions, _) = extract(
            &[record],
            AnnotationFormat::Bed,
            &FeatureOptions {
                name_field: Some(b"absent".to_vec()),
                ..options(FeatureType::Span)
            },
        )
        .expect("extracted");
        assert_eq!(regions[0].key.name(), b"fallback");
    }

    #[test]
    fn region_keys_order_by_annotation_input_order() {
        let mut keys = [
            RegionKey::new(2, b"zzz".to_vec()),
            RegionKey::new(0, b"mmm".to_vec()),
            RegionKey::new(1, b"aaa".to_vec()),
        ];
        keys.sort();
        assert_eq!(
            keys.iter().map(RegionKey::ordinal).collect::<Vec<_>>(),
            [0, 1, 2]
        );
    }

    #[test]
    fn window_sizes_parse_in_decimal_si_units() {
        assert_eq!(parse_window_size("1000").expect("valid"), 1_000);
        assert_eq!(parse_window_size("50M").expect("valid"), 50_000_000);
        assert_eq!(parse_window_size("2k").expect("valid"), 2_000);
        assert_eq!(parse_window_size("1G").expect("valid"), 1_000_000_000);
        assert_eq!(parse_window_size(" 100 ").expect("valid"), 100);

        for bad in ["", "0", "0M", "abc", "10X", "-5", "1.5M"] {
            assert!(parse_window_size(bad).is_err(), "{bad:?} must be rejected");
        }
        assert!(parse_window_size("99999999999999999999G").is_err());
    }

    #[test]
    fn windows_tile_each_reference_and_truncate_the_last() {
        let references = vec![(b"chr1".to_vec(), 250u32), (b"chr2".to_vec(), 100)];
        let windows = generate_windows(&references, 100);

        let described: Vec<(String, i64, i64)> = windows
            .iter()
            .map(|region| {
                (
                    String::from_utf8_lossy(region.key.name()).into_owned(),
                    region.envelope.start,
                    region.envelope.end,
                )
            })
            .collect();
        assert_eq!(
            described,
            [
                ("chr1:1-100".to_string(), 0, 100),
                ("chr1:101-200".to_string(), 100, 200),
                ("chr1:201-250".to_string(), 200, 250),
                ("chr2:1-100".to_string(), 0, 100),
            ]
        );
        // Ordinals are dense and in reference order.
        assert_eq!(
            windows.iter().map(|r| r.source_ordinal).collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
    }

    #[test]
    fn a_window_larger_than_the_reference_yields_one_window() {
        let windows = generate_windows(&[(b"chrM".to_vec(), 16_569)], 50_000_000);
        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].envelope, Interval::new(0, 16_569));
    }

    #[test]
    fn a_zero_length_reference_yields_no_windows() {
        assert!(generate_windows(&[(b"empty".to_vec(), 0)], 100).is_empty());
    }

    #[test]
    fn feature_types_round_trip_their_option_values() {
        for feature in [
            FeatureType::Span,
            FeatureType::Exon,
            FeatureType::Intron,
            FeatureType::Cds,
            FeatureType::Utr,
            FeatureType::FiveUtr,
            FeatureType::ThreeUtr,
        ] {
            assert_eq!(FeatureType::parse(feature.as_str()), Some(feature));
        }
        assert_eq!(FeatureType::parse("nonesuch"), None);
        assert!(FeatureType::Cds.needs_coding_bounds());
        assert!(!FeatureType::Exon.needs_coding_bounds());
        assert!(FeatureType::FiveUtr.needs_strand());
        assert!(!FeatureType::Utr.needs_strand());
        assert!(!FeatureType::Span.needs_blocks());
    }
}
