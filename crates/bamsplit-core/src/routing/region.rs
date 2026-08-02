// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Routing by annotated region.
//!
//! # Feature and assignment are independent
//!
//! `--feature` decides *which segments of an annotation record* participate.
//! `--assignment` decides *how an alignment is matched against them*. The two
//! compose freely, so `--feature exon --assignment overlap` compares a read's
//! aligned blocks against the union of a transcript's exons.
//!
//! ```text
//! transcript   ▓▓▓▓▓        ▓▓▓▓▓▓▓          ▓▓▓▓        (--feature exon)
//! read A        ══                                        overlap ✓  contained ✓
//! read B            ══════════                            overlap ✓  contained ✗
//! read C     ══                                           overlap ✗
//! read D        ══───────────══                           spliced: both blocks inside
//! ```
//!
//! Read D is the case that matters for RNA-seq: its `N` gap is not required to
//! be exonic, so a correctly spliced read *is* contained even though the
//! intron it skips is not part of the feature.
//!
//! # The only duplicating mode
//!
//! `overlap` emits a record to every region it touches, and is the sole reason
//! [`Route::Many`] exists. Every other mode emits at most once, and the
//! manifest's conservation equation switches accordingly.
//!
//! # Determinism
//!
//! Candidate regions arrive from the index in annotation input order, and every
//! tie-break ends in that order, so the same input always produces the same
//! assignment — including which region wins a tie.

use std::cell::RefCell;
use std::sync::atomic::{AtomicU64, Ordering};

use smallvec::SmallVec;

use crate::annotation::features::{LogicalRegion, RegionKey};
use crate::annotation::interval_index::IntervalIndex;
use crate::bam::cigar::{AlignmentGeometry, DeletionPolicy};
use crate::bam::header::BamHeader;
use crate::bam::raw_record::RawRecord;
use crate::error::RoutingError;
use crate::interval::{Interval, contains_all_sorted, overlap_len_sorted, total_len};
use crate::routing::{DropReason, Route, Router};

/// How an alignment is matched against a region's segments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AssignmentMode {
    /// By the alignment's leftmost reference position. At most one output.
    #[default]
    Start,
    /// By the midpoint of the alignment's reference span. At most one output.
    Midpoint,
    /// Only when the whole alignment sits inside the feature. At most one.
    Contained,
    /// To the region sharing the most reference bases. Exactly one, when any
    /// region qualifies.
    BestOverlap,
    /// To **every** region with a qualifying overlap. May duplicate.
    Overlap,
}

impl AssignmentMode {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Start => "start",
            Self::Midpoint => "midpoint",
            Self::Contained => "contained",
            Self::BestOverlap => "best-overlap",
            Self::Overlap => "overlap",
        }
    }

    /// Parses an option value.
    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "start" => Self::Start,
            "midpoint" => Self::Midpoint,
            "contained" => Self::Contained,
            "best-overlap" => Self::BestOverlap,
            "overlap" => Self::Overlap,
            _ => return None,
        })
    }

    /// Whether this mode can send one record to several outputs.
    #[must_use]
    pub const fn may_duplicate(self) -> bool {
        matches!(self, Self::Overlap)
    }
}

impl std::fmt::Display for AssignmentMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// How region routing behaves.
#[derive(Debug, Clone, Copy)]
pub struct RegionRouterOptions {
    /// How an alignment is matched.
    pub assignment: AssignmentMode,
    /// Whether to compare against aligned blocks or the outer span.
    pub geometry: AlignmentGeometry,
    /// Whether a deletion counts as aligned overlap.
    pub deletions: DeletionPolicy,
}

impl Default for RegionRouterOptions {
    fn default() -> Self {
        Self {
            assignment: AssignmentMode::Start,
            geometry: AlignmentGeometry::Blocks,
            deletions: DeletionPolicy::ExcludeFromBlocks,
        }
    }
}

/// What routing observed, for the manifest.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RegionRoutingStats {
    /// Records where more than one region qualified and a tie-break decided.
    pub ambiguous_assignments: u64,
    /// Records that matched no region.
    pub unmatched_records: u64,
    /// Records that reached at least one region.
    pub matched_records: u64,
    /// Emissions beyond the first, only possible in `overlap` mode.
    pub duplicate_emissions: u64,
    /// The largest number of regions one record reached.
    pub max_matches_for_one_record: u64,
}

thread_local! {
    /// Per-thread scratch, so routing a record allocates nothing.
    ///
    /// The indexed engine routes from several threads; each gets its own.
    static SCRATCH: RefCell<Scratch> = const { RefCell::new(Scratch::new()) };
}

struct Scratch {
    geometry: Vec<Interval>,
    candidates: Vec<usize>,
    qualified: Vec<Candidate>,
}

impl Scratch {
    const fn new() -> Self {
        Self {
            geometry: Vec::new(),
            candidates: Vec::new(),
            qualified: Vec::new(),
        }
    }
}

/// A region that passed the mode's test, with the numbers a tie-break needs.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    region: usize,
    ordinal: u64,
    /// Reference bases shared with the alignment. Zero for the point modes,
    /// which do not measure overlap.
    overlap: i64,
    /// Total bases the region's segments cover; the "smaller feature" rule.
    span: i64,
}

/// Routes records to annotated regions.
pub struct RegionRouter {
    index: IntervalIndex,
    options: RegionRouterOptions,
    ambiguous: AtomicU64,
    unmatched: AtomicU64,
    matched: AtomicU64,
    duplicates: AtomicU64,
    max_matches: AtomicU64,
}

impl std::fmt::Debug for RegionRouter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RegionRouter")
            .field("regions", &self.index.region_count())
            .field("assignment", &self.options.assignment)
            .field("geometry", &self.options.geometry)
            .finish_non_exhaustive()
    }
}

impl RegionRouter {
    /// Creates a router over an interval index.
    #[must_use]
    pub fn new(index: IntervalIndex, options: RegionRouterOptions) -> Self {
        Self {
            index,
            options,
            ambiguous: AtomicU64::new(0),
            unmatched: AtomicU64::new(0),
            matched: AtomicU64::new(0),
            duplicates: AtomicU64::new(0),
            max_matches: AtomicU64::new(0),
        }
    }

    /// The index being routed against.
    #[must_use]
    pub const fn index(&self) -> &IntervalIndex {
        &self.index
    }

    /// The counters accumulated so far.
    #[must_use]
    pub fn stats(&self) -> RegionRoutingStats {
        RegionRoutingStats {
            ambiguous_assignments: self.ambiguous.load(Ordering::Relaxed),
            unmatched_records: self.unmatched.load(Ordering::Relaxed),
            matched_records: self.matched.load(Ordering::Relaxed),
            duplicate_emissions: self.duplicates.load(Ordering::Relaxed),
            max_matches_for_one_record: self.max_matches.load(Ordering::Relaxed),
        }
    }

    /// Fills `out` with the alignment geometry a record contributes.
    ///
    /// Returns the envelope covering it, or [`None`] when the record has no
    /// position at all.
    fn geometry_of(
        &self,
        record: &RawRecord<'_>,
        out: &mut Vec<Interval>,
    ) -> Result<Option<Interval>, RoutingError> {
        out.clear();
        let Some(start) = record.alignment_start().map_err(Box::new)? else {
            return Ok(None);
        };
        let start = i64::from(start);
        let ops = record.cigar_ops().map_err(Box::new)?;

        match self.options.geometry {
            AlignmentGeometry::Span => {
                let span =
                    ops.reference_interval(start)
                        .map_err(|source| RoutingError::Geometry {
                            source: Box::new(source),
                            location: crate::bam::raw_record::RecordLocation::UNKNOWN,
                        })?;
                out.push(span);
            }
            AlignmentGeometry::Blocks => {
                ops.for_each_aligned_block(start, self.options.deletions, |block| {
                    out.push(block);
                })
                .map_err(|source| RoutingError::Geometry {
                    source: Box::new(source),
                    location: crate::bam::raw_record::RecordLocation::UNKNOWN,
                })?;
                if out.is_empty() {
                    // A CIGAR with no aligned bases — all insertion, or absent
                    // entirely — still sits somewhere. Falling back to the
                    // one-base reference interval keeps such a record routable
                    // instead of silently unmatched.
                    let span =
                        ops.reference_interval(start)
                            .map_err(|source| RoutingError::Geometry {
                                source: Box::new(source),
                                location: crate::bam::raw_record::RecordLocation::UNKNOWN,
                            })?;
                    out.push(span);
                }
            }
        }

        Ok(match (out.first(), out.last()) {
            (Some(first), Some(last)) => Some(Interval::new(first.start, last.end)),
            _ => None,
        })
    }

    /// Records the outcome of one routing decision.
    fn observe(&self, emissions: usize, tied: bool) {
        if emissions == 0 {
            self.unmatched.fetch_add(1, Ordering::Relaxed);
            return;
        }
        self.matched.fetch_add(1, Ordering::Relaxed);
        if emissions > 1 {
            self.duplicates
                .fetch_add(emissions as u64 - 1, Ordering::Relaxed);
        }
        self.max_matches
            .fetch_max(emissions as u64, Ordering::Relaxed);
        if tied {
            self.ambiguous.fetch_add(1, Ordering::Relaxed);
        }
    }
}

/// Orders two candidates for a unique-assignment mode.
///
/// Smaller feature first, then earlier annotation order — the documented rule,
/// and the reason "which region does an ambiguous read belong to" has one
/// answer rather than whichever the hash map happened to yield.
fn by_smallest_then_earliest(left: &Candidate, right: &Candidate) -> std::cmp::Ordering {
    left.span
        .cmp(&right.span)
        .then(left.ordinal.cmp(&right.ordinal))
}

/// Orders two candidates for `best-overlap`.
///
/// Greatest overlap wins. The specification's next criterion — greatest overlap
/// *fraction* relative to the alignment's reference-consuming length — cannot
/// break a tie here, because that denominator is the same record for every
/// candidate; it is therefore folded into the overlap comparison rather than
/// pretended to be a separate step. What actually decides a tie is the smaller
/// feature, then annotation order.
fn by_best_overlap(left: &Candidate, right: &Candidate) -> std::cmp::Ordering {
    right
        .overlap
        .cmp(&left.overlap)
        .then_with(|| by_smallest_then_earliest(left, right))
}

impl Router for RegionRouter {
    type Key = RegionKey;

    fn route(
        &self,
        header: &BamHeader,
        record: &RawRecord<'_>,
    ) -> Result<Route<Self::Key>, RoutingError> {
        let Some(reference_id) = record.reference_sequence_id().map_err(Box::new)? else {
            // An unplaced record has no coordinate, so no region can contain it.
            self.observe(0, false);
            return Ok(Route::Drop(DropReason::Unmatched));
        };
        let Some(chrom) = usize::try_from(reference_id)
            .ok()
            .and_then(|id| header.reference_name(id))
        else {
            return Err(RoutingError::InvalidReferenceId {
                id: reference_id,
                reference_count: header.reference_count(),
                location: crate::bam::raw_record::RecordLocation::UNKNOWN,
            });
        };

        SCRATCH.with(|scratch| {
            let scratch = &mut *scratch.borrow_mut();
            let Scratch {
                geometry,
                candidates,
                qualified,
            } = scratch;

            let Some(envelope) = self.geometry_of(record, geometry)? else {
                self.observe(0, false);
                return Ok(Route::Drop(DropReason::Unmatched));
            };

            self.index.query(chrom, envelope, candidates);
            qualified.clear();

            for &index in candidates.iter() {
                let Some(region) = self.index.region(index) else {
                    continue;
                };
                if region.segments.is_empty() {
                    continue;
                }
                if let Some(candidate) = qualify(region, index, geometry, envelope, self.options) {
                    qualified.push(candidate);
                }
            }

            if qualified.is_empty() {
                self.observe(0, false);
                return Ok(Route::Drop(DropReason::Unmatched));
            }

            let tied = qualified.len() > 1;
            match self.options.assignment {
                AssignmentMode::Overlap => {
                    // Already in annotation input order, which fixes the order
                    // outputs are created in.
                    let keys: SmallVec<[RegionKey; 2]> = qualified
                        .iter()
                        .filter_map(|candidate| self.index.region(candidate.region))
                        .map(|region| region.key.clone())
                        .collect();
                    self.observe(keys.len(), false);
                    Ok(if keys.len() == 1 {
                        // A single match is still a single match; `One` keeps
                        // the manifest's duplicate accounting honest.
                        Route::One(
                            keys.into_iter().next().unwrap_or_else(|| {
                                unreachable!("length was just checked to be one")
                            }),
                        )
                    } else {
                        Route::Many(keys)
                    })
                }
                AssignmentMode::BestOverlap => {
                    let best = qualified
                        .iter()
                        .min_by(|left, right| by_best_overlap(left, right))
                        .copied();
                    Ok(self.finish_unique(best, tied))
                }
                AssignmentMode::Start | AssignmentMode::Midpoint | AssignmentMode::Contained => {
                    let best = qualified
                        .iter()
                        .min_by(|left, right| by_smallest_then_earliest(left, right))
                        .copied();
                    Ok(self.finish_unique(best, tied))
                }
            }
        })
    }

    fn mode(&self) -> &'static str {
        "region"
    }

    fn may_duplicate(&self) -> bool {
        self.options.assignment.may_duplicate()
    }

    fn declared_keys(&self, _header: &BamHeader) -> Vec<Self::Key> {
        self.index
            .regions()
            .iter()
            .map(|region| region.key.clone())
            .collect()
    }

    fn is_grouped_by_coordinate(&self) -> bool {
        // Regions can overlap and nest, so even a coordinate-sorted BAM visits
        // their keys in an interleaved order. Promising otherwise would make the
        // streaming engine finalize an output it later has to reopen.
        false
    }
}

impl RegionRouter {
    fn finish_unique(&self, best: Option<Candidate>, tied: bool) -> Route<RegionKey> {
        let Some(region) = best.and_then(|candidate| self.index.region(candidate.region)) else {
            self.observe(0, false);
            return Route::Drop(DropReason::Unmatched);
        };
        self.observe(1, tied);
        Route::One(region.key.clone())
    }
}

/// Whether a region qualifies for this alignment under the chosen mode.
fn qualify(
    region: &LogicalRegion,
    index: usize,
    geometry: &[Interval],
    envelope: Interval,
    options: RegionRouterOptions,
) -> Option<Candidate> {
    let span = total_len(&region.segments);
    let base = Candidate {
        region: index,
        ordinal: region.source_ordinal,
        overlap: 0,
        span,
    };

    match options.assignment {
        AssignmentMode::Start => {
            let point = envelope.start;
            region
                .segments
                .iter()
                .any(|segment| segment.contains(point))
                .then_some(base)
        }
        AssignmentMode::Midpoint => {
            // The midpoint of the alignment's outer span, not of its blocks: a
            // spliced read's midpoint can legitimately fall inside its intron,
            // and that is the position the mode names.
            let point = envelope.midpoint();
            region
                .segments
                .iter()
                .any(|segment| segment.contains(point))
                .then_some(base)
        }
        AssignmentMode::Contained => {
            // Every reference-consuming aligned block must be inside the
            // feature. The gaps between blocks — skipped `N` regions — are not
            // required to be, which is what makes a correctly spliced read
            // contained in the transcript it came from.
            contains_all_sorted(&region.segments, geometry).then_some(base)
        }
        AssignmentMode::BestOverlap | AssignmentMode::Overlap => {
            let overlap = overlap_len_sorted(geometry, &region.segments);
            (overlap > 0).then_some(Candidate { overlap, ..base })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::features::{RegionKey as Key, RegionMetadata};
    use crate::bam::header::ReferenceSequence;
    use crate::routing::RoutingKey as _;

    fn header() -> BamHeader {
        BamHeader::from_parts(
            b"@HD\tVN:1.6\tSO:coordinate\n".to_vec(),
            vec![
                ReferenceSequence {
                    name: b"chr1".to_vec(),
                    length: 1_000_000,
                },
                ReferenceSequence {
                    name: b"chr2".to_vec(),
                    length: 1_000_000,
                },
            ],
        )
        .expect("valid header")
    }

    fn region(ordinal: u64, name: &str, chrom: &[u8], segments: &[(i64, i64)]) -> LogicalRegion {
        let segments: Vec<Interval> = segments
            .iter()
            .map(|(start, end)| Interval::new(*start, *end))
            .collect();
        let envelope = Interval::new(
            segments.first().map_or(0, |first| first.start),
            segments.last().map_or(0, |last| last.end),
        );
        LogicalRegion {
            key: Key::new(ordinal, name.as_bytes().to_vec()),
            chrom: chrom.to_vec(),
            segments,
            envelope,
            source_ordinal: ordinal,
            metadata: RegionMetadata {
                original_name: name.as_bytes().to_vec(),
                generated_name: false,
                duplicate_suffix: 0,
                span_as_exon_fallback: false,
                strand: None,
                derived_segment_count: 0,
            },
        }
    }

    /// A record on `chr1` with the given CIGAR, as `(length, op_code)`.
    fn record(position: i32, cigar: &[(u32, u32)]) -> Vec<u8> {
        let read_length: usize = cigar
            .iter()
            .filter(|(_, kind)| matches!(kind, 0 | 1 | 4 | 7 | 8))
            .map(|(length, _)| *length as usize)
            .sum();
        let mut body = Vec::new();
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&position.to_le_bytes());
        body.extend_from_slice(&[2, 60]);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&u16::try_from(cigar.len()).expect("short").to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&i32::try_from(read_length).expect("short").to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(b"r\0");
        for (length, kind) in cigar {
            body.extend_from_slice(&((length << 4) | kind).to_le_bytes());
        }
        body.extend(std::iter::repeat_n(0u8, read_length.div_ceil(2)));
        body.extend(std::iter::repeat_n(0xffu8, read_length));
        body
    }

    fn unplaced() -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&[2, 0]);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&4u16.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(b"r\0");
        body
    }

    fn router(regions: Vec<LogicalRegion>, assignment: AssignmentMode) -> RegionRouter {
        RegionRouter::new(
            IntervalIndex::build(regions).expect("built"),
            RegionRouterOptions {
                assignment,
                ..RegionRouterOptions::default()
            },
        )
    }

    fn names(route: &Route<RegionKey>) -> Vec<String> {
        route
            .keys()
            .iter()
            .map(|key| String::from_utf8_lossy(key.logical()).into_owned())
            .collect()
    }

    fn route_of(router: &RegionRouter, body: &[u8]) -> Route<RegionKey> {
        router
            .route(&header(), &RawRecord::new(body).expect("valid"))
            .expect("routable")
    }

    #[test]
    fn start_assigns_by_the_leftmost_position() {
        let router = router(
            vec![
                region(0, "a", b"chr1", &[(100, 200)]),
                region(1, "b", b"chr1", &[(300, 400)]),
            ],
            AssignmentMode::Start,
        );
        assert_eq!(names(&route_of(&router, &record(150, &[(10, 0)]))), ["a"]);
        assert_eq!(names(&route_of(&router, &record(350, &[(10, 0)]))), ["b"]);
        // A read starting before the region does not qualify even if it reaches
        // into it.
        assert_eq!(
            route_of(&router, &record(90, &[(100, 0)])),
            Route::Drop(DropReason::Unmatched)
        );
    }

    #[test]
    fn start_prefers_the_smallest_containing_feature_then_input_order() {
        let router = router(
            vec![
                region(0, "wide", b"chr1", &[(0, 1000)]),
                region(1, "narrow", b"chr1", &[(100, 200)]),
                region(2, "same-size", b"chr1", &[(100, 200)]),
            ],
            AssignmentMode::Start,
        );
        // The narrow feature wins over the wide one.
        assert_eq!(
            names(&route_of(&router, &record(150, &[(10, 0)]))),
            ["narrow"]
        );
        // Between two equally small features, the earlier one wins.
        assert_eq!(router.stats().ambiguous_assignments, 1);
    }

    #[test]
    fn midpoint_uses_the_span_not_the_blocks() {
        let router = router(
            vec![region(0, "intronic", b"chr1", &[(1000, 2000)])],
            AssignmentMode::Midpoint,
        );
        // 100M 1000N 100M starting at 500: span 500..1700, midpoint 1100, which
        // falls inside the intron and inside the region.
        assert_eq!(
            names(&route_of(
                &router,
                &record(500, &[(100, 0), (1000, 3), (100, 0)])
            )),
            ["intronic"]
        );
    }

    #[test]
    fn contained_requires_every_aligned_block_inside() {
        let exons = vec![region(
            0,
            "transcript",
            b"chr1",
            &[(1000, 1100), (2000, 2100)],
        )];
        let router = router(exons.clone(), AssignmentMode::Contained);

        // A spliced read whose two blocks land in the two exons is contained,
        // even though the intron it skips is not part of the feature.
        let spliced = record(1050, &[(50, 0), (900, 3), (50, 0)]);
        assert_eq!(names(&route_of(&router, &spliced)), ["transcript"]);

        // A read that runs off the end of an exon is not contained.
        let overhang = record(1050, &[(100, 0)]);
        assert_eq!(
            route_of(&router, &overhang),
            Route::Drop(DropReason::Unmatched)
        );

        // An intron-retention read spans the gap contiguously and is not
        // contained either.
        let retained = record(1050, &[(1000, 0)]);
        assert_eq!(
            route_of(&router, &retained),
            Route::Drop(DropReason::Unmatched)
        );
    }

    #[test]
    fn best_overlap_picks_the_largest_shared_span() {
        let router = router(
            vec![
                region(0, "small", b"chr1", &[(1000, 1020)]),
                region(1, "large", b"chr1", &[(1000, 1200)]),
            ],
            AssignmentMode::BestOverlap,
        );
        // 100M from 1000: 20 bases with `small`, 100 with `large`.
        assert_eq!(
            names(&route_of(&router, &record(1000, &[(100, 0)]))),
            ["large"]
        );
    }

    #[test]
    fn best_overlap_breaks_ties_by_smaller_feature_then_input_order() {
        let router = router(
            vec![
                region(0, "wide", b"chr1", &[(0, 10_000)]),
                region(1, "narrow", b"chr1", &[(1000, 1100)]),
            ],
            AssignmentMode::BestOverlap,
        );
        // A 100-base read inside both overlaps each by exactly 100.
        assert_eq!(
            names(&route_of(&router, &record(1000, &[(100, 0)]))),
            ["narrow"]
        );
        assert_eq!(router.stats().ambiguous_assignments, 1);
    }

    #[test]
    fn overlap_emits_to_every_matching_region() {
        let router = router(
            vec![
                region(0, "a", b"chr1", &[(1000, 1100)]),
                region(1, "b", b"chr1", &[(1050, 1200)]),
                region(2, "c", b"chr1", &[(5000, 5100)]),
            ],
            AssignmentMode::Overlap,
        );
        let route = route_of(&router, &record(1040, &[(100, 0)]));
        assert_eq!(names(&route), ["a", "b"]);
        assert_eq!(route.emission_count(), 2);
        assert!(router.may_duplicate());

        let stats = router.stats();
        assert_eq!(stats.duplicate_emissions, 1);
        assert_eq!(stats.max_matches_for_one_record, 2);
    }

    #[test]
    fn overlap_with_one_match_is_reported_as_a_single_emission() {
        let router = router(
            vec![region(0, "a", b"chr1", &[(1000, 1100)])],
            AssignmentMode::Overlap,
        );
        let route = route_of(&router, &record(1000, &[(50, 0)]));
        assert!(matches!(route, Route::One(_)));
        assert_eq!(router.stats().duplicate_emissions, 0);
    }

    #[test]
    fn a_spliced_read_does_not_match_the_intron_it_skips() {
        let router = router(
            vec![region(0, "intron", b"chr1", &[(1100, 1900)])],
            AssignmentMode::Overlap,
        );
        // 100M 800N 100M from 1000: the blocks are 1000..1100 and 1900..2000,
        // neither of which touches the intron feature.
        let spliced = record(1000, &[(100, 0), (800, 3), (100, 0)]);
        assert_eq!(
            route_of(&router, &spliced),
            Route::Drop(DropReason::Unmatched)
        );

        // With `--alignment-geometry span` the same read *does* match, because
        // the outer span covers the intron.
        let span_router = RegionRouter::new(
            IntervalIndex::build(vec![region(0, "intron", b"chr1", &[(1100, 1900)])])
                .expect("built"),
            RegionRouterOptions {
                assignment: AssignmentMode::Overlap,
                geometry: AlignmentGeometry::Span,
                ..RegionRouterOptions::default()
            },
        );
        assert_eq!(names(&route_of(&span_router, &spliced)), ["intron"]);
    }

    #[test]
    fn an_intron_retention_read_matches_the_intron() {
        let router = router(
            vec![region(0, "intron", b"chr1", &[(1100, 1900)])],
            AssignmentMode::Overlap,
        );
        // A contiguous 1000M read across the same span aligns through it.
        assert_eq!(
            names(&route_of(&router, &record(1000, &[(1000, 0)]))),
            ["intron"]
        );
    }

    #[test]
    fn a_deletion_does_not_contribute_overlap() {
        let router = router(
            vec![region(0, "gap", b"chr1", &[(1100, 1120)])],
            AssignmentMode::Overlap,
        );
        // 100M 50D 100M from 1000: the deletion covers 1100..1150, which is
        // exactly the feature, but deletions are not aligned bases.
        let deleted = record(1000, &[(100, 0), (50, 2), (100, 0)]);
        assert_eq!(
            route_of(&router, &deleted),
            Route::Drop(DropReason::Unmatched)
        );

        // Including deletions in blocks flips the answer, which is why the
        // policy is explicit.
        let inclusive = RegionRouter::new(
            IntervalIndex::build(vec![region(0, "gap", b"chr1", &[(1100, 1120)])]).expect("built"),
            RegionRouterOptions {
                assignment: AssignmentMode::Overlap,
                deletions: DeletionPolicy::IncludeInBlocks,
                ..RegionRouterOptions::default()
            },
        );
        assert_eq!(names(&route_of(&inclusive, &deleted)), ["gap"]);
    }

    #[test]
    fn regions_on_another_reference_never_match() {
        let router = router(
            vec![region(0, "elsewhere", b"chr2", &[(1000, 2000)])],
            AssignmentMode::Overlap,
        );
        assert_eq!(
            route_of(&router, &record(1500, &[(10, 0)])),
            Route::Drop(DropReason::Unmatched)
        );
    }

    #[test]
    fn an_unplaced_record_is_unmatched_rather_than_an_error() {
        let router = router(
            vec![region(0, "a", b"chr1", &[(0, 1_000_000)])],
            AssignmentMode::Overlap,
        );
        assert_eq!(
            route_of(&router, &unplaced()),
            Route::Drop(DropReason::Unmatched)
        );
        assert_eq!(router.stats().unmatched_records, 1);
    }

    #[test]
    fn an_empty_region_never_matches_but_is_still_declared() {
        let mut empty = region(1, "empty", b"chr1", &[]);
        empty.envelope = Interval::empty_at(1500);
        let router = router(
            vec![region(0, "real", b"chr1", &[(1000, 2000)]), empty],
            AssignmentMode::Overlap,
        );
        assert_eq!(
            names(&route_of(&router, &record(1500, &[(10, 0)]))),
            ["real"]
        );

        let declared: Vec<String> = router
            .declared_keys(&header())
            .iter()
            .map(|key| String::from_utf8_lossy(key.logical()).into_owned())
            .collect();
        assert_eq!(declared, ["real", "empty"], "--emit-empty needs both");
    }

    #[test]
    fn a_record_with_no_aligned_bases_still_routes_by_position() {
        let router = router(
            vec![region(0, "a", b"chr1", &[(1000, 2000)])],
            AssignmentMode::Start,
        );
        // 10I consumes no reference; the record still sits at 1500.
        assert_eq!(names(&route_of(&router, &record(1500, &[(10, 1)]))), ["a"]);
    }

    #[test]
    fn unique_modes_never_duplicate() {
        for assignment in [
            AssignmentMode::Start,
            AssignmentMode::Midpoint,
            AssignmentMode::Contained,
            AssignmentMode::BestOverlap,
        ] {
            let router = router(
                vec![
                    region(0, "a", b"chr1", &[(0, 10_000)]),
                    region(1, "b", b"chr1", &[(0, 10_000)]),
                    region(2, "c", b"chr1", &[(0, 10_000)]),
                ],
                assignment,
            );
            let route = route_of(&router, &record(1000, &[(100, 0)]));
            assert_eq!(route.emission_count(), 1, "{assignment}");
            assert!(!router.may_duplicate(), "{assignment}");
            assert_eq!(router.stats().duplicate_emissions, 0, "{assignment}");
        }
    }

    #[test]
    fn routing_is_reproducible() {
        let build = || {
            router(
                vec![
                    region(0, "a", b"chr1", &[(1000, 1500)]),
                    region(1, "b", b"chr1", &[(1200, 1700)]),
                    region(2, "c", b"chr1", &[(1000, 1500)]),
                ],
                AssignmentMode::BestOverlap,
            )
        };
        let first = build();
        let second = build();
        for position in (900..1800).step_by(37) {
            let body = record(position, &[(100, 0)]);
            assert_eq!(
                names(&route_of(&first, &body)),
                names(&route_of(&second, &body)),
                "position {position}"
            );
        }
    }

    #[test]
    fn assignment_modes_round_trip_their_option_values() {
        for mode in [
            AssignmentMode::Start,
            AssignmentMode::Midpoint,
            AssignmentMode::Contained,
            AssignmentMode::BestOverlap,
            AssignmentMode::Overlap,
        ] {
            assert_eq!(AssignmentMode::parse(mode.as_str()), Some(mode));
        }
        assert_eq!(AssignmentMode::parse("nonesuch"), None);
        assert!(AssignmentMode::Overlap.may_duplicate());
        assert!(!AssignmentMode::BestOverlap.may_duplicate());
    }
}
