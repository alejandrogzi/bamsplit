// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! A per-reference interval index over logical regions.
//!
//! # The shape
//!
//! Per reference, regions are sorted by envelope start and each entry carries a
//! running **prefix maximum** of the envelope ends seen so far:
//!
//! ```text
//! entry     0        1        2        3        4
//! start    100      150      400      900     1200
//! end      300      800      500     1000     1400
//! max_end  300      800      800     1000     1400
//!                    ▲
//!                    └── entry 1 reaches further right than entry 2
//! ```
//!
//! A query for `[950, 960)` binary-searches for the first entry starting at or
//! after 960 — entry 4 — then walks *backwards* while `max_end > 950`. That
//! stops at entry 2 (`max_end = 800 <= 950`), so entries 3 and 2 are examined
//! and the first two are never touched.
//!
//! # Why not a tree
//!
//! A centered interval tree has the same asymptotics and better worst-case
//! behaviour, but this structure is two flat `Vec`s: it builds with one sort,
//! has no pointer chasing, and scans contiguous memory. For the query pattern
//! here — millions of short queries against a static set — that wins in
//! practice.
//!
//! The honest weakness: the backwards walk stops at the first entry whose prefix
//! maximum falls behind the query, so **one region with a very wide envelope
//! poisons every query to its right**. A gene spanning a whole chromosome would
//! force every later query to walk back to it. Real annotations do not look like
//! that; a file that does would be better served by a tree, and swapping the
//! implementation is a local change behind [`IntervalIndex::query`].
//!
//! # Envelopes are a filter, not an answer
//!
//! An envelope covers a transcript's introns as well as its exons, so a hit here
//! only means "worth checking". The caller verifies against the actual segments.

use std::collections::HashMap;

use crate::annotation::features::LogicalRegion;
use crate::error::AnnotationError;
use crate::interval::Interval;

/// One indexed region, flattened for cache-friendly scanning.
#[derive(Debug, Clone, Copy)]
struct Entry {
    start: i64,
    end: i64,
    /// The largest `end` among this entry and every earlier one.
    max_end: i64,
    /// Index into [`IntervalIndex::regions`].
    region: u32,
}

/// The entries of one reference sequence, sorted by envelope start.
#[derive(Debug, Default)]
struct ReferenceIndex {
    entries: Vec<Entry>,
}

/// A static interval index over every logical region.
#[derive(Debug)]
pub struct IntervalIndex {
    regions: Vec<LogicalRegion>,
    by_reference: HashMap<Vec<u8>, ReferenceIndex>,
    covered_bases: u64,
}

impl IntervalIndex {
    /// Builds the index.
    ///
    /// Regions with no segments are indexed too — `--emit-empty` still needs
    /// them to exist — but with an empty envelope, so they match nothing.
    ///
    /// # Errors
    ///
    /// Returns [`AnnotationError::InvalidCoordinates`] if a region's envelope is
    /// reversed, which would silently break the binary search.
    pub fn build(regions: Vec<LogicalRegion>) -> Result<Self, AnnotationError> {
        let mut by_reference: HashMap<Vec<u8>, ReferenceIndex> = HashMap::new();
        let mut covered_bases = 0u64;

        for (index, region) in regions.iter().enumerate() {
            covered_bases += region.covered_bases().max(0) as u64;
            if region.segments.is_empty() {
                continue;
            }
            if region.envelope.end < region.envelope.start {
                return Err(AnnotationError::InvalidCoordinates {
                    name: String::from_utf8_lossy(region.key.name()).into_owned(),
                    chrom: String::from_utf8_lossy(&region.chrom).into_owned(),
                    ordinal: region.source_ordinal,
                    start: region.envelope.start.max(0) as u64,
                    end: region.envelope.end.max(0) as u64,
                    reason: "the envelope is reversed",
                });
            }
            by_reference
                .entry(region.chrom.clone())
                .or_default()
                .entries
                .push(Entry {
                    start: region.envelope.start,
                    end: region.envelope.end,
                    max_end: region.envelope.end,
                    region: u32::try_from(index).unwrap_or(u32::MAX),
                });
        }

        for reference in by_reference.values_mut() {
            // Sort by start, then end, then region index: a total order, so the
            // index is byte-identical across runs.
            reference.entries.sort_unstable_by(|left, right| {
                left.start
                    .cmp(&right.start)
                    .then(left.end.cmp(&right.end))
                    .then(left.region.cmp(&right.region))
            });
            let mut running = i64::MIN;
            for entry in &mut reference.entries {
                running = running.max(entry.end);
                entry.max_end = running;
            }
        }

        Ok(Self {
            regions,
            by_reference,
            covered_bases,
        })
    }

    /// How many regions are indexed, including empty ones.
    #[must_use]
    pub fn region_count(&self) -> usize {
        self.regions.len()
    }

    /// Whether the index holds no regions.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }

    /// Every region, in annotation input order.
    #[must_use]
    pub fn regions(&self) -> &[LogicalRegion] {
        &self.regions
    }

    /// One region by its index.
    #[must_use]
    pub fn region(&self, index: usize) -> Option<&LogicalRegion> {
        self.regions.get(index)
    }

    /// The reference names the annotation mentions.
    pub fn reference_names(&self) -> impl Iterator<Item = &[u8]> {
        self.by_reference.keys().map(Vec::as_slice)
    }

    /// The total number of reference bases every region covers.
    ///
    /// Used by the planner to decide whether an annotation is dense enough that
    /// one sequential pass beats many index queries.
    #[must_use]
    pub const fn covered_bases(&self) -> u64 {
        self.covered_bases
    }

    /// Collects the regions on `chrom` whose envelope intersects `query`.
    ///
    /// `out` is cleared first and filled with indices into [`Self::regions`], in
    /// annotation input order so every downstream tie-break is deterministic.
    pub fn query(&self, chrom: &[u8], query: Interval, out: &mut Vec<usize>) {
        out.clear();
        if query.is_empty() {
            return;
        }
        let Some(reference) = self.by_reference.get(chrom) else {
            return;
        };

        let entries = &reference.entries;
        // The first entry that starts at or after the query's end cannot
        // overlap, and neither can anything after it.
        let upper = entries.partition_point(|entry| entry.start < query.end);

        // Walk back while some earlier entry could still reach into the query.
        for entry in entries[..upper].iter().rev() {
            if entry.max_end <= query.start {
                break;
            }
            if entry.end > query.start {
                out.push(entry.region as usize);
            }
        }
        out.sort_unstable();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::features::{RegionKey, RegionMetadata};

    fn region(ordinal: u64, chrom: &[u8], segments: &[(i64, i64)]) -> LogicalRegion {
        let segments: Vec<Interval> = segments
            .iter()
            .map(|(start, end)| Interval::new(*start, *end))
            .collect();
        let envelope = segments.iter().fold(
            Interval::empty_at(segments.first().map_or(0, |first| first.start)),
            |envelope, segment| {
                if envelope.is_empty() {
                    *segment
                } else {
                    envelope.envelope(*segment)
                }
            },
        );
        let name = format!("r{ordinal}").into_bytes();
        LogicalRegion {
            key: RegionKey::new(ordinal, name.clone()),
            chrom: chrom.to_vec(),
            segments,
            envelope,
            source_ordinal: ordinal,
            metadata: RegionMetadata {
                original_name: name,
                generated_name: false,
                duplicate_suffix: 0,
                span_as_exon_fallback: false,
                strand: None,
                derived_segment_count: 0,
            },
        }
    }

    /// The obvious, obviously-correct implementation, for differential testing.
    fn slow_query(index: &IntervalIndex, chrom: &[u8], query: Interval) -> Vec<usize> {
        index
            .regions()
            .iter()
            .enumerate()
            .filter(|(_, region)| {
                region.chrom == chrom
                    && !region.segments.is_empty()
                    && region.envelope.intersects(query)
            })
            .map(|(index, _)| index)
            .collect()
    }

    #[test]
    fn an_empty_index_answers_nothing() {
        let index = IntervalIndex::build(Vec::new()).expect("built");
        assert!(index.is_empty());
        assert_eq!(index.region_count(), 0);
        let mut out = Vec::new();
        index.query(b"chr1", Interval::new(0, 100), &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn a_single_region_is_found_and_bounded() {
        let index = IntervalIndex::build(vec![region(0, b"chr1", &[(100, 200)])]).expect("built");
        let mut out = Vec::new();

        for (query, expected) in [
            (Interval::new(150, 160), vec![0]),
            (Interval::new(0, 101), vec![0]),
            (Interval::new(199, 500), vec![0]),
            (Interval::new(0, 100), vec![]),
            (Interval::new(200, 300), vec![]),
        ] {
            index.query(b"chr1", query, &mut out);
            assert_eq!(out, expected, "{query}");
        }
    }

    #[test]
    fn an_unknown_reference_answers_nothing() {
        let index = IntervalIndex::build(vec![region(0, b"chr1", &[(100, 200)])]).expect("built");
        let mut out = Vec::new();
        index.query(b"chrZ", Interval::new(100, 200), &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn nested_and_identical_envelopes_are_all_found() {
        let index = IntervalIndex::build(vec![
            region(0, b"chr1", &[(100, 1000)]),
            region(1, b"chr1", &[(200, 300)]),
            region(2, b"chr1", &[(200, 300)]),
            region(3, b"chr1", &[(250, 260)]),
        ])
        .expect("built");
        let mut out = Vec::new();
        index.query(b"chr1", Interval::new(255, 256), &mut out);
        assert_eq!(out, [0, 1, 2, 3]);
    }

    #[test]
    fn results_are_in_annotation_input_order() {
        // Deliberately built with envelopes in the opposite order to the
        // ordinals, so a sort by position would give a different answer.
        let index = IntervalIndex::build(vec![
            region(0, b"chr1", &[(900, 1000)]),
            region(1, b"chr1", &[(500, 600)]),
            region(2, b"chr1", &[(100, 200)]),
        ])
        .expect("built");
        let mut out = Vec::new();
        index.query(b"chr1", Interval::new(0, 2000), &mut out);
        assert_eq!(out, [0, 1, 2]);
    }

    #[test]
    fn an_empty_region_is_kept_but_never_matches() {
        let mut empty = region(1, b"chr1", &[]);
        empty.envelope = Interval::empty_at(150);
        let index =
            IntervalIndex::build(vec![region(0, b"chr1", &[(100, 200)]), empty]).expect("built");
        assert_eq!(index.region_count(), 2);
        let mut out = Vec::new();
        index.query(b"chr1", Interval::new(140, 160), &mut out);
        assert_eq!(out, [0]);
    }

    #[test]
    fn an_empty_query_matches_nothing() {
        let index = IntervalIndex::build(vec![region(0, b"chr1", &[(100, 200)])]).expect("built");
        let mut out = Vec::new();
        index.query(b"chr1", Interval::new(150, 150), &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn covered_bases_sum_every_segment() {
        let index = IntervalIndex::build(vec![
            region(0, b"chr1", &[(0, 100), (200, 250)]),
            region(1, b"chr2", &[(0, 10)]),
        ])
        .expect("built");
        assert_eq!(index.covered_bases(), 160);
    }

    #[test]
    fn reference_names_cover_every_annotated_contig() {
        let index = IntervalIndex::build(vec![
            region(0, b"chr1", &[(0, 10)]),
            region(1, b"chr2", &[(0, 10)]),
            region(2, b"chr1", &[(20, 30)]),
        ])
        .expect("built");
        let mut names: Vec<Vec<u8>> = index.reference_names().map(<[u8]>::to_vec).collect();
        names.sort();
        assert_eq!(names, [b"chr1".to_vec(), b"chr2".to_vec()]);
    }

    #[test]
    fn the_fast_query_agrees_with_a_linear_scan() {
        // A deterministic pseudo-random annotation: nested, overlapping, and
        // disjoint regions across three references.
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let references: [&[u8]; 3] = [b"chr1", b"chr2", b"chr3"];
        let regions: Vec<LogicalRegion> = (0..400u64)
            .map(|ordinal| {
                let chrom = references[(next() % 3) as usize];
                let start = (next() % 10_000) as i64;
                let width = 1 + (next() % 900) as i64;
                region(ordinal, chrom, &[(start, start + width)])
            })
            .collect();

        let index = IntervalIndex::build(regions).expect("built");
        let mut out = Vec::new();
        for _ in 0..3_000 {
            let chrom = references[(next() % 3) as usize];
            let start = (next() % 11_000) as i64;
            let width = 1 + (next() % 500) as i64;
            let query = Interval::new(start, start + width);

            index.query(chrom, query, &mut out);
            let mut expected = slow_query(&index, chrom, query);
            expected.sort_unstable();
            assert_eq!(out, expected, "{} {query}", String::from_utf8_lossy(chrom));
        }
    }

    #[test]
    fn a_reversed_envelope_is_rejected() {
        let mut broken = region(0, b"chr1", &[(100, 200)]);
        broken.envelope = Interval::new(200, 100);
        let error = IntervalIndex::build(vec![broken]).expect_err("must reject");
        assert!(
            matches!(error, AnnotationError::InvalidCoordinates { .. }),
            "{error}"
        );
    }

    #[test]
    fn a_wide_region_does_not_hide_narrow_ones() {
        // The documented weakness: a whole-chromosome region forces the walk
        // back. Correctness must not suffer even though speed does.
        let mut regions = vec![region(0, b"chr1", &[(0, 1_000_000)])];
        regions.extend((1..50u64).map(|ordinal| {
            let start = (ordinal as i64) * 1_000;
            region(ordinal, b"chr1", &[(start, start + 100)])
        }));
        let index = IntervalIndex::build(regions).expect("built");

        let mut out = Vec::new();
        let query = Interval::new(40_000, 40_050);
        index.query(b"chr1", query, &mut out);
        let mut expected = slow_query(&index, b"chr1", query);
        expected.sort_unstable();
        assert_eq!(out, expected);
        assert!(out.contains(&0), "the wide region must still be reported");
    }
}
