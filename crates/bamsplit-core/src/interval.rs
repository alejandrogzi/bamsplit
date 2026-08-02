// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! A single half-open genomic interval type, shared by every subsystem.
//!
//! BAM positions, BED coordinates, and the `genepred` model all use 0-based
//! half-open intervals, so `bamsplit` does too: `[start, end)`. Converting to
//! the 1-based inclusive coordinates that SAM text and index queries want
//! happens exactly once, at the boundary, via
//! [`one_based_inclusive`](Interval::one_based_inclusive).
//!
//! Coordinates are `i64` so that arithmetic on BAM's signed 32-bit positions
//! cannot overflow, and so an unset position can be represented without an
//! extra `Option` layer inside hot loops.
//!
//! # Example
//!
//! ```
//! use bamsplit_core::interval::Interval;
//!
//! let exon = Interval::new(100, 200);
//! let read = Interval::new(150, 250);
//!
//! assert_eq!(exon.overlap_len(read), 50);
//! assert!(exon.intersects(read));
//! assert!(!exon.contains_interval(read));
//! assert_eq!(exon.one_based_inclusive(), (101, 200));
//! ```

use std::cmp::Ordering;

/// A half-open interval `[start, end)` on one reference sequence.
///
/// An interval with `end <= start` is *empty*: it overlaps nothing and contains
/// nothing. Empty intervals are permitted so that degenerate annotation records
/// can be represented and reported rather than silently dropped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Interval {
    /// 0-based inclusive start.
    pub start: i64,
    /// 0-based exclusive end.
    pub end: i64,
}

impl Interval {
    /// Creates an interval.
    #[must_use]
    pub const fn new(start: i64, end: i64) -> Self {
        Self { start, end }
    }

    /// An empty interval anchored at `position`.
    #[must_use]
    pub const fn empty_at(position: i64) -> Self {
        Self {
            start: position,
            end: position,
        }
    }

    /// The number of positions covered, saturating at zero.
    #[must_use]
    pub const fn len(&self) -> i64 {
        if self.end > self.start {
            self.end - self.start
        } else {
            0
        }
    }

    /// Whether the interval covers no positions.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.end <= self.start
    }

    /// Whether `position` falls inside the interval.
    #[must_use]
    pub const fn contains(&self, position: i64) -> bool {
        position >= self.start && position < self.end
    }

    /// Whether `other` is entirely inside this interval.
    ///
    /// An empty `other` is contained when its anchor is inside, or when it sits
    /// exactly at this interval's end and this interval is itself empty.
    #[must_use]
    pub const fn contains_interval(&self, other: Self) -> bool {
        if other.is_empty() {
            return self.contains(other.start) || (self.is_empty() && self.start == other.start);
        }
        other.start >= self.start && other.end <= self.end
    }

    /// Whether the two intervals share at least one position.
    #[must_use]
    pub const fn intersects(&self, other: Self) -> bool {
        self.start < other.end && other.start < self.end && !self.is_empty() && !other.is_empty()
    }

    /// The number of positions the two intervals share.
    #[must_use]
    pub const fn overlap_len(&self, other: Self) -> i64 {
        let start = if self.start > other.start {
            self.start
        } else {
            other.start
        };
        let end = if self.end < other.end {
            self.end
        } else {
            other.end
        };
        if end > start { end - start } else { 0 }
    }

    /// The intersection, or [`None`] when the intervals are disjoint.
    #[must_use]
    pub fn intersection(&self, other: Self) -> Option<Self> {
        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        (end > start).then_some(Self { start, end })
    }

    /// The smallest interval containing both.
    #[must_use]
    pub fn envelope(&self, other: Self) -> Self {
        Self {
            start: self.start.min(other.start),
            end: self.end.max(other.end),
        }
    }

    /// Whether the two intervals touch or overlap, so they could be merged.
    #[must_use]
    pub const fn is_adjacent_or_overlapping(&self, other: Self) -> bool {
        self.start <= other.end && other.start <= self.end
    }

    /// The interval as 1-based inclusive `(start, end)`, the convention SAM
    /// text and binning-index queries use.
    ///
    /// For an empty interval both coordinates collapse onto the anchor, which
    /// keeps `start <= end` and avoids constructing an invalid query.
    #[must_use]
    pub const fn one_based_inclusive(&self) -> (i64, i64) {
        if self.is_empty() {
            (self.start + 1, self.start + 1)
        } else {
            (self.start + 1, self.end)
        }
    }

    /// The midpoint, rounding down.
    ///
    /// Defined for empty intervals as the anchor itself.
    #[must_use]
    pub const fn midpoint(&self) -> i64 {
        if self.is_empty() {
            self.start
        } else {
            self.start + (self.end - self.start) / 2
        }
    }

    /// Deterministic ordering: by start, then by end.
    #[must_use]
    pub fn cmp_by_position(&self, other: &Self) -> Ordering {
        self.start
            .cmp(&other.start)
            .then_with(|| self.end.cmp(&other.end))
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_by_position(other)
    }
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {})", self.start, self.end)
    }
}

/// Sorts and merges a set of intervals in place, collapsing overlaps and
/// abutments.
///
/// Empty intervals are discarded, because they can neither be overlapped nor
/// contain anything and would otherwise create spurious merge boundaries.
pub fn merge_in_place(intervals: &mut Vec<Interval>) {
    intervals.retain(|interval| !interval.is_empty());
    if intervals.len() < 2 {
        return;
    }
    intervals.sort_unstable();

    let mut write = 0usize;
    for read in 1..intervals.len() {
        let current = intervals[read];
        if current.start <= intervals[write].end {
            intervals[write].end = intervals[write].end.max(current.end);
        } else {
            write += 1;
            intervals[write] = current;
        }
    }
    intervals.truncate(write + 1);
}

/// The total number of positions covered by a **sorted, merged** interval set.
#[must_use]
pub fn total_len(intervals: &[Interval]) -> i64 {
    intervals.iter().map(Interval::len).sum()
}

/// The number of positions shared between two **sorted, merged** interval sets.
///
/// Runs in `O(n + m)` with a two-pointer sweep rather than the `O(n * m)` naive
/// product, which matters because this is called once per candidate region per
/// BAM record in `best-overlap` mode.
#[must_use]
pub fn overlap_len_sorted(left: &[Interval], right: &[Interval]) -> i64 {
    let (mut i, mut j) = (0usize, 0usize);
    let mut total = 0i64;
    while i < left.len() && j < right.len() {
        total += left[i].overlap_len(right[j]);
        if left[i].end < right[j].end {
            i += 1;
        } else {
            j += 1;
        }
    }
    total
}

/// Whether every interval in `parts` lies inside the **sorted, merged** set
/// `whole`.
#[must_use]
pub fn contains_all_sorted(whole: &[Interval], parts: &[Interval]) -> bool {
    let mut cursor = 0usize;
    for part in parts {
        if part.is_empty() {
            continue;
        }
        while cursor < whole.len() && whole[cursor].end <= part.start {
            cursor += 1;
        }
        match whole.get(cursor) {
            Some(candidate) if candidate.contains_interval(*part) => {}
            _ => return false,
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn length_and_emptiness() {
        assert_eq!(Interval::new(10, 20).len(), 10);
        assert_eq!(Interval::new(10, 10).len(), 0);
        assert_eq!(Interval::new(20, 10).len(), 0);
        assert!(Interval::new(10, 10).is_empty());
        assert!(Interval::new(20, 10).is_empty());
    }

    #[test]
    fn containment_and_intersection() {
        let outer = Interval::new(0, 100);
        assert!(outer.contains(0));
        assert!(!outer.contains(100));
        assert!(outer.contains_interval(Interval::new(10, 90)));
        assert!(outer.contains_interval(Interval::new(0, 100)));
        assert!(!outer.contains_interval(Interval::new(0, 101)));
        assert!(outer.intersects(Interval::new(99, 200)));
        assert!(!outer.intersects(Interval::new(100, 200)));
        assert!(!outer.intersects(Interval::new(50, 50)));
    }

    #[test]
    fn overlap_arithmetic() {
        assert_eq!(Interval::new(0, 10).overlap_len(Interval::new(5, 15)), 5);
        assert_eq!(Interval::new(0, 10).overlap_len(Interval::new(10, 15)), 0);
        assert_eq!(
            Interval::new(0, 10).intersection(Interval::new(5, 15)),
            Some(Interval::new(5, 10))
        );
        assert_eq!(
            Interval::new(0, 10).intersection(Interval::new(10, 15)),
            None
        );
    }

    #[test]
    fn coordinate_conversion() {
        assert_eq!(Interval::new(0, 10).one_based_inclusive(), (1, 10));
        assert_eq!(Interval::new(99, 100).one_based_inclusive(), (100, 100));
        assert_eq!(Interval::new(5, 5).one_based_inclusive(), (6, 6));
    }

    #[test]
    fn midpoints_round_down() {
        assert_eq!(Interval::new(0, 10).midpoint(), 5);
        assert_eq!(Interval::new(0, 11).midpoint(), 5);
        assert_eq!(Interval::new(7, 7).midpoint(), 7);
    }

    #[test]
    fn merging_collapses_overlaps_and_abutments() {
        let mut intervals = vec![
            Interval::new(30, 40),
            Interval::new(0, 10),
            Interval::new(10, 20),
            Interval::new(5, 8),
            Interval::new(100, 100),
        ];
        merge_in_place(&mut intervals);
        assert_eq!(intervals, vec![Interval::new(0, 20), Interval::new(30, 40)]);
    }

    #[test]
    fn merging_is_idempotent() {
        let mut intervals = vec![Interval::new(0, 10), Interval::new(20, 30)];
        let expected = intervals.clone();
        merge_in_place(&mut intervals);
        assert_eq!(intervals, expected);
        merge_in_place(&mut intervals);
        assert_eq!(intervals, expected);
    }

    #[test]
    fn sorted_overlap_matches_the_naive_product() {
        let left = vec![
            Interval::new(0, 10),
            Interval::new(20, 30),
            Interval::new(40, 50),
        ];
        let right = vec![Interval::new(5, 25), Interval::new(45, 60)];
        let naive: i64 = left
            .iter()
            .flat_map(|l| right.iter().map(move |r| l.overlap_len(*r)))
            .sum();
        assert_eq!(overlap_len_sorted(&left, &right), naive);
        assert_eq!(overlap_len_sorted(&left, &right), 5 + 5 + 5);
    }

    #[test]
    fn sorted_overlap_handles_empty_sets() {
        assert_eq!(overlap_len_sorted(&[], &[Interval::new(0, 10)]), 0);
        assert_eq!(overlap_len_sorted(&[Interval::new(0, 10)], &[]), 0);
    }

    #[test]
    fn containment_across_a_merged_set() {
        let whole = vec![Interval::new(0, 100), Interval::new(200, 300)];
        assert!(contains_all_sorted(
            &whole,
            &[Interval::new(10, 20), Interval::new(250, 260)]
        ));
        assert!(!contains_all_sorted(
            &whole,
            &[Interval::new(10, 20), Interval::new(150, 160)]
        ));
        // A block straddling the gap is not contained.
        assert!(!contains_all_sorted(&whole, &[Interval::new(90, 210)]));
        assert!(contains_all_sorted(&whole, &[]));
    }

    #[test]
    fn total_length_sums_a_merged_set() {
        assert_eq!(
            total_len(&[Interval::new(0, 10), Interval::new(20, 25)]),
            15
        );
        assert_eq!(total_len(&[]), 0);
    }

    #[test]
    fn ordering_is_by_start_then_end() {
        let mut intervals = vec![
            Interval::new(10, 20),
            Interval::new(0, 100),
            Interval::new(0, 50),
        ];
        intervals.sort();
        assert_eq!(
            intervals,
            vec![
                Interval::new(0, 50),
                Interval::new(0, 100),
                Interval::new(10, 20)
            ]
        );
    }
}
