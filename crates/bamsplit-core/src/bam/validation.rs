// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Cross-checks between records and the header, plus the counters that
//! `bamsplit inspect --full` and the manifest report.
//!
//! Everything here is *observational*: it never mutates a record and never
//! decides where one goes. Routers call [`check_reference_ids`] because an
//! out-of-range `ref_id` is an input error they must reject; the rest is used
//! by inspection and by the conservation checks.

use crate::bam::header::BamHeader;
use crate::bam::raw_record::{RawRecord, RecordLocation};
use crate::error::BamRecordError;

/// Verifies that a record's reference identifiers address the header
/// dictionary.
///
/// A `ref_id` of `-1` is legal (the record is unplaced) and is not an error.
/// Any other value must be a valid index; anything else means the BAM and its
/// header disagree, which `bamsplit` treats as malformed input rather than
/// silently routing the record somewhere arbitrary.
///
/// # Errors
///
/// Returns [`BamRecordError::InvalidReferenceSequenceId`] when `ref_id` or
/// `next_ref_id` is out of range, located at `location`.
pub fn check_reference_ids(
    header: &BamHeader,
    record: &RawRecord<'_>,
    location: RecordLocation,
) -> Result<(), BamRecordError> {
    let count = header.reference_count();
    for (kind, id) in [
        ("reference", record.reference_sequence_id()),
        ("mate reference", record.mate_reference_sequence_id()),
    ] {
        let id = id.map_err(|error| error.relocate(location))?;
        if let Some(id) = id
            && !header.is_valid_reference_id(id)
        {
            return Err(BamRecordError::InvalidReferenceSequenceId {
                kind,
                id,
                reference_count: count,
                location,
            });
        }
    }
    Ok(())
}

/// Where a record sits in the placed/unplaced/mapped taxonomy.
///
/// The distinction that matters for `bamsplit chrom` is between *placed
/// unmapped* and *unplaced*: a record can carry the `UNMAPPED` flag and still
/// have a perfectly good `ref_id` and position, because coordinate-sorted BAMs
/// park an unmapped read next to its mapped mate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Placement {
    /// Has a reference id and is not flagged unmapped.
    Mapped,
    /// Has a reference id but is flagged unmapped.
    PlacedUnmapped,
    /// Has no reference id (`ref_id == -1`).
    Unplaced,
}

impl Placement {
    /// Classifies a record.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the fixed core cannot be read.
    pub fn of(record: &RawRecord<'_>) -> Result<Self, BamRecordError> {
        Ok(match record.reference_sequence_id()? {
            None => Self::Unplaced,
            Some(_) if record.is_unmapped()? => Self::PlacedUnmapped,
            Some(_) => Self::Mapped,
        })
    }

    /// A short label for manifests and reports.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Mapped => "mapped",
            Self::PlacedUnmapped => "placed-unmapped",
            Self::Unplaced => "unplaced-unmapped",
        }
    }
}

/// Running counts over a stream of records.
///
/// Kept as plain `u64` fields rather than a map so that the hot loop is a
/// handful of increments.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RecordCounters {
    /// Every record seen.
    pub total: u64,
    /// Neither secondary nor supplementary.
    pub primary: u64,
    /// `SECONDARY` set.
    pub secondary: u64,
    /// `SUPPLEMENTARY` set.
    pub supplementary: u64,
    /// Placed and not flagged unmapped.
    pub mapped: u64,
    /// Placed but flagged unmapped.
    pub placed_unmapped: u64,
    /// Not placed.
    pub unplaced_unmapped: u64,
    /// `DUPLICATE` set.
    pub duplicate: u64,
    /// `QC_FAIL` set.
    pub qc_fail: u64,
}

impl RecordCounters {
    /// Folds one record into the counters.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the fixed core cannot be read.
    pub fn observe(&mut self, record: &RawRecord<'_>) -> Result<(), BamRecordError> {
        use crate::bam::flags;

        let flags_value = record.flags()?;
        self.total += 1;
        if flags_value & flags::SECONDARY != 0 {
            self.secondary += 1;
        } else if flags_value & flags::SUPPLEMENTARY != 0 {
            self.supplementary += 1;
        } else {
            self.primary += 1;
        }
        if flags_value & flags::DUPLICATE != 0 {
            self.duplicate += 1;
        }
        if flags_value & flags::QC_FAIL != 0 {
            self.qc_fail += 1;
        }
        match Placement::of(record)? {
            Placement::Mapped => self.mapped += 1,
            Placement::PlacedUnmapped => self.placed_unmapped += 1,
            Placement::Unplaced => self.unplaced_unmapped += 1,
        }
        Ok(())
    }

    /// Adds another set of counters.
    pub fn merge(&mut self, other: &Self) {
        self.total += other.total;
        self.primary += other.primary;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.mapped += other.mapped;
        self.placed_unmapped += other.placed_unmapped;
        self.unplaced_unmapped += other.unplaced_unmapped;
        self.duplicate += other.duplicate;
        self.qc_fail += other.qc_fail;
    }

    /// Whether the placement categories sum to the total.
    ///
    /// A `false` here means the counters were corrupted, which the manifest
    /// turns into a validation failure.
    #[must_use]
    pub const fn is_consistent(&self) -> bool {
        self.mapped + self.placed_unmapped + self.unplaced_unmapped == self.total
            && self.primary + self.secondary + self.supplementary == self.total
    }
}

/// Tracks whether a stream really is in coordinate order.
///
/// The stream engine relies on grouping, and grouping is only guaranteed when
/// the input is genuinely sorted — a header that *claims* `SO:coordinate` is
/// not enough. This tracker is one comparison per record.
///
/// Unplaced records are expected to arrive in one run at the very end of a
/// coordinate-sorted BAM; a placed record appearing *after* an unplaced one is
/// therefore a violation too.
#[derive(Debug, Clone, Copy, Default)]
pub struct CoordinateOrderTracker {
    previous_reference_id: Option<i32>,
    previous_position: i32,
    seen_unplaced: bool,
    violations: u64,
    first_violation: Option<OrderViolation>,
}

/// The details of the first ordering violation seen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrderViolation {
    /// Ordinal of the offending record.
    pub record_number: u64,
    /// Its reference id, or `-1` when unplaced.
    pub reference_id: i32,
    /// Its 0-based position, or `-1` when absent.
    pub position: i32,
    /// The previous record's reference id, or `-1`.
    pub previous_reference_id: i32,
    /// The previous record's 0-based position, or `-1`.
    pub previous_position: i32,
}

impl CoordinateOrderTracker {
    /// Creates a tracker.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            previous_reference_id: None,
            previous_position: -1,
            seen_unplaced: false,
            violations: 0,
            first_violation: None,
        }
    }

    /// Folds one record in, returning whether it was in order.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the fixed core cannot be read.
    pub fn observe(
        &mut self,
        record: &RawRecord<'_>,
        record_number: u64,
    ) -> Result<bool, BamRecordError> {
        let reference_id = record.reference_sequence_id()?;
        let position = record.alignment_start()?.unwrap_or(-1);

        let in_order = match (self.previous_reference_id, reference_id) {
            // The first record is always in order.
            (None, _) if self.violations == 0 && !self.seen_unplaced => true,
            (_, None) => true,
            (Some(_), Some(_)) if self.seen_unplaced => false,
            (Some(previous_id), Some(current_id)) => match current_id.cmp(&previous_id) {
                std::cmp::Ordering::Less => false,
                std::cmp::Ordering::Greater => true,
                std::cmp::Ordering::Equal => position >= self.previous_position,
            },
            (None, Some(_)) => !self.seen_unplaced,
        };

        if !in_order {
            self.violations += 1;
            if self.first_violation.is_none() {
                self.first_violation = Some(OrderViolation {
                    record_number,
                    reference_id: reference_id.unwrap_or(-1),
                    position,
                    previous_reference_id: self.previous_reference_id.unwrap_or(-1),
                    previous_position: self.previous_position,
                });
            }
        }

        match reference_id {
            Some(id) => {
                self.previous_reference_id = Some(id);
                self.previous_position = position;
            }
            None => self.seen_unplaced = true,
        }
        Ok(in_order)
    }

    /// How many out-of-order records were seen.
    #[must_use]
    pub const fn violations(&self) -> u64 {
        self.violations
    }

    /// The first violation, if any.
    #[must_use]
    pub const fn first_violation(&self) -> Option<OrderViolation> {
        self.first_violation
    }

    /// Whether the stream was fully in coordinate order.
    #[must_use]
    pub const fn is_sorted(&self) -> bool {
        self.violations == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::header::ReferenceSequence;

    fn header(reference_count: usize) -> BamHeader {
        let references = (0..reference_count)
            .map(|index| ReferenceSequence {
                name: format!("chr{}", index + 1).into_bytes(),
                length: 1_000_000,
            })
            .collect();
        BamHeader::from_parts(b"@HD\tVN:1.6\tSO:coordinate\n".to_vec(), references)
            .expect("valid header")
    }

    fn record(reference_id: i32, position: i32, flags: u16) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&reference_id.to_le_bytes());
        body.extend_from_slice(&position.to_le_bytes());
        body.push(2);
        body.push(60);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&1u16.to_le_bytes());
        body.extend_from_slice(&flags.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(b"r\0");
        body.extend_from_slice(&(10u32 << 4).to_le_bytes()); // 10M (op code 0)
        body
    }

    #[test]
    fn accepts_valid_and_absent_reference_ids() {
        let header = header(3);
        for id in [-1, 0, 1, 2] {
            let body = record(id, 0, 0);
            let parsed = RawRecord::new(&body).expect("valid");
            check_reference_ids(&header, &parsed, RecordLocation::new(1, 0)).expect("in range");
        }
    }

    #[test]
    fn rejects_an_out_of_range_reference_id() {
        let header = header(2);
        let body = record(5, 0, 0);
        let parsed = RawRecord::new(&body).expect("valid");
        let error =
            check_reference_ids(&header, &parsed, RecordLocation::new(9, 0)).expect_err("reject");
        assert!(
            matches!(
                error,
                BamRecordError::InvalidReferenceSequenceId {
                    kind: "reference",
                    id: 5,
                    reference_count: 2,
                    ..
                }
            ),
            "{error}"
        );
        assert_eq!(error.location().record_number, 9);
    }

    #[test]
    fn rejects_an_out_of_range_mate_reference_id() {
        let header = header(1);
        let mut body = record(0, 0, 0);
        body[20..24].clone_from_slice(&7i32.to_le_bytes());
        let parsed = RawRecord::new(&body).expect("valid");
        let error =
            check_reference_ids(&header, &parsed, RecordLocation::new(1, 0)).expect_err("reject");
        assert!(
            matches!(
                error,
                BamRecordError::InvalidReferenceSequenceId {
                    kind: "mate reference",
                    ..
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn classifies_placement() {
        let mapped = record(0, 100, 0);
        assert_eq!(
            Placement::of(&RawRecord::new(&mapped).expect("valid")).expect("valid"),
            Placement::Mapped
        );

        let placed_unmapped = record(0, 100, 0x4);
        assert_eq!(
            Placement::of(&RawRecord::new(&placed_unmapped).expect("valid")).expect("valid"),
            Placement::PlacedUnmapped
        );

        let unplaced = record(-1, -1, 0x4);
        assert_eq!(
            Placement::of(&RawRecord::new(&unplaced).expect("valid")).expect("valid"),
            Placement::Unplaced
        );
    }

    #[test]
    fn counters_stay_consistent() {
        let mut counters = RecordCounters::default();
        for body in [
            record(0, 1, 0),
            record(0, 2, 0x100),
            record(0, 3, 0x800),
            record(0, 4, 0x4),
            record(-1, -1, 0x4),
            record(0, 5, 0x400 | 0x200),
        ] {
            counters
                .observe(&RawRecord::new(&body).expect("valid"))
                .expect("valid");
        }
        assert_eq!(counters.total, 6);
        assert_eq!(counters.primary, 4);
        assert_eq!(counters.secondary, 1);
        assert_eq!(counters.supplementary, 1);
        assert_eq!(counters.mapped, 4);
        assert_eq!(counters.placed_unmapped, 1);
        assert_eq!(counters.unplaced_unmapped, 1);
        assert_eq!(counters.duplicate, 1);
        assert_eq!(counters.qc_fail, 1);
        assert!(counters.is_consistent());
    }

    #[test]
    fn merging_counters_sums_every_field() {
        let mut left = RecordCounters {
            total: 2,
            primary: 2,
            mapped: 2,
            ..RecordCounters::default()
        };
        let right = RecordCounters {
            total: 3,
            primary: 3,
            unplaced_unmapped: 3,
            ..RecordCounters::default()
        };
        left.merge(&right);
        assert_eq!(left.total, 5);
        assert_eq!(left.mapped, 2);
        assert_eq!(left.unplaced_unmapped, 3);
        assert!(left.is_consistent());
    }

    #[test]
    fn tracker_accepts_a_sorted_stream_with_trailing_unplaced_records() {
        let mut tracker = CoordinateOrderTracker::new();
        let bodies = [
            record(0, 10, 0),
            record(0, 10, 0),
            record(0, 200, 0),
            record(1, 5, 0),
            record(-1, -1, 0x4),
            record(-1, -1, 0x4),
        ];
        for (index, body) in bodies.iter().enumerate() {
            let parsed = RawRecord::new(body).expect("valid");
            assert!(
                tracker.observe(&parsed, index as u64 + 1).expect("valid"),
                "record {index} should be in order"
            );
        }
        assert!(tracker.is_sorted());
    }

    #[test]
    fn tracker_flags_a_backwards_position() {
        let mut tracker = CoordinateOrderTracker::new();
        for (index, body) in [record(0, 100, 0), record(0, 50, 0)].iter().enumerate() {
            let parsed = RawRecord::new(body).expect("valid");
            let _ = tracker.observe(&parsed, index as u64 + 1);
        }
        assert!(!tracker.is_sorted());
        let violation = tracker.first_violation().expect("recorded");
        assert_eq!(violation.record_number, 2);
        assert_eq!(violation.position, 50);
        assert_eq!(violation.previous_position, 100);
    }

    #[test]
    fn tracker_flags_a_backwards_reference() {
        let mut tracker = CoordinateOrderTracker::new();
        for (index, body) in [record(1, 10, 0), record(0, 10, 0)].iter().enumerate() {
            let parsed = RawRecord::new(body).expect("valid");
            let _ = tracker.observe(&parsed, index as u64 + 1);
        }
        assert_eq!(tracker.violations(), 1);
        let violation = tracker.first_violation().expect("recorded");
        assert_eq!(violation.reference_id, 0);
        assert_eq!(violation.previous_reference_id, 1);
    }

    #[test]
    fn tracker_flags_a_placed_record_after_an_unplaced_one() {
        let mut tracker = CoordinateOrderTracker::new();
        for (index, body) in [record(0, 10, 0), record(-1, -1, 0x4), record(0, 20, 0)]
            .iter()
            .enumerate()
        {
            let parsed = RawRecord::new(body).expect("valid");
            let _ = tracker.observe(&parsed, index as u64 + 1);
        }
        assert_eq!(tracker.violations(), 1);
        assert_eq!(
            tracker.first_violation().expect("recorded").record_number,
            3
        );
    }
}
