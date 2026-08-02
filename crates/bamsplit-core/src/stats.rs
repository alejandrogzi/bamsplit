// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The counters a run accumulates, and the digest that proves it was lossless.
//!
//! # The digest
//!
//! Counting records is not enough to show a split was lossless: two records
//! could be swapped, or one truncated, without changing a count. So every
//! output also carries an order-sensitive `xxh3_64` digest over the *record
//! bodies* it wrote — the exact bytes, excluding the `block_size` framing that
//! `bamsplit` recomputes.
//!
//! Because the digest is over bodies only, it is invariant to compression
//! level, block boundaries, and BGZF framing. That is what makes it directly
//! comparable with a digest computed from `samtools view -b`, which is how the
//! differential test-suite proves equivalence.
//!
//! The digest is *not* a cryptographic checksum. It detects the accidental
//! corruption a splitter can cause; it is not a defence against a deliberate
//! forgery.

use serde::{Deserialize, Serialize};

use crate::bam::raw_record::RawRecord;
use crate::error::BamRecordError;

/// An order-sensitive digest over a sequence of record bodies.
#[derive(Clone)]
pub struct RecordDigest {
    hasher: xxhash_rust::xxh3::Xxh3,
    records: u64,
}

impl std::fmt::Debug for RecordDigest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RecordDigest")
            .field("records", &self.records)
            .field("digest", &self.to_hex())
            .finish_non_exhaustive()
    }
}

impl Default for RecordDigest {
    fn default() -> Self {
        Self::new()
    }
}

impl RecordDigest {
    /// Creates an empty digest.
    #[must_use]
    pub fn new() -> Self {
        Self {
            hasher: xxhash_rust::xxh3::Xxh3::new(),
            records: 0,
        }
    }

    /// Folds one record body in.
    ///
    /// The body length is mixed in as well, so `["ab", "c"]` and `["a", "bc"]`
    /// produce different digests.
    pub fn update(&mut self, body: &[u8]) {
        self.hasher.update(&(body.len() as u64).to_le_bytes());
        self.hasher.update(body);
        self.records += 1;
    }

    /// How many records have been folded in.
    #[must_use]
    pub const fn record_count(&self) -> u64 {
        self.records
    }

    /// The digest so far.
    #[must_use]
    pub fn finish(&self) -> u64 {
        self.hasher.digest()
    }

    /// The digest rendered as 16 lower-case hex digits, for the manifest.
    #[must_use]
    pub fn to_hex(&self) -> String {
        format!("{:016x}", self.finish())
    }
}

/// A genomic coordinate, as recorded for the manifest.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Coordinate {
    /// The reference sequence id.
    pub reference_id: i32,
    /// The 0-based leftmost position.
    pub position: i32,
}

/// Per-output counters.
///
/// Every field is populated as records are written, so finalizing an output
/// costs nothing beyond reading these out.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct OutputStats {
    /// Records written.
    pub record_count: u64,
    /// Placed and not flagged unmapped.
    pub mapped_count: u64,
    /// Placed but flagged unmapped.
    pub placed_unmapped_count: u64,
    /// Not placed.
    pub unplaced_unmapped_count: u64,
    /// Neither secondary nor supplementary.
    pub primary_count: u64,
    /// `SECONDARY` set.
    pub secondary_count: u64,
    /// `SUPPLEMENTARY` set.
    pub supplementary_count: u64,
    /// `DUPLICATE` set.
    pub duplicate_count: u64,
    /// `QC_FAIL` set.
    pub qc_fail_count: u64,
    /// Total record-body bytes, before compression.
    pub uncompressed_bytes: u64,
    /// Size of the finished BAM on disk.
    pub compressed_bytes: u64,
    /// The first placed coordinate written.
    pub first_coordinate: Option<Coordinate>,
    /// The last placed coordinate written.
    pub last_coordinate: Option<Coordinate>,
    /// Whether every placed record was written in non-decreasing coordinate
    /// order.
    pub coordinate_sorted: bool,
    /// The record-body digest, as 16 hex digits.
    pub raw_record_digest: String,
}

impl OutputStats {
    /// Creates counters for an output that has not been written to.
    ///
    /// `coordinate_sorted` starts `true`: an empty or single-record output is
    /// trivially sorted, and the flag is only ever cleared by an observation
    /// that contradicts it.
    #[must_use]
    pub fn new() -> Self {
        Self {
            coordinate_sorted: true,
            ..Self::default()
        }
    }

    /// Whether the placement counts sum to the record count.
    #[must_use]
    pub const fn is_consistent(&self) -> bool {
        self.mapped_count + self.placed_unmapped_count + self.unplaced_unmapped_count
            == self.record_count
            && self.primary_count + self.secondary_count + self.supplementary_count
                == self.record_count
    }
}

/// Accumulates [`OutputStats`] while records are written.
#[derive(Debug)]
pub struct OutputStatsBuilder {
    stats: OutputStats,
    digest: RecordDigest,
    previous: Option<Coordinate>,
}

impl Default for OutputStatsBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl OutputStatsBuilder {
    /// Creates a builder.
    #[must_use]
    pub fn new() -> Self {
        Self {
            stats: OutputStats::new(),
            digest: RecordDigest::new(),
            previous: None,
        }
    }

    /// The counters accumulated so far.
    #[must_use]
    pub const fn stats(&self) -> &OutputStats {
        &self.stats
    }

    /// The digest accumulated so far.
    #[must_use]
    pub const fn digest(&self) -> &RecordDigest {
        &self.digest
    }

    /// How many records have been observed.
    #[must_use]
    pub const fn record_count(&self) -> u64 {
        self.stats.record_count
    }

    /// Folds one record in.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the fixed core cannot be read, which
    /// cannot happen for a record that came from
    /// [`RawRecord::new`](crate::bam::RawRecord::new).
    pub fn observe(&mut self, record: &RawRecord<'_>) -> Result<(), BamRecordError> {
        use crate::bam::flags;

        let body = record.raw_bytes();
        self.digest.update(body);
        self.stats.record_count += 1;
        self.stats.uncompressed_bytes += body.len() as u64;

        let flag_bits = record.flags()?;
        if flag_bits & flags::SECONDARY != 0 {
            self.stats.secondary_count += 1;
        } else if flag_bits & flags::SUPPLEMENTARY != 0 {
            self.stats.supplementary_count += 1;
        } else {
            self.stats.primary_count += 1;
        }
        if flag_bits & flags::DUPLICATE != 0 {
            self.stats.duplicate_count += 1;
        }
        if flag_bits & flags::QC_FAIL != 0 {
            self.stats.qc_fail_count += 1;
        }

        match record.reference_sequence_id()? {
            None => self.stats.unplaced_unmapped_count += 1,
            Some(reference_id) => {
                if flag_bits & flags::UNMAPPED != 0 {
                    self.stats.placed_unmapped_count += 1;
                } else {
                    self.stats.mapped_count += 1;
                }
                let coordinate = Coordinate {
                    reference_id,
                    position: record.alignment_start()?.unwrap_or(-1),
                };
                if self.stats.first_coordinate.is_none() {
                    self.stats.first_coordinate = Some(coordinate);
                }
                if let Some(previous) = self.previous
                    && (coordinate.reference_id, coordinate.position)
                        < (previous.reference_id, previous.position)
                {
                    self.stats.coordinate_sorted = false;
                }
                self.stats.last_coordinate = Some(coordinate);
                self.previous = Some(coordinate);
            }
        }
        Ok(())
    }

    /// Finalizes the counters, attaching the on-disk size.
    #[must_use]
    pub fn finish(mut self, compressed_bytes: u64) -> OutputStats {
        self.stats.compressed_bytes = compressed_bytes;
        self.stats.raw_record_digest = self.digest.to_hex();
        self.stats
    }
}

/// Run-wide counters, used for the conservation checks in the manifest.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct RunStats {
    /// Records read from the input.
    pub input_records: u64,
    /// Records emitted to at least one output.
    pub unique_emitted_records: u64,
    /// The sum of every output's record count.
    pub total_output_emissions: u64,
    /// Records deliberately discarded by a routing policy.
    pub dropped_records: u64,
    /// Records that matched no output key.
    pub unmatched_records: u64,
    /// Emissions beyond the first for a record, only possible in `overlap`
    /// mode.
    pub duplicate_emissions: u64,
    /// The largest number of outputs one record went to.
    pub max_emissions_for_one_record: u64,
    /// Temporary bytes written by the spool engine.
    pub temporary_bytes: u64,
    /// Records whose assignment was ambiguous and resolved by tie-breaking.
    pub ambiguous_assignments: u64,
}

impl RunStats {
    /// Folds in one record's emission count.
    pub fn observe_emissions(&mut self, emissions: u64) {
        self.input_records += 1;
        if emissions == 0 {
            self.unmatched_records += 1;
            return;
        }
        self.unique_emitted_records += 1;
        self.total_output_emissions += emissions;
        self.duplicate_emissions += emissions - 1;
        self.max_emissions_for_one_record = self.max_emissions_for_one_record.max(emissions);
    }

    /// Folds in one deliberately dropped record.
    pub fn observe_drop(&mut self) {
        self.input_records += 1;
        self.dropped_records += 1;
    }

    /// Merges counters from another shard of the same run.
    pub fn merge(&mut self, other: &Self) {
        self.input_records += other.input_records;
        self.unique_emitted_records += other.unique_emitted_records;
        self.total_output_emissions += other.total_output_emissions;
        self.dropped_records += other.dropped_records;
        self.unmatched_records += other.unmatched_records;
        self.duplicate_emissions += other.duplicate_emissions;
        self.max_emissions_for_one_record = self
            .max_emissions_for_one_record
            .max(other.max_emissions_for_one_record);
        self.temporary_bytes += other.temporary_bytes;
        self.ambiguous_assignments += other.ambiguous_assignments;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(reference_id: i32, position: i32, flags: u16) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&reference_id.to_le_bytes());
        body.extend_from_slice(&position.to_le_bytes());
        body.push(2);
        body.push(60);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&flags.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(b"r\0");
        body
    }

    #[test]
    fn the_digest_is_order_sensitive_and_length_aware() {
        let mut forward = RecordDigest::new();
        forward.update(b"ab");
        forward.update(b"c");

        let mut reversed = RecordDigest::new();
        reversed.update(b"c");
        reversed.update(b"ab");

        let mut split_differently = RecordDigest::new();
        split_differently.update(b"a");
        split_differently.update(b"bc");

        assert_ne!(forward.finish(), reversed.finish());
        assert_ne!(forward.finish(), split_differently.finish());
        assert_eq!(forward.record_count(), 2);
        assert_eq!(forward.to_hex().len(), 16);
    }

    #[test]
    fn the_digest_is_stable_across_runs() {
        let build = || {
            let mut digest = RecordDigest::new();
            for index in 0..100u32 {
                digest.update(&index.to_le_bytes());
            }
            digest.finish()
        };
        assert_eq!(build(), build());
    }

    #[test]
    fn an_empty_digest_is_well_defined() {
        let digest = RecordDigest::new();
        assert_eq!(digest.record_count(), 0);
        assert_eq!(digest.to_hex().len(), 16);
    }

    #[test]
    fn output_counters_classify_every_record() {
        let mut builder = OutputStatsBuilder::new();
        for body in [
            record(0, 10, 0),
            record(0, 20, 0x100),
            record(0, 30, 0x800),
            record(0, 40, 0x4),
            record(-1, -1, 0x4),
            record(0, 50, 0x400 | 0x200),
        ] {
            builder
                .observe(&RawRecord::new(&body).expect("valid"))
                .expect("observed");
        }
        let stats = builder.finish(1234);
        assert_eq!(stats.record_count, 6);
        assert_eq!(stats.mapped_count, 4);
        assert_eq!(stats.placed_unmapped_count, 1);
        assert_eq!(stats.unplaced_unmapped_count, 1);
        assert_eq!(stats.primary_count, 4);
        assert_eq!(stats.secondary_count, 1);
        assert_eq!(stats.supplementary_count, 1);
        assert_eq!(stats.duplicate_count, 1);
        assert_eq!(stats.qc_fail_count, 1);
        assert_eq!(stats.compressed_bytes, 1234);
        assert!(stats.coordinate_sorted);
        assert_eq!(
            stats.first_coordinate,
            Some(Coordinate {
                reference_id: 0,
                position: 10
            })
        );
        assert_eq!(
            stats.last_coordinate,
            Some(Coordinate {
                reference_id: 0,
                position: 50
            })
        );
        assert!(stats.is_consistent());
        assert_eq!(stats.raw_record_digest.len(), 16);
    }

    #[test]
    fn unsorted_output_is_detected() {
        let mut builder = OutputStatsBuilder::new();
        for body in [record(0, 100, 0), record(0, 50, 0)] {
            builder
                .observe(&RawRecord::new(&body).expect("valid"))
                .expect("observed");
        }
        assert!(!builder.finish(0).coordinate_sorted);
    }

    #[test]
    fn unplaced_records_do_not_disturb_sortedness() {
        let mut builder = OutputStatsBuilder::new();
        for body in [record(0, 100, 0), record(-1, -1, 0x4), record(0, 200, 0)] {
            builder
                .observe(&RawRecord::new(&body).expect("valid"))
                .expect("observed");
        }
        let stats = builder.finish(0);
        assert!(stats.coordinate_sorted);
        assert_eq!(stats.unplaced_unmapped_count, 1);
    }

    #[test]
    fn an_empty_output_is_trivially_sorted_and_consistent() {
        let stats = OutputStatsBuilder::new().finish(28);
        assert_eq!(stats.record_count, 0);
        assert!(stats.coordinate_sorted);
        assert!(stats.is_consistent());
        assert!(stats.first_coordinate.is_none());
    }

    #[test]
    fn run_counters_track_emissions_and_duplicates() {
        let mut stats = RunStats::default();
        stats.observe_emissions(1);
        stats.observe_emissions(3);
        stats.observe_emissions(0);
        stats.observe_drop();

        assert_eq!(stats.input_records, 4);
        assert_eq!(stats.unique_emitted_records, 2);
        assert_eq!(stats.total_output_emissions, 4);
        assert_eq!(stats.duplicate_emissions, 2);
        assert_eq!(stats.unmatched_records, 1);
        assert_eq!(stats.dropped_records, 1);
        assert_eq!(stats.max_emissions_for_one_record, 3);
    }

    #[test]
    fn merging_run_counters_sums_and_maximizes() {
        let mut left = RunStats {
            input_records: 2,
            max_emissions_for_one_record: 5,
            ..RunStats::default()
        };
        let right = RunStats {
            input_records: 3,
            max_emissions_for_one_record: 2,
            temporary_bytes: 99,
            ..RunStats::default()
        };
        left.merge(&right);
        assert_eq!(left.input_records, 5);
        assert_eq!(left.max_emissions_for_one_record, 5);
        assert_eq!(left.temporary_bytes, 99);
    }
}
