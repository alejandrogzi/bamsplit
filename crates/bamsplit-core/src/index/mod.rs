// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Streaming construction of BAI and CSI indexes for the outputs `bamsplit`
//! writes.
//!
//! # Why streaming
//!
//! An index is a map from genomic interval to BGZF virtual-offset ranges. Both
//! ends of that map are available while writing:
//!
//! ```text
//! start = writer.virtual_position()
//! write(record)
//! end   = writer.virtual_position()
//! index.add(reference, start_pos, end_pos, mapped, start..end)
//! ```
//!
//! so no second pass over the finished BAM is needed. That matters: a second
//! pass would re-inflate and re-parse every byte just written, roughly doubling
//! the cost of `--index auto`.
//!
//! The one constraint this imposes is that the BGZF writer must be able to
//! report its position. `noodles`' single-threaded `bgzf::io::Writer` can;
//! `MultithreadedWriter` cannot. [`crate::output`] therefore selects the
//! writer based on whether an index was requested, and the trade-off is
//! documented in `docs/performance.md`.
//!
//! # Choosing between BAI and CSI
//!
//! Decided **before writing**, from the reference dictionary:
//!
//! | `--index` | longest reference ≤ 2^29-1 | longest reference > 2^29-1 |
//! | --- | --- | --- |
//! | `auto` | BAI | CSI at the smallest sufficient depth |
//! | `bai` | BAI | error — BAI cannot address those coordinates |
//! | `csi` | CSI | CSI |
//! | `none` | no index | no index |
//!
//! Deciding up front keeps the choice deterministic and means a run never
//! discovers halfway through that it picked the wrong format.
//!
//! # When no index is written
//!
//! A binning index describes a coordinate-sorted file. If the records routed to
//! an output are not in coordinate order — an unsorted input, or a
//! query-name-sorted one — an index would be actively misleading. The builder
//! watches the reference ids and positions it is fed and, on the first
//! violation, abandons the index and records why. `auto` reports the reason and
//! writes no index; an explicit `--index bai`/`--index csi` turns it into an
//! error, because the user asked for something that cannot be produced.

pub mod bai;
pub mod csi;

use std::path::{Path, PathBuf};

use noodles_bgzf::VirtualPosition;
use noodles_core::Position;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::error::IndexError;

/// What the user asked for on the command line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum IndexMode {
    /// Pick BAI or CSI from the reference dictionary; skip when not indexable.
    #[default]
    Auto,
    /// Force BAI; fail if the coordinates do not fit.
    Bai,
    /// Force CSI.
    Csi,
    /// Write no index.
    None,
}

impl IndexMode {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Bai => "bai",
            Self::Csi => "csi",
            Self::None => "none",
        }
    }

    /// Whether the user pinned a specific format.
    #[must_use]
    pub const fn is_explicit(self) -> bool {
        matches!(self, Self::Bai | Self::Csi)
    }
}

impl std::fmt::Display for IndexMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The format actually selected, with its binning parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    /// A BAI index: `min_shift = 14`, `depth = 5`.
    Bai,
    /// A CSI index at the given depth.
    Csi {
        /// The minimum shift.
        min_shift: u8,
        /// The R-tree depth.
        depth: u8,
    },
}

impl IndexKind {
    /// The filename extension, without the dot.
    #[must_use]
    pub const fn extension(self) -> &'static str {
        match self {
            Self::Bai => bai::EXTENSION,
            Self::Csi { .. } => csi::EXTENSION,
        }
    }

    /// The largest 1-based position this index can address.
    #[must_use]
    pub fn max_position(self) -> u64 {
        match self {
            Self::Bai => bai::MAX_POSITION,
            Self::Csi { min_shift, depth } => csi::addressable(min_shift, depth),
        }
    }
}

impl std::fmt::Display for IndexKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bai => f.write_str("bai"),
            Self::Csi { depth, .. } => write!(f, "csi(depth={depth})"),
        }
    }
}

/// The outcome of index planning, including the reason, for the manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexPlan {
    /// The chosen format, or [`None`] when no index will be written.
    pub kind: Option<IndexKind>,
    /// A one-line explanation, always populated.
    pub reason: String,
}

/// Chooses an index format from the mode and the reference dictionary.
///
/// # Errors
///
/// Returns [`IndexError::BaiLimitExceeded`] when `--index bai` is forced but a
/// reference is longer than BAI can address.
pub fn plan(
    mode: IndexMode,
    max_reference_length: u64,
    longest_reference_name: &str,
    output_hint: &Path,
) -> Result<IndexPlan, IndexError> {
    Ok(match mode {
        IndexMode::None => IndexPlan {
            kind: None,
            reason: "`--index none` was requested".to_string(),
        },
        IndexMode::Bai => {
            if !bai::can_represent(max_reference_length) {
                return Err(IndexError::BaiLimitExceeded {
                    path: output_hint.to_path_buf(),
                    reference: longest_reference_name.to_string(),
                    position: max_reference_length,
                    limit: bai::MAX_POSITION,
                });
            }
            IndexPlan {
                kind: Some(IndexKind::Bai),
                reason: "`--index bai` was requested and every coordinate fits".to_string(),
            }
        }
        IndexMode::Csi => {
            let depth = csi::depth_for(max_reference_length);
            IndexPlan {
                kind: Some(IndexKind::Csi {
                    min_shift: csi::MIN_SHIFT,
                    depth,
                }),
                reason: format!(
                    "`--index csi` was requested; depth {depth} addresses {} bases",
                    csi::addressable(csi::MIN_SHIFT, depth)
                ),
            }
        }
        IndexMode::Auto => {
            if bai::can_represent(max_reference_length) {
                IndexPlan {
                    kind: Some(IndexKind::Bai),
                    reason: format!(
                        "the longest reference is {max_reference_length} bases, within BAI's \
                         {} limit",
                        bai::MAX_POSITION
                    ),
                }
            } else {
                let depth = csi::depth_for(max_reference_length);
                IndexPlan {
                    kind: Some(IndexKind::Csi {
                        min_shift: csi::MIN_SHIFT,
                        depth,
                    }),
                    reason: format!(
                        "reference {longest_reference_name:?} is {max_reference_length} bases, \
                         beyond BAI's {} limit, so CSI at depth {depth} was selected",
                        bai::MAX_POSITION
                    ),
                }
            }
        }
    })
}

/// The alignment context one record contributes to an index.
///
/// `None` means the record is unplaced: it has no reference and only bumps the
/// index's unplaced-unmapped counter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentContext {
    /// The reference the record is placed on.
    pub reference_id: usize,
    /// 1-based inclusive start.
    pub start: u64,
    /// 1-based inclusive end.
    pub end: u64,
    /// Whether the record is mapped (the `UNMAPPED` flag is clear).
    pub mapped: bool,
}

/// Why an index was abandoned partway through.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AbandonReason {
    /// A record referenced an earlier reference than one already indexed.
    OutOfOrderReference {
        /// The reference id that went backwards.
        reference_id: usize,
        /// The highest reference id seen before it.
        previous_reference_id: usize,
    },
    /// A record started before the previous record on the same reference.
    OutOfOrderPosition {
        /// The reference id.
        reference_id: usize,
        /// The offending 1-based start.
        start: u64,
        /// The previous 1-based start.
        previous_start: u64,
    },
    /// A coordinate exceeded what the chosen format can address.
    CoordinateOutOfRange {
        /// The offending 1-based position.
        position: u64,
        /// The format's ceiling.
        limit: u64,
    },
    /// `noodles` rejected the record.
    Rejected {
        /// The rejection message.
        message: String,
    },
}

impl std::fmt::Display for AbandonReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::OutOfOrderReference {
                reference_id,
                previous_reference_id,
            } => write!(
                f,
                "reference {reference_id} follows reference {previous_reference_id}, so the \
                 output is not coordinate-sorted"
            ),
            Self::OutOfOrderPosition {
                reference_id,
                start,
                previous_start,
            } => write!(
                f,
                "position {start} follows {previous_start} on reference {reference_id}, so the \
                 output is not coordinate-sorted"
            ),
            Self::CoordinateOutOfRange { position, limit } => write!(
                f,
                "position {position} exceeds the index's addressable maximum of {limit}"
            ),
            Self::Rejected { message } => write!(f, "the index rejected a record: {message}"),
        }
    }
}

/// Accumulates an index while records are written.
///
/// The builder is fed one call per record, in write order. It is deliberately
/// tolerant: rather than failing the whole run when an output turns out not to
/// be coordinate-sorted, it abandons *that output's* index and remembers why,
/// so the caller can decide whether that is fatal.
#[derive(Debug)]
pub struct IndexBuilder {
    kind: IndexKind,
    reference_count: usize,
    inner: Inner,
    previous_reference_id: Option<usize>,
    previous_start: u64,
    abandoned: Option<AbandonReason>,
    records: u64,
}

#[derive(Debug)]
enum Inner {
    Bai(Box<bai::Indexer>),
    Csi(Box<csi::Indexer>),
    Abandoned,
}

impl IndexBuilder {
    /// Creates a builder for `kind` over a dictionary of `reference_count`
    /// references.
    #[must_use]
    pub fn new(kind: IndexKind, reference_count: usize) -> Self {
        let inner = match kind {
            IndexKind::Bai => Inner::Bai(Box::new(bai::Indexer::new(bai::MIN_SHIFT, bai::DEPTH))),
            IndexKind::Csi { min_shift, depth } => {
                Inner::Csi(Box::new(csi::Indexer::new(min_shift, depth)))
            }
        };
        Self {
            kind,
            reference_count,
            inner,
            previous_reference_id: None,
            previous_start: 0,
            abandoned: None,
            records: 0,
        }
    }

    /// The format being built.
    #[must_use]
    pub const fn kind(&self) -> IndexKind {
        self.kind
    }

    /// Whether the index has been abandoned.
    #[must_use]
    pub const fn is_abandoned(&self) -> bool {
        self.abandoned.is_some()
    }

    /// Why the index was abandoned, if it was.
    #[must_use]
    pub const fn abandon_reason(&self) -> Option<&AbandonReason> {
        self.abandoned.as_ref()
    }

    /// How many records have been offered, including after abandonment.
    #[must_use]
    pub const fn record_count(&self) -> u64 {
        self.records
    }

    /// Adds one record.
    ///
    /// `chunk` is the record's `[start, end)` virtual-offset range in the
    /// output stream. Once the index has been abandoned this is a cheap no-op,
    /// so the write loop does not need to branch.
    ///
    /// # Errors
    ///
    /// Never fails: an unusable record abandons the index instead. The
    /// [`Result`] is kept so the signature does not change if a future format
    /// needs to report a hard failure.
    pub fn add(
        &mut self,
        context: Option<AlignmentContext>,
        chunk_start: u64,
        chunk_end: u64,
    ) -> Result<(), IndexError> {
        self.records += 1;
        if matches!(self.inner, Inner::Abandoned) {
            return Ok(());
        }

        let chunk = Chunk::new(
            VirtualPosition::from(chunk_start),
            VirtualPosition::from(chunk_end),
        );

        let Some(context) = context else {
            self.push(None, chunk);
            return Ok(());
        };

        let limit = self.kind.max_position();
        if context.end > limit || context.start > limit {
            self.abandon(AbandonReason::CoordinateOutOfRange {
                position: context.end.max(context.start),
                limit,
            });
            return Ok(());
        }

        match self.previous_reference_id {
            Some(previous) if context.reference_id < previous => {
                self.abandon(AbandonReason::OutOfOrderReference {
                    reference_id: context.reference_id,
                    previous_reference_id: previous,
                });
                return Ok(());
            }
            Some(previous)
                if context.reference_id == previous && context.start < self.previous_start =>
            {
                self.abandon(AbandonReason::OutOfOrderPosition {
                    reference_id: context.reference_id,
                    start: context.start,
                    previous_start: self.previous_start,
                });
                return Ok(());
            }
            _ => {}
        }

        // `Position` is 1-based and non-zero; `AlignmentContext` guarantees
        // both, but a hostile record could still present a zero, so convert
        // fallibly rather than unwrapping.
        let (Some(start), Some(end)) = (
            usize::try_from(context.start).ok().and_then(Position::new),
            usize::try_from(context.end).ok().and_then(Position::new),
        ) else {
            self.abandon(AbandonReason::Rejected {
                message: format!(
                    "the interval {}..={} is not a valid 1-based range",
                    context.start, context.end
                ),
            });
            return Ok(());
        };

        self.previous_reference_id = Some(context.reference_id);
        self.previous_start = context.start;
        self.push(
            Some((context.reference_id, start, end, context.mapped)),
            chunk,
        );
        Ok(())
    }

    fn push(&mut self, context: Option<(usize, Position, Position, bool)>, chunk: Chunk) {
        let result = match &mut self.inner {
            Inner::Bai(indexer) => indexer.add_record(context, chunk),
            Inner::Csi(indexer) => indexer.add_record(context, chunk),
            Inner::Abandoned => return,
        };
        if let Err(error) = result {
            self.abandon(AbandonReason::Rejected {
                message: error.to_string(),
            });
        }
    }

    fn abandon(&mut self, reason: AbandonReason) {
        if self.abandoned.is_none() {
            self.abandoned = Some(reason);
        }
        // Drop the partially built index so its memory is released as soon as
        // it is known to be useless. Outputs with millions of bins make this
        // worth doing eagerly.
        self.inner = Inner::Abandoned;
    }

    /// Finishes the index and writes it to `path`.
    ///
    /// Returns [`None`] without touching the filesystem when the index was
    /// abandoned.
    ///
    /// # Errors
    ///
    /// Returns [`IndexError::Io`] if the index cannot be written.
    pub fn finish(self, path: &Path) -> Result<Option<FinishedIndex>, IndexError> {
        match self.inner {
            Inner::Abandoned => Ok(None),
            Inner::Bai(indexer) => {
                let index = indexer.build(self.reference_count);
                bai::write(path, &index)?;
                Ok(Some(FinishedIndex {
                    kind: self.kind,
                    path: path.to_path_buf(),
                }))
            }
            Inner::Csi(indexer) => {
                let index = indexer.build(self.reference_count);
                csi::write(path, &index)?;
                Ok(Some(FinishedIndex {
                    kind: self.kind,
                    path: path.to_path_buf(),
                }))
            }
        }
    }
}

/// A written index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinishedIndex {
    /// The format.
    pub kind: IndexKind,
    /// Where it was written.
    pub path: PathBuf,
}

#[cfg(test)]
#[allow(clippy::items_after_statements, clippy::unnecessary_wraps)]
mod tests {
    use super::*;

    fn context(reference_id: usize, start: u64, end: u64) -> Option<AlignmentContext> {
        Some(AlignmentContext {
            reference_id,
            start,
            end,
            mapped: true,
        })
    }

    #[test]
    fn auto_prefers_bai_for_ordinary_genomes() {
        let plan =
            plan(IndexMode::Auto, 248_956_422, "chr1", Path::new("out.bam")).expect("planned");
        assert_eq!(plan.kind, Some(IndexKind::Bai));
        assert!(plan.reason.contains("within BAI"), "{}", plan.reason);
    }

    #[test]
    fn auto_falls_back_to_csi_beyond_the_bai_limit() {
        let plan = plan(
            IndexMode::Auto,
            1_000_000_000,
            "chrHuge",
            Path::new("out.bam"),
        )
        .expect("planned");
        assert_eq!(
            plan.kind,
            Some(IndexKind::Csi {
                min_shift: 14,
                depth: 6
            })
        );
        assert!(plan.reason.contains("beyond BAI"), "{}", plan.reason);
    }

    #[test]
    fn forced_bai_is_rejected_beyond_its_limit() {
        let error = plan(
            IndexMode::Bai,
            bai::MAX_POSITION + 1,
            "chrHuge",
            Path::new("out.bam"),
        )
        .expect_err("must reject");
        assert!(
            matches!(error, IndexError::BaiLimitExceeded { .. }),
            "{error}"
        );
    }

    #[test]
    fn forced_bai_is_accepted_exactly_at_its_limit() {
        let plan = plan(
            IndexMode::Bai,
            bai::MAX_POSITION,
            "chrBig",
            Path::new("out.bam"),
        )
        .expect("planned");
        assert_eq!(plan.kind, Some(IndexKind::Bai));
    }

    #[test]
    fn none_writes_nothing_and_says_so() {
        let plan = plan(IndexMode::None, 1, "chr1", Path::new("out.bam")).expect("planned");
        assert!(plan.kind.is_none());
        assert!(plan.reason.contains("none"), "{}", plan.reason);
    }

    #[test]
    fn a_sorted_stream_builds_and_writes() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("out.bam.bai");

        use noodles_csi::BinningIndex as _;

        let mut builder = IndexBuilder::new(IndexKind::Bai, 2);
        builder.add(context(0, 1, 100), 0, 200).expect("added");
        builder.add(context(0, 500, 600), 200, 400).expect("added");
        builder.add(context(1, 1, 100), 400, 600).expect("added");
        builder.add(None, 600, 800).expect("added");
        assert!(!builder.is_abandoned());
        assert_eq!(builder.record_count(), 4);

        let finished = builder.finish(&path).expect("written").expect("some");
        assert_eq!(finished.kind, IndexKind::Bai);
        assert!(path.exists());

        let read_back = bai::read(&path).expect("readable");
        assert_eq!(read_back.min_shift(), 14);
        assert_eq!(read_back.depth(), 5);
        assert_eq!(read_back.unplaced_unmapped_record_count(), Some(1));
    }

    #[test]
    fn a_backwards_reference_abandons_the_index() {
        let mut builder = IndexBuilder::new(IndexKind::Bai, 2);
        builder.add(context(1, 1, 100), 0, 100).expect("added");
        builder.add(context(0, 1, 100), 100, 200).expect("added");
        assert!(builder.is_abandoned());
        assert!(matches!(
            builder.abandon_reason(),
            Some(AbandonReason::OutOfOrderReference { .. })
        ));

        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("out.bam.bai");
        assert!(builder.finish(&path).expect("no error").is_none());
        assert!(!path.exists(), "an abandoned index must not be written");
    }

    #[test]
    fn a_backwards_position_abandons_the_index() {
        let mut builder = IndexBuilder::new(IndexKind::Bai, 1);
        builder.add(context(0, 500, 600), 0, 100).expect("added");
        builder.add(context(0, 100, 200), 100, 200).expect("added");
        assert!(matches!(
            builder.abandon_reason(),
            Some(AbandonReason::OutOfOrderPosition {
                start: 100,
                previous_start: 500,
                ..
            })
        ));
    }

    #[test]
    fn equal_positions_are_in_order() {
        let mut builder = IndexBuilder::new(IndexKind::Bai, 1);
        builder.add(context(0, 100, 200), 0, 100).expect("added");
        builder.add(context(0, 100, 200), 100, 200).expect("added");
        assert!(!builder.is_abandoned());
    }

    #[test]
    fn a_coordinate_beyond_the_format_limit_abandons_the_index() {
        let mut builder = IndexBuilder::new(IndexKind::Bai, 1);
        builder
            .add(context(0, 1, bai::MAX_POSITION + 1), 0, 100)
            .expect("added");
        assert!(matches!(
            builder.abandon_reason(),
            Some(AbandonReason::CoordinateOutOfRange { .. })
        ));
    }

    #[test]
    fn csi_accepts_coordinates_bai_cannot() {
        use noodles_csi::BinningIndex as _;

        let mut builder = IndexBuilder::new(
            IndexKind::Csi {
                min_shift: 14,
                depth: 6,
            },
            1,
        );
        builder
            .add(context(0, 1, bai::MAX_POSITION + 1_000), 0, 100)
            .expect("added");
        assert!(!builder.is_abandoned());

        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("out.bam.csi");
        let finished = builder.finish(&path).expect("written").expect("some");
        assert_eq!(finished.kind.extension(), "csi");

        let read_back = csi::read(&path).expect("readable");
        assert_eq!(read_back.depth(), 6);
    }

    #[test]
    fn a_zero_start_is_rejected_rather_than_panicking() {
        let mut builder = IndexBuilder::new(IndexKind::Bai, 1);
        builder
            .add(
                Some(AlignmentContext {
                    reference_id: 0,
                    start: 0,
                    end: 0,
                    mapped: true,
                }),
                0,
                100,
            )
            .expect("added");
        assert!(matches!(
            builder.abandon_reason(),
            Some(AbandonReason::Rejected { .. })
        ));
    }

    #[test]
    fn an_output_with_no_records_still_writes_a_valid_index() {
        use noodles_csi::BinningIndex;

        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("empty.bam.bai");
        let builder = IndexBuilder::new(IndexKind::Bai, 3);
        assert!(builder.finish(&path).expect("written").is_some());

        let index = bai::read(&path).expect("readable");
        assert_eq!(BinningIndex::reference_sequences(&index).count(), 3);
    }

    #[test]
    fn adding_after_abandonment_is_a_no_op() {
        let mut builder = IndexBuilder::new(IndexKind::Bai, 1);
        builder.add(context(0, 500, 600), 0, 100).expect("added");
        builder.add(context(0, 1, 2), 100, 200).expect("added");
        let reason = builder.abandon_reason().cloned();
        builder.add(context(0, 1, 2), 200, 300).expect("added");
        assert_eq!(builder.abandon_reason().cloned(), reason);
        assert_eq!(builder.record_count(), 3);
    }

    #[test]
    fn index_modes_round_trip_their_option_values() {
        for (mode, text) in [
            (IndexMode::Auto, "auto"),
            (IndexMode::Bai, "bai"),
            (IndexMode::Csi, "csi"),
            (IndexMode::None, "none"),
        ] {
            assert_eq!(mode.to_string(), text);
        }
        assert!(IndexMode::Bai.is_explicit());
        assert!(!IndexMode::Auto.is_explicit());
    }
}
