// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! CIGAR decoding and the reusable alignment-geometry view.
//!
//! A BAM CIGAR is a packed `u32` array: `(length << 4) | op_code`. This module
//! turns it into the two geometric views region routing needs.
//!
//! # Reference span versus aligned blocks
//!
//! ```text
//! CIGAR   10M  5I  20M  1000N  15M  4D  10M  20S
//!         └────────────────────────────────────┘  read
//!         └──┘    └──┘ └────┘ └──┘└─┘└──┘         reference span (55 + 1000 + 4)
//!         ▓▓▓▓    ▓▓▓▓        ▓▓▓▓    ▓▓▓▓         aligned blocks (M/=/X only)
//! ```
//!
//! * **Reference span** — every operation that consumes reference bases:
//!   `M`, `D`, `N`, `=`, `X`. This is the outer footprint of the alignment.
//! * **Aligned blocks** — maximal runs of `M`, `=`, and `X`. These are the
//!   positions where a base of the read is actually placed on the reference.
//!
//! `N` separates blocks: it is the spliced-out intron of an RNA-seq alignment
//! and must never count as exonic overlap.
//!
//! `D` is the interesting case, and `bamsplit` makes a deliberate choice:
//!
//! > **A deletion contributes to the reference span but not to aligned-base
//! > overlap.**
//!
//! So `20M 4D 10M` yields two blocks, `[p, p+20)` and `[p+24, p+34)`, and a
//! reference span of 34. This keeps `best-overlap` scoring honest — a read is
//! credited only for bases it actually aligns — and it makes `contained` reject
//! a read whose deletion falls outside the selected feature only when the
//! flanking aligned bases do too. The alternative (folding `D` into the
//! surrounding block) is available as
//! [`DeletionPolicy::IncludeInBlocks`] for callers that prefer htslib-style
//! block merging; nothing in the CLI selects it today.
//!
//! `I`, `S`, `H`, and `P` consume no reference bases and never appear in either
//! view.
//!
//! # Example
//!
//! ```
//! use bamsplit_core::bam::cigar::{CigarOps, DeletionPolicy, Kind};
//! use bamsplit_core::interval::Interval;
//!
//! // 10M 5N 10M, packed little-endian.
//! let mut packed = Vec::new();
//! for (length, kind) in [(10u32, 0u32), (5, 3), (10, 0)] {
//!     packed.extend_from_slice(&((length << 4) | kind).to_le_bytes());
//! }
//!
//! let ops = CigarOps::new(&packed);
//! let span = ops.reference_span()?;
//! assert_eq!(span.total(), 25);
//! assert_eq!(span.aligned, 20);
//! assert_eq!(span.skipped, 5);
//!
//! let blocks = ops.aligned_blocks(100, DeletionPolicy::ExcludeFromBlocks)?;
//! assert_eq!(blocks, vec![Interval::new(100, 110), Interval::new(115, 125)]);
//! # Ok::<_, bamsplit_core::error::CigarError>(())
//! ```

use crate::error::CigarError;
use crate::interval::Interval;

/// A CIGAR operation kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Kind {
    /// `M`: alignment match or mismatch.
    Match,
    /// `I`: insertion into the reference.
    Insertion,
    /// `D`: deletion from the reference.
    Deletion,
    /// `N`: skipped region of the reference (an intron, for RNA-seq).
    Skip,
    /// `S`: soft clip.
    SoftClip,
    /// `H`: hard clip.
    HardClip,
    /// `P`: padding.
    Pad,
    /// `=`: sequence match.
    SequenceMatch,
    /// `X`: sequence mismatch.
    SequenceMismatch,
}

impl Kind {
    /// Decodes the low four bits of a packed operation.
    ///
    /// # Errors
    ///
    /// Returns [`CigarError::UnknownOperation`] for codes outside `0..=8`.
    pub const fn from_code(code: u32) -> Result<Self, CigarError> {
        Ok(match code {
            0 => Self::Match,
            1 => Self::Insertion,
            2 => Self::Deletion,
            3 => Self::Skip,
            4 => Self::SoftClip,
            5 => Self::HardClip,
            6 => Self::Pad,
            7 => Self::SequenceMatch,
            8 => Self::SequenceMismatch,
            other => return Err(CigarError::UnknownOperation { code: other }),
        })
    }

    /// The SAM character for this kind.
    #[must_use]
    pub const fn as_char(self) -> char {
        match self {
            Self::Match => 'M',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
            Self::Skip => 'N',
            Self::SoftClip => 'S',
            Self::HardClip => 'H',
            Self::Pad => 'P',
            Self::SequenceMatch => '=',
            Self::SequenceMismatch => 'X',
        }
    }

    /// Whether the operation advances along the reference.
    #[must_use]
    pub const fn consumes_reference(self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Deletion
                | Self::Skip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }

    /// Whether the operation consumes bases of the read.
    #[must_use]
    pub const fn consumes_read(self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Insertion
                | Self::SoftClip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }

    /// Whether the operation places read bases on reference bases.
    ///
    /// This is exactly `M`, `=`, and `X` — the operations that make up an
    /// aligned block.
    #[must_use]
    pub const fn is_aligned_match(self) -> bool {
        matches!(
            self,
            Self::Match | Self::SequenceMatch | Self::SequenceMismatch
        )
    }
}

/// One decoded CIGAR operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Op {
    kind: Kind,
    length: u32,
}

impl Op {
    /// Creates an operation.
    #[must_use]
    pub const fn new(kind: Kind, length: u32) -> Self {
        Self { kind, length }
    }

    /// Decodes a packed `u32`.
    ///
    /// # Errors
    ///
    /// Returns [`CigarError::UnknownOperation`] for an operation code outside
    /// `0..=8`.
    pub const fn decode(packed: u32) -> Result<Self, CigarError> {
        match Kind::from_code(packed & 0xf) {
            Ok(kind) => Ok(Self {
                kind,
                length: packed >> 4,
            }),
            Err(error) => Err(error),
        }
    }

    /// The operation kind.
    #[must_use]
    pub const fn kind(&self) -> Kind {
        self.kind
    }

    /// The operation length.
    #[must_use]
    pub const fn length(&self) -> u32 {
        self.length
    }
}

impl std::fmt::Display for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.length, self.kind.as_char())
    }
}

/// How deletions are treated when building aligned blocks.
///
/// See the module documentation for why `bamsplit` defaults to excluding them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DeletionPolicy {
    /// A `D` run breaks the block; its reference bases belong to no block.
    ///
    /// This is the `bamsplit` default and the behaviour used by
    /// `--assignment best-overlap`.
    #[default]
    ExcludeFromBlocks,
    /// A `D` run is folded into the surrounding block, htslib-style.
    IncludeInBlocks,
}

/// A breakdown of the reference bases an alignment covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ReferenceSpan {
    /// Bases covered by `M`, `=`, and `X`.
    pub aligned: u32,
    /// Bases covered by `D`.
    pub deleted: u32,
    /// Bases covered by `N`.
    pub skipped: u32,
}

impl ReferenceSpan {
    /// The full reference footprint: aligned + deleted + skipped.
    #[must_use]
    pub const fn total(&self) -> u32 {
        // Each component is bounded by the sum of CIGAR op lengths, which
        // `CigarOps::reference_span` has already checked against `u32::MAX`.
        self.aligned
            .saturating_add(self.deleted)
            .saturating_add(self.skipped)
    }
}

/// An iterator over the packed operations of a CIGAR.
///
/// The iterator borrows the packed bytes and decodes lazily, so a router that
/// only needs the reference span never materializes a `Vec<Op>`.
#[derive(Debug, Clone, Copy)]
pub struct CigarOps<'a> {
    packed: &'a [u8],
    offset: usize,
}

impl<'a> CigarOps<'a> {
    /// Wraps a packed CIGAR.
    ///
    /// A trailing partial operation (fewer than four bytes) is ignored by
    /// iteration; [`RawRecord::new`](crate::bam::raw_record::RawRecord::new)
    /// rejects such records before they reach this type.
    #[must_use]
    pub const fn new(packed: &'a [u8]) -> Self {
        Self { packed, offset: 0 }
    }

    /// The number of complete operations.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.packed.len() / 4
    }

    /// Whether the CIGAR is absent.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.packed.len() < 4
    }

    /// The packed bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &'a [u8] {
        self.packed
    }

    /// Computes the reference footprint in one pass.
    ///
    /// # Errors
    ///
    /// Returns [`CigarError::UnknownOperation`] for an unknown operation code
    /// and [`CigarError::SpanOverflow`] if the total exceeds `u32::MAX`.
    pub fn reference_span(&self) -> Result<ReferenceSpan, CigarError> {
        let mut span = ReferenceSpan::default();
        for op in *self {
            let op = op?;
            let slot = match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => &mut span.aligned,
                Kind::Deletion => &mut span.deleted,
                Kind::Skip => &mut span.skipped,
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => continue,
            };
            *slot = slot
                .checked_add(op.length())
                .ok_or(CigarError::SpanOverflow)?;
        }
        span.aligned
            .checked_add(span.deleted)
            .and_then(|partial| partial.checked_add(span.skipped))
            .ok_or(CigarError::SpanOverflow)?;
        Ok(span)
    }

    /// The number of read bases the CIGAR accounts for, excluding hard clips.
    ///
    /// Used by validation to cross-check `l_seq`.
    ///
    /// # Errors
    ///
    /// Returns [`CigarError`] for an unknown operation or an overflowing sum.
    pub fn read_length(&self) -> Result<u32, CigarError> {
        let mut total = 0u32;
        for op in *self {
            let op = op?;
            if op.kind().consumes_read() {
                total = total
                    .checked_add(op.length())
                    .ok_or(CigarError::SpanOverflow)?;
            }
        }
        Ok(total)
    }

    /// Builds the aligned blocks of an alignment starting at `start`.
    ///
    /// `start` is the 0-based leftmost reference position. The returned
    /// intervals are 0-based half-open, sorted, and non-overlapping.
    ///
    /// # Errors
    ///
    /// Returns [`CigarError`] for an unknown operation or an overflowing span.
    pub fn aligned_blocks(
        &self,
        start: i64,
        policy: DeletionPolicy,
    ) -> Result<Vec<Interval>, CigarError> {
        let mut blocks = Vec::new();
        self.for_each_aligned_block(start, policy, |block| blocks.push(block))?;
        Ok(blocks)
    }

    /// Streams the aligned blocks into `sink` without allocating a `Vec`.
    ///
    /// The hot path in region routing uses this to reuse one scratch buffer
    /// across every record.
    ///
    /// # Errors
    ///
    /// Returns [`CigarError`] for an unknown operation or an overflowing span.
    pub fn for_each_aligned_block<F>(
        &self,
        start: i64,
        policy: DeletionPolicy,
        mut sink: F,
    ) -> Result<(), CigarError>
    where
        F: FnMut(Interval),
    {
        let mut cursor = start;
        let mut open: Option<Interval> = None;

        for op in *self {
            let op = op?;
            let length = i64::from(op.length());
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let end = cursor.checked_add(length).ok_or(CigarError::SpanOverflow)?;
                    match &mut open {
                        Some(block) => block.end = end,
                        None => open = Some(Interval::new(cursor, end)),
                    }
                    cursor = end;
                }
                Kind::Deletion => {
                    let end = cursor.checked_add(length).ok_or(CigarError::SpanOverflow)?;
                    match policy {
                        DeletionPolicy::IncludeInBlocks => {
                            if let Some(block) = &mut open {
                                block.end = end;
                            }
                            // A leading `D` with no preceding aligned base
                            // opens no block: there is nothing aligned yet.
                        }
                        DeletionPolicy::ExcludeFromBlocks => {
                            if let Some(block) = open.take() {
                                sink(block);
                            }
                        }
                    }
                    cursor = end;
                }
                Kind::Skip => {
                    if let Some(block) = open.take() {
                        sink(block);
                    }
                    cursor = cursor.checked_add(length).ok_or(CigarError::SpanOverflow)?;
                }
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
            }
        }

        if let Some(block) = open.take() {
            sink(block);
        }
        Ok(())
    }

    /// The full reference footprint as one interval, `[start, start + span)`.
    ///
    /// A CIGAR that consumes no reference bases still yields a one-base
    /// interval, matching `bam_endpos` in htslib and
    /// [`RawRecord::alignment_end`](crate::bam::raw_record::RawRecord::alignment_end).
    ///
    /// # Errors
    ///
    /// Returns [`CigarError`] for an unknown operation or an overflowing span.
    pub fn reference_interval(&self, start: i64) -> Result<Interval, CigarError> {
        let span = self.reference_span()?;
        let width = i64::from(span.total().max(1));
        let end = start.checked_add(width).ok_or(CigarError::SpanOverflow)?;
        Ok(Interval::new(start, end))
    }

    /// Renders the CIGAR in SAM text form, for diagnostics.
    ///
    /// # Errors
    ///
    /// Returns [`CigarError::UnknownOperation`] for an unknown operation code.
    pub fn to_sam_string(&self) -> Result<String, CigarError> {
        use std::fmt::Write as _;
        let mut out = String::new();
        for op in *self {
            let op = op?;
            let _ = write!(out, "{op}");
        }
        if out.is_empty() {
            out.push('*');
        }
        Ok(out)
    }
}

// `CigarOps` is `Copy` so a router can hand the same view to several helpers
// without re-slicing the record; iterating a copy is intentional here.
#[allow(clippy::copy_iterator)]
impl Iterator for CigarOps<'_> {
    type Item = Result<Op, CigarError>;

    fn next(&mut self) -> Option<Self::Item> {
        let rest = self.packed.get(self.offset..)?;
        let chunk: &[u8; 4] = rest.get(..4)?.try_into().ok()?;
        self.offset += 4;
        Some(Op::decode(u32::from_le_bytes(*chunk)))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = (self.packed.len().saturating_sub(self.offset)) / 4;
        (remaining, Some(remaining))
    }
}

impl std::iter::ExactSizeIterator for CigarOps<'_> {}

/// Which geometry a region router compares against annotation segments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlignmentGeometry {
    /// Use the aligned blocks (`M`/`=`/`X` runs). The default.
    #[default]
    Blocks,
    /// Use the outer reference span, ignoring introns and deletions.
    Span,
}

impl AlignmentGeometry {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Blocks => "blocks",
            Self::Span => "span",
        }
    }
}

impl std::fmt::Display for AlignmentGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pack(ops: &[(u32, u32)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (length, kind) in ops {
            out.extend_from_slice(&((length << 4) | kind).to_le_bytes());
        }
        out
    }

    const M: u32 = 0;
    const I: u32 = 1;
    const D: u32 = 2;
    const N: u32 = 3;
    const S: u32 = 4;
    const H: u32 = 5;
    const P: u32 = 6;
    const EQ: u32 = 7;
    const X: u32 = 8;

    #[test]
    fn decodes_every_operation_kind() {
        for (code, expected, character) in [
            (M, Kind::Match, 'M'),
            (I, Kind::Insertion, 'I'),
            (D, Kind::Deletion, 'D'),
            (N, Kind::Skip, 'N'),
            (S, Kind::SoftClip, 'S'),
            (H, Kind::HardClip, 'H'),
            (P, Kind::Pad, 'P'),
            (EQ, Kind::SequenceMatch, '='),
            (X, Kind::SequenceMismatch, 'X'),
        ] {
            let op = Op::decode((7 << 4) | code).expect("known code");
            assert_eq!(op.kind(), expected);
            assert_eq!(op.length(), 7);
            assert_eq!(expected.as_char(), character);
        }
    }

    #[test]
    fn rejects_unknown_operation_codes() {
        for code in 9..16u32 {
            let error = Op::decode(code).expect_err("must reject");
            assert!(
                matches!(error, CigarError::UnknownOperation { .. }),
                "{error}"
            );
        }
    }

    #[test]
    fn reference_consumption_matches_the_specification() {
        assert!(Kind::Match.consumes_reference());
        assert!(Kind::Deletion.consumes_reference());
        assert!(Kind::Skip.consumes_reference());
        assert!(Kind::SequenceMatch.consumes_reference());
        assert!(Kind::SequenceMismatch.consumes_reference());
        for kind in [Kind::Insertion, Kind::SoftClip, Kind::HardClip, Kind::Pad] {
            assert!(!kind.consumes_reference(), "{kind:?}");
        }
    }

    #[test]
    fn read_consumption_matches_the_specification() {
        for kind in [
            Kind::Match,
            Kind::Insertion,
            Kind::SoftClip,
            Kind::SequenceMatch,
            Kind::SequenceMismatch,
        ] {
            assert!(kind.consumes_read(), "{kind:?}");
        }
        for kind in [Kind::Deletion, Kind::Skip, Kind::HardClip, Kind::Pad] {
            assert!(!kind.consumes_read(), "{kind:?}");
        }
    }

    #[test]
    fn splits_a_spliced_alignment_into_exonic_blocks() {
        let packed = pack(&[(10, S), (20, M), (1000, N), (30, M), (5, S)]);
        let ops = CigarOps::new(&packed);

        let span = ops.reference_span().expect("valid");
        assert_eq!(span.aligned, 50);
        assert_eq!(span.skipped, 1000);
        assert_eq!(span.deleted, 0);
        assert_eq!(span.total(), 1050);

        let blocks = ops
            .aligned_blocks(1_000, DeletionPolicy::ExcludeFromBlocks)
            .expect("valid");
        assert_eq!(
            blocks,
            vec![Interval::new(1_000, 1_020), Interval::new(2_020, 2_050)]
        );
        assert_eq!(
            ops.reference_interval(1_000).expect("valid"),
            Interval::new(1_000, 2_050)
        );
    }

    #[test]
    fn deletions_break_blocks_by_default_but_stay_in_the_span() {
        let packed = pack(&[(20, M), (4, D), (10, M)]);
        let ops = CigarOps::new(&packed);

        let span = ops.reference_span().expect("valid");
        assert_eq!(span.aligned, 30);
        assert_eq!(span.deleted, 4);
        assert_eq!(span.total(), 34);

        let excluded = ops
            .aligned_blocks(0, DeletionPolicy::ExcludeFromBlocks)
            .expect("valid");
        assert_eq!(excluded, vec![Interval::new(0, 20), Interval::new(24, 34)]);

        let included = ops
            .aligned_blocks(0, DeletionPolicy::IncludeInBlocks)
            .expect("valid");
        assert_eq!(included, vec![Interval::new(0, 34)]);
    }

    #[test]
    fn insertions_and_clips_do_not_advance_the_reference() {
        let packed = pack(&[(5, H), (10, S), (20, M), (7, I), (20, M), (10, S), (5, H)]);
        let ops = CigarOps::new(&packed);
        let blocks = ops
            .aligned_blocks(100, DeletionPolicy::ExcludeFromBlocks)
            .expect("valid");
        // The insertion consumes no reference, so the two `M` runs are
        // contiguous and merge into one block.
        assert_eq!(blocks, vec![Interval::new(100, 140)]);
        assert_eq!(ops.reference_span().expect("valid").total(), 40);
        assert_eq!(ops.read_length().expect("valid"), 10 + 20 + 7 + 20 + 10);
    }

    #[test]
    fn padding_is_ignored_entirely() {
        let packed = pack(&[(10, M), (5, P), (10, M)]);
        let ops = CigarOps::new(&packed);
        assert_eq!(
            ops.aligned_blocks(0, DeletionPolicy::ExcludeFromBlocks)
                .expect("valid"),
            vec![Interval::new(0, 20)]
        );
    }

    #[test]
    fn sequence_match_and_mismatch_behave_like_match() {
        let packed = pack(&[(10, EQ), (5, X), (5, N), (10, EQ)]);
        let ops = CigarOps::new(&packed);
        assert_eq!(
            ops.aligned_blocks(0, DeletionPolicy::ExcludeFromBlocks)
                .expect("valid"),
            vec![Interval::new(0, 15), Interval::new(20, 30)]
        );
    }

    #[test]
    fn a_leading_deletion_opens_no_block() {
        let packed = pack(&[(5, D), (10, M)]);
        let ops = CigarOps::new(&packed);
        for policy in [
            DeletionPolicy::ExcludeFromBlocks,
            DeletionPolicy::IncludeInBlocks,
        ] {
            assert_eq!(
                ops.aligned_blocks(0, policy).expect("valid"),
                vec![Interval::new(5, 15)],
                "{policy:?}"
            );
        }
    }

    #[test]
    fn an_empty_cigar_yields_no_blocks_but_a_one_base_interval() {
        let ops = CigarOps::new(&[]);
        assert!(ops.is_empty());
        assert_eq!(ops.len(), 0);
        assert_eq!(
            ops.aligned_blocks(42, DeletionPolicy::ExcludeFromBlocks)
                .expect("valid"),
            Vec::<Interval>::new()
        );
        assert_eq!(
            ops.reference_interval(42).expect("valid"),
            Interval::new(42, 43)
        );
        assert_eq!(ops.to_sam_string().expect("valid"), "*");
    }

    #[test]
    fn detects_span_overflow() {
        let packed = pack(&[((1 << 28) - 1, M); 20]);
        let ops = CigarOps::new(&packed);
        let error = ops.reference_span().expect_err("must overflow");
        assert!(matches!(error, CigarError::SpanOverflow), "{error}");
    }

    #[test]
    fn iteration_reports_an_exact_size() {
        let packed = pack(&[(1, M), (2, I), (3, D)]);
        let ops = CigarOps::new(&packed);
        assert_eq!(ops.len(), 3);
        assert_eq!(ops.size_hint(), (3, Some(3)));
        assert_eq!(ops.count(), 3);
    }

    #[test]
    fn renders_sam_text() {
        let packed = pack(&[(10, S), (20, M), (1000, N), (30, M)]);
        assert_eq!(
            CigarOps::new(&packed).to_sam_string().expect("valid"),
            "10S20M1000N30M"
        );
    }

    #[test]
    fn intron_retention_alignment_produces_one_block() {
        // A read that stays aligned across what the annotation calls an intron.
        let packed = pack(&[(200, M)]);
        let ops = CigarOps::new(&packed);
        assert_eq!(
            ops.aligned_blocks(1_000, DeletionPolicy::ExcludeFromBlocks)
                .expect("valid"),
            vec![Interval::new(1_000, 1_200)]
        );
    }

    #[test]
    fn streaming_and_collecting_agree() {
        let packed = pack(&[(10, M), (5, N), (10, M), (3, D), (7, M)]);
        let ops = CigarOps::new(&packed);
        let collected = ops
            .aligned_blocks(0, DeletionPolicy::ExcludeFromBlocks)
            .expect("valid");
        let mut streamed = Vec::new();
        ops.for_each_aligned_block(0, DeletionPolicy::ExcludeFromBlocks, |block| {
            streamed.push(block);
        })
        .expect("valid");
        assert_eq!(collected, streamed);
    }

    #[test]
    fn geometry_renders_its_option_value() {
        assert_eq!(AlignmentGeometry::Blocks.to_string(), "blocks");
        assert_eq!(AlignmentGeometry::Span.to_string(), "span");
        assert_eq!(AlignmentGeometry::default(), AlignmentGeometry::Blocks);
    }
}
