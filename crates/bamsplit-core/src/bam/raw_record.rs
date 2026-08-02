// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! A minimal, bounds-checked reader over decompressed BAM records.
//!
//! This is the layer that makes `bamsplit chrom` fast: a record is never
//! decoded into a large owned Rust value and re-serialized. Instead the reader
//! copies the record body into a reusable buffer once, hands out a borrowed
//! [`RawRecord`] that parses only the fields the active router asks for, and
//! the writer emits the *same bytes* back out.
//!
//! # Layout
//!
//! On disk a BAM record is
//!
//! ```text
//! block_size: u32 LE      <- framing, not part of the body
//! ├─ 0..4    ref_id       i32 LE   (-1 => unplaced)
//! ├─ 4..8    pos          i32 LE   (-1 => no position; otherwise 0-based)
//! ├─ 8       l_read_name  u8       (includes the NUL terminator)
//! ├─ 9       mapq         u8
//! ├─ 10..12  bin          u16 LE   (advisory; bamsplit never reads it)
//! ├─ 12..14  n_cigar_op   u16 LE
//! ├─ 14..16  flag         u16 LE
//! ├─ 16..20  l_seq        i32 LE
//! ├─ 20..24  next_ref_id  i32 LE
//! ├─ 24..28  next_pos     i32 LE
//! ├─ 28..32  tlen         i32 LE
//! ├─ read_name  l_read_name bytes, NUL-terminated
//! ├─ cigar      n_cigar_op * 4 bytes
//! ├─ seq        (l_seq + 1) / 2 bytes
//! ├─ qual       l_seq bytes
//! └─ data       the remainder
//! ```
//!
//! [`RawRecord::raw_bytes`] returns the **body** — everything after
//! `block_size`. That is exactly what [`crate::output`] writes back out, behind
//! a freshly computed `block_size` prefix, so the body is preserved verbatim.
//!
//! # Safety and robustness
//!
//! * There is no `unsafe` in this module and no unchecked indexing on input.
//! * Multi-byte scalars are rebuilt with `from_le_bytes`, so no unaligned
//!   pointer dereference can occur.
//! * Every length is validated before use, and every offset computation uses
//!   `checked_add`/`checked_mul`.
//! * [`RawRecord::new`] fully validates the record shape once, so the field
//!   accessors cannot read out of bounds afterwards.
//!
//! # Example
//!
//! ```
//! use bamsplit_core::bam::raw_record::{RawRecord, RawRecordReader};
//!
//! // A minimal unmapped record: ref_id = -1, pos = -1, name "r0".
//! let mut body = Vec::new();
//! body.extend_from_slice(&(-1i32).to_le_bytes()); // ref_id
//! body.extend_from_slice(&(-1i32).to_le_bytes()); // pos
//! body.push(3);                                   // l_read_name ("r0\0")
//! body.push(255);                                 // mapq
//! body.extend_from_slice(&4680u16.to_le_bytes()); // bin
//! body.extend_from_slice(&0u16.to_le_bytes());    // n_cigar_op
//! body.extend_from_slice(&4u16.to_le_bytes());    // flag = UNMAPPED
//! body.extend_from_slice(&0i32.to_le_bytes());    // l_seq
//! body.extend_from_slice(&(-1i32).to_le_bytes()); // next_ref_id
//! body.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
//! body.extend_from_slice(&0i32.to_le_bytes());    // tlen
//! body.extend_from_slice(b"r0\0");                // read_name
//!
//! let record = RawRecord::new(&body)?;
//! assert_eq!(record.reference_sequence_id()?, None);
//! assert_eq!(record.qname()?, b"r0");
//! assert!(record.flags()? & 0x4 != 0);
//!
//! // Framed as a stream, the reader round-trips the body byte for byte.
//! let mut framed = (body.len() as u32).to_le_bytes().to_vec();
//! framed.extend_from_slice(&body);
//! let mut reader = RawRecordReader::new(&framed[..]);
//! let read_back = reader.read_record()?.expect("one record");
//! assert_eq!(read_back.raw_bytes(), &body[..]);
//! # Ok::<_, bamsplit_core::error::BamRecordError>(())
//! ```

use std::io::Read;

use crate::bam::cigar::{CigarOps, ReferenceSpan};
use crate::bam::tags::{RawTagValue, Tag, TagReader, find_tag};
use crate::error::{BamRecordError, CigarError, TagError};

/// The size of the fixed BAM record core, in bytes.
pub const FIXED_CORE_SIZE: usize = 32;

/// The default upper bound on a single record body, in bytes.
///
/// The SAM specification caps `block_size` at `i32::MAX`, but a legitimate
/// record is orders of magnitude smaller. Rejecting absurd values early turns a
/// corrupt length prefix into a clear error instead of a multi-gigabyte
/// allocation.
pub const DEFAULT_MAX_RECORD_SIZE: usize = 64 * 1024 * 1024;

/// The sentinel `ref_id`/`pos` value meaning "absent".
const MISSING_I32: i32 = -1;

const REFERENCE_SEQUENCE_ID: std::ops::Range<usize> = 0..4;
const ALIGNMENT_START: std::ops::Range<usize> = 4..8;
const NAME_LENGTH: usize = 8;
const MAPPING_QUALITY: usize = 9;
const CIGAR_OP_COUNT: std::ops::Range<usize> = 12..14;
const FLAGS: std::ops::Range<usize> = 14..16;
const SEQUENCE_LENGTH: std::ops::Range<usize> = 16..20;
const MATE_REFERENCE_SEQUENCE_ID: std::ops::Range<usize> = 20..24;
const MATE_ALIGNMENT_START: std::ops::Range<usize> = 24..28;
const TEMPLATE_LENGTH: std::ops::Range<usize> = 28..32;

/// Where a record began, for diagnostics.
///
/// Both coordinates are recorded because either can be the useful one: the
/// ordinal is what a user counts to in `samtools view`, and the virtual offset
/// is what a BGZF-aware tool seeks to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RecordLocation {
    /// 1-based ordinal of the record within the input stream.
    ///
    /// Zero means "not known", which happens when a [`RawRecord`] accessor
    /// fails outside a reader loop; engines call
    /// [`BamRecordError::relocate`] to fill it in.
    pub record_number: u64,
    /// The BGZF virtual offset the record's `block_size` prefix started at.
    pub virtual_offset: u64,
}

impl RecordLocation {
    /// A location that carries no information.
    pub const UNKNOWN: Self = Self {
        record_number: 0,
        virtual_offset: 0,
    };

    /// Creates a location.
    #[must_use]
    pub const fn new(record_number: u64, virtual_offset: u64) -> Self {
        Self {
            record_number,
            virtual_offset,
        }
    }

    /// Whether this location carries no information.
    #[must_use]
    pub const fn is_unknown(&self) -> bool {
        self.record_number == 0 && self.virtual_offset == 0
    }
}

impl std::fmt::Display for RecordLocation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_unknown() {
            return f.write_str("an unknown position");
        }
        write!(
            f,
            "record {} (virtual offset {}:{})",
            self.record_number,
            self.virtual_offset >> 16,
            self.virtual_offset & 0xffff
        )
    }
}

impl BamRecordError {
    /// Fills in an unknown [`RecordLocation`] after the fact.
    ///
    /// [`RawRecord`] holds only a byte slice, exactly as specified, so its
    /// accessors cannot know where the record came from. Readers and engines,
    /// which do know, call this so the user-facing message names the record.
    #[must_use]
    pub fn relocate(mut self, location: RecordLocation) -> Self {
        let slot: &mut RecordLocation = match &mut self {
            Self::NegativeBlockSize { location, .. }
            | Self::BlockSizeTooSmall { location, .. }
            | Self::BlockSizeTooLarge { location, .. }
            | Self::TruncatedRecord { location, .. }
            | Self::TruncatedBlockSize { location, .. }
            | Self::FieldOutOfBounds { location, .. }
            | Self::LayoutOverflow { location, .. }
            | Self::EmptyReadName { location }
            | Self::UnterminatedReadName { location }
            | Self::NegativeSequenceLength { location, .. }
            | Self::InvalidReferenceSequenceId { location, .. }
            | Self::InvalidPosition { location, .. }
            | Self::Tag { location, .. }
            | Self::Cigar { location, .. }
            | Self::Io { location, .. } => location,
        };
        if slot.is_unknown() {
            *slot = location;
        }
        self
    }
}

/// A borrowed, bounds-checked view over one BAM record body.
///
/// The body excludes the 4-byte `block_size` prefix; see the module
/// documentation for the exact layout.
///
/// Construct one with [`RawRecord::new`], which validates the record shape
/// once. After that the fixed-core accessors are infallible in practice — they
/// still return [`Result`] so that the API stays stable if a future field needs
/// deeper validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RawRecord<'a> {
    bytes: &'a [u8],
}

impl<'a> RawRecord<'a> {
    /// Validates and wraps a record body.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the 32-byte fixed
    /// core, if `l_seq` is negative, if the read name is empty or not
    /// NUL-terminated, or if the variable-length sections do not fit inside the
    /// body.
    ///
    /// The returned errors carry [`RecordLocation::UNKNOWN`]; callers that know
    /// where the record came from should apply [`BamRecordError::relocate`].
    pub fn new(bytes: &'a [u8]) -> Result<Self, BamRecordError> {
        let record = Self { bytes };
        record.validate_layout()?;
        Ok(record)
    }

    /// Wraps a body that a caller has already validated.
    ///
    /// This exists for the spool engine, which re-reads bodies it validated on
    /// the way in and does not need to pay for a second shape check. Passing an
    /// unvalidated slice is not *unsound* — every accessor still bounds-checks
    /// — but accessors may then return errors the caller did not expect.
    #[must_use]
    pub const fn from_validated_bytes(bytes: &'a [u8]) -> Self {
        Self { bytes }
    }

    /// The exact record body, ready to be written back out unchanged.
    #[must_use]
    pub const fn raw_bytes(&self) -> &'a [u8] {
        self.bytes
    }

    /// The length of the record body, which is also its `block_size`.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Whether the body is empty. Always `false` for a validated record.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    fn validate_layout(&self) -> Result<(), BamRecordError> {
        if self.bytes.len() < FIXED_CORE_SIZE {
            return Err(BamRecordError::BlockSizeTooSmall {
                block_size: self.bytes.len(),
                minimum: FIXED_CORE_SIZE,
                location: RecordLocation::UNKNOWN,
            });
        }

        let name_len = usize::from(self.bytes[NAME_LENGTH]);
        if name_len == 0 {
            return Err(BamRecordError::EmptyReadName {
                location: RecordLocation::UNKNOWN,
            });
        }

        let l_seq_raw = self.read_i32(SEQUENCE_LENGTH);
        if l_seq_raw < 0 {
            return Err(BamRecordError::NegativeSequenceLength {
                l_seq: l_seq_raw,
                location: RecordLocation::UNKNOWN,
            });
        }
        let l_seq = l_seq_raw as usize;

        let cigar_ops = usize::from(self.read_u16(CIGAR_OP_COUNT));
        let cigar_len = cigar_ops
            .checked_mul(4)
            .ok_or(BamRecordError::LayoutOverflow {
                context: "n_cigar_op * 4",
                location: RecordLocation::UNKNOWN,
            })?;

        let seq_len = l_seq.div_ceil(2);

        let mut end = FIXED_CORE_SIZE;
        for (field, width) in [
            ("read_name", name_len),
            ("cigar", cigar_len),
            ("seq", seq_len),
            ("qual", l_seq),
        ] {
            let next = end
                .checked_add(width)
                .ok_or(BamRecordError::LayoutOverflow {
                    context: "variable-length section offset",
                    location: RecordLocation::UNKNOWN,
                })?;
            if next > self.bytes.len() {
                return Err(BamRecordError::FieldOutOfBounds {
                    field: match field {
                        "read_name" => "read_name",
                        "cigar" => "cigar",
                        "seq" => "seq",
                        _ => "qual",
                    },
                    offset: end,
                    needed: width,
                    body_len: self.bytes.len(),
                    location: RecordLocation::UNKNOWN,
                });
            }
            end = next;
        }

        // `name_len` counts the NUL terminator, so the last name byte must be
        // zero. `name_len >= 1` was checked above, so this index is in bounds.
        if self.bytes[FIXED_CORE_SIZE + name_len - 1] != 0 {
            return Err(BamRecordError::UnterminatedReadName {
                location: RecordLocation::UNKNOWN,
            });
        }

        Ok(())
    }

    fn read_i32(&self, range: std::ops::Range<usize>) -> i32 {
        // Callers only pass fixed-core ranges, and `validate_layout` has
        // established that the body is at least `FIXED_CORE_SIZE` bytes.
        let slice = &self.bytes[range];
        i32::from_le_bytes([slice[0], slice[1], slice[2], slice[3]])
    }

    fn read_u16(&self, range: std::ops::Range<usize>) -> u16 {
        let slice = &self.bytes[range];
        u16::from_le_bytes([slice[0], slice[1]])
    }

    fn require_core(&self) -> Result<(), BamRecordError> {
        if self.bytes.len() < FIXED_CORE_SIZE {
            return Err(BamRecordError::BlockSizeTooSmall {
                block_size: self.bytes.len(),
                minimum: FIXED_CORE_SIZE,
                location: RecordLocation::UNKNOWN,
            });
        }
        Ok(())
    }

    /// The reference sequence id, or [`None`] for an unplaced record.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core,
    /// or if the id is negative but not the `-1` sentinel.
    pub fn reference_sequence_id(&self) -> Result<Option<i32>, BamRecordError> {
        self.require_core()?;
        decode_reference_id(self.read_i32(REFERENCE_SEQUENCE_ID), "reference")
    }

    /// The mate's reference sequence id, or [`None`] when absent.
    ///
    /// # Errors
    ///
    /// As [`reference_sequence_id`](Self::reference_sequence_id).
    pub fn mate_reference_sequence_id(&self) -> Result<Option<i32>, BamRecordError> {
        self.require_core()?;
        decode_reference_id(self.read_i32(MATE_REFERENCE_SEQUENCE_ID), "mate reference")
    }

    /// The 0-based leftmost alignment position, or [`None`] when absent.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core,
    /// or if the position is negative but not the `-1` sentinel.
    pub fn alignment_start(&self) -> Result<Option<i32>, BamRecordError> {
        self.require_core()?;
        decode_position(self.read_i32(ALIGNMENT_START), "alignment")
    }

    /// The mate's 0-based leftmost alignment position, or [`None`].
    ///
    /// # Errors
    ///
    /// As [`alignment_start`](Self::alignment_start).
    pub fn mate_alignment_start(&self) -> Result<Option<i32>, BamRecordError> {
        self.require_core()?;
        decode_position(self.read_i32(MATE_ALIGNMENT_START), "mate alignment")
    }

    /// The SAM flags.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core.
    pub fn flags(&self) -> Result<u16, BamRecordError> {
        self.require_core()?;
        Ok(self.read_u16(FLAGS))
    }

    /// The mapping quality. `255` means "unavailable".
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core.
    pub fn mapping_quality(&self) -> Result<u8, BamRecordError> {
        self.require_core()?;
        Ok(self.bytes[MAPPING_QUALITY])
    }

    /// The template length (`TLEN`).
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core.
    pub fn template_length(&self) -> Result<i32, BamRecordError> {
        self.require_core()?;
        Ok(self.read_i32(TEMPLATE_LENGTH))
    }

    /// The number of bases in the read.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core or
    /// `l_seq` is negative.
    pub fn sequence_length(&self) -> Result<usize, BamRecordError> {
        self.require_core()?;
        let raw = self.read_i32(SEQUENCE_LENGTH);
        if raw < 0 {
            return Err(BamRecordError::NegativeSequenceLength {
                l_seq: raw,
                location: RecordLocation::UNKNOWN,
            });
        }
        Ok(raw as usize)
    }

    /// The read name, without its NUL terminator.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the name is empty, unterminated, or
    /// overruns the body.
    pub fn qname(&self) -> Result<&'a [u8], BamRecordError> {
        self.require_core()?;
        let name_len = usize::from(self.bytes[NAME_LENGTH]);
        if name_len == 0 {
            return Err(BamRecordError::EmptyReadName {
                location: RecordLocation::UNKNOWN,
            });
        }
        let end = FIXED_CORE_SIZE
            .checked_add(name_len)
            .ok_or(BamRecordError::LayoutOverflow {
                context: "read_name offset",
                location: RecordLocation::UNKNOWN,
            })?;
        if end > self.bytes.len() {
            return Err(BamRecordError::FieldOutOfBounds {
                field: "read_name",
                offset: FIXED_CORE_SIZE,
                needed: name_len,
                body_len: self.bytes.len(),
                location: RecordLocation::UNKNOWN,
            });
        }
        let raw = &self.bytes[FIXED_CORE_SIZE..end];
        // `name_len >= 1`, so `raw` is non-empty.
        let (last, head) = raw.split_last().unwrap_or((&0, &[]));
        if *last != 0 {
            return Err(BamRecordError::UnterminatedReadName {
                location: RecordLocation::UNKNOWN,
            });
        }
        Ok(head)
    }

    /// The raw packed CIGAR operations, as they appear in the record.
    ///
    /// This does **not** resolve the long-CIGAR `CG` convention; use
    /// [`cigar_ops`](Self::cigar_ops) for that.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the CIGAR section overruns the body.
    pub fn packed_cigar(&self) -> Result<&'a [u8], BamRecordError> {
        self.require_core()?;
        let (start, len) = self.cigar_bounds()?;
        Ok(&self.bytes[start..start + len])
    }

    fn cigar_bounds(&self) -> Result<(usize, usize), BamRecordError> {
        let name_len = usize::from(self.bytes[NAME_LENGTH]);
        let ops = usize::from(self.read_u16(CIGAR_OP_COUNT));
        let len = ops.checked_mul(4).ok_or(BamRecordError::LayoutOverflow {
            context: "n_cigar_op * 4",
            location: RecordLocation::UNKNOWN,
        })?;
        let start =
            FIXED_CORE_SIZE
                .checked_add(name_len)
                .ok_or(BamRecordError::LayoutOverflow {
                    context: "cigar offset",
                    location: RecordLocation::UNKNOWN,
                })?;
        let end = start
            .checked_add(len)
            .ok_or(BamRecordError::LayoutOverflow {
                context: "cigar offset",
                location: RecordLocation::UNKNOWN,
            })?;
        if end > self.bytes.len() {
            return Err(BamRecordError::FieldOutOfBounds {
                field: "cigar",
                offset: start,
                needed: len,
                body_len: self.bytes.len(),
                location: RecordLocation::UNKNOWN,
            });
        }
        Ok((start, len))
    }

    /// The auxiliary data section.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the preceding sections overrun the body.
    pub fn data(&self) -> Result<&'a [u8], BamRecordError> {
        self.require_core()?;
        let name_len = usize::from(self.bytes[NAME_LENGTH]);
        let l_seq = self.sequence_length()?;
        let (cigar_start, cigar_len) = self.cigar_bounds()?;
        let _ = name_len;

        let mut offset = cigar_start + cigar_len;
        for (field, width) in [("seq", l_seq.div_ceil(2)), ("qual", l_seq)] {
            let next = offset
                .checked_add(width)
                .ok_or(BamRecordError::LayoutOverflow {
                    context: "auxiliary data offset",
                    location: RecordLocation::UNKNOWN,
                })?;
            if next > self.bytes.len() {
                return Err(BamRecordError::FieldOutOfBounds {
                    field: if field == "seq" { "seq" } else { "qual" },
                    offset,
                    needed: width,
                    body_len: self.bytes.len(),
                    location: RecordLocation::UNKNOWN,
                });
            }
            offset = next;
        }
        Ok(&self.bytes[offset..])
    }

    /// Looks up one auxiliary tag.
    ///
    /// Decoding stops at the first match, so this is cheap for the common case
    /// of routing on a single tag near the front of the section.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError::Tag`] if the auxiliary section is malformed
    /// before the requested tag is found.
    pub fn tag(&self, tag: Tag) -> Result<Option<RawTagValue<'a>>, BamRecordError> {
        let data = self.data()?;
        find_tag(data, tag).map_err(|source| BamRecordError::Tag {
            source: Box::new(source),
            location: RecordLocation::UNKNOWN,
        })
    }

    /// Iterates every auxiliary field.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the preceding sections overrun the body.
    /// Per-field errors surface lazily from the returned reader.
    pub fn tags(&self) -> Result<TagReader<'a>, BamRecordError> {
        Ok(TagReader::new(self.data()?))
    }

    /// The effective CIGAR operations, resolving the long-CIGAR `CG`
    /// convention when it applies.
    ///
    /// Per SAM §4.2.2, a record with more than 65 535 CIGAR operations stores a
    /// two-operation placeholder (`<l_seq>S<ref_span>N`) inline and the real
    /// operations in a `CG:B:I` tag. This method transparently returns the
    /// real operations in that case.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the CIGAR or auxiliary sections are
    /// malformed, or if the `CG` convention is signalled but the tag is absent
    /// or has the wrong type.
    pub fn cigar_ops(&self) -> Result<CigarOps<'a>, BamRecordError> {
        let packed = self.packed_cigar()?;
        if !self.uses_long_cigar_placeholder(packed)? {
            return Ok(CigarOps::new(packed));
        }

        let Some(value) = self.tag(*b"CG")? else {
            return Err(BamRecordError::Cigar {
                source: Box::new(CigarError::MissingLongCigarTag),
                location: RecordLocation::UNKNOWN,
            });
        };
        let Some(array) = value.as_array() else {
            return Err(BamRecordError::Cigar {
                source: Box::new(CigarError::MalformedLongCigar {
                    reason: "the `CG` tag must be a `B:I` array of packed CIGAR operations",
                }),
                location: RecordLocation::UNKNOWN,
            });
        };
        if array.subtype() != crate::bam::tags::ArraySubtype::UInt32 {
            return Err(BamRecordError::Cigar {
                source: Box::new(CigarError::MalformedLongCigar {
                    reason: "the `CG` tag must use the `I` (u32) subtype",
                }),
                location: RecordLocation::UNKNOWN,
            });
        }
        Ok(CigarOps::new(array.payload()))
    }

    /// Whether the inline CIGAR is the two-operation long-CIGAR placeholder.
    fn uses_long_cigar_placeholder(&self, packed: &[u8]) -> Result<bool, BamRecordError> {
        // op 0 must be `<l_seq>S`, op 1 must be an `N`.
        const SOFT_CLIP: u32 = 4;
        const SKIP: u32 = 3;

        if packed.len() != 8 {
            return Ok(false);
        }
        let first = u32::from_le_bytes([packed[0], packed[1], packed[2], packed[3]]);
        let second = u32::from_le_bytes([packed[4], packed[5], packed[6], packed[7]]);
        if first & 0xf != SOFT_CLIP || second & 0xf != SKIP {
            return Ok(false);
        }
        let l_seq = self.sequence_length()?;
        Ok(u64::from(first >> 4) == l_seq as u64)
    }

    /// The reference span of the alignment, in bases.
    ///
    /// Only `M`, `D`, `N`, `=`, and `X` consume reference bases.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the CIGAR is malformed.
    pub fn reference_span(&self) -> Result<ReferenceSpan, BamRecordError> {
        let ops = self.cigar_ops()?;
        ops.reference_span()
            .map_err(|source| BamRecordError::Cigar {
                source: Box::new(source),
                location: RecordLocation::UNKNOWN,
            })
    }

    /// The 0-based, **exclusive** end of the alignment on the reference.
    ///
    /// Together with [`alignment_start`](Self::alignment_start) this forms the
    /// half-open interval `[start, end)`. A record whose CIGAR consumes no
    /// reference bases still occupies one position, matching `bam_endpos` in
    /// htslib; that keeps `end > start` for every placed record so downstream
    /// interval logic never sees an empty span.
    ///
    /// Returns [`None`] when the record has no position.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the CIGAR is malformed or the span
    /// overflows `i64`.
    pub fn alignment_end(&self) -> Result<Option<i64>, BamRecordError> {
        let Some(start) = self.alignment_start()? else {
            return Ok(None);
        };
        let span = self.reference_span()?;
        let effective = i64::from(span.total().max(1));
        i64::from(start)
            .checked_add(effective)
            .ok_or(BamRecordError::Cigar {
                source: Box::new(CigarError::SpanOverflow),
                location: RecordLocation::UNKNOWN,
            })
            .map(Some)
    }

    /// Whether the `UNMAPPED` (`0x4`) flag is set.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core.
    pub fn is_unmapped(&self) -> Result<bool, BamRecordError> {
        Ok(self.flags()? & 0x4 != 0)
    }

    /// Whether the `SECONDARY` (`0x100`) flag is set.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core.
    pub fn is_secondary(&self) -> Result<bool, BamRecordError> {
        Ok(self.flags()? & 0x100 != 0)
    }

    /// Whether the `SUPPLEMENTARY` (`0x800`) flag is set.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core.
    pub fn is_supplementary(&self) -> Result<bool, BamRecordError> {
        Ok(self.flags()? & 0x800 != 0)
    }

    /// Whether the record is neither secondary nor supplementary.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] if the body is shorter than the fixed core.
    pub fn is_primary(&self) -> Result<bool, BamRecordError> {
        let flags = self.flags()?;
        Ok(flags & 0x100 == 0 && flags & 0x800 == 0)
    }

    /// Fully validates the record, including its auxiliary section.
    ///
    /// [`RawRecord::new`] checks the record *shape*; this additionally walks
    /// every auxiliary field. It is used by `bamsplit inspect --full` and by
    /// the differential test-suite, not on the hot path.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] for any structural problem.
    pub fn validate_deeply(&self, reject_duplicate_tags: bool) -> Result<(), BamRecordError> {
        self.validate_layout()?;
        self.qname()?;
        for op in self.cigar_ops()? {
            op.map_err(|source| BamRecordError::Cigar {
                source: Box::new(source),
                location: RecordLocation::UNKNOWN,
            })?;
        }
        crate::bam::tags::validate_section(self.data()?, reject_duplicate_tags).map_err(
            |source: TagError| BamRecordError::Tag {
                source: Box::new(source),
                location: RecordLocation::UNKNOWN,
            },
        )?;
        Ok(())
    }
}

fn decode_reference_id(raw: i32, kind: &'static str) -> Result<Option<i32>, BamRecordError> {
    match raw {
        MISSING_I32 => Ok(None),
        value if value < 0 => Err(BamRecordError::InvalidReferenceSequenceId {
            kind,
            id: value,
            reference_count: 0,
            location: RecordLocation::UNKNOWN,
        }),
        value => Ok(Some(value)),
    }
}

fn decode_position(raw: i32, kind: &'static str) -> Result<Option<i32>, BamRecordError> {
    match raw {
        MISSING_I32 => Ok(None),
        value if value < 0 => Err(BamRecordError::InvalidPosition {
            kind,
            position: value,
            location: RecordLocation::UNKNOWN,
        }),
        value => Ok(Some(value)),
    }
}

/// Something that can report a BGZF virtual position.
///
/// Implemented for the BGZF readers `bamsplit` uses, so the record reader can
/// tag every failure with the offset it happened at without depending on which
/// decompressor is underneath.
pub trait VirtualPositionSource {
    /// The current virtual position, as a packed `u64`.
    fn virtual_position(&self) -> u64;
}

impl<R: std::io::Read> VirtualPositionSource for noodles_bgzf::io::Reader<R> {
    fn virtual_position(&self) -> u64 {
        u64::from(noodles_bgzf::io::Reader::virtual_position(self))
    }
}

impl<R> VirtualPositionSource for noodles_bgzf::io::MultithreadedReader<R> {
    fn virtual_position(&self) -> u64 {
        u64::from(noodles_bgzf::io::MultithreadedReader::virtual_position(
            self,
        ))
    }
}

/// A reader that reports a virtual position of zero.
///
/// Wrap a plain [`Read`] in this when offsets are not meaningful, for example
/// when replaying an uncompressed spool file.
#[derive(Debug)]
pub struct Unpositioned<R>(pub R);

impl<R: Read> Read for Unpositioned<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.0.read(buf)
    }
}

impl<R> VirtualPositionSource for Unpositioned<R> {
    fn virtual_position(&self) -> u64 {
        0
    }
}

impl VirtualPositionSource for &[u8] {
    fn virtual_position(&self) -> u64 {
        0
    }
}

/// A streaming reader over framed BAM records.
///
/// The reader owns one growable buffer and reuses it for every record, so a
/// full pass over a BAM performs a constant number of heap allocations rather
/// than one per record.
///
/// The buffer is exposed through [`RawRecordReader::read_record`], which
/// returns a [`RawRecord`] borrowing it. Because the borrow ties up the reader,
/// callers that need to keep a record across iterations should copy
/// [`RawRecord::raw_bytes`].
#[derive(Debug)]
pub struct RawRecordReader<R> {
    inner: R,
    buffer: Vec<u8>,
    max_record_size: usize,
    record_number: u64,
    last_location: RecordLocation,
}

impl<R> RawRecordReader<R> {
    /// Creates a reader with the default record-size limit.
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buffer: Vec::with_capacity(8 * 1024),
            max_record_size: DEFAULT_MAX_RECORD_SIZE,
            record_number: 0,
            last_location: RecordLocation::UNKNOWN,
        }
    }

    /// Overrides the maximum accepted `block_size`.
    ///
    /// Values above `i32::MAX` are clamped, because `block_size` is a signed
    /// 32-bit field on disk and a larger limit could never trigger.
    #[must_use]
    pub fn with_max_record_size(mut self, limit: usize) -> Self {
        self.max_record_size = limit.min(i32::MAX as usize);
        self
    }

    /// The number of records read so far.
    #[must_use]
    pub const fn record_count(&self) -> u64 {
        self.record_number
    }

    /// Where the most recently read record began.
    #[must_use]
    pub const fn last_location(&self) -> RecordLocation {
        self.last_location
    }

    /// Borrows the underlying reader.
    pub const fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Mutably borrows the underlying reader.
    pub const fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R: Read + VirtualPositionSource> RawRecordReader<R> {
    /// Reads the next record into the internal buffer.
    ///
    /// Returns [`None`] at a clean end of stream.
    ///
    /// # Errors
    ///
    /// Returns [`BamRecordError`] for a negative, undersized, or oversized
    /// `block_size`, for a truncated record, or for an I/O failure. Every error
    /// carries the record ordinal and the BGZF virtual offset the record began
    /// at.
    pub fn read_record(&mut self) -> Result<Option<RawRecord<'_>>, BamRecordError> {
        let virtual_offset = self.inner.virtual_position();
        let location = RecordLocation::new(self.record_number + 1, virtual_offset);

        let mut prefix = [0u8; 4];
        match read_full(&mut self.inner, &mut prefix) {
            Ok(0) => return Ok(None),
            Ok(4) => {}
            Ok(partial) => {
                return Err(BamRecordError::TruncatedBlockSize {
                    actual: partial,
                    location,
                });
            }
            Err(source) => return Err(BamRecordError::Io { source, location }),
        }

        let block_size = i32::from_le_bytes(prefix);
        if block_size < 0 {
            return Err(BamRecordError::NegativeBlockSize {
                block_size,
                location,
            });
        }
        let block_size = block_size as usize;
        if block_size < FIXED_CORE_SIZE {
            return Err(BamRecordError::BlockSizeTooSmall {
                block_size,
                minimum: FIXED_CORE_SIZE,
                location,
            });
        }
        if block_size > self.max_record_size {
            return Err(BamRecordError::BlockSizeTooLarge {
                block_size,
                limit: self.max_record_size,
                location,
            });
        }

        self.buffer.clear();
        self.buffer.resize(block_size, 0);
        match read_full(&mut self.inner, &mut self.buffer) {
            Ok(read) if read == block_size => {}
            Ok(read) => {
                return Err(BamRecordError::TruncatedRecord {
                    expected: block_size,
                    actual: read,
                    location,
                });
            }
            Err(source) => return Err(BamRecordError::Io { source, location }),
        }

        self.record_number += 1;
        self.last_location = location;

        RawRecord::new(&self.buffer)
            .map(Some)
            .map_err(|error| error.relocate(location))
    }
}

/// Reads until `buf` is full, the stream ends, or an error occurs.
///
/// Returns how many bytes were read. Unlike [`Read::read_exact`], a short read
/// is reported rather than turned into an error, so the caller can distinguish
/// "clean EOF" from "truncated record".
fn read_full<R: Read>(reader: &mut R, buf: &mut [u8]) -> std::io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(filled)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a minimal but valid record body.
    pub(crate) fn build_record(
        reference_id: i32,
        position: i32,
        name: &[u8],
        flags: u16,
        cigar: &[(u8, u32)],
        sequence_length: usize,
        data: &[u8],
    ) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&reference_id.to_le_bytes());
        body.extend_from_slice(&position.to_le_bytes());
        body.push(u8::try_from(name.len() + 1).expect("short name"));
        body.push(60);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(
            &u16::try_from(cigar.len())
                .expect("short cigar")
                .to_le_bytes(),
        );
        body.extend_from_slice(&flags.to_le_bytes());
        body.extend_from_slice(
            &i32::try_from(sequence_length)
                .expect("short sequence")
                .to_le_bytes(),
        );
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(name);
        body.push(0);
        for (kind, len) in cigar {
            body.extend_from_slice(&((len << 4) | u32::from(*kind)).to_le_bytes());
        }
        body.extend(std::iter::repeat_n(0u8, sequence_length.div_ceil(2)));
        body.extend(std::iter::repeat_n(0xffu8, sequence_length));
        body.extend_from_slice(data);
        body
    }

    fn frame(body: &[u8]) -> Vec<u8> {
        let mut out = (body.len() as u32).to_le_bytes().to_vec();
        out.extend_from_slice(body);
        out
    }

    #[test]
    fn parses_the_fixed_core() {
        let body = build_record(
            2,
            99,
            b"read1",
            0x63,
            &[(0, 100)],
            100,
            b"NMi\x01\x00\x00\x00",
        );
        let record = RawRecord::new(&body).expect("valid");

        assert_eq!(record.reference_sequence_id().expect("valid"), Some(2));
        assert_eq!(record.alignment_start().expect("valid"), Some(99));
        assert_eq!(record.mate_reference_sequence_id().expect("valid"), None);
        assert_eq!(record.mate_alignment_start().expect("valid"), None);
        assert_eq!(record.flags().expect("valid"), 0x63);
        assert_eq!(record.mapping_quality().expect("valid"), 60);
        assert_eq!(record.qname().expect("valid"), b"read1");
        assert_eq!(record.sequence_length().expect("valid"), 100);
        assert_eq!(record.template_length().expect("valid"), 0);
        assert_eq!(record.alignment_end().expect("valid"), Some(199));
        assert_eq!(record.len(), body.len());
        assert_eq!(record.raw_bytes(), &body[..]);
    }

    #[test]
    fn unplaced_records_report_none() {
        let body = build_record(-1, -1, b"r", 0x4, &[], 0, b"");
        let record = RawRecord::new(&body).expect("valid");
        assert_eq!(record.reference_sequence_id().expect("valid"), None);
        assert_eq!(record.alignment_start().expect("valid"), None);
        assert_eq!(record.alignment_end().expect("valid"), None);
        assert!(record.is_unmapped().expect("valid"));
    }

    #[test]
    fn placed_record_with_no_reference_consuming_cigar_spans_one_base() {
        // 10I consumes no reference bases; htslib's `bam_endpos` still yields
        // pos + 1, and so do we.
        let body = build_record(0, 50, b"r", 0, &[(1, 10)], 10, b"");
        let record = RawRecord::new(&body).expect("valid");
        assert_eq!(record.alignment_end().expect("valid"), Some(51));
    }

    #[test]
    fn rejects_a_body_shorter_than_the_fixed_core() {
        let error = RawRecord::new(&[0u8; 31]).expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::BlockSizeTooSmall { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_zero_length_read_name() {
        let mut body = build_record(0, 0, b"r", 0, &[], 0, b"");
        body[NAME_LENGTH] = 0;
        let error = RawRecord::new(&body).expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::EmptyReadName { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_an_unterminated_read_name() {
        let mut body = build_record(0, 0, b"read", 0, &[], 0, b"");
        let terminator = FIXED_CORE_SIZE + 4;
        body[terminator] = b'X';
        let error = RawRecord::new(&body).expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::UnterminatedReadName { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_negative_sequence_length() {
        let mut body = build_record(0, 0, b"r", 0, &[], 0, b"");
        body[SEQUENCE_LENGTH].clone_from_slice(&(-5i32).to_le_bytes());
        let error = RawRecord::new(&body).expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::NegativeSequenceLength { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_cigar_that_overruns_the_body() {
        let mut body = build_record(0, 0, b"r", 0, &[(0, 1)], 0, b"");
        body[CIGAR_OP_COUNT].clone_from_slice(&1000u16.to_le_bytes());
        let error = RawRecord::new(&body).expect_err("must reject");
        assert!(
            matches!(
                error,
                BamRecordError::FieldOutOfBounds { field: "cigar", .. }
            ),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_sequence_that_overruns_the_body() {
        let mut body = build_record(0, 0, b"r", 0, &[], 4, b"");
        body[SEQUENCE_LENGTH].clone_from_slice(&100_000i32.to_le_bytes());
        let error = RawRecord::new(&body).expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::FieldOutOfBounds { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_an_out_of_range_negative_reference_id() {
        let mut body = build_record(0, 0, b"r", 0, &[], 0, b"");
        body[REFERENCE_SEQUENCE_ID].clone_from_slice(&(-7i32).to_le_bytes());
        let record = RawRecord::new(&body).expect("shape is valid");
        let error = record.reference_sequence_id().expect_err("must reject");
        assert!(
            matches!(
                error,
                BamRecordError::InvalidReferenceSequenceId { id: -7, .. }
            ),
            "{error}"
        );
    }

    #[test]
    fn rejects_an_out_of_range_negative_position() {
        let mut body = build_record(0, 0, b"r", 0, &[], 0, b"");
        body[ALIGNMENT_START].clone_from_slice(&(-9i32).to_le_bytes());
        let record = RawRecord::new(&body).expect("shape is valid");
        let error = record.alignment_start().expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::InvalidPosition { position: -9, .. }),
            "{error}"
        );
    }

    #[test]
    fn finds_auxiliary_tags() {
        let body = build_record(0, 0, b"r", 0, &[], 0, b"NMi\x07\x00\x00\x00RGZgroup\x00");
        let record = RawRecord::new(&body).expect("valid");
        let nm = record.tag(*b"NM").expect("valid").expect("present");
        assert_eq!(nm.as_integer(), Some(7));
        let rg = record.tag(*b"RG").expect("valid").expect("present");
        assert_eq!(rg.as_bytes(), Some(&b"group"[..]));
        assert!(record.tag(*b"ZZ").expect("valid").is_none());
    }

    #[test]
    fn reads_a_framed_stream_and_preserves_bodies() {
        let bodies = [
            build_record(0, 10, b"a", 0, &[(0, 5)], 5, b""),
            build_record(0, 20, b"bb", 0, &[(0, 5)], 5, b"NMi\x00\x00\x00\x00"),
            build_record(-1, -1, b"ccc", 0x4, &[], 0, b""),
        ];
        let mut stream = Vec::new();
        for body in &bodies {
            stream.extend_from_slice(&frame(body));
        }

        let mut reader = RawRecordReader::new(&stream[..]);
        let mut seen = Vec::new();
        while let Some(record) = reader.read_record().expect("valid") {
            seen.push(record.raw_bytes().to_vec());
        }
        assert_eq!(seen.len(), 3);
        for (expected, actual) in bodies.iter().zip(&seen) {
            assert_eq!(expected, actual);
        }
        assert_eq!(reader.record_count(), 3);
    }

    #[test]
    fn reports_a_clean_end_of_stream() {
        let mut reader = RawRecordReader::new(&[][..]);
        assert!(reader.read_record().expect("valid").is_none());
    }

    #[test]
    fn rejects_a_negative_block_size() {
        let mut stream = (-1i32).to_le_bytes().to_vec();
        stream.extend_from_slice(&[0u8; 32]);
        let mut reader = RawRecordReader::new(&stream[..]);
        let error = reader.read_record().expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::NegativeBlockSize { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_an_undersized_block_size() {
        let mut stream = 8u32.to_le_bytes().to_vec();
        stream.extend_from_slice(&[0u8; 8]);
        let mut reader = RawRecordReader::new(&stream[..]);
        let error = reader.read_record().expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::BlockSizeTooSmall { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_an_oversized_block_size() {
        let mut stream = 5_000_000u32.to_le_bytes().to_vec();
        stream.extend_from_slice(&[0u8; 64]);
        let mut reader = RawRecordReader::new(&stream[..]).with_max_record_size(1024);
        let error = reader.read_record().expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::BlockSizeTooLarge { limit: 1024, .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_truncated_record_body() {
        let body = build_record(0, 0, b"r", 0, &[], 0, b"");
        let mut stream = frame(&body);
        stream.truncate(stream.len() - 4);
        let mut reader = RawRecordReader::new(&stream[..]);
        let error = reader.read_record().expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::TruncatedRecord { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_a_truncated_length_prefix() {
        let stream = [1u8, 2, 3];
        let mut reader = RawRecordReader::new(&stream[..]);
        let error = reader.read_record().expect_err("must reject");
        assert!(
            matches!(error, BamRecordError::TruncatedBlockSize { actual: 3, .. }),
            "{error}"
        );
    }

    #[test]
    fn errors_carry_the_record_ordinal() {
        let good = build_record(0, 0, b"r", 0, &[], 0, b"");
        let mut stream = frame(&good);
        stream.extend_from_slice(&(-1i32).to_le_bytes());
        stream.extend_from_slice(&[0u8; 32]);

        let mut reader = RawRecordReader::new(&stream[..]);
        assert!(reader.read_record().expect("valid").is_some());
        let error = reader.read_record().expect_err("must reject");
        assert_eq!(error.location().record_number, 2);
    }

    #[test]
    fn relocate_only_fills_unknown_locations() {
        let error = BamRecordError::EmptyReadName {
            location: RecordLocation::UNKNOWN,
        };
        let located = error.relocate(RecordLocation::new(7, 42));
        assert_eq!(located.location().record_number, 7);

        let already = BamRecordError::EmptyReadName {
            location: RecordLocation::new(1, 1),
        };
        let unchanged = already.relocate(RecordLocation::new(7, 42));
        assert_eq!(unchanged.location().record_number, 1);
    }

    #[test]
    fn resolves_the_long_cigar_convention() {
        // Placeholder: 8S 100N, real CIGAR in `CG:B:I` as 8M.
        let mut data = Vec::new();
        data.extend_from_slice(b"CGB");
        data.push(b'I');
        data.extend_from_slice(&1u32.to_le_bytes());
        data.extend_from_slice(&(8u32 << 4).to_le_bytes()); // 8M (op code 0)

        let body = build_record(0, 0, b"r", 0, &[(4, 8), (3, 100)], 8, &data);
        let record = RawRecord::new(&body).expect("valid");
        let ops: Vec<_> = record
            .cigar_ops()
            .expect("valid")
            .collect::<Result<_, _>>()
            .expect("valid ops");
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0].length(), 8);
        assert_eq!(record.reference_span().expect("valid").total(), 8);
    }

    #[test]
    fn does_not_misread_a_genuine_two_op_cigar_as_a_long_cigar() {
        // 8S100N is a real (if unusual) CIGAR when there is no `CG` tag and
        // l_seq matches; without the tag we must report a clear error rather
        // than silently mis-parsing.
        let body = build_record(0, 0, b"r", 0, &[(4, 8), (3, 100)], 8, b"");
        let record = RawRecord::new(&body).expect("valid");
        let error = record.cigar_ops().expect_err("must reject");
        assert!(matches!(error, BamRecordError::Cigar { .. }), "{error}");

        // A two-op CIGAR whose soft clip does not match `l_seq` is ordinary.
        let body = build_record(0, 0, b"r", 0, &[(4, 3), (3, 100)], 8, b"");
        let record = RawRecord::new(&body).expect("valid");
        assert_eq!(record.reference_span().expect("valid").total(), 100);
    }

    #[test]
    fn rejects_a_long_cigar_tag_with_the_wrong_subtype() {
        let mut data = Vec::new();
        data.extend_from_slice(b"CGB");
        data.push(b'S');
        data.extend_from_slice(&1u32.to_le_bytes());
        data.extend_from_slice(&1u16.to_le_bytes());

        let body = build_record(0, 0, b"r", 0, &[(4, 8), (3, 100)], 8, &data);
        let record = RawRecord::new(&body).expect("valid");
        let error = record.cigar_ops().expect_err("must reject");
        assert!(matches!(error, BamRecordError::Cigar { .. }), "{error}");
    }

    #[test]
    fn deep_validation_catches_malformed_tags() {
        let body = build_record(0, 0, b"r", 0, &[], 0, b"XXq\x00");
        let record = RawRecord::new(&body).expect("shape is valid");
        let error = record.validate_deeply(false).expect_err("must reject");
        assert!(matches!(error, BamRecordError::Tag { .. }), "{error}");
    }
}
