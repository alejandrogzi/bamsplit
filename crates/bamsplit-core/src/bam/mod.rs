// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The BAM format layer: headers, raw records, CIGARs, and auxiliary tags.
//!
//! Nothing here opens files or makes routing decisions. The layer's job is to
//! turn bytes into borrowed, bounds-checked views as cheaply as possible, and
//! to reject malformed input with an error that names the offending record.

pub mod cigar;
pub mod header;
pub mod raw_record;
pub mod tags;
pub mod validation;

pub use cigar::{AlignmentGeometry, CigarOps, DeletionPolicy, Kind, Op, ReferenceSpan};
pub use header::{
    BamHeader, MAGIC_NUMBER, ReadGroupField, ReadGroupInfo, ReferenceSequence, SortOrder,
};
pub use raw_record::{
    FIXED_CORE_SIZE, RawRecord, RawRecordReader, RecordLocation, VirtualPositionSource,
};
pub use tags::{ArraySubtype, RawTagArray, RawTagValue, Tag, TagReader};

/// SAM flag bits, named so call sites read as prose.
pub mod flags {
    /// `0x1`: the template has multiple segments.
    pub const SEGMENTED: u16 = 0x1;
    /// `0x2`: every segment is properly aligned.
    pub const PROPERLY_SEGMENTED: u16 = 0x2;
    /// `0x4`: the segment is unmapped.
    pub const UNMAPPED: u16 = 0x4;
    /// `0x8`: the next segment is unmapped.
    pub const MATE_UNMAPPED: u16 = 0x8;
    /// `0x10`: the segment is reverse-complemented.
    pub const REVERSE_COMPLEMENTED: u16 = 0x10;
    /// `0x20`: the next segment is reverse-complemented.
    pub const MATE_REVERSE_COMPLEMENTED: u16 = 0x20;
    /// `0x40`: the first segment of the template.
    pub const FIRST_SEGMENT: u16 = 0x40;
    /// `0x80`: the last segment of the template.
    pub const LAST_SEGMENT: u16 = 0x80;
    /// `0x100`: a secondary alignment.
    pub const SECONDARY: u16 = 0x100;
    /// `0x200`: the record failed quality control.
    pub const QC_FAIL: u16 = 0x200;
    /// `0x400`: an optical or PCR duplicate.
    pub const DUPLICATE: u16 = 0x400;
    /// `0x800`: a supplementary alignment.
    pub const SUPPLEMENTARY: u16 = 0x800;
}
