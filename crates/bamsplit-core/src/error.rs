// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Typed error hierarchy for `bamsplit`.
//!
//! Every fallible operation in [`bamsplit_core`](crate) returns a typed error.
//! The crate never panics because of malformed input, invalid configuration,
//! integer overflow, or filesystem failure; see the crate-level documentation
//! for the (short) list of documented internal invariants that are allowed to
//! panic.
//!
//! Errors are layered:
//!
//! * leaf errors describe one subsystem ([`BamRecordError`], [`HeaderError`],
//!   [`RoutingError`], [`OutputError`], [`IndexError`], [`AnnotationError`],
//!   [`SpoolError`], [`ManifestError`], [`ConfigError`]);
//! * [`EngineError`] aggregates whatever an execution engine can hit;
//! * [`Error`] is the crate-wide type returned by the high-level entry points
//!   and carries an [`ExitCode`] so the CLI never has to re-classify failures.

use std::path::PathBuf;

use crate::bam::raw_record::RecordLocation;

/// Process exit codes, as documented in the `bamsplit` manual.
///
/// The numeric values are part of the public interface: shell pipelines and
/// workflow managers branch on them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum ExitCode {
    /// Everything completed successfully.
    Success = 0,
    /// A general runtime failure that does not fit a more specific code.
    Runtime = 1,
    /// The invocation itself was invalid (bad flag, bad value, bad combination).
    InvalidArguments = 2,
    /// The input BAM, index, or annotation is invalid or malformed.
    InvalidInput = 3,
    /// An output path already exists, or two logical keys collide on disk.
    OutputConflict = 4,
    /// A post-condition failed: record conservation, digest, or index checks.
    ValidationFailure = 5,
    /// Execution was interrupted (SIGINT/SIGTERM) and rolled back.
    Interrupted = 6,
}

impl ExitCode {
    /// Returns the numeric exit status.
    #[must_use]
    pub const fn code(self) -> i32 {
        self as i32
    }
}

impl std::fmt::Display for ExitCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.code())
    }
}

/// Anything that can classify itself into an [`ExitCode`].
pub trait Classify {
    /// The exit code this error should terminate the process with.
    fn exit_code(&self) -> ExitCode;
}

// ---------------------------------------------------------------------------
// Raw BAM record layer
// ---------------------------------------------------------------------------

/// Failures raised while parsing the binary body of a BAM record.
///
/// Every variant carries enough context ([`RecordLocation`]) to point a user at
/// the offending record: its ordinal in the stream and the BGZF virtual offset
/// it started at.
#[derive(Debug, thiserror::Error)]
pub enum BamRecordError {
    /// The 4-byte little-endian `block_size` prefix was negative.
    #[error("negative BAM record block_size {block_size} at {location}")]
    NegativeBlockSize {
        /// The decoded (signed) block size.
        block_size: i32,
        /// Where the record began.
        location: RecordLocation,
    },

    /// `block_size` was below the 32-byte fixed core.
    #[error(
        "BAM record block_size {block_size} is smaller than the {minimum}-byte fixed core at {location}"
    )]
    BlockSizeTooSmall {
        /// The decoded block size.
        block_size: usize,
        /// The smallest legal block size.
        minimum: usize,
        /// Where the record began.
        location: RecordLocation,
    },

    /// `block_size` exceeded the configured safety limit.
    #[error(
        "BAM record block_size {block_size} exceeds the safety limit of {limit} bytes at {location}; \
         raise it with `RawRecordReader::with_max_record_size` if this input is legitimate"
    )]
    BlockSizeTooLarge {
        /// The decoded block size.
        block_size: usize,
        /// The configured limit.
        limit: usize,
        /// Where the record began.
        location: RecordLocation,
    },

    /// The stream ended in the middle of a record.
    #[error(
        "truncated BAM record at {location}: expected {expected} bytes of record body, found {actual}"
    )]
    TruncatedRecord {
        /// How many body bytes the header promised.
        expected: usize,
        /// How many were actually available.
        actual: usize,
        /// Where the record began.
        location: RecordLocation,
    },

    /// The stream ended in the middle of the 4-byte `block_size` prefix.
    #[error("truncated BAM record length prefix at {location}: only {actual} of 4 bytes available")]
    TruncatedBlockSize {
        /// How many prefix bytes were available.
        actual: usize,
        /// Where the record began.
        location: RecordLocation,
    },

    /// A variable-length section declared a size that does not fit the record.
    #[error(
        "BAM record field `{field}` overruns the record body at {location}: \
         needs {needed} bytes from offset {offset}, body is {body_len} bytes"
    )]
    FieldOutOfBounds {
        /// Which logical field overran.
        field: &'static str,
        /// Offset the field starts at, relative to the record body.
        offset: usize,
        /// How many bytes the field claimed.
        needed: usize,
        /// Total body length.
        body_len: usize,
        /// Where the record began.
        location: RecordLocation,
    },

    /// An arithmetic step in the layout computation overflowed.
    #[error("integer overflow computing BAM record layout (`{context}`) at {location}")]
    LayoutOverflow {
        /// Which computation overflowed.
        context: &'static str,
        /// Where the record began.
        location: RecordLocation,
    },

    /// `l_read_name` was zero, so the NUL terminator cannot exist.
    #[error("BAM record has zero-length read name at {location}")]
    EmptyReadName {
        /// Where the record began.
        location: RecordLocation,
    },

    /// The read name was not NUL-terminated.
    #[error("BAM record read name is not NUL-terminated at {location}")]
    UnterminatedReadName {
        /// Where the record began.
        location: RecordLocation,
    },

    /// `l_seq` was negative.
    #[error("BAM record has negative sequence length {l_seq} at {location}")]
    NegativeSequenceLength {
        /// The decoded value.
        l_seq: i32,
        /// Where the record began.
        location: RecordLocation,
    },

    /// A reference identifier was neither `-1` nor a valid dictionary index.
    #[error(
        "BAM record references {kind} sequence id {id}, but the header declares {reference_count} \
         reference sequences, at {location}"
    )]
    InvalidReferenceSequenceId {
        /// `"reference"` or `"mate reference"`.
        kind: &'static str,
        /// The offending identifier.
        id: i32,
        /// How many references the header declares.
        reference_count: usize,
        /// Where the record began.
        location: RecordLocation,
    },

    /// An alignment start position was below `-1`.
    #[error("BAM record has invalid {kind} position {position} at {location}")]
    InvalidPosition {
        /// `"alignment"` or `"mate alignment"`.
        kind: &'static str,
        /// The offending 0-based position.
        position: i32,
        /// Where the record began.
        location: RecordLocation,
    },

    /// An auxiliary data field could not be parsed.
    #[error("malformed auxiliary data at {location}: {source}")]
    Tag {
        /// The underlying tag error.
        #[source]
        source: Box<TagError>,
        /// Where the record began.
        location: RecordLocation,
    },

    /// A CIGAR string could not be interpreted.
    #[error("malformed CIGAR at {location}: {source}")]
    Cigar {
        /// The underlying CIGAR error.
        #[source]
        source: Box<CigarError>,
        /// Where the record began.
        location: RecordLocation,
    },

    /// An I/O failure while filling the record buffer.
    #[error("i/o error reading BAM record at {location}")]
    Io {
        /// The underlying error.
        #[source]
        source: std::io::Error,
        /// Where the record began.
        location: RecordLocation,
    },
}

impl BamRecordError {
    /// Returns the location the failure was detected at.
    #[must_use]
    pub fn location(&self) -> RecordLocation {
        match self {
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
            | Self::Io { location, .. } => *location,
        }
    }
}

/// Failures raised while decoding BAM auxiliary (tag) data.
#[derive(Debug, thiserror::Error)]
pub enum TagError {
    /// The auxiliary section ended in the middle of a field.
    #[error("truncated auxiliary field `{tag}`: needed {needed} bytes, {available} remain")]
    Truncated {
        /// The two-character tag, rendered for humans.
        tag: String,
        /// Bytes required.
        needed: usize,
        /// Bytes available.
        available: usize,
    },

    /// The auxiliary section ended before a complete tag/type prefix.
    #[error("truncated auxiliary field header at offset {offset}: {available} bytes remain")]
    TruncatedHeader {
        /// Offset into the auxiliary section.
        offset: usize,
        /// Bytes available.
        available: usize,
    },

    /// The one-byte value type code is not one of `AcCsSiIfZHB`.
    #[error("unknown auxiliary value type {type_code:?} (0x{byte:02x}) for tag `{tag}`")]
    UnknownValueType {
        /// The tag whose value could not be typed.
        tag: String,
        /// The type code rendered as a character, when printable.
        type_code: char,
        /// The raw byte.
        byte: u8,
    },

    /// A `B` array declared a subtype that is not one of `cCsSiIf`.
    #[error("unknown auxiliary array subtype {subtype:?} (0x{byte:02x}) for tag `{tag}`")]
    UnknownArraySubtype {
        /// The tag whose array could not be typed.
        tag: String,
        /// The subtype rendered as a character, when printable.
        subtype: char,
        /// The raw byte.
        byte: u8,
    },

    /// `count * element_size` overflowed, or the payload does not fit.
    #[error(
        "auxiliary array `{tag}` declares {count} elements of {element_size} bytes, \
         which overflows or exceeds the {available} remaining bytes"
    )]
    ArrayLengthOverflow {
        /// The tag.
        tag: String,
        /// Declared element count.
        count: u32,
        /// Size of one element in bytes.
        element_size: usize,
        /// Bytes available.
        available: usize,
    },

    /// A `Z` or `H` string was not NUL-terminated before the section ended.
    #[error("unterminated `{type_code}` auxiliary string for tag `{tag}`")]
    UnterminatedString {
        /// The tag.
        tag: String,
        /// `'Z'` or `'H'`.
        type_code: char,
    },

    /// A `H` (hex) value had an odd length or a non-hex digit.
    #[error("malformed hex auxiliary value for tag `{tag}`: {reason}")]
    MalformedHex {
        /// The tag.
        tag: String,
        /// Why it was rejected.
        reason: &'static str,
    },

    /// The same tag appeared twice in one record.
    #[error("duplicate auxiliary tag `{tag}` in one record")]
    DuplicateTag {
        /// The tag.
        tag: String,
    },
}

/// Failures raised while interpreting a BAM CIGAR.
#[derive(Debug, thiserror::Error)]
pub enum CigarError {
    /// A CIGAR operation code was outside `0..=8`.
    #[error("unknown CIGAR operation code {code}")]
    UnknownOperation {
        /// The 4-bit operation code.
        code: u32,
    },

    /// The reference span computation overflowed `i64`.
    #[error("integer overflow computing CIGAR reference span")]
    SpanOverflow,

    /// The `CG` long-CIGAR placeholder was present but the `CG:B:I` tag was not.
    #[error("record uses the long-CIGAR `CG` convention but has no `CG:B:I` tag")]
    MissingLongCigarTag,

    /// The `CG:B:I` tag was present but its contents were unusable.
    #[error("malformed long-CIGAR `CG:B:I` tag: {reason}")]
    MalformedLongCigar {
        /// Why it was rejected.
        reason: &'static str,
    },
}

// ---------------------------------------------------------------------------
// Header layer
// ---------------------------------------------------------------------------

/// Failures raised while reading or validating a BAM header.
#[derive(Debug, thiserror::Error)]
pub enum HeaderError {
    /// The file did not start with `BAM\1`.
    #[error("not a BAM file: expected magic {expected:02x?}, found {found:02x?}")]
    BadMagic {
        /// `b"BAM\x01"`.
        expected: [u8; 4],
        /// What was actually read.
        found: [u8; 4],
    },

    /// The stream ended before the header was complete.
    #[error("truncated BAM header while reading {section}")]
    Truncated {
        /// Which part of the header was being read.
        section: &'static str,
    },

    /// A declared length was negative or absurd.
    #[error("BAM header declares an invalid {field} of {value}")]
    InvalidLength {
        /// Which length field.
        field: &'static str,
        /// The decoded value.
        value: i64,
    },

    /// The textual SAM header could not be parsed.
    #[error("malformed SAM header text: {message}")]
    MalformedText {
        /// Parser diagnostics.
        message: String,
    },

    /// The binary reference dictionary and the `@SQ` lines disagree.
    #[error(
        "BAM header inconsistency: binary dictionary has {binary} reference sequences \
         but the text header declares {text} `@SQ` lines"
    )]
    ReferenceCountMismatch {
        /// Count from the binary dictionary.
        binary: usize,
        /// Count from the `@SQ` lines.
        text: usize,
    },

    /// A binary reference name does not match the corresponding `@SQ` `SN`.
    #[error(
        "BAM header inconsistency at reference {index}: binary name {binary:?} \
         does not match `@SQ SN:{text:?}`"
    )]
    ReferenceNameMismatch {
        /// Dictionary index.
        index: usize,
        /// Name from the binary dictionary.
        binary: String,
        /// Name from the text header.
        text: String,
    },

    /// A binary reference length does not match the corresponding `@SQ` `LN`.
    #[error(
        "BAM header inconsistency at reference {index} ({name:?}): binary length {binary} \
         does not match `@SQ LN:{text}`"
    )]
    ReferenceLengthMismatch {
        /// Dictionary index.
        index: usize,
        /// Reference name.
        name: String,
        /// Length from the binary dictionary.
        binary: u32,
        /// Length from the text header.
        text: u32,
    },

    /// Two references share a name.
    #[error("duplicate reference sequence name {name:?} at indices {first} and {second}")]
    DuplicateReferenceName {
        /// The repeated name.
        name: String,
        /// First occurrence.
        first: usize,
        /// Second occurrence.
        second: usize,
    },

    /// A reference name was empty or contained a forbidden byte.
    #[error("invalid reference sequence name at index {index}: {reason}")]
    InvalidReferenceName {
        /// Dictionary index.
        index: usize,
        /// Why it was rejected.
        reason: String,
    },

    /// A reference declared a zero or negative length.
    #[error("reference sequence {name:?} at index {index} declares an invalid length {length}")]
    InvalidReferenceLength {
        /// Reference name.
        name: String,
        /// Dictionary index.
        index: usize,
        /// The declared length.
        length: i64,
    },

    /// Two `@RG` lines share an `ID`.
    #[error("duplicate read group ID {id:?}")]
    DuplicateReadGroupId {
        /// The repeated identifier.
        id: String,
    },

    /// Two `@PG` lines share an `ID`.
    #[error("duplicate program record ID {id:?}")]
    DuplicateProgramId {
        /// The repeated identifier.
        id: String,
    },

    /// A `@PG` `PP` pointed at an identifier that does not exist.
    #[error("program record {id:?} declares a previous program {previous:?} that does not exist")]
    DanglingProgramLink {
        /// The referring program.
        id: String,
        /// The missing target.
        previous: String,
    },

    /// The `@HD SO` value is not one of the SAM-specified sort orders.
    #[error("unrecognized `@HD SO:{value}` sort order")]
    UnknownSortOrder {
        /// The raw value.
        value: String,
    },

    /// Re-serializing the header failed.
    #[error("failed to serialize BAM header")]
    Serialize {
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },

    /// An I/O failure while reading the header.
    #[error("i/o error reading BAM header")]
    Io {
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------------
// Routing layer
// ---------------------------------------------------------------------------

/// Failures raised by a [`Router`](crate::routing::Router).
#[derive(Debug, thiserror::Error)]
pub enum RoutingError {
    /// The record could not be parsed well enough to route it.
    #[error(transparent)]
    Record(#[from] Box<BamRecordError>),

    /// A record carried a reference id outside the header dictionary.
    #[error("record at {location} has reference sequence id {id}, outside 0..{reference_count}")]
    InvalidReferenceId {
        /// The offending identifier.
        id: i32,
        /// Size of the dictionary.
        reference_count: usize,
        /// Where the record began.
        location: RecordLocation,
    },

    /// A record has no `RG` tag and the policy is `error`.
    #[error("record at {location} has no `RG` tag (policy: error)")]
    MissingReadGroup {
        /// Where the record began.
        location: RecordLocation,
    },

    /// A record's `RG` is not declared in the header and the policy is `error`.
    #[error("record at {location} references unknown read group {id:?} (policy: error)")]
    UnknownReadGroup {
        /// The unknown read-group identifier.
        id: String,
        /// Where the record began.
        location: RecordLocation,
    },

    /// The requested `@RG` sub-field is absent from an otherwise valid group.
    #[error("read group {id:?} has no `{field}` field (policy: error)")]
    MissingReadGroupField {
        /// The read group.
        id: String,
        /// The missing SAM tag, e.g. `SM`.
        field: &'static str,
    },

    /// A record lacks the routing tag and the policy is `error`.
    #[error("record at {location} has no `{tag}` tag (policy: error)")]
    MissingTag {
        /// The requested tag.
        tag: String,
        /// Where the record began.
        location: RecordLocation,
    },

    /// The routing tag held a `B` array and array tags were not enabled.
    #[error(
        "tag `{tag}` at {location} holds a `B` array; array-valued tags are rejected by default \
         because their filename and equality semantics are surprising — pass `--allow-array-tags` \
         to opt in"
    )]
    ArrayValuedTag {
        /// The requested tag.
        tag: String,
        /// Where the record began.
        location: RecordLocation,
    },

    /// An unplaced record was seen and the policy is `error`.
    #[error("record at {location} is unplaced (refID == -1) and `--unplaced error` was requested")]
    UnplacedRecord {
        /// Where the record began.
        location: RecordLocation,
    },

    /// A record did not match any annotation region and the policy is `error`.
    #[error("record at {location} matched no annotation region (policy: error)")]
    UnmatchedRecord {
        /// Where the record began.
        location: RecordLocation,
    },

    /// The number of distinct output keys exceeded `--max-outputs`.
    #[error(
        "output cardinality limit exceeded: {discovered} distinct keys discovered but \
         `--max-outputs` is {limit}. For high-cardinality keys prefer \
         `bamsplit shard --key tag:{hint}`, or raise the limit with `--allow-high-cardinality`"
    )]
    CardinalityExceeded {
        /// How many keys had been discovered when the limit tripped.
        discovered: usize,
        /// The configured limit.
        limit: usize,
        /// A tag name to suggest in the hint.
        hint: String,
    },

    /// A CIGAR could not be interpreted while computing region geometry.
    #[error("cannot compute alignment geometry for record at {location}")]
    Geometry {
        /// The underlying CIGAR error.
        #[source]
        source: Box<CigarError>,
        /// Where the record began.
        location: RecordLocation,
    },
}

// ---------------------------------------------------------------------------
// Output layer
// ---------------------------------------------------------------------------

/// Failures raised by the transactional output manager.
#[derive(Debug, thiserror::Error)]
pub enum OutputError {
    /// The final path already exists and `--force` was not given.
    #[error("output {path} already exists; pass `--force` to replace it")]
    AlreadyExists {
        /// The conflicting path.
        path: PathBuf,
    },

    /// Two distinct logical keys encoded to the same filename.
    #[error(
        "output filename collision: logical keys {first:?} and {second:?} both encode to {encoded:?}"
    )]
    EncodedKeyCollision {
        /// First logical key.
        first: String,
        /// Second logical key.
        second: String,
        /// The shared encoded name.
        encoded: String,
    },

    /// A logical key produced an empty encoded name.
    #[error("logical key {key:?} encodes to an empty filename")]
    EmptyEncodedKey {
        /// The offending key.
        key: String,
    },

    /// The output directory could not be created or is not a directory.
    #[error("cannot use output directory {path}: {reason}")]
    BadOutputDirectory {
        /// The directory.
        path: PathBuf,
        /// Why it is unusable.
        reason: String,
    },

    /// A filename template referenced an unknown placeholder.
    #[error("unknown placeholder `{{{placeholder}}}` in filename template {template:?}")]
    UnknownTemplatePlaceholder {
        /// The unknown placeholder name.
        placeholder: String,
        /// The template as given.
        template: String,
    },

    /// A filename template was structurally invalid.
    #[error("invalid filename template {template:?}: {reason}")]
    InvalidTemplate {
        /// The template as given.
        template: String,
        /// Why it was rejected.
        reason: String,
    },

    /// The final rename from the temporary path failed.
    #[error("failed to commit {temporary} to {final_path}")]
    Commit {
        /// The temporary path.
        temporary: PathBuf,
        /// The intended final path.
        final_path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },

    /// A post-write validation check failed.
    #[error("output validation failed for {path}: {reason}")]
    Validation {
        /// The output.
        path: PathBuf,
        /// What was wrong.
        reason: String,
    },

    /// Too many outputs would have to be open simultaneously.
    #[error(
        "the requested split needs more than {limit} concurrently open outputs; \
         raise `--max-open-files` or use `--engine spool`"
    )]
    TooManyOpenFiles {
        /// The configured limit.
        limit: usize,
    },

    /// An I/O failure attributable to a specific path.
    #[error("i/o error on {path}")]
    Io {
        /// The path involved.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------------
// Index layer
// ---------------------------------------------------------------------------

/// Failures raised while building or writing an output index.
#[derive(Debug, thiserror::Error)]
pub enum IndexError {
    /// BAI was forced but a coordinate exceeds what BAI can address.
    #[error(
        "cannot build a BAI index for {path}: reference {reference:?} needs position {position}, \
         but BAI addresses at most {limit}. Use `--index csi` or `--index auto`"
    )]
    BaiLimitExceeded {
        /// The output being indexed.
        path: PathBuf,
        /// The reference that overflowed.
        reference: String,
        /// The offending 1-based position.
        position: u64,
        /// The BAI addressable maximum.
        limit: u64,
    },

    /// An index was requested for output that is not coordinate-sorted.
    #[error("cannot index {path}: {reason}")]
    NotIndexable {
        /// The output.
        path: PathBuf,
        /// Why it cannot be indexed.
        reason: String,
    },

    /// A virtual offset pair was non-monotonic.
    #[error("non-monotonic virtual offsets while indexing {path}: chunk start {start} > end {end}")]
    NonMonotonicChunk {
        /// The output.
        path: PathBuf,
        /// Chunk start.
        start: u64,
        /// Chunk end.
        end: u64,
    },

    /// The input index could not be read.
    #[error("failed to read index {path}")]
    ReadInput {
        /// The index path.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },

    /// The input index does not describe the input BAM.
    #[error(
        "index {path} describes {index_references} reference sequences but the BAM header \
         declares {header_references}"
    )]
    IndexHeaderMismatch {
        /// The index path.
        path: PathBuf,
        /// References in the index.
        index_references: usize,
        /// References in the header.
        header_references: usize,
    },

    /// An I/O failure while writing the index.
    #[error("i/o error writing index {path}")]
    Io {
        /// The index path.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------------
// Annotation layer
// ---------------------------------------------------------------------------

/// Failures raised by the `genepred` annotation adapter.
#[derive(Debug, thiserror::Error)]
pub enum AnnotationError {
    /// Neither the extension nor the contents identified the format.
    #[error(
        "cannot determine the annotation format of {path}: {reason}. \
         Pass `--format bed|gtf|gff` to disambiguate"
    )]
    AmbiguousFormat {
        /// The annotation path.
        path: PathBuf,
        /// What the detector saw.
        reason: String,
    },

    /// The BED column count could not be resolved to a supported width.
    #[error(
        "cannot determine the BED width of {path}: {reason}. \
         Pass `--type 3|4|5|6|8|9|12` to disambiguate"
    )]
    AmbiguousBedType {
        /// The annotation path.
        path: PathBuf,
        /// What the detector saw.
        reason: String,
    },

    /// A BED width was requested or detected that `bamsplit` does not support.
    #[error("unsupported BED width {width}; supported widths are 3, 4, 5, 6, 8, 9, and 12")]
    UnsupportedBedType {
        /// The offending width.
        width: usize,
    },

    /// Rows in one BED file disagreed about the column count.
    #[error("mixed BED widths in {path}: row 1 has {expected} columns but row {line} has {actual}")]
    MixedBedWidths {
        /// The annotation path.
        path: PathBuf,
        /// Width chosen from the first data row.
        expected: usize,
        /// Width of the offending row.
        actual: usize,
        /// 1-based line number of the offending row.
        line: usize,
    },

    /// `--type` was given for a non-BED input.
    #[error("`--type` applies to BED input only, but {path} was detected as {format}")]
    BedTypeOnNonBed {
        /// The annotation path.
        path: PathBuf,
        /// The detected format.
        format: &'static str,
    },

    /// The annotation file contained no usable records.
    #[error("annotation {path} contains no records")]
    Empty {
        /// The annotation path.
        path: PathBuf,
    },

    /// A record was rejected by the `genepred` reader.
    #[error("failed to parse annotation {path}: {message}")]
    Parse {
        /// The annotation path.
        path: PathBuf,
        /// The `genepred` diagnostic.
        message: String,
    },

    /// A record's coordinates are not a valid half-open interval.
    #[error(
        "annotation record {name:?} on {chrom} at ordinal {ordinal} has invalid coordinates \
         {start}..{end}: {reason}"
    )]
    InvalidCoordinates {
        /// Logical name, if any.
        name: String,
        /// Reference name.
        chrom: String,
        /// 0-based input ordinal.
        ordinal: u64,
        /// Start.
        start: u64,
        /// End.
        end: u64,
        /// Why it was rejected.
        reason: &'static str,
    },

    /// Block structure was required for the requested feature but is absent.
    #[error(
        "annotation record {name:?} (ordinal {ordinal}) has no exon/block structure, which \
         `--feature {feature}` requires{hint}"
    )]
    MissingBlocks {
        /// Logical name, if any.
        name: String,
        /// 0-based input ordinal.
        ordinal: u64,
        /// The requested feature.
        feature: &'static str,
        /// A trailing hint, possibly empty.
        hint: &'static str,
    },

    /// Blocks were unsorted, overlapping, or zero-length.
    #[error("annotation record {name:?} (ordinal {ordinal}) has invalid block structure: {reason}")]
    InvalidBlocks {
        /// Logical name, if any.
        name: String,
        /// 0-based input ordinal.
        ordinal: u64,
        /// Why it was rejected.
        reason: String,
    },

    /// Coding boundaries were required but absent or nonsensical.
    #[error(
        "annotation record {name:?} (ordinal {ordinal}) has no usable coding boundaries, which \
         `--feature {feature}` requires; use `--missing-feature skip|empty` to tolerate this"
    )]
    MissingCoding {
        /// Logical name, if any.
        name: String,
        /// 0-based input ordinal.
        ordinal: u64,
        /// The requested feature.
        feature: &'static str,
    },

    /// Strand was required for a strand-aware feature but absent.
    #[error(
        "annotation record {name:?} (ordinal {ordinal}) has no strand, which `--feature {feature}` \
         requires; use `--missing-feature skip|empty` to tolerate this"
    )]
    MissingStrand {
        /// Logical name, if any.
        name: String,
        /// 0-based input ordinal.
        ordinal: u64,
        /// The requested feature.
        feature: &'static str,
    },

    /// A `--window-size` value could not be parsed or was zero.
    #[error("invalid window size {value:?}: {reason}")]
    InvalidWindowSize {
        /// The raw value.
        value: String,
        /// Why it was rejected.
        reason: &'static str,
    },

    /// No annotation reference name matched any BAM reference name.
    #[error(
        "no annotation reference name matches any BAM reference sequence; \
         the annotation uses names like {annotation_examples} and the BAM uses {bam_examples}"
    )]
    NoReferenceOverlap {
        /// A few names from the annotation.
        annotation_examples: String,
        /// A few names from the BAM.
        bam_examples: String,
    },

    /// An I/O failure while reading the annotation.
    #[error("i/o error reading annotation {path}")]
    Io {
        /// The annotation path.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------------
// Spool layer
// ---------------------------------------------------------------------------

/// Failures raised by the spool-and-finalize engine's temporary storage.
#[derive(Debug, thiserror::Error)]
pub enum SpoolError {
    /// The spool file did not start with the expected magic.
    #[error("spool file {path} has bad magic {found:02x?}, expected {expected:02x?}")]
    BadMagic {
        /// The spool path.
        path: PathBuf,
        /// `b"BSPL"`.
        expected: [u8; 4],
        /// What was read.
        found: [u8; 4],
    },

    /// The spool format version is not understood.
    #[error(
        "spool file {path} declares format version {version}, but only {supported} is supported"
    )]
    UnsupportedVersion {
        /// The spool path.
        path: PathBuf,
        /// Version read from the file.
        version: u16,
        /// Version this build writes.
        supported: u16,
    },

    /// The spool file ended before its footer.
    #[error("spool file {path} is truncated: {reason}")]
    Truncated {
        /// The spool path.
        path: PathBuf,
        /// What was missing.
        reason: String,
    },

    /// The footer checksum did not match the payload.
    #[error(
        "spool file {path} failed its checksum: expected {expected:#018x}, computed {actual:#018x}"
    )]
    ChecksumMismatch {
        /// The spool path.
        path: PathBuf,
        /// The stored checksum.
        expected: u64,
        /// The recomputed checksum.
        actual: u64,
    },

    /// The footer record count did not match the number of records read back.
    #[error("spool file {path} declares {expected} records but contains {actual}")]
    RecordCountMismatch {
        /// The spool path.
        path: PathBuf,
        /// The stored count.
        expected: u64,
        /// The count observed on replay.
        actual: u64,
    },

    /// A length prefix was implausible.
    #[error(
        "spool file {path} declares a {length}-byte record at offset {offset}, exceeding the {limit}-byte limit"
    )]
    RecordTooLarge {
        /// The spool path.
        path: PathBuf,
        /// The declared length.
        length: u64,
        /// Byte offset of the prefix.
        offset: u64,
        /// The configured limit.
        limit: u64,
    },

    /// The temporary directory could not be created or is unusable.
    #[error("cannot use temporary directory {path}: {reason}")]
    BadTempDirectory {
        /// The directory.
        path: PathBuf,
        /// Why it is unusable.
        reason: String,
    },

    /// An I/O failure attributable to a spool path.
    #[error("i/o error on spool file {path}")]
    Io {
        /// The spool path.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------------
// Manifest layer
// ---------------------------------------------------------------------------

/// Failures raised while producing or validating the run manifest.
#[derive(Debug, thiserror::Error)]
pub enum ManifestError {
    /// Record conservation did not hold for a non-duplicating routing mode.
    #[error(
        "record conservation failed: {input_records} input records != \
         {unique_emitted} emitted + {dropped} dropped + {unmatched} unmatched \
         (difference {difference})"
    )]
    ConservationFailed {
        /// Records read from the input.
        input_records: u64,
        /// Records emitted at least once.
        unique_emitted: u64,
        /// Records deliberately dropped.
        dropped: u64,
        /// Records that matched nothing.
        unmatched: u64,
        /// Signed difference, for quick diagnosis.
        difference: i128,
    },

    /// Emission accounting did not hold for `overlap` mode.
    #[error(
        "emission conservation failed: {total_emissions} total emissions != \
         {unique_emitted} unique + {duplicates} duplicates (difference {difference})"
    )]
    EmissionConservationFailed {
        /// Sum of per-output record counts.
        total_emissions: u64,
        /// Records emitted at least once.
        unique_emitted: u64,
        /// Extra copies beyond the first.
        duplicates: u64,
        /// Signed difference.
        difference: i128,
    },

    /// A per-output record count disagreed with the sum of its categories.
    #[error(
        "per-output accounting failed for {key:?}: record_count {record_count} != \
         mapped {mapped} + placed-unmapped {placed_unmapped} + unplaced-unmapped {unplaced_unmapped}"
    )]
    PerOutputAccountingFailed {
        /// The logical key.
        key: String,
        /// Declared total.
        record_count: u64,
        /// Mapped records.
        mapped: u64,
        /// Placed but unmapped records.
        placed_unmapped: u64,
        /// Unplaced records.
        unplaced_unmapped: u64,
    },

    /// The manifest could not be serialized.
    #[error("failed to serialize the manifest as {format}")]
    Serialize {
        /// `"json"` or `"tsv"`.
        format: &'static str,
        /// The underlying error.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    /// An I/O failure while writing the manifest.
    #[error("i/o error writing manifest {path}")]
    Io {
        /// The manifest path.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------------
// Configuration layer
// ---------------------------------------------------------------------------

/// Failures raised while validating a run configuration.
///
/// These are surfaced with [`ExitCode::InvalidArguments`]. The CLI performs its
/// own `clap`-level validation first; this type covers the semantic checks that
/// belong to the library, so the Rust API rejects the same nonsense the CLI
/// does.
#[derive(Debug, thiserror::Error)]
pub enum ConfigError {
    /// A numeric option was outside its permitted range.
    #[error("`{option}` must be {constraint}, got {value}")]
    OutOfRange {
        /// The option name as the user writes it.
        option: &'static str,
        /// A human description of the valid range.
        constraint: &'static str,
        /// The offending value.
        value: String,
    },

    /// Two options cannot be combined.
    #[error("`{first}` cannot be combined with `{second}`: {reason}")]
    Conflict {
        /// First option.
        first: &'static str,
        /// Second option.
        second: &'static str,
        /// Why they conflict.
        reason: &'static str,
    },

    /// Exactly one of a set of options is required.
    #[error("exactly one of {options} is required")]
    ExactlyOneRequired {
        /// A rendered list such as ``"`--tag`, `--field`"``.
        options: &'static str,
    },

    /// A requested engine cannot run this job.
    #[error("`--engine {engine}` is not usable here: {reason}")]
    IncompatibleEngine {
        /// The requested engine.
        engine: &'static str,
        /// Why it cannot be used.
        reason: String,
    },

    /// A requested I/O backend cannot be used.
    #[error("`--io {backend}` is not usable here: {reason}")]
    IncompatibleIoBackend {
        /// The requested backend.
        backend: &'static str,
        /// Why it cannot be used.
        reason: String,
    },

    /// A reference name requested via `--include`/`--exclude` does not exist.
    #[error(
        "unknown reference sequence name{plural} {names} requested via `{option}`; \
         pass `--ignore-missing-references` to skip them instead"
    )]
    UnknownReferenceNames {
        /// `""` or `"s"`.
        plural: &'static str,
        /// A rendered, quoted list.
        names: String,
        /// Which option they came from.
        option: &'static str,
    },

    /// `--include` and `--exclude` left nothing to write.
    #[error("the reference selection is empty: {reason}")]
    EmptySelection {
        /// Why nothing was selected.
        reason: String,
    },

    /// A `--reference-list` file could not be read.
    #[error("cannot read reference list {path}")]
    ReferenceList {
        /// The list path.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },

    /// The input path is not usable.
    #[error("cannot use input {path}: {reason}")]
    BadInput {
        /// The input path.
        path: PathBuf,
        /// Why it is unusable.
        reason: String,
    },

    /// A tag name was not exactly two alphanumeric characters.
    #[error("invalid tag {tag:?}: a BAM tag is exactly two characters, `[A-Za-z][A-Za-z0-9]`")]
    InvalidTag {
        /// The offending value.
        tag: String,
    },
}

// ---------------------------------------------------------------------------
// Engine layer
// ---------------------------------------------------------------------------

/// Failures raised by an execution engine.
#[derive(Debug, thiserror::Error)]
pub enum EngineError {
    /// A record could not be parsed.
    #[error(transparent)]
    Record(#[from] Box<BamRecordError>),

    /// The header could not be read, validated, or rewritten.
    #[error(transparent)]
    Header(#[from] Box<HeaderError>),

    /// Routing failed.
    #[error(transparent)]
    Routing(#[from] Box<RoutingError>),

    /// Writing or committing an output failed.
    #[error(transparent)]
    Output(#[from] Box<OutputError>),

    /// Index construction failed.
    #[error(transparent)]
    Index(#[from] Box<IndexError>),

    /// Spool storage failed.
    #[error(transparent)]
    Spool(#[from] Box<SpoolError>),

    /// The configuration was rejected.
    #[error(transparent)]
    Config(#[from] Box<ConfigError>),

    /// The annotation could not be loaded.
    #[error(transparent)]
    Annotation(#[from] Box<AnnotationError>),

    /// A post-run invariant failed.
    #[error(transparent)]
    Manifest(#[from] Box<ManifestError>),

    /// The stream engine detected that the input is not actually sorted.
    #[error(
        "input is not coordinate-sorted: record {record_number} at ({reference_id}, {position}) \
         follows ({previous_reference_id}, {previous_position})"
    )]
    SortOrderViolation {
        /// Ordinal of the offending record.
        record_number: u64,
        /// Its reference id.
        reference_id: i32,
        /// Its 0-based position.
        position: i32,
        /// The previous record's reference id.
        previous_reference_id: i32,
        /// The previous record's 0-based position.
        previous_position: i32,
    },

    /// The stream engine saw a routing key it had already finalized.
    #[error(
        "routing key {key:?} reappeared after being finalized at record {record_number}; \
         the input is not grouped by routing key. Use `--engine spool`"
    )]
    UngroupedInput {
        /// The key that came back.
        key: String,
        /// Ordinal of the offending record.
        record_number: u64,
    },

    /// The run was interrupted by a signal.
    #[error(
        "execution interrupted after {records_processed} records; partial outputs were removed"
    )]
    Interrupted {
        /// How far the run got.
        records_processed: u64,
    },

    /// A worker thread panicked or disappeared.
    #[error("worker thread failed: {reason}")]
    WorkerFailure {
        /// What happened.
        reason: String,
    },

    /// An I/O failure that is not attributable to one path.
    #[error("i/o error during execution")]
    Io {
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------------
// Crate-wide error
// ---------------------------------------------------------------------------

/// The crate-wide error type returned by the high-level entry points.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Configuration was rejected.
    #[error(transparent)]
    Config(#[from] ConfigError),

    /// The BAM header was invalid.
    #[error(transparent)]
    Header(#[from] HeaderError),

    /// A BAM record was invalid.
    #[error(transparent)]
    Record(#[from] BamRecordError),

    /// Routing failed.
    #[error(transparent)]
    Routing(#[from] RoutingError),

    /// An execution engine failed.
    #[error(transparent)]
    Engine(#[from] EngineError),

    /// Output management failed.
    #[error(transparent)]
    Output(#[from] OutputError),

    /// Index construction failed.
    #[error(transparent)]
    Index(#[from] IndexError),

    /// The annotation could not be loaded.
    #[error(transparent)]
    Annotation(#[from] AnnotationError),

    /// Spool storage failed.
    #[error(transparent)]
    Spool(#[from] SpoolError),

    /// The manifest failed to serialize or to validate.
    #[error(transparent)]
    Manifest(#[from] ManifestError),

    /// An I/O failure that is not attributable to a subsystem.
    #[error("i/o error")]
    Io(#[from] std::io::Error),
}

/// A convenient result alias for the crate-wide [`Error`].
pub type Result<T, E = Error> = std::result::Result<T, E>;

impl Classify for ConfigError {
    fn exit_code(&self) -> ExitCode {
        ExitCode::InvalidArguments
    }
}

impl Classify for HeaderError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::Io { .. } | Self::Serialize { .. } => ExitCode::Runtime,
            _ => ExitCode::InvalidInput,
        }
    }
}

impl Classify for BamRecordError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::Io { .. } => ExitCode::Runtime,
            _ => ExitCode::InvalidInput,
        }
    }
}

impl Classify for TagError {
    fn exit_code(&self) -> ExitCode {
        ExitCode::InvalidInput
    }
}

impl Classify for CigarError {
    fn exit_code(&self) -> ExitCode {
        ExitCode::InvalidInput
    }
}

impl Classify for RoutingError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::CardinalityExceeded { .. } => ExitCode::OutputConflict,
            Self::Record(inner) => inner.exit_code(),
            _ => ExitCode::InvalidInput,
        }
    }
}

impl Classify for OutputError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::AlreadyExists { .. }
            | Self::EncodedKeyCollision { .. }
            | Self::TooManyOpenFiles { .. } => ExitCode::OutputConflict,
            Self::Validation { .. } => ExitCode::ValidationFailure,
            Self::EmptyEncodedKey { .. }
            | Self::UnknownTemplatePlaceholder { .. }
            | Self::InvalidTemplate { .. } => ExitCode::InvalidArguments,
            Self::BadOutputDirectory { .. } | Self::Commit { .. } | Self::Io { .. } => {
                ExitCode::Runtime
            }
        }
    }
}

impl Classify for IndexError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::BaiLimitExceeded { .. } | Self::NotIndexable { .. } => ExitCode::InvalidArguments,
            Self::NonMonotonicChunk { .. } => ExitCode::ValidationFailure,
            Self::ReadInput { .. } | Self::IndexHeaderMismatch { .. } => ExitCode::InvalidInput,
            Self::Io { .. } => ExitCode::Runtime,
        }
    }
}

impl Classify for AnnotationError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::UnsupportedBedType { .. }
            | Self::BedTypeOnNonBed { .. }
            | Self::InvalidWindowSize { .. } => ExitCode::InvalidArguments,
            Self::Io { .. } => ExitCode::Runtime,
            _ => ExitCode::InvalidInput,
        }
    }
}

impl Classify for SpoolError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::BadTempDirectory { .. } | Self::Io { .. } => ExitCode::Runtime,
            _ => ExitCode::ValidationFailure,
        }
    }
}

impl Classify for ManifestError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::Serialize { .. } | Self::Io { .. } => ExitCode::Runtime,
            _ => ExitCode::ValidationFailure,
        }
    }
}

impl Classify for EngineError {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::Record(inner) => inner.exit_code(),
            Self::Header(inner) => inner.exit_code(),
            Self::Routing(inner) => inner.exit_code(),
            Self::Output(inner) => inner.exit_code(),
            Self::Index(inner) => inner.exit_code(),
            Self::Spool(inner) => inner.exit_code(),
            Self::Config(inner) => inner.exit_code(),
            Self::Annotation(inner) => inner.exit_code(),
            Self::Manifest(inner) => inner.exit_code(),
            Self::SortOrderViolation { .. } | Self::UngroupedInput { .. } => ExitCode::InvalidInput,
            Self::Interrupted { .. } => ExitCode::Interrupted,
            Self::WorkerFailure { .. } | Self::Io { .. } => ExitCode::Runtime,
        }
    }
}

impl Classify for Error {
    fn exit_code(&self) -> ExitCode {
        match self {
            Self::Config(inner) => inner.exit_code(),
            Self::Header(inner) => inner.exit_code(),
            Self::Record(inner) => inner.exit_code(),
            Self::Routing(inner) => inner.exit_code(),
            Self::Engine(inner) => inner.exit_code(),
            Self::Output(inner) => inner.exit_code(),
            Self::Index(inner) => inner.exit_code(),
            Self::Annotation(inner) => inner.exit_code(),
            Self::Spool(inner) => inner.exit_code(),
            Self::Manifest(inner) => inner.exit_code(),
            Self::Io(_) => ExitCode::Runtime,
        }
    }
}

macro_rules! boxed_from {
    ($($outer:ident :: $variant:ident <- $inner:ty),* $(,)?) => {
        $(
            impl From<$inner> for $outer {
                fn from(value: $inner) -> Self {
                    $outer::$variant(Box::new(value))
                }
            }
        )*
    };
}

boxed_from! {
    EngineError::Record <- BamRecordError,
    EngineError::Header <- HeaderError,
    EngineError::Routing <- RoutingError,
    EngineError::Output <- OutputError,
    EngineError::Index <- IndexError,
    EngineError::Spool <- SpoolError,
    EngineError::Config <- ConfigError,
    EngineError::Annotation <- AnnotationError,
    EngineError::Manifest <- ManifestError,
    RoutingError::Record <- BamRecordError,
}

/// Renders a chain of [`std::error::Error`] sources as `a: b: c`.
///
/// The CLI uses this so a nested failure prints one actionable line instead of
/// a bare top-level message.
#[must_use]
pub fn render_chain(error: &dyn std::error::Error) -> String {
    let mut rendered = error.to_string();
    let mut source = error.source();
    while let Some(current) = source {
        let text = current.to_string();
        // `#[error(transparent)]` variants re-print their child verbatim; do
        // not repeat identical links in the chain.
        if !rendered.ends_with(&text) {
            rendered.push_str(": ");
            rendered.push_str(&text);
        }
        source = current.source();
    }
    rendered
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exit_codes_match_the_documented_table() {
        assert_eq!(ExitCode::Success.code(), 0);
        assert_eq!(ExitCode::Runtime.code(), 1);
        assert_eq!(ExitCode::InvalidArguments.code(), 2);
        assert_eq!(ExitCode::InvalidInput.code(), 3);
        assert_eq!(ExitCode::OutputConflict.code(), 4);
        assert_eq!(ExitCode::ValidationFailure.code(), 5);
        assert_eq!(ExitCode::Interrupted.code(), 6);
    }

    #[test]
    fn malformed_input_classifies_as_invalid_input() {
        let error = Error::from(HeaderError::BadMagic {
            expected: *b"BAM\x01",
            found: *b"SAM\n",
        });
        assert_eq!(error.exit_code(), ExitCode::InvalidInput);
    }

    #[test]
    fn output_collisions_classify_as_output_conflict() {
        let error = Error::from(OutputError::AlreadyExists {
            path: PathBuf::from("out/chr1.bam"),
        });
        assert_eq!(error.exit_code(), ExitCode::OutputConflict);
    }

    #[test]
    fn conservation_failure_classifies_as_validation_failure() {
        let error = Error::from(ManifestError::ConservationFailed {
            input_records: 10,
            unique_emitted: 9,
            dropped: 0,
            unmatched: 0,
            difference: 1,
        });
        assert_eq!(error.exit_code(), ExitCode::ValidationFailure);
    }

    #[test]
    fn interruption_classifies_as_interrupted() {
        let error = Error::from(EngineError::Interrupted {
            records_processed: 1234,
        });
        assert_eq!(error.exit_code(), ExitCode::Interrupted);
    }

    #[test]
    fn engine_errors_inherit_the_inner_classification() {
        let error = Error::from(EngineError::from(OutputError::AlreadyExists {
            path: PathBuf::from("out/chr1.bam"),
        }));
        assert_eq!(error.exit_code(), ExitCode::OutputConflict);
    }

    #[test]
    fn render_chain_does_not_duplicate_transparent_links() {
        let error = Error::from(EngineError::from(ConfigError::InvalidTag {
            tag: "XYZ".to_string(),
        }));
        let rendered = render_chain(&error);
        assert_eq!(rendered.matches("invalid tag").count(), 1, "{rendered}");
    }
}
