// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! BAM header parsing, validation, amendment, and serialization.
//!
//! # What is preserved
//!
//! A BAM header is two things stacked on top of each other:
//!
//! ```text
//! magic     "BAM\1"
//! l_text    i32 LE
//! text      the SAM header, verbatim (@HD/@SQ/@RG/@PG/@CO)
//! n_ref     i32 LE
//! refs      n_ref * { l_name: i32 LE, name (NUL-terminated), l_ref: i32 LE }
//! ```
//!
//! The **binary dictionary** is authoritative for reference identifiers: a
//! record's `ref_id` indexes it. The **text** is what humans and most SAM
//! tooling read. `bamsplit` keeps both, byte-for-byte where it can, and
//! cross-checks them.
//!
//! # Why every output keeps the full dictionary
//!
//! `bamsplit` never reduces an output header to the one chromosome it holds.
//! Doing so would force it to rewrite every record's `ref_id` *and*
//! `next_ref_id`, which would:
//!
//! * break raw-record transfer, the whole point of the fast path;
//! * silently invalidate mate coordinates for cross-chromosome pairs;
//! * confuse downstream tools that expect the original dictionary.
//!
//! Keeping the dictionary costs a few kilobytes per output and keeps every
//! record byte-identical to its input. See `docs/bam-semantics.md`.

use std::collections::HashMap;
use std::io::{Read, Write};

use crate::error::HeaderError;

/// The BAM magic number.
pub const MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";

/// An upper bound on the SAM header text, to keep a corrupt `l_text` from
/// triggering a huge allocation. Real headers with hundreds of thousands of
/// contigs stay well under this.
const MAX_HEADER_TEXT: usize = 1 << 30;

/// An upper bound on `n_ref`, for the same reason.
const MAX_REFERENCE_COUNT: usize = 1 << 26;

/// An upper bound on a single reference name, for the same reason.
const MAX_REFERENCE_NAME: usize = 1 << 20;

/// The declared sort order of a BAM.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SortOrder {
    /// No `@HD SO` field, or `SO:unknown`.
    #[default]
    Unknown,
    /// `SO:unsorted`.
    Unsorted,
    /// `SO:queryname`.
    QueryName,
    /// `SO:coordinate`.
    Coordinate,
}

impl SortOrder {
    /// Parses a `SO` value.
    ///
    /// # Errors
    ///
    /// Returns [`HeaderError::UnknownSortOrder`] for a value outside the SAM
    /// specification.
    pub fn parse(value: &[u8]) -> Result<Self, HeaderError> {
        Ok(match value {
            b"unknown" => Self::Unknown,
            b"unsorted" => Self::Unsorted,
            b"queryname" => Self::QueryName,
            b"coordinate" => Self::Coordinate,
            other => {
                return Err(HeaderError::UnknownSortOrder {
                    value: String::from_utf8_lossy(other).into_owned(),
                });
            }
        })
    }

    /// The `SO` value as it appears in the header.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::Unsorted => "unsorted",
            Self::QueryName => "queryname",
            Self::Coordinate => "coordinate",
        }
    }

    /// Whether the header claims coordinate sorting.
    #[must_use]
    pub const fn is_coordinate(self) -> bool {
        matches!(self, Self::Coordinate)
    }
}

impl std::fmt::Display for SortOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One entry of the binary reference dictionary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSequence {
    /// The reference name, without its NUL terminator.
    pub name: Vec<u8>,
    /// The reference length in bases.
    pub length: u32,
}

impl ReferenceSequence {
    /// The name rendered lossily, for diagnostics and manifests.
    #[must_use]
    pub fn display_name(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(&self.name)
    }
}

/// Read-group metadata resolved from the `@RG` lines.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadGroupInfo {
    /// `ID`.
    pub id: Vec<u8>,
    /// `SM`, the sample.
    pub sample: Option<Vec<u8>>,
    /// `LB`, the library.
    pub library: Option<Vec<u8>>,
    /// `PL`, the sequencing platform.
    pub platform: Option<Vec<u8>>,
    /// `PU`, the platform unit.
    pub platform_unit: Option<Vec<u8>>,
    /// `CN`, the sequencing center.
    pub sequencing_center: Option<Vec<u8>>,
}

impl ReadGroupInfo {
    /// Looks up one metadata field by its logical name.
    #[must_use]
    pub fn field(&self, field: ReadGroupField) -> Option<&[u8]> {
        match field {
            ReadGroupField::ReadGroup => Some(&self.id),
            ReadGroupField::Sample => self.sample.as_deref(),
            ReadGroupField::Library => self.library.as_deref(),
            ReadGroupField::Platform => self.platform.as_deref(),
            ReadGroupField::PlatformUnit => self.platform_unit.as_deref(),
            ReadGroupField::SequencingCenter => self.sequencing_center.as_deref(),
        }
    }
}

/// A `@RG` sub-field that can drive routing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ReadGroupField {
    /// The read-group `ID` itself.
    ReadGroup,
    /// `SM`.
    Sample,
    /// `LB`.
    Library,
    /// `PL`.
    Platform,
    /// `PU`.
    PlatformUnit,
    /// `CN`.
    SequencingCenter,
}

impl ReadGroupField {
    /// The option value a user writes, e.g. `platform-unit`.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ReadGroup => "read-group",
            Self::Sample => "sample",
            Self::Library => "library",
            Self::Platform => "platform",
            Self::PlatformUnit => "platform-unit",
            Self::SequencingCenter => "sequencing-center",
        }
    }

    /// The SAM tag this field comes from, e.g. `PU`.
    #[must_use]
    pub const fn sam_tag(self) -> &'static str {
        match self {
            Self::ReadGroup => "ID",
            Self::Sample => "SM",
            Self::Library => "LB",
            Self::Platform => "PL",
            Self::PlatformUnit => "PU",
            Self::SequencingCenter => "CN",
        }
    }
}

impl std::fmt::Display for ReadGroupField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A parsed, validated BAM header.
#[derive(Debug, Clone)]
pub struct BamHeader {
    sam: noodles_sam::Header,
    references: Vec<ReferenceSequence>,
    reference_ids: HashMap<Vec<u8>, usize>,
    read_groups: Vec<ReadGroupInfo>,
    read_group_ids: HashMap<Vec<u8>, usize>,
    sort_order: SortOrder,
    original_text: Vec<u8>,
}

impl BamHeader {
    /// Reads and validates a BAM header from a decompressed stream.
    ///
    /// The reader is left positioned at the first record.
    ///
    /// # Errors
    ///
    /// Returns [`HeaderError`] for a bad magic number, a truncated header, an
    /// implausible declared length, malformed header text, or any validation
    /// failure listed on [`BamHeader::validate`].
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self, HeaderError> {
        let mut magic = [0u8; 4];
        read_exact(reader, &mut magic, "magic number")?;
        if magic != MAGIC_NUMBER {
            return Err(HeaderError::BadMagic {
                expected: MAGIC_NUMBER,
                found: magic,
            });
        }

        let l_text = read_i32(reader, "l_text")?;
        let text_len = plausible_length(l_text, MAX_HEADER_TEXT, "l_text")?;
        let mut text = vec![0u8; text_len];
        read_exact(reader, &mut text, "header text")?;

        let n_ref = read_i32(reader, "n_ref")?;
        let reference_count = plausible_length(n_ref, MAX_REFERENCE_COUNT, "n_ref")?;

        let mut references = Vec::with_capacity(reference_count.min(4096));
        for index in 0..reference_count {
            let l_name = read_i32(reader, "l_name")?;
            let name_len = plausible_length(l_name, MAX_REFERENCE_NAME, "l_name")?;
            if name_len == 0 {
                return Err(HeaderError::InvalidReferenceName {
                    index,
                    reason: "l_name must include the NUL terminator, so it cannot be zero"
                        .to_string(),
                });
            }
            let mut name = vec![0u8; name_len];
            read_exact(reader, &mut name, "reference name")?;
            if name.pop() != Some(0) {
                return Err(HeaderError::InvalidReferenceName {
                    index,
                    reason: "the reference name is not NUL-terminated".to_string(),
                });
            }
            let l_ref = read_i32(reader, "l_ref")?;
            if l_ref < 0 {
                return Err(HeaderError::InvalidReferenceLength {
                    name: String::from_utf8_lossy(&name).into_owned(),
                    index,
                    length: i64::from(l_ref),
                });
            }
            references.push(ReferenceSequence {
                name,
                length: l_ref as u32,
            });
        }

        Self::from_parts(text, references)
    }

    /// Builds a header from an already-read text block and dictionary.
    ///
    /// # Errors
    ///
    /// Returns [`HeaderError`] if the text cannot be parsed or validation
    /// fails.
    pub fn from_parts(
        text: Vec<u8>,
        references: Vec<ReferenceSequence>,
    ) -> Result<Self, HeaderError> {
        let sam = parse_sam_header(&text)?;
        let sort_order = read_sort_order(&sam)?;
        let read_groups = collect_read_groups(&sam);

        let mut reference_ids = HashMap::with_capacity(references.len());
        for (index, reference) in references.iter().enumerate() {
            if let Some(first) = reference_ids.insert(reference.name.clone(), index) {
                return Err(HeaderError::DuplicateReferenceName {
                    name: reference.display_name().into_owned(),
                    first,
                    second: index,
                });
            }
        }

        let mut read_group_ids = HashMap::with_capacity(read_groups.len());
        for (index, group) in read_groups.iter().enumerate() {
            if read_group_ids.insert(group.id.clone(), index).is_some() {
                return Err(HeaderError::DuplicateReadGroupId {
                    id: String::from_utf8_lossy(&group.id).into_owned(),
                });
            }
        }

        let header = Self {
            sam,
            references,
            reference_ids,
            read_groups,
            read_group_ids,
            sort_order,
            original_text: text,
        };
        header.validate()?;
        Ok(header)
    }

    /// Validates the header's internal consistency.
    ///
    /// Checks, in order:
    ///
    /// * reference names are non-empty and free of tabs and newlines;
    /// * reference names are unique (already enforced by construction);
    /// * the `@SQ` lines, when present, agree with the binary dictionary in
    ///   count, name, and length;
    /// * read-group identifiers are unique;
    /// * program identifiers are unique and every `PP` resolves;
    /// * the declared sort order is one the specification defines.
    ///
    /// The `@SQ` cross-check is skipped when the text carries no `@SQ` lines at
    /// all, which is legal and common — the binary dictionary is authoritative.
    ///
    /// # Errors
    ///
    /// Returns the first [`HeaderError`] found.
    pub fn validate(&self) -> Result<(), HeaderError> {
        for (index, reference) in self.references.iter().enumerate() {
            if reference.name.is_empty() {
                return Err(HeaderError::InvalidReferenceName {
                    index,
                    reason: "the reference name is empty".to_string(),
                });
            }
            if let Some(byte) = reference
                .name
                .iter()
                .find(|byte| matches!(byte, b'\t' | b'\n' | b'\r' | 0))
            {
                return Err(HeaderError::InvalidReferenceName {
                    index,
                    reason: format!("the reference name contains the forbidden byte 0x{byte:02x}"),
                });
            }
            if reference.length == 0 {
                return Err(HeaderError::InvalidReferenceLength {
                    name: reference.display_name().into_owned(),
                    index,
                    length: 0,
                });
            }
        }

        let text_references = self.sam.reference_sequences();
        if !text_references.is_empty() {
            if text_references.len() != self.references.len() {
                return Err(HeaderError::ReferenceCountMismatch {
                    binary: self.references.len(),
                    text: text_references.len(),
                });
            }
            for (index, ((text_name, text_map), binary)) in
                text_references.iter().zip(&self.references).enumerate()
            {
                if text_name.as_slice() != binary.name.as_slice() {
                    return Err(HeaderError::ReferenceNameMismatch {
                        index,
                        binary: binary.display_name().into_owned(),
                        text: String::from_utf8_lossy(text_name.as_slice()).into_owned(),
                    });
                }
                let text_length = u32::try_from(text_map.length().get()).unwrap_or(u32::MAX);
                if text_length != binary.length {
                    return Err(HeaderError::ReferenceLengthMismatch {
                        index,
                        name: binary.display_name().into_owned(),
                        binary: binary.length,
                        text: text_length,
                    });
                }
            }
        }

        self.validate_program_records()?;
        Ok(())
    }

    fn validate_program_records(&self) -> Result<(), HeaderError> {
        use noodles_sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID;

        let programs = self.sam.programs().as_ref();
        // `IndexMap` keys are unique by construction, and the noodles parser
        // rejects a duplicate `@PG ID` before we get here; this loop exists so
        // a header assembled through the Rust API is held to the same standard.
        let mut seen: HashMap<&[u8], usize> = HashMap::with_capacity(programs.len());
        for (index, id) in programs.keys().enumerate() {
            if seen.insert(id.as_slice(), index).is_some() {
                return Err(HeaderError::DuplicateProgramId {
                    id: String::from_utf8_lossy(id.as_slice()).into_owned(),
                });
            }
        }
        for (id, program) in programs {
            if let Some(previous) = program.other_fields().get(&PREVIOUS_PROGRAM_ID)
                && !programs.contains_key(previous.as_slice())
            {
                return Err(HeaderError::DanglingProgramLink {
                    id: String::from_utf8_lossy(id.as_slice()).into_owned(),
                    previous: String::from_utf8_lossy(previous.as_slice()).into_owned(),
                });
            }
        }
        Ok(())
    }

    /// The parsed SAM header.
    #[must_use]
    pub const fn sam(&self) -> &noodles_sam::Header {
        &self.sam
    }

    /// The binary reference dictionary, in `ref_id` order.
    #[must_use]
    pub fn references(&self) -> &[ReferenceSequence] {
        &self.references
    }

    /// The number of reference sequences.
    #[must_use]
    pub fn reference_count(&self) -> usize {
        self.references.len()
    }

    /// The name of reference `id`, or [`None`] when out of range.
    #[must_use]
    pub fn reference_name(&self, id: usize) -> Option<&[u8]> {
        self.references.get(id).map(|reference| &reference.name[..])
    }

    /// The length of reference `id`, or [`None`] when out of range.
    #[must_use]
    pub fn reference_length(&self, id: usize) -> Option<u32> {
        self.references.get(id).map(|reference| reference.length)
    }

    /// The `ref_id` of a reference name.
    #[must_use]
    pub fn reference_id(&self, name: &[u8]) -> Option<usize> {
        self.reference_ids.get(name).copied()
    }

    /// Whether `id` addresses a reference in the dictionary.
    #[must_use]
    pub fn is_valid_reference_id(&self, id: i32) -> bool {
        usize::try_from(id).is_ok_and(|id| id < self.references.len())
    }

    /// The resolved `@RG` metadata, in header order.
    #[must_use]
    pub fn read_groups(&self) -> &[ReadGroupInfo] {
        &self.read_groups
    }

    /// Looks up read-group metadata by `ID`.
    #[must_use]
    pub fn read_group(&self, id: &[u8]) -> Option<&ReadGroupInfo> {
        self.read_group_ids
            .get(id)
            .and_then(|index| self.read_groups.get(*index))
    }

    /// The declared sort order.
    #[must_use]
    pub const fn sort_order(&self) -> SortOrder {
        self.sort_order
    }

    /// The original header text, exactly as it appeared on disk.
    #[must_use]
    pub fn original_text(&self) -> &[u8] {
        &self.original_text
    }

    /// The largest reference length in the dictionary.
    ///
    /// Used to decide between BAI and CSI: BAI addresses at most
    /// [`crate::index::bai::MAX_POSITION`] bases.
    #[must_use]
    pub fn max_reference_length(&self) -> u32 {
        self.references
            .iter()
            .map(|reference| reference.length)
            .max()
            .unwrap_or(0)
    }

    /// A stable digest of the header, for the manifest.
    ///
    /// Covers the text and the binary dictionary, so two runs over the same
    /// input produce the same value and a re-headered input produces a
    /// different one.
    #[must_use]
    pub fn checksum(&self) -> u64 {
        use xxhash_rust::xxh3::Xxh3;

        let mut hasher = Xxh3::new();
        hasher.update(&self.original_text);
        hasher.update(&(self.references.len() as u64).to_le_bytes());
        for reference in &self.references {
            hasher.update(&(reference.name.len() as u64).to_le_bytes());
            hasher.update(&reference.name);
            hasher.update(&reference.length.to_le_bytes());
        }
        hasher.digest()
    }

    /// Overrides the declared sort order, updating both the parsed header and
    /// the text that will be written.
    ///
    /// Used when an output's ordering differs from the input's — for example
    /// when the spool engine emits records in input order from an unsorted BAM.
    pub fn set_sort_order(&mut self, sort_order: SortOrder) {
        use noodles_sam::header::record::value::map::{self, Map, header as hd};

        self.sort_order = sort_order;
        let entry = self
            .sam
            .header_mut()
            .get_or_insert_with(|| Map::<map::Header>::new(hd::Version::new(1, 6)));
        entry
            .other_fields_mut()
            .insert(hd::tag::SORT_ORDER, bstring(sort_order.as_str().as_bytes()));
    }

    /// Appends an `@PG` record describing this invocation.
    ///
    /// The identifier is `bamsplit`; if that is taken, `bamsplit.1`,
    /// `bamsplit.2`, and so on are tried in order, so the result is
    /// deterministic for a given input header.
    ///
    /// `PP` is linked to the last leaf program in the existing chain, which
    /// preserves the chain without rewriting any existing `@PG` record. When
    /// the header has no programs, no `PP` is emitted.
    ///
    /// Returns the identifier that was used.
    ///
    /// # Errors
    ///
    /// Returns [`HeaderError`] if the program record cannot be built.
    pub fn append_program_record(
        &mut self,
        version: &str,
        command_line: &str,
    ) -> Result<Vec<u8>, HeaderError> {
        use noodles_sam::header::record::value::{Map, map::Program, map::program::tag};

        let id = self.next_program_id();
        let previous = self.leaf_program_id();

        let mut builder = Map::<Program>::builder()
            .insert(tag::NAME, "bamsplit")
            .insert(tag::VERSION, version)
            .insert(tag::COMMAND_LINE, command_line);
        if let Some(previous) = &previous {
            builder = builder.insert(tag::PREVIOUS_PROGRAM_ID, bstring(previous));
        }
        let program = builder
            .build()
            .map_err(|error| HeaderError::MalformedText {
                message: format!("cannot build the bamsplit @PG record: {error}"),
            })?;

        self.sam
            .programs_mut()
            .as_mut()
            .insert(bstring(&id), program);
        Ok(id)
    }

    fn next_program_id(&self) -> Vec<u8> {
        let programs = self.sam.programs().as_ref();
        if !programs.contains_key(&b"bamsplit"[..]) {
            return b"bamsplit".to_vec();
        }
        // Bounded by the number of existing programs plus one, so this always
        // terminates.
        for suffix in 1..=programs.len().saturating_add(1) {
            let candidate = format!("bamsplit.{suffix}").into_bytes();
            if !programs.contains_key(candidate.as_slice()) {
                return candidate;
            }
        }
        // Unreachable: the loop tries more candidates than there are programs.
        format!("bamsplit.{}", programs.len() + 1).into_bytes()
    }

    /// The last `@PG` record that no other record names as its `PP`.
    fn leaf_program_id(&self) -> Option<Vec<u8>> {
        use noodles_sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID;

        let programs = self.sam.programs().as_ref();
        if programs.is_empty() {
            return None;
        }
        let referenced: std::collections::HashSet<&[u8]> = programs
            .values()
            .filter_map(|program| program.other_fields().get(&PREVIOUS_PROGRAM_ID))
            .map(|value| value.as_slice())
            .collect();
        programs
            .keys()
            .rev()
            .find(|id| !referenced.contains(id.as_slice()))
            .map(|id| id.as_slice().to_vec())
    }

    /// Serializes the header text as it would appear in a BAM.
    ///
    /// # Errors
    ///
    /// Returns [`HeaderError::Serialize`] if `noodles` cannot render the
    /// header.
    pub fn to_text(&self) -> Result<Vec<u8>, HeaderError> {
        let mut writer = noodles_sam::io::Writer::new(Vec::new());
        writer
            .write_header(&self.sam)
            .map_err(|source| HeaderError::Serialize { source })?;
        Ok(writer.into_inner())
    }

    /// Writes the complete binary BAM header.
    ///
    /// # Errors
    ///
    /// Returns [`HeaderError`] if serialization or the underlying write fails.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<(), HeaderError> {
        let text = self.to_text()?;
        let l_text = i32::try_from(text.len()).map_err(|_| HeaderError::InvalidLength {
            field: "l_text",
            value: text.len() as i64,
        })?;
        let n_ref =
            i32::try_from(self.references.len()).map_err(|_| HeaderError::InvalidLength {
                field: "n_ref",
                value: self.references.len() as i64,
            })?;

        let io = |source| HeaderError::Io { source };
        writer.write_all(&MAGIC_NUMBER).map_err(io)?;
        writer.write_all(&l_text.to_le_bytes()).map_err(io)?;
        writer.write_all(&text).map_err(io)?;
        writer.write_all(&n_ref.to_le_bytes()).map_err(io)?;

        for reference in &self.references {
            let l_name = i32::try_from(reference.name.len() + 1).map_err(|_| {
                HeaderError::InvalidLength {
                    field: "l_name",
                    value: reference.name.len() as i64 + 1,
                }
            })?;
            writer.write_all(&l_name.to_le_bytes()).map_err(io)?;
            writer.write_all(&reference.name).map_err(io)?;
            writer.write_all(&[0]).map_err(io)?;
            let l_ref =
                i32::try_from(reference.length).map_err(|_| HeaderError::InvalidLength {
                    field: "l_ref",
                    value: i64::from(reference.length),
                })?;
            writer.write_all(&l_ref.to_le_bytes()).map_err(io)?;
        }
        Ok(())
    }

    /// The serialized size of the binary header, in bytes.
    ///
    /// # Errors
    ///
    /// Returns [`HeaderError::Serialize`] if the text cannot be rendered.
    pub fn serialized_len(&self) -> Result<usize, HeaderError> {
        let text = self.to_text()?;
        let dictionary: usize = self
            .references
            .iter()
            .map(|reference| 4 + reference.name.len() + 1 + 4)
            .sum();
        Ok(4 + 4 + text.len() + 4 + dictionary)
    }
}

fn parse_sam_header(text: &[u8]) -> Result<noodles_sam::Header, HeaderError> {
    let mut parser = noodles_sam::header::Parser::default();
    for (index, line) in text.split(|byte| *byte == b'\n').enumerate() {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        parser
            .parse_partial(line)
            .map_err(|error| HeaderError::MalformedText {
                message: format!("line {}: {error}", index + 1),
            })?;
    }
    Ok(parser.finish())
}

fn read_sort_order(header: &noodles_sam::Header) -> Result<SortOrder, HeaderError> {
    use noodles_sam::header::record::value::map::header::tag::SORT_ORDER;

    let Some(value) = header
        .header()
        .and_then(|entry| entry.other_fields().get(&SORT_ORDER))
    else {
        return Ok(SortOrder::Unknown);
    };
    SortOrder::parse(value.as_slice())
}

fn collect_read_groups(header: &noodles_sam::Header) -> Vec<ReadGroupInfo> {
    use noodles_sam::header::record::value::map::read_group::tag;

    header
        .read_groups()
        .iter()
        .map(|(id, map)| {
            let get = |key| {
                map.other_fields()
                    .get(&key)
                    .map(|value| value.as_slice().to_vec())
            };
            ReadGroupInfo {
                id: id.as_slice().to_vec(),
                sample: get(tag::SAMPLE),
                library: get(tag::LIBRARY),
                platform: get(tag::PLATFORM),
                platform_unit: get(tag::PLATFORM_UNIT),
                sequencing_center: get(tag::SEQUENCING_CENTER),
            }
        })
        .collect()
}

/// `noodles-sam` stores every header value as a [`bstr::BString`]; this keeps
/// the conversion in one place.
fn bstring(bytes: &[u8]) -> bstr::BString {
    bstr::BString::from(bytes.to_vec())
}

fn read_exact<R: Read>(
    reader: &mut R,
    buf: &mut [u8],
    section: &'static str,
) -> Result<(), HeaderError> {
    match reader.read_exact(buf) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => {
            Err(HeaderError::Truncated { section })
        }
        Err(source) => Err(HeaderError::Io { source }),
    }
}

fn read_i32<R: Read>(reader: &mut R, section: &'static str) -> Result<i32, HeaderError> {
    let mut buf = [0u8; 4];
    read_exact(reader, &mut buf, section)?;
    Ok(i32::from_le_bytes(buf))
}

fn plausible_length(value: i32, limit: usize, field: &'static str) -> Result<usize, HeaderError> {
    if value < 0 || value as usize > limit {
        return Err(HeaderError::InvalidLength {
            field,
            value: i64::from(value),
        });
    }
    Ok(value as usize)
}
