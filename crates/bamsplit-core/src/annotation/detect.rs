// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Working out what an annotation file actually is.
//!
//! # Why the extension is not enough
//!
//! `genes.bed` containing GTF is common enough that trusting the extension is a
//! bug waiting to happen: `genepred` would parse the GTF's `seqname`,
//! `source`, and `feature` columns as `chrom`, `start`, and `end`, fail on the
//! non-numeric start, and report a parse error that says nothing about the real
//! problem.
//!
//! So detection is layered, and **content wins over extension**:
//!
//! 1. an explicit `--format`;
//! 2. the extension, after stripping any compression suffix;
//! 3. the contents of the first non-empty, non-comment records;
//! 4. otherwise a clear error naming what was seen.
//!
//! Steps 2 and 3 both run. When they disagree, step 3 decides and the reason
//! string records the disagreement, so a surprising result can be traced.
//!
//! # Compressed input
//!
//! Compression is detected from **magic bytes**, not the extension, so a
//! gzip-compressed file named `.bed` is still inspected properly. The first
//! 64 KiB are decompressed for content inspection; the whole file is never
//! materialized just to be looked at.
//!
//! # BED width
//!
//! With `--type auto` the first data row's column count picks the highest
//! supported width consistent with it, and every later row must agree. Columns
//! beyond that width become additional fields rather than an error.

use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;

use crate::error::AnnotationError;

/// How many bytes of a file are decompressed for content inspection.
///
/// Large enough for a GFF3's `##` pragma block plus a few real records, small
/// enough that inspecting a 20 GiB annotation is instant.
pub const INSPECTION_BYTES: usize = 64 * 1024;

/// How many data rows are checked for a consistent column count.
///
/// Checking every row would mean reading the whole file twice. A few thousand
/// catches the realistic mistake — a file concatenated from two sources — and
/// `genepred` rejects anything that slips past.
pub const CONSISTENCY_ROWS: usize = 4_096;

/// The annotation formats `bamsplit` accepts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationFormat {
    /// BED, of some width.
    Bed,
    /// GTF, also known as GFF2.
    Gtf,
    /// GFF3.
    Gff,
}

impl AnnotationFormat {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Bed => "bed",
            Self::Gtf => "gtf",
            Self::Gff => "gff",
        }
    }

    /// Parses an explicit `--format` value.
    ///
    /// # Errors
    ///
    /// Returns [`AnnotationError::UnsupportedBedType`] never; an unrecognized
    /// value yields [`None`] so the caller can raise a configuration error with
    /// its own option name.
    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "bed" => Self::Bed,
            "gtf" | "gff2" => Self::Gtf,
            "gff" | "gff3" => Self::Gff,
            _ => return None,
        })
    }

    /// Whether this format carries per-transcript block structure.
    #[must_use]
    pub const fn is_gxf(self) -> bool {
        matches!(self, Self::Gtf | Self::Gff)
    }
}

impl std::fmt::Display for AnnotationFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The supported BED widths.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum BedType {
    /// `chrom start end`
    Bed3,
    /// … `name`
    Bed4,
    /// … `score`
    Bed5,
    /// … `strand`
    Bed6,
    /// … `thickStart thickEnd`
    Bed8,
    /// … `itemRgb`
    Bed9,
    /// … `blockCount blockSizes blockStarts`
    Bed12,
}

/// Every supported width, ascending.
pub const BED_TYPES: [BedType; 7] = [
    BedType::Bed3,
    BedType::Bed4,
    BedType::Bed5,
    BedType::Bed6,
    BedType::Bed8,
    BedType::Bed9,
    BedType::Bed12,
];

impl BedType {
    /// The number of standard columns.
    #[must_use]
    pub const fn width(self) -> usize {
        match self {
            Self::Bed3 => 3,
            Self::Bed4 => 4,
            Self::Bed5 => 5,
            Self::Bed6 => 6,
            Self::Bed8 => 8,
            Self::Bed9 => 9,
            Self::Bed12 => 12,
        }
    }

    /// Whether this width carries `thickStart`/`thickEnd`, which `cds` and the
    /// UTR features need.
    #[must_use]
    pub const fn has_coding_bounds(self) -> bool {
        self.width() >= 8
    }

    /// Whether this width carries block structure.
    #[must_use]
    pub const fn has_blocks(self) -> bool {
        matches!(self, Self::Bed12)
    }

    /// Whether this width carries a strand.
    #[must_use]
    pub const fn has_strand(self) -> bool {
        self.width() >= 6
    }

    /// Parses an explicit `--type` value.
    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        BED_TYPES
            .into_iter()
            .find(|candidate| candidate.width().to_string() == value)
    }

    /// The widest supported type no wider than `columns`.
    #[must_use]
    pub fn widest_fitting(columns: usize) -> Option<Self> {
        BED_TYPES
            .into_iter()
            .rev()
            .find(|candidate| candidate.width() <= columns)
    }
}

impl std::fmt::Display for BedType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.width())
    }
}

/// The compression a file actually uses, from its magic bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    /// Plain text.
    None,
    /// gzip.
    Gzip,
    /// BGZF — gzip with a `BC` extra field, so gzip decoders read it.
    Bgzf,
    /// Zstandard.
    Zstd,
    /// bzip2.
    Bzip2,
}

impl Compression {
    /// A short label for the manifest.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Gzip => "gzip",
            Self::Bgzf => "bgzf",
            Self::Zstd => "zstd",
            Self::Bzip2 => "bzip2",
        }
    }

    /// Whether the bytes are compressed at all.
    #[must_use]
    pub const fn is_compressed(self) -> bool {
        !matches!(self, Self::None)
    }

    /// The filename suffix `genepred` recognizes for this codec.
    ///
    /// [`Compression::Bgzf`] maps to `.gz`, because `genepred` has no `.bgz`
    /// mapping but BGZF is gzip-compatible.
    #[must_use]
    pub const fn genepred_suffix(self) -> Option<&'static str> {
        match self {
            Self::None => None,
            Self::Gzip | Self::Bgzf => Some("gz"),
            Self::Zstd => Some("zst"),
            Self::Bzip2 => Some("bz2"),
        }
    }

    /// Detects the codec from a file's first bytes.
    #[must_use]
    pub fn from_magic(head: &[u8]) -> Self {
        match head {
            [0x1f, 0x8b, ..] => {
                // BGZF is gzip with FEXTRA and a `BC` subfield at offset 12.
                if head.len() >= 16 && head[3] & 0x04 != 0 && head[12] == b'B' && head[13] == b'C' {
                    Self::Bgzf
                } else {
                    Self::Gzip
                }
            }
            [0x28, 0xb5, 0x2f, 0xfd, ..] => Self::Zstd,
            [0x42, 0x5a, 0x68, ..] => Self::Bzip2,
            _ => Self::None,
        }
    }
}

impl std::fmt::Display for Compression {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// What detection concluded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Detection {
    /// The format to parse as.
    pub format: AnnotationFormat,
    /// The BED width, for BED input only.
    pub bed_type: Option<BedType>,
    /// How the file is compressed.
    pub compression: Compression,
    /// Columns beyond the chosen BED width.
    pub trailing_columns: usize,
    /// How the conclusion was reached.
    pub reason: String,
}

/// Reads and decompresses the first [`INSPECTION_BYTES`] of a file.
///
/// Returns the decompressed head and the codec that was used.
///
/// # Errors
///
/// Returns [`AnnotationError::Io`] if the file cannot be opened or read. A
/// truncated compressed stream is *not* an error here: whatever decoded before
/// the stream ran out is enough to identify a format, and a genuinely broken
/// file will fail loudly when `genepred` reads it.
pub fn read_head(path: &Path) -> Result<(Vec<u8>, Compression), AnnotationError> {
    let io = |source| AnnotationError::Io {
        path: path.to_path_buf(),
        source,
    };

    let mut file = std::fs::File::open(path).map_err(io)?;
    let mut magic = [0u8; 16];
    let magic_len = read_up_to(&mut file, &mut magic).map_err(io)?;
    let compression = Compression::from_magic(&magic[..magic_len]);

    let mut raw = magic[..magic_len].to_vec();
    // Compressed formats need more input than output to be worth decoding, and
    // an uncompressed file needs no more than the inspection budget.
    let budget = if compression.is_compressed() {
        INSPECTION_BYTES * 4
    } else {
        INSPECTION_BYTES
    };
    let mut rest = vec![0u8; budget.saturating_sub(raw.len())];
    let read = read_up_to(&mut file, &mut rest).map_err(io)?;
    raw.extend_from_slice(&rest[..read]);

    let head = decompress_head(&raw, compression);
    Ok((head, compression))
}

/// Decodes as much of `raw` as the codec can manage, tolerating truncation.
fn decompress_head(raw: &[u8], compression: Compression) -> Vec<u8> {
    let mut out = Vec::with_capacity(INSPECTION_BYTES);
    let taken = match compression {
        Compression::None => {
            out.extend_from_slice(&raw[..raw.len().min(INSPECTION_BYTES)]);
            return out;
        }
        Compression::Gzip | Compression::Bgzf => {
            // `MultiGzDecoder` handles both a plain gzip member and the
            // concatenated members BGZF is made of.
            let mut decoder = flate2::read::MultiGzDecoder::new(raw);
            decoder
                .by_ref()
                .take(INSPECTION_BYTES as u64)
                .read_to_end(&mut out)
        }
        Compression::Zstd => match zstd::stream::read::Decoder::new(raw) {
            Ok(mut decoder) => decoder
                .by_ref()
                .take(INSPECTION_BYTES as u64)
                .read_to_end(&mut out),
            Err(error) => Err(error),
        },
        Compression::Bzip2 => {
            let mut decoder = bzip2::read::MultiBzDecoder::new(raw);
            decoder
                .by_ref()
                .take(INSPECTION_BYTES as u64)
                .read_to_end(&mut out)
        }
    };
    // A truncated tail is expected: the input was cut mid-stream on purpose.
    let _ = taken;
    out
}

fn read_up_to<R: Read>(reader: &mut R, buf: &mut [u8]) -> std::io::Result<usize> {
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

/// Whether a line carries data rather than metadata.
fn is_data_line(line: &[u8]) -> bool {
    !(line.is_empty()
        || line.starts_with(b"#")
        || line.starts_with(b"track")
        || line.starts_with(b"browser"))
}

/// The data rows of a decompressed head, as tab-separated column slices.
///
/// The final line is dropped when the head was truncated mid-line, so a
/// half-read row cannot be mistaken for a short one.
fn data_rows(head: &[u8], truncated: bool) -> Vec<Vec<&[u8]>> {
    let mut lines: Vec<&[u8]> = head
        .split(|byte| *byte == b'\n')
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
        .collect();
    if truncated {
        lines.pop();
    }
    lines
        .into_iter()
        .filter(|line| is_data_line(line))
        .map(|line| line.split(|byte| *byte == b'\t').collect())
        .collect()
}

fn is_integer(field: &[u8]) -> bool {
    !field.is_empty() && field.iter().all(u8::is_ascii_digit)
}

/// Classifies a single data row.
fn classify_row(columns: &[&[u8]]) -> Option<AnnotationFormat> {
    // GTF/GFF: nine columns, numeric start and end in 4 and 5, a strand in 7.
    if columns.len() >= 9
        && is_integer(columns[3])
        && is_integer(columns[4])
        && matches!(columns[6], b"+" | b"-" | b"." | b"?")
    {
        let attributes = columns[8];
        // GTF writes `key "value";`, GFF3 writes `key=value;`. A GTF attribute
        // block can contain `=` inside a quoted value, so the quote test comes
        // first.
        return Some(if attributes.contains(&b'"') {
            AnnotationFormat::Gtf
        } else if attributes.contains(&b'=') {
            AnnotationFormat::Gff
        } else {
            AnnotationFormat::Gtf
        });
    }
    // BED: at least three columns, numeric start and end in 2 and 3.
    if columns.len() >= 3 && is_integer(columns[1]) && is_integer(columns[2]) {
        return Some(AnnotationFormat::Bed);
    }
    None
}

/// The format the file extension advertises, ignoring any compression suffix.
fn format_from_extension(path: &Path) -> Option<AnnotationFormat> {
    let mut name = path.file_name()?.to_string_lossy().into_owned();
    name.make_ascii_lowercase();
    for suffix in [".gz", ".bgz", ".zst", ".zstd", ".bz2", ".bzip2"] {
        if let Some(stripped) = name.strip_suffix(suffix) {
            name = stripped.to_string();
            break;
        }
    }
    // `name` was lower-cased above, so a plain suffix test is already
    // case-insensitive; `Path::extension` is avoided because a `.gff3` after a
    // stripped `.gz` is not what it would report.
    let suffix = name.rsplit_once('.').map(|(_, suffix)| suffix);
    match suffix {
        Some("bed") => Some(AnnotationFormat::Bed),
        Some("gtf") => Some(AnnotationFormat::Gtf),
        Some("gff" | "gff3") => Some(AnnotationFormat::Gff),
        _ => None,
    }
}

/// Detects the format and, for BED, the width.
///
/// # Errors
///
/// Returns [`AnnotationError::Empty`] for a file with no data rows,
/// [`AnnotationError::AmbiguousFormat`] when neither the extension nor the
/// contents identify it, [`AnnotationError::AmbiguousBedType`] when the width
/// cannot be resolved, [`AnnotationError::MixedBedWidths`] when rows disagree,
/// [`AnnotationError::BedTypeOnNonBed`] when `--type` is given for GTF/GFF, and
/// [`AnnotationError::Io`] if the file cannot be read.
pub fn detect(
    path: &Path,
    requested_format: Option<AnnotationFormat>,
    requested_bed_type: Option<BedType>,
) -> Result<Detection, AnnotationError> {
    let (head, compression) = read_head(path)?;
    let truncated = head.len() >= INSPECTION_BYTES;
    let rows = data_rows(&head, truncated);

    if rows.is_empty() {
        return Err(AnnotationError::Empty {
            path: path.to_path_buf(),
        });
    }

    let from_extension = format_from_extension(path);
    let from_content = rows.iter().find_map(|row| classify_row(row));

    let (format, mut reason) = match (requested_format, from_content, from_extension) {
        (Some(explicit), _, _) => (
            explicit,
            format!("`--format {explicit}` was requested explicitly"),
        ),
        (None, Some(content), Some(extension)) if content == extension => (
            content,
            format!("the extension and the contents both say {content}"),
        ),
        (None, Some(content), Some(extension)) => (
            content,
            format!("the extension says {extension} but the contents are {content}; contents win"),
        ),
        (None, Some(content), None) => (
            content,
            format!("the contents are {content}; the extension says nothing"),
        ),
        (None, None, Some(extension)) => (
            extension,
            format!("no data row could be classified, so the {extension} extension was trusted"),
        ),
        (None, None, None) => {
            let sample = rows.first().map_or_else(
                || "no data rows".to_string(),
                |row| format!("{} column(s)", row.len()),
            );
            return Err(AnnotationError::AmbiguousFormat {
                path: path.to_path_buf(),
                reason: format!(
                    "the extension is not recognized and the first data row has {sample}, which \
                     matches neither BED (numeric columns 2 and 3) nor GTF/GFF (nine columns with \
                     numeric 4 and 5 and a strand in 7)"
                ),
            });
        }
    };

    if format.is_gxf() {
        if requested_bed_type.is_some() {
            return Err(AnnotationError::BedTypeOnNonBed {
                path: path.to_path_buf(),
                format: format.as_str(),
            });
        }
        return Ok(Detection {
            format,
            bed_type: None,
            compression,
            trailing_columns: 0,
            reason,
        });
    }

    let (bed_type, columns) = resolve_bed_type(path, &rows, requested_bed_type, &mut reason)?;
    Ok(Detection {
        format,
        bed_type: Some(bed_type),
        compression,
        trailing_columns: columns.saturating_sub(bed_type.width()),
        reason,
    })
}

/// Picks a BED width and checks that every inspected row agrees.
fn resolve_bed_type(
    path: &Path,
    rows: &[Vec<&[u8]>],
    requested: Option<BedType>,
    reason: &mut String,
) -> Result<(BedType, usize), AnnotationError> {
    let Some(first) = rows.first() else {
        return Err(AnnotationError::Empty {
            path: path.to_path_buf(),
        });
    };
    let columns = first.len();

    // A row's column count must not change partway through the file. Checking a
    // bounded prefix catches the realistic mistake without a second full pass.
    for (index, row) in rows.iter().take(CONSISTENCY_ROWS).enumerate() {
        if row.len() != columns {
            return Err(AnnotationError::MixedBedWidths {
                path: path.to_path_buf(),
                expected: columns,
                actual: row.len(),
                line: index + 1,
            });
        }
    }

    if let Some(requested) = requested {
        if requested.width() > columns {
            return Err(AnnotationError::AmbiguousBedType {
                path: path.to_path_buf(),
                reason: format!(
                    "`--type {requested}` needs {} columns but the first data row has {columns}",
                    requested.width()
                ),
            });
        }
        let _ = write!(reason, "; `--type {requested}` was requested explicitly");
        return Ok((requested, columns));
    }

    // The widest supported type the row can support, stepping down when the
    // extra columns do not actually parse as the fields that width requires.
    let mut candidate =
        BedType::widest_fitting(columns).ok_or_else(|| AnnotationError::AmbiguousBedType {
            path: path.to_path_buf(),
            reason: format!("the first data row has {columns} column(s); BED needs at least 3"),
        })?;
    while !row_supports(first, candidate) {
        let Some(narrower) = BED_TYPES
            .into_iter()
            .rev()
            .find(|next| *next < candidate && row_supports(first, *next))
        else {
            return Err(AnnotationError::AmbiguousBedType {
                path: path.to_path_buf(),
                reason: format!(
                    "the first data row has {columns} column(s) but none of the supported widths \
                     fit its contents"
                ),
            });
        };
        candidate = narrower;
    }

    let _ = write!(reason, "; {columns} column(s) resolved to BED{candidate}");
    Ok((candidate, columns))
}

/// Whether a row's fields actually parse as the given width requires.
fn row_supports(row: &[&[u8]], bed_type: BedType) -> bool {
    let width = bed_type.width();
    if row.len() < width {
        return false;
    }
    // Columns 2 and 3 are the start and end for every width.
    if !is_integer(row[1]) || !is_integer(row[2]) {
        return false;
    }
    if width >= 5 && !is_integer(row[4]) {
        return false;
    }
    if width >= 6 && !matches!(row[5], b"+" | b"-" | b"." | b"?") {
        return false;
    }
    if width >= 8 && (!is_integer(row[6]) || !is_integer(row[7])) {
        return false;
    }
    if width >= 12 && !is_integer(row[9]) {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write(directory: &tempfile::TempDir, name: &str, contents: &[u8]) -> std::path::PathBuf {
        let path = directory.path().join(name);
        std::fs::write(&path, contents).expect("written");
        path
    }

    fn gzip(contents: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(contents).expect("compressed");
        encoder.finish().expect("finished")
    }

    const BED12: &[u8] =
        b"chr1\t100\t500\tgene1\t0\t+\t150\t450\t0,0,0\t2\t100,100\t0,300\nchr1\t600\t900\tgene2\t0\t-\t650\t850\t0,0,0\t1\t300\t0\n";
    const BED6: &[u8] = b"chr1\t100\t500\tgene1\t0\t+\nchr1\t600\t900\tgene2\t0\t-\n";
    const BED3: &[u8] = b"chr1\t100\t500\nchr1\t600\t900\n";
    const GTF: &[u8] = b"##description: test\nchr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\nchr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n";
    const GFF3: &[u8] = b"##gff-version 3\nchr1\ttest\texon\t101\t200\t.\t+\t.\tID=e1;Parent=t1\nchr1\ttest\texon\t301\t400\t.\t+\t.\tID=e2;Parent=t1\n";

    #[test]
    fn detects_bed_widths_from_content() {
        let directory = tempfile::tempdir().expect("temp dir");
        for (contents, expected, trailing) in [
            (BED3, BedType::Bed3, 0),
            (BED6, BedType::Bed6, 0),
            (BED12, BedType::Bed12, 0),
        ] {
            let path = write(&directory, "a.bed", contents);
            let detection = detect(&path, None, None).expect("detected");
            assert_eq!(detection.format, AnnotationFormat::Bed);
            assert_eq!(detection.bed_type, Some(expected));
            assert_eq!(detection.trailing_columns, trailing);
        }
    }

    #[test]
    fn trailing_columns_are_reported_not_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        // 14 columns: BED12 plus two extras.
        let contents = b"chr1\t100\t500\tg\t0\t+\t150\t450\t0,0,0\t1\t400\t0\textra1\textra2\n";
        let path = write(&directory, "a.bed", contents);
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.bed_type, Some(BedType::Bed12));
        assert_eq!(detection.trailing_columns, 2);

        // 7 columns: BED6 plus one extra, because BED8 needs 8.
        let path = write(&directory, "b.bed", b"chr1\t100\t500\tg\t0\t+\tnote\n");
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.bed_type, Some(BedType::Bed6));
        assert_eq!(detection.trailing_columns, 1);
    }

    #[test]
    fn a_width_whose_fields_do_not_parse_steps_down() {
        let directory = tempfile::tempdir().expect("temp dir");
        // Nine columns, but column 7 is not an integer, so BED9 and BED8 are
        // out and BED6 is the widest that fits.
        let path = write(
            &directory,
            "a.bed",
            b"chr1\t100\t500\tg\t0\t+\tnotanumber\tx\ty\n",
        );
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.bed_type, Some(BedType::Bed6));
        assert_eq!(detection.trailing_columns, 3);
    }

    #[test]
    fn detects_gtf_and_gff3_from_content() {
        let directory = tempfile::tempdir().expect("temp dir");
        let gtf = write(&directory, "a.gtf", GTF);
        assert_eq!(
            detect(&gtf, None, None).expect("detected").format,
            AnnotationFormat::Gtf
        );
        let gff = write(&directory, "a.gff3", GFF3);
        assert_eq!(
            detect(&gff, None, None).expect("detected").format,
            AnnotationFormat::Gff
        );
    }

    #[test]
    fn content_beats_a_lying_extension() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "genes.bed", GTF);
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.format, AnnotationFormat::Gtf);
        assert!(
            detection.reason.contains("contents win"),
            "{}",
            detection.reason
        );

        let path = write(&directory, "genes.gtf", BED12);
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.format, AnnotationFormat::Bed);
    }

    #[test]
    fn an_explicit_format_wins_over_everything() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "genes.bed", BED12);
        let detection = detect(&path, Some(AnnotationFormat::Gtf), None).expect("detected");
        assert_eq!(detection.format, AnnotationFormat::Gtf);
        assert!(detection.reason.contains("explicitly"));
    }

    #[test]
    fn compression_is_detected_from_magic_not_extension() {
        let directory = tempfile::tempdir().expect("temp dir");
        // A gzip file that does not admit it in its name.
        let path = write(&directory, "genes.bed", &gzip(BED12));
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.compression, Compression::Gzip);
        assert_eq!(detection.format, AnnotationFormat::Bed);
        assert_eq!(detection.bed_type, Some(BedType::Bed12));
    }

    #[test]
    fn a_properly_named_compressed_file_is_detected_too() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "genes.gtf.gz", &gzip(GTF));
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.compression, Compression::Gzip);
        assert_eq!(detection.format, AnnotationFormat::Gtf);
    }

    #[test]
    fn magic_bytes_classify_every_codec() {
        assert_eq!(Compression::from_magic(b"chr1\t1"), Compression::None);
        assert_eq!(
            Compression::from_magic(&[0x1f, 0x8b, 0x08, 0x00]),
            Compression::Gzip
        );
        assert_eq!(
            Compression::from_magic(&[
                0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0
            ]),
            Compression::Bgzf
        );
        assert_eq!(
            Compression::from_magic(&[0x28, 0xb5, 0x2f, 0xfd]),
            Compression::Zstd
        );
        assert_eq!(Compression::from_magic(b"BZh9"), Compression::Bzip2);
    }

    #[test]
    fn comments_and_track_lines_are_skipped() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut contents = b"# a comment\ntrack name=x\nbrowser position chr1\n\n".to_vec();
        contents.extend_from_slice(BED6);
        let path = write(&directory, "a.bed", &contents);
        let detection = detect(&path, None, None).expect("detected");
        assert_eq!(detection.bed_type, Some(BedType::Bed6));
    }

    #[test]
    fn mixed_bed_widths_are_rejected_with_a_line_number() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(
            &directory,
            "a.bed",
            b"chr1\t100\t500\tg1\t0\t+\nchr1\t600\t900\n",
        );
        let error = detect(&path, None, None).expect_err("must reject");
        assert!(
            matches!(
                error,
                AnnotationError::MixedBedWidths {
                    expected: 6,
                    actual: 3,
                    line: 2,
                    ..
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn an_empty_or_comment_only_file_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        for contents in [&b""[..], b"# nothing here\n\n"] {
            let path = write(&directory, "a.bed", contents);
            let error = detect(&path, None, None).expect_err("must reject");
            assert!(matches!(error, AnnotationError::Empty { .. }), "{error}");
        }
    }

    #[test]
    fn unclassifiable_content_with_no_extension_hint_is_ambiguous() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "mystery.dat", b"alpha\tbeta\tgamma\n");
        let error = detect(&path, None, None).expect_err("must reject");
        assert!(
            matches!(error, AnnotationError::AmbiguousFormat { .. }),
            "{error}"
        );
        assert!(error.to_string().contains("--format"), "{error}");
    }

    #[test]
    fn a_bed_type_on_gtf_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.gtf", GTF);
        let error = detect(&path, None, Some(BedType::Bed12)).expect_err("must reject");
        assert!(
            matches!(error, AnnotationError::BedTypeOnNonBed { .. }),
            "{error}"
        );
    }

    #[test]
    fn an_explicit_type_wider_than_the_data_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.bed", BED3);
        let error = detect(&path, None, Some(BedType::Bed12)).expect_err("must reject");
        assert!(
            matches!(error, AnnotationError::AmbiguousBedType { .. }),
            "{error}"
        );
    }

    #[test]
    fn an_explicit_type_narrower_than_the_data_keeps_the_rest_as_extras() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.bed", BED12);
        let detection = detect(&path, None, Some(BedType::Bed6)).expect("detected");
        assert_eq!(detection.bed_type, Some(BedType::Bed6));
        assert_eq!(detection.trailing_columns, 6);
    }

    #[test]
    fn bed_type_capabilities_match_their_widths() {
        assert!(!BedType::Bed6.has_coding_bounds());
        assert!(BedType::Bed8.has_coding_bounds());
        assert!(BedType::Bed12.has_blocks());
        assert!(!BedType::Bed9.has_blocks());
        assert!(BedType::Bed6.has_strand());
        assert!(!BedType::Bed5.has_strand());
        assert_eq!(BedType::widest_fitting(11), Some(BedType::Bed9));
        assert_eq!(BedType::widest_fitting(2), None);
        assert_eq!(BedType::parse("12"), Some(BedType::Bed12));
        assert_eq!(BedType::parse("7"), None);
    }

    #[test]
    fn a_missing_file_is_an_io_error() {
        let error = detect(Path::new("/nonexistent/a.bed"), None, None).expect_err("must fail");
        assert!(matches!(error, AnnotationError::Io { .. }), "{error}");
    }
}
