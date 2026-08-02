// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Loading an annotation through `genepred`, with a real typed dispatch.
//!
//! # No parser of our own
//!
//! `bamsplit` contains no BED, GTF, or GFF parser and should never grow one.
//! Those formats have enough dialects that a second implementation is a second
//! set of bugs, and `genepred` already normalizes all of them into one
//! `GenePred` shape.
//!
//! # Typed dispatch, not type erasure
//!
//! `genepred`'s reader is generic over the BED width — `Reader<Bed3>` through
//! `Reader<Bed12>` are distinct types with distinct field sets. The detected
//! width therefore selects a concrete reader in a `match`. There is no
//! `Box<dyn Any>`, no transmute, and no "parse as BED12 and hope": a BED6 file
//! read as BED12 would be a parse error, and a BED12 file read as BED6 would
//! silently discard the block structure that `--feature exon` depends on.
//!
//! # Normalizing compressed input
//!
//! `genepred` picks its decompressor from the **file extension**, and its
//! GTF/GFF path only accepts a path, not a reader. So a compressed file whose
//! name does not advertise its codec — a gzip stream called `genes.bed`, or a
//! `.bgz` that `genepred` has no mapping for — is decompressed into a temporary
//! file with a name `genepred` understands. A correctly named compressed file is
//! streamed directly and never copied.

use std::path::{Path, PathBuf};

use genepred::{Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, Bed12, GenePred, Gff, Gtf, Reader};

use crate::annotation::detect::{AnnotationFormat, BedType, Compression, Detection, detect};
use crate::error::AnnotationError;

/// How an annotation is loaded.
#[derive(Debug, Clone, Default)]
pub struct LoadOptions {
    /// Override the detected format.
    pub format: Option<AnnotationFormat>,
    /// Override the detected BED width.
    pub bed_type: Option<BedType>,
    /// Prefer a memory mapping, where the codec allows one.
    pub use_mmap: bool,
}

/// A loaded annotation.
#[derive(Debug)]
pub struct LoadedAnnotations {
    /// Every record, in input order.
    pub records: Vec<GenePred>,
    /// What detection concluded.
    pub detection: Detection,
}

impl LoadedAnnotations {
    /// How many records were read.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the annotation is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

/// Reads an annotation into memory.
///
/// The whole annotation is materialized because the interval index needs it: a
/// per-reference sorted structure cannot be built from a stream. Annotations are
/// small next to BAMs — a full GENCODE GTF is a few hundred thousand records —
/// so this is bounded and predictable.
///
/// # Errors
///
/// Returns [`AnnotationError`] for a detection failure, a `genepred` parse
/// error, an empty file, or an I/O failure.
pub fn load(path: &Path, options: &LoadOptions) -> Result<LoadedAnnotations, AnnotationError> {
    let detection = detect(path, options.format, options.bed_type)?;
    let source = NormalizedSource::prepare(path, detection.compression)?;
    let records = read_records(path, source.path(), &detection, options)?;

    if records.is_empty() {
        return Err(AnnotationError::Empty {
            path: path.to_path_buf(),
        });
    }
    Ok(LoadedAnnotations { records, detection })
}

/// Reads every record with the reader the detection selected.
fn read_records(
    original: &Path,
    source: &Path,
    detection: &Detection,
    options: &LoadOptions,
) -> Result<Vec<GenePred>, AnnotationError> {
    // A mapping cannot be combined with decompression, and `genepred` rejects
    // the pair explicitly, so honour `use_mmap` only for plain input.
    let use_mmap = options.use_mmap && !detection.compression.is_compressed();

    match detection.format {
        AnnotationFormat::Gtf => collect(original, gxf_reader::<Gtf>(original, source, use_mmap)?),
        AnnotationFormat::Gff => collect(original, gxf_reader::<Gff>(original, source, use_mmap)?),
        AnnotationFormat::Bed => {
            let extra = detection.trailing_columns;
            // The one place the width actually matters. Each arm instantiates a
            // different concrete reader with a different field set.
            match detection.bed_type.unwrap_or(BedType::Bed3) {
                BedType::Bed3 => collect(
                    original,
                    bed_reader::<Bed3>(original, source, extra, use_mmap)?,
                ),
                BedType::Bed4 => collect(
                    original,
                    bed_reader::<Bed4>(original, source, extra, use_mmap)?,
                ),
                BedType::Bed5 => collect(
                    original,
                    bed_reader::<Bed5>(original, source, extra, use_mmap)?,
                ),
                BedType::Bed6 => collect(
                    original,
                    bed_reader::<Bed6>(original, source, extra, use_mmap)?,
                ),
                BedType::Bed8 => collect(
                    original,
                    bed_reader::<Bed8>(original, source, extra, use_mmap)?,
                ),
                BedType::Bed9 => collect(
                    original,
                    bed_reader::<Bed9>(original, source, extra, use_mmap)?,
                ),
                BedType::Bed12 => collect(
                    original,
                    bed_reader::<Bed12>(original, source, extra, use_mmap)?,
                ),
            }
        }
    }
}

fn bed_reader<F>(
    original: &Path,
    source: &Path,
    additional_fields: usize,
    use_mmap: bool,
) -> Result<Reader<F>, AnnotationError>
where
    F: genepred::BedFormat + Into<GenePred>,
{
    let options = genepred::ReaderOptions::new().additional_fields(additional_fields);
    let reader = if use_mmap {
        Reader::<F>::from_mmap_with_custom_fields(source, options)
    } else {
        Reader::<F>::from_path_with_custom_fields(source, options)
    };
    reader.map_err(|error| parse_error(original, &error))
}

/// Builds an aggregating GTF or GFF reader.
///
/// `genepred` exposes these as inherent methods on the concrete `Reader<Gtf>`
/// and `Reader<Gff>` rather than through a trait, so the two are bridged here
/// with a small sealed trait rather than duplicating the call site.
fn gxf_reader<F>(
    original: &Path,
    source: &Path,
    use_mmap: bool,
) -> Result<Reader<F>, AnnotationError>
where
    F: GxfReadable,
{
    F::open(source, use_mmap).map_err(|error| parse_error(original, &error))
}

/// The two aggregating formats, bridged so one call site can serve both.
///
/// `genepred` exposes GTF and GFF aggregation as inherent methods on the
/// concrete `Reader<Gtf>` and `Reader<Gff>` rather than through a trait, so
/// this supplies the trait it does not.
pub trait GxfReadable: genepred::BedFormat + Into<GenePred> + Sized {
    /// Opens an aggregating reader over `path`.
    ///
    /// # Errors
    ///
    /// Returns whatever `genepred` reports.
    fn open(path: &Path, use_mmap: bool) -> genepred::ReaderResult<Reader<Self>>;
}

impl GxfReadable for Gtf {
    fn open(path: &Path, use_mmap: bool) -> genepred::ReaderResult<Reader<Self>> {
        if use_mmap {
            Reader::<Gtf>::from_mmap_with_options(path, genepred::ReaderOptions::new())
        } else {
            Reader::<Gtf>::from_gxf(path)
        }
    }
}

impl GxfReadable for Gff {
    fn open(path: &Path, use_mmap: bool) -> genepred::ReaderResult<Reader<Self>> {
        if use_mmap {
            Reader::<Gff>::from_mmap_with_options(path, genepred::ReaderOptions::new())
        } else {
            Reader::<Gff>::from_gxf(path)
        }
    }
}

fn collect<F>(original: &Path, mut reader: Reader<F>) -> Result<Vec<GenePred>, AnnotationError>
where
    F: genepred::BedFormat + Into<GenePred>,
{
    let mut records = Vec::new();
    for record in reader.records() {
        records.push(record.map_err(|error| parse_error(original, &error))?);
    }
    Ok(records)
}

fn parse_error(path: &Path, error: &genepred::reader::ReaderError) -> AnnotationError {
    AnnotationError::Parse {
        path: path.to_path_buf(),
        message: error.to_string(),
    }
}

/// A path `genepred` can decompress, plus the temporary file backing it when one
/// had to be created.
struct NormalizedSource {
    path: PathBuf,
    // Held only for its `Drop`: the temporary file must outlive the reader.
    _temporary: Option<tempfile::TempPath>,
}

impl NormalizedSource {
    fn prepare(path: &Path, compression: Compression) -> Result<Self, AnnotationError> {
        let Some(suffix) = compression.genepred_suffix() else {
            return Ok(Self {
                path: path.to_path_buf(),
                _temporary: None,
            });
        };
        if advertises(path, suffix) {
            return Ok(Self {
                path: path.to_path_buf(),
                _temporary: None,
            });
        }

        // The name does not tell `genepred` how to decode this, and its GTF/GFF
        // path takes a path rather than a reader, so give it one it understands.
        let io = |source| AnnotationError::Io {
            path: path.to_path_buf(),
            source,
        };
        let stem = path.file_name().map_or_else(
            || "annotation".to_string(),
            |name| name.to_string_lossy().into_owned(),
        );
        let temporary = tempfile::Builder::new()
            .prefix("bamsplit-annotation-")
            .suffix(&format!("-{stem}.{suffix}"))
            .tempfile()
            .map_err(io)?;
        std::fs::copy(path, temporary.path()).map_err(io)?;

        Ok(Self {
            path: temporary.path().to_path_buf(),
            _temporary: Some(temporary.into_temp_path()),
        })
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

/// Whether a filename already ends with the suffix `genepred` needs.
fn advertises(path: &Path, suffix: &str) -> bool {
    let Some(name) = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
    else {
        return false;
    };
    // The name is lower-cased first, so a plain suffix test is already
    // case-insensitive.
    let name = name.to_ascii_lowercase();
    let actual = name.rsplit_once('.').map(|(_, suffix)| suffix);
    match suffix {
        "gz" => actual == Some("gz"),
        "zst" => matches!(actual, Some("zst" | "zstd")),
        "bz2" => matches!(actual, Some("bz2" | "bzip2")),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    const BED12: &[u8] = b"chr1\t100\t500\tgene1\t0\t+\t150\t450\t0,0,0\t2\t100,100\t0,300\n\
                           chr1\t600\t900\tgene2\t0\t-\t650\t850\t0,0,0\t1\t300\t0\n";
    const BED6: &[u8] = b"chr1\t100\t500\tgene1\t0\t+\nchr2\t600\t900\tgene2\t0\t-\n";
    const GTF: &[u8] = b"chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
                         chr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
                         chr1\ttest\tCDS\t151\t180\t.\t+\t0\tgene_id \"g1\"; transcript_id \"t1\";\n";
    const GFF3: &[u8] = b"##gff-version 3\n\
                          chr1\ttest\texon\t101\t200\t.\t+\t.\tID=e1;Parent=t1\n\
                          chr1\ttest\texon\t301\t400\t.\t+\t.\tID=e2;Parent=t1\n";

    fn write(directory: &tempfile::TempDir, name: &str, contents: &[u8]) -> PathBuf {
        let path = directory.path().join(name);
        std::fs::write(&path, contents).expect("written");
        path
    }

    fn gzip(contents: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(contents).expect("compressed");
        encoder.finish().expect("finished")
    }

    fn bzip2(contents: &[u8]) -> Vec<u8> {
        let mut encoder = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::fast());
        encoder.write_all(contents).expect("compressed");
        encoder.finish().expect("finished")
    }

    #[test]
    fn loads_bed12_with_block_structure() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.bed", BED12);
        let loaded = load(&path, &LoadOptions::default()).expect("loaded");

        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded.detection.bed_type, Some(BedType::Bed12));
        let first = &loaded.records[0];
        assert_eq!(first.chrom, b"chr1");
        assert_eq!(first.start, 100);
        assert_eq!(first.end, 500);
        assert_eq!(first.name.as_deref(), Some(&b"gene1"[..]));
        assert_eq!(first.exons(), vec![(100, 200), (400, 500)]);
        assert_eq!(first.thick_start, Some(150));
        assert_eq!(first.thick_end, Some(450));
    }

    #[test]
    fn loads_every_supported_bed_width() {
        let directory = tempfile::tempdir().expect("temp dir");
        let rows: [(&str, &[u8], BedType); 7] = [
            ("b3.bed", b"chr1\t100\t200\n", BedType::Bed3),
            ("b4.bed", b"chr1\t100\t200\tn\n", BedType::Bed4),
            ("b5.bed", b"chr1\t100\t200\tn\t0\n", BedType::Bed5),
            ("b6.bed", b"chr1\t100\t200\tn\t0\t+\n", BedType::Bed6),
            (
                "b8.bed",
                b"chr1\t100\t200\tn\t0\t+\t120\t180\n",
                BedType::Bed8,
            ),
            (
                "b9.bed",
                b"chr1\t100\t200\tn\t0\t+\t120\t180\t0,0,0\n",
                BedType::Bed9,
            ),
            (
                "b12.bed",
                b"chr1\t100\t200\tn\t0\t+\t120\t180\t0,0,0\t1\t100\t0\n",
                BedType::Bed12,
            ),
        ];
        for (name, contents, expected) in rows {
            let path = write(&directory, name, contents);
            let loaded = load(&path, &LoadOptions::default())
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            assert_eq!(loaded.detection.bed_type, Some(expected), "{name}");
            assert_eq!(loaded.len(), 1, "{name}");
            assert_eq!(loaded.records[0].start, 100);
        }
    }

    #[test]
    fn narrower_widths_do_not_invent_fields() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.bed", BED6);
        let loaded = load(&path, &LoadOptions::default()).expect("loaded");
        assert!(loaded.records[0].thick_start.is_none());
        assert!(loaded.records[0].block_starts.is_none());
        assert!(loaded.records[0].strand.is_some());
    }

    #[test]
    fn an_explicit_narrower_type_discards_the_extra_fields() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.bed", BED12);
        let loaded = load(
            &path,
            &LoadOptions {
                bed_type: Some(BedType::Bed6),
                ..LoadOptions::default()
            },
        )
        .expect("loaded");
        assert_eq!(loaded.detection.bed_type, Some(BedType::Bed6));
        // Read as BED6, the block columns are extras rather than structure.
        assert!(loaded.records[0].block_starts.is_none());
    }

    #[test]
    fn aggregates_gtf_into_transcripts() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.gtf", GTF);
        let loaded = load(&path, &LoadOptions::default()).expect("loaded");
        assert_eq!(loaded.detection.format, AnnotationFormat::Gtf);
        assert_eq!(loaded.len(), 1, "two exons of one transcript aggregate");
        let record = &loaded.records[0];
        assert_eq!(record.chrom, b"chr1");
        // GTF is 1-based inclusive; genepred normalizes to 0-based half-open.
        assert_eq!(record.start, 100);
        assert_eq!(record.end, 400);
        assert_eq!(record.exons(), vec![(100, 200), (300, 400)]);
    }

    #[test]
    fn aggregates_gff3_into_transcripts() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.gff3", GFF3);
        let loaded = load(&path, &LoadOptions::default()).expect("loaded");
        assert_eq!(loaded.detection.format, AnnotationFormat::Gff);
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded.records[0].exons(), vec![(100, 200), (300, 400)]);
    }

    #[test]
    fn reads_gzip_regardless_of_the_extension() {
        let directory = tempfile::tempdir().expect("temp dir");
        // Correctly named: streamed directly.
        let named = write(&directory, "a.bed.gz", &gzip(BED12));
        let loaded = load(&named, &LoadOptions::default()).expect("loaded");
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded.detection.compression, Compression::Gzip);

        // Misleadingly named: normalized through a temporary copy.
        let misnamed = write(&directory, "b.bed", &gzip(BED12));
        let loaded = load(&misnamed, &LoadOptions::default()).expect("loaded");
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded.detection.compression, Compression::Gzip);
    }

    #[test]
    fn reads_bzip2() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.bed.bz2", &bzip2(BED12));
        let loaded = load(&path, &LoadOptions::default()).expect("loaded");
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded.detection.compression, Compression::Bzip2);
    }

    #[test]
    fn reads_a_compressed_gtf() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.gtf.gz", &gzip(GTF));
        let loaded = load(&path, &LoadOptions::default()).expect("loaded");
        assert_eq!(loaded.detection.format, AnnotationFormat::Gtf);
        assert_eq!(loaded.len(), 1);
    }

    #[test]
    fn a_malformed_row_is_a_parse_error_naming_the_original_path() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(
            &directory,
            "a.bed",
            b"chr1\t100\t200\nchr1\tNOTANUMBER\t300\n",
        );
        let error = load(&path, &LoadOptions::default()).expect_err("must reject");
        match &error {
            AnnotationError::Parse {
                path: reported,
                message,
            } => {
                assert_eq!(reported, &path);
                assert!(!message.is_empty());
            }
            other => panic!("expected a parse error, got {other}"),
        }
    }

    #[test]
    fn an_empty_annotation_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(&directory, "a.bed", b"# only a comment\n");
        let error = load(&path, &LoadOptions::default()).expect_err("must reject");
        assert!(matches!(error, AnnotationError::Empty { .. }), "{error}");
    }

    #[test]
    fn mmap_is_used_only_for_uncompressed_input() {
        let directory = tempfile::tempdir().expect("temp dir");
        let plain = write(&directory, "a.bed", BED12);
        let loaded = load(
            &plain,
            &LoadOptions {
                use_mmap: true,
                ..LoadOptions::default()
            },
        )
        .expect("loaded");
        assert_eq!(loaded.len(), 2);

        // Asking for a mapping over a compressed file must not fail; the loader
        // silently falls back, because `genepred` cannot do both at once.
        let compressed = write(&directory, "b.bed.gz", &gzip(BED12));
        let loaded = load(
            &compressed,
            &LoadOptions {
                use_mmap: true,
                ..LoadOptions::default()
            },
        )
        .expect("loaded");
        assert_eq!(loaded.len(), 2);
    }

    #[test]
    fn record_order_matches_input_order() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = write(
            &directory,
            "a.bed",
            b"chr1\t300\t400\tthird\nchr1\t100\t200\tfirst\nchr1\t200\t300\tsecond\n",
        );
        let loaded = load(&path, &LoadOptions::default()).expect("loaded");
        let names: Vec<&[u8]> = loaded
            .records
            .iter()
            .filter_map(|record| record.name.as_deref())
            .collect();
        assert_eq!(names, [&b"third"[..], b"first", b"second"]);
    }
}
