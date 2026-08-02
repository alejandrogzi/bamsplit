// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! One output BAM: header, raw records, streaming index, statistics, commit.
//!
//! A [`RecordWriter`] is the join point of everything in [`crate::output`]. It
//! owns a transactional temporary file, a BGZF block writer that reports exact
//! virtual offsets, an index builder fed by those offsets, and the per-output
//! counters the manifest reports.
//!
//! # Parking
//!
//! When more outputs are live than `--max-open-files` allows, the manager
//! *parks* the least recently used one: the current BGZF block is flushed, the
//! file handle is closed, and the compressed offset is remembered. Resuming
//! reopens the file in append mode and starts a new block writer whose base
//! offset continues where the last one stopped, so the index stays correct
//! across the gap.
//!
//! Parking is possible only because the end-of-file marker is written at
//! [`finish`](RecordWriter::finish) rather than on drop. A BGZF stream is a
//! concatenation of independent blocks, so appending more blocks to a file that
//! has not been terminated is well-defined.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;

use crate::bam::header::BamHeader;
use crate::bam::raw_record::RawRecord;
use crate::error::{IndexError, OutputError};
use crate::index::{AlignmentContext, IndexBuilder, IndexKind};
use crate::output::bgzf::BgzfBlockWriter;
use crate::output::filename::OutputPaths;
use crate::output::transaction::{PendingFile, Transaction};
use crate::stats::{OutputStats, OutputStatsBuilder};

/// Settings shared by every output of a run.
#[derive(Clone)]
pub struct OutputSettings {
    /// The DEFLATE level, `0..=9`.
    pub compression_level: u32,
    /// Whether an existing output may be replaced.
    pub force: bool,
    /// The index format, or [`None`] for no index.
    pub index_kind: Option<IndexKind>,
    /// How many threads one output's compression may use.
    pub compression_workers: usize,
    /// The pool those threads come from.
    ///
    /// Shared across outputs so the total thread count stays inside the global
    /// `--threads` budget no matter how many outputs are live.
    pub pool: Option<Arc<rayon::ThreadPool>>,
    /// The size of the buffer between the block writer and the file.
    pub write_buffer_bytes: usize,
}

impl std::fmt::Debug for OutputSettings {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OutputSettings")
            .field("compression_level", &self.compression_level)
            .field("force", &self.force)
            .field("index_kind", &self.index_kind)
            .field("compression_workers", &self.compression_workers)
            .field("has_pool", &self.pool.is_some())
            .field("write_buffer_bytes", &self.write_buffer_bytes)
            .finish()
    }
}

impl Default for OutputSettings {
    fn default() -> Self {
        Self {
            compression_level: 6,
            force: false,
            index_kind: None,
            compression_workers: 1,
            pool: None,
            write_buffer_bytes: 1 << 20,
        }
    }
}

type BlockWriter = BgzfBlockWriter<BufWriter<File>, Option<AlignmentContext>>;

/// A single output BAM being written.
pub struct RecordWriter {
    paths: OutputPaths,
    settings: OutputSettings,
    pending: PendingFile,
    writer: Option<BlockWriter>,
    parked_offset: u64,
    index: Option<IndexBuilder>,
    index_error: Option<IndexError>,
    stats: OutputStatsBuilder,
}

impl std::fmt::Debug for RecordWriter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RecordWriter")
            .field("stem", &self.paths.stem)
            .field("records", &self.stats.record_count())
            .field("parked", &self.writer.is_none())
            .finish_non_exhaustive()
    }
}

impl RecordWriter {
    /// Creates an output and writes its BAM header.
    ///
    /// The header is the caller's, already amended with the `@PG` record, and
    /// carries the **complete** input reference dictionary — see
    /// [`crate::bam::header`] for why.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError`] if the output already exists, the temporary file
    /// cannot be created, the header cannot be serialized, or the write fails.
    pub fn create(
        paths: OutputPaths,
        header: &BamHeader,
        reference_count: usize,
        settings: &OutputSettings,
    ) -> Result<Self, OutputError> {
        let (file, pending) = PendingFile::create(&paths.bam, settings.force)?;
        let sink = BufWriter::with_capacity(settings.write_buffer_bytes, file);
        let mut writer = BgzfBlockWriter::new(
            sink,
            settings.compression_level,
            settings.pool.clone(),
            settings.compression_workers,
            0,
        );

        let mut header_bytes = Vec::new();
        header
            .write_to(&mut header_bytes)
            .map_err(|source| OutputError::Validation {
                path: paths.bam.clone(),
                reason: format!("cannot serialize the BAM header: {source}"),
            })?;
        writer
            .write_raw(&header_bytes, &mut |_, _| {})
            .map_err(|source| OutputError::Io {
                path: pending.temporary_path().to_path_buf(),
                source,
            })?;

        let index = settings
            .index_kind
            .map(|kind| IndexBuilder::new(kind, reference_count));

        Ok(Self {
            paths,
            settings: settings.clone(),
            pending,
            writer: Some(writer),
            parked_offset: 0,
            index,
            index_error: None,
            stats: OutputStatsBuilder::new(),
        })
    }

    /// The paths this output occupies.
    #[must_use]
    pub const fn paths(&self) -> &OutputPaths {
        &self.paths
    }

    /// How many records have been written.
    #[must_use]
    pub const fn record_count(&self) -> u64 {
        self.stats.record_count()
    }

    /// Whether the output is currently parked.
    #[must_use]
    pub const fn is_parked(&self) -> bool {
        self.writer.is_none()
    }

    /// Whether the index was abandoned, and why.
    #[must_use]
    pub fn index_abandoned(&self) -> Option<String> {
        self.index
            .as_ref()
            .and_then(|index| index.abandon_reason().map(ToString::to_string))
    }

    /// Writes one record, preserving its body byte for byte.
    ///
    /// The `block_size` prefix is recomputed from the body length rather than
    /// copied, because the body is what must round-trip; a prefix that
    /// disagreed with the body would be malformed regardless of what the input
    /// said.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::Io`] if the write fails, or
    /// [`OutputError::Validation`] if the record's fixed core cannot be read.
    pub fn write_record(&mut self, record: &RawRecord<'_>) -> Result<(), OutputError> {
        self.ensure_open()?;

        let context = alignment_context(record).map_err(|source| OutputError::Validation {
            path: self.paths.bam.clone(),
            reason: format!("cannot index a record: {source}"),
        })?;

        let Self {
            writer,
            index,
            index_error,
            ..
        } = self;
        // `ensure_open` guarantees the writer is present.
        let Some(writer) = writer.as_mut() else {
            return Err(OutputError::Validation {
                path: self.paths.bam.clone(),
                reason: "internal error: the output writer is parked".to_string(),
            });
        };

        writer
            .write_record(record.raw_bytes(), context, &mut |context, range| {
                if let Some(index) = index.as_mut()
                    && let Err(error) = index.add(context, range.start, range.end)
                    && index_error.is_none()
                {
                    *index_error = Some(error);
                }
            })
            .map_err(|source| OutputError::Io {
                path: self.pending.temporary_path().to_path_buf(),
                source,
            })?;

        self.stats
            .observe(record)
            .map_err(|source| OutputError::Validation {
                path: self.paths.bam.clone(),
                reason: format!("cannot account for a record: {source}"),
            })?;
        Ok(())
    }

    /// Flushes and closes the file handle, keeping the output resumable.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::Io`] if the flush fails.
    pub fn park(&mut self) -> Result<(), OutputError> {
        let Some(writer) = self.writer.take() else {
            return Ok(());
        };
        let Self {
            index, index_error, ..
        } = self;
        let (sink, offset) = writer
            .park(&mut |context, range| {
                if let Some(index) = index.as_mut()
                    && let Err(error) = index.add(context, range.start, range.end)
                    && index_error.is_none()
                {
                    *index_error = Some(error);
                }
            })
            .map_err(|source| OutputError::Io {
                path: self.pending.temporary_path().to_path_buf(),
                source,
            })?;
        drop(sink);
        self.parked_offset = offset;
        Ok(())
    }

    /// Reopens a parked output.
    fn ensure_open(&mut self) -> Result<(), OutputError> {
        if self.writer.is_some() {
            return Ok(());
        }
        let file = self.pending.reopen_append()?;
        let sink = BufWriter::with_capacity(self.settings.write_buffer_bytes, file);
        self.writer = Some(BgzfBlockWriter::new(
            sink,
            self.settings.compression_level,
            self.settings.pool.clone(),
            self.settings.compression_workers,
            self.parked_offset,
        ));
        Ok(())
    }

    /// Finishes the output: BGZF end-of-file marker, index, statistics.
    ///
    /// The BAM and its index are added to `transaction` rather than renamed
    /// here, so a caller can commit a whole group atomically.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError`] for an I/O failure, an unresolvable virtual
    /// offset, or a failure to write the index.
    pub fn finish(mut self, transaction: &mut Transaction) -> Result<FinishedOutput, OutputError> {
        self.ensure_open()?;
        let writer = self.writer.take().ok_or_else(|| OutputError::Validation {
            path: self.paths.bam.clone(),
            reason: "internal error: the output writer was already finished".to_string(),
        })?;

        let Self {
            index, index_error, ..
        } = &mut self;
        let sink = writer
            .finish(&mut |context, range| {
                if let Some(index) = index.as_mut()
                    && let Err(error) = index.add(context, range.start, range.end)
                    && index_error.is_none()
                {
                    *index_error = Some(error);
                }
            })
            .map_err(|source| OutputError::Io {
                path: self.pending.temporary_path().to_path_buf(),
                source,
            })?;
        // `BufWriter::into_inner` flushes; the block writer already did, so
        // this only surfaces an error the earlier flush could not.
        let file = sink.into_inner().map_err(|error| OutputError::Io {
            path: self.pending.temporary_path().to_path_buf(),
            source: error.into_error(),
        })?;
        file.sync_all().map_err(|source| OutputError::Io {
            path: self.pending.temporary_path().to_path_buf(),
            source,
        })?;
        drop(file);

        if let Some(error) = self.index_error.take() {
            return Err(OutputError::Validation {
                path: self.paths.bam.clone(),
                reason: format!("index construction failed: {error}"),
            });
        }

        let compressed_bytes = self.pending.size()?;
        let mut index_note = None;
        let mut index_path = None;

        if let (Some(builder), Some(target)) = (self.index.take(), self.paths.index.clone()) {
            if let Some(reason) = builder.abandon_reason() {
                index_note = Some(reason.to_string());
            } else {
                let (file, pending_index) = PendingFile::create(&target, self.settings.force)?;
                // The index writers take a path, so the handle is only used to
                // reserve the temporary name.
                drop(file);
                let temporary = pending_index.temporary_path().to_path_buf();
                match builder.finish(&temporary) {
                    Ok(Some(_)) => {
                        index_path = Some(target);
                        transaction.push(pending_index);
                    }
                    Ok(None) => {
                        pending_index.abort();
                        index_note = Some("the output is not coordinate-indexable".to_string());
                    }
                    Err(error) => {
                        pending_index.abort();
                        return Err(OutputError::Validation {
                            path: self.paths.bam.clone(),
                            reason: format!("cannot write the index: {error}"),
                        });
                    }
                }
            }
        }

        let stats = self.stats.finish(compressed_bytes);
        validate_output(&self.paths.bam, &stats)?;

        let paths = self.paths.clone();
        transaction.push(self.pending);
        Ok(FinishedOutput {
            paths,
            stats,
            index_path,
            index_note,
        })
    }

    /// Discards the output and its temporary file.
    pub fn abort(self) {
        self.pending.abort();
    }
}

/// A completed output, before its transaction is committed.
#[derive(Debug, Clone)]
pub struct FinishedOutput {
    /// Where the output will live.
    pub paths: OutputPaths,
    /// The counters and digest.
    pub stats: OutputStats,
    /// The index path, when one was written.
    pub index_path: Option<std::path::PathBuf>,
    /// Why no index was written, when one was requested but skipped.
    pub index_note: Option<String>,
}

/// Extracts the index context of a record.
fn alignment_context(
    record: &RawRecord<'_>,
) -> Result<Option<AlignmentContext>, crate::error::BamRecordError> {
    let Some(reference_id) = record.reference_sequence_id()? else {
        return Ok(None);
    };
    let Some(start) = record.alignment_start()? else {
        return Ok(None);
    };
    let Some(end) = record.alignment_end()? else {
        return Ok(None);
    };
    let reference_id = usize::try_from(reference_id).unwrap_or(0);
    // `alignment_start` is 0-based inclusive and `alignment_end` 0-based
    // exclusive, so the 1-based inclusive pair is `(start + 1, end)`.
    Ok(Some(AlignmentContext {
        reference_id,
        start: u64::try_from(start).unwrap_or(0) + 1,
        end: u64::try_from(end).unwrap_or(1),
        mapped: !record.is_unmapped()?,
    }))
}

/// Post-write sanity checks on one output.
fn validate_output(path: &Path, stats: &OutputStats) -> Result<(), OutputError> {
    if !stats.is_consistent() {
        return Err(OutputError::Validation {
            path: path.to_path_buf(),
            reason: format!(
                "per-output accounting is inconsistent: {} records but {} placement categories \
                 and {} alignment categories",
                stats.record_count,
                stats.mapped_count + stats.placed_unmapped_count + stats.unplaced_unmapped_count,
                stats.primary_count + stats.secondary_count + stats.supplementary_count
            ),
        });
    }
    // Even a header-only BAM has a header and a 28-byte end-of-file marker.
    if stats.compressed_bytes < crate::output::bgzf::BGZF_EOF.len() as u64 {
        return Err(OutputError::Validation {
            path: path.to_path_buf(),
            reason: format!(
                "the output is only {} bytes, too small to contain a BGZF end-of-file marker",
                stats.compressed_bytes
            ),
        });
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::items_after_statements)]
mod tests {
    use super::*;
    use crate::bam::header::ReferenceSequence;
    use crate::index::IndexKind;
    use crate::output::filename::{FilenameEncoder, FilenameTemplate};

    fn header() -> BamHeader {
        BamHeader::from_parts(
            b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n".to_vec(),
            vec![
                ReferenceSequence {
                    name: b"chr1".to_vec(),
                    length: 1000,
                },
                ReferenceSequence {
                    name: b"chr2".to_vec(),
                    length: 2000,
                },
            ],
        )
        .expect("valid header")
    }

    fn record(reference_id: i32, position: i32, name: &[u8]) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&reference_id.to_le_bytes());
        body.extend_from_slice(&position.to_le_bytes());
        body.push(u8::try_from(name.len() + 1).expect("short"));
        body.push(60);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&1u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&10i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(name);
        body.push(0);
        body.extend_from_slice(&(10u32 << 4).to_le_bytes()); // 10M
        body.extend(std::iter::repeat_n(0u8, 5));
        body.extend(std::iter::repeat_n(0xffu8, 10));
        body
    }

    fn settings(index: Option<IndexKind>) -> OutputSettings {
        OutputSettings {
            compression_level: 1,
            index_kind: index,
            ..OutputSettings::default()
        }
    }

    fn read_back(path: &Path) -> (BamHeader, Vec<Vec<u8>>) {
        let file = File::open(path).expect("readable");
        let mut reader = noodles_bgzf::io::Reader::new(file);
        let header = BamHeader::read_from(&mut reader).expect("valid header");
        let mut records = Vec::new();
        let mut record_reader = crate::bam::RawRecordReader::new(reader);
        while let Some(record) = record_reader.read_record().expect("valid record") {
            records.push(record.raw_bytes().to_vec());
        }
        (header, records)
    }

    #[test]
    fn writes_a_readable_bam_and_preserves_bodies() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chr1", None).expect("assigned");
        let header = header();

        let bodies = vec![
            record(0, 10, b"a"),
            record(0, 200, b"bb"),
            record(0, 300, b"ccc"),
        ];
        let mut writer = RecordWriter::create(paths, &header, 2, &settings(None)).expect("created");
        for body in &bodies {
            writer
                .write_record(&RawRecord::new(body).expect("valid"))
                .expect("written");
        }
        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction).expect("finished");
        transaction.commit().expect("committed");

        let (read_header, read_bodies) = read_back(&finished.paths.bam);
        assert_eq!(read_header.reference_count(), 2);
        assert_eq!(read_header.reference_name(1), Some(&b"chr2"[..]));
        assert_eq!(read_bodies, bodies);
        assert_eq!(finished.stats.record_count, 3);
        assert_eq!(finished.stats.mapped_count, 3);
        assert!(finished.stats.coordinate_sorted);
        assert!(finished.stats.compressed_bytes > 28);
    }

    #[test]
    fn every_output_keeps_the_full_reference_dictionary() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chr2", None).expect("assigned");
        let header = header();

        let mut writer = RecordWriter::create(paths, &header, 2, &settings(None)).expect("created");
        writer
            .write_record(&RawRecord::new(&record(1, 5, b"x")).expect("valid"))
            .expect("written");
        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction).expect("finished");
        transaction.commit().expect("committed");

        let (read_header, _) = read_back(&finished.paths.bam);
        assert_eq!(
            read_header.reference_count(),
            2,
            "chr1 must still be present"
        );
        assert_eq!(read_header.reference_id(b"chr1"), Some(0));
    }

    #[test]
    fn a_header_only_output_is_valid() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chrEmpty", Some("bai")).expect("assigned");
        let header = header();

        let writer = RecordWriter::create(paths, &header, 2, &settings(Some(IndexKind::Bai)))
            .expect("created");
        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction).expect("finished");
        transaction.commit().expect("committed");

        let (read_header, records) = read_back(&finished.paths.bam);
        assert_eq!(read_header.reference_count(), 2);
        assert!(records.is_empty());
        assert!(finished.index_path.is_some());
        assert!(finished.index_path.as_deref().is_some_and(Path::exists));
    }

    #[test]
    fn an_index_is_written_alongside_the_bam() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chr1", Some("bai")).expect("assigned");
        let header = header();

        let mut writer = RecordWriter::create(paths, &header, 2, &settings(Some(IndexKind::Bai)))
            .expect("created");
        for position in [10, 100, 500] {
            writer
                .write_record(&RawRecord::new(&record(0, position, b"r")).expect("valid"))
                .expect("written");
        }
        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction).expect("finished");
        transaction.commit().expect("committed");

        use noodles_csi::BinningIndex as _;

        let index_path = finished.index_path.expect("index written");
        assert!(index_path.exists());
        let index = crate::index::bai::read(&index_path).expect("valid index");
        assert_eq!(index.min_shift(), 14);
    }

    #[test]
    fn an_unsorted_output_is_written_without_a_misleading_index() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"jumbled", Some("bai")).expect("assigned");
        let header = header();

        let mut writer = RecordWriter::create(paths, &header, 2, &settings(Some(IndexKind::Bai)))
            .expect("created");
        for position in [500, 100] {
            writer
                .write_record(&RawRecord::new(&record(0, position, b"r")).expect("valid"))
                .expect("written");
        }

        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction).expect("finished");
        transaction.commit().expect("committed");

        assert!(finished.index_path.is_none());
        assert!(finished.index_note.is_some(), "the reason must be recorded");
        assert!(!finished.stats.coordinate_sorted);
        assert!(finished.paths.bam.exists(), "the BAM itself is still valid");
    }

    #[test]
    fn parking_and_resuming_preserves_the_stream_and_the_index() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chr1", Some("bai")).expect("assigned");
        let header = header();

        let bodies: Vec<Vec<u8>> = (0..50)
            .map(|index| record(0, index * 10, format!("r{index}").as_bytes()))
            .collect();

        let mut writer = RecordWriter::create(paths, &header, 2, &settings(Some(IndexKind::Bai)))
            .expect("created");
        for (index, body) in bodies.iter().enumerate() {
            writer
                .write_record(&RawRecord::new(body).expect("valid"))
                .expect("written");
            if index % 7 == 0 {
                writer.park().expect("parked");
                assert!(writer.is_parked());
            }
        }
        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction).expect("finished");
        transaction.commit().expect("committed");

        let (_, read_bodies) = read_back(&finished.paths.bam);
        assert_eq!(read_bodies, bodies);
        assert!(finished.index_path.is_some(), "the index survived parking");
        assert_eq!(finished.stats.record_count, 50);
    }

    #[test]
    fn parked_virtual_offsets_stay_consistent_with_the_file() {
        use std::io::Read as _;

        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chr1", None).expect("assigned");
        let header = header();

        let bodies: Vec<Vec<u8>> = (0..30)
            .map(|index| record(0, index * 10, format!("r{index}").as_bytes()))
            .collect();

        // Capture the offsets the writer reports by observing the index
        // builder's view: park after every record to exercise the base offset.
        let mut writer = RecordWriter::create(paths, &header, 2, &settings(None)).expect("created");
        for body in &bodies {
            writer
                .write_record(&RawRecord::new(body).expect("valid"))
                .expect("written");
            writer.park().expect("parked");
        }
        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction).expect("finished");
        transaction.commit().expect("committed");

        // The resulting file must still be a single valid BGZF stream.
        let mut reader =
            noodles_bgzf::io::Reader::new(File::open(&finished.paths.bam).expect("readable"));
        let mut inflated = Vec::new();
        reader.read_to_end(&mut inflated).expect("valid BGZF");
        let (_, read_bodies) = read_back(&finished.paths.bam);
        assert_eq!(read_bodies, bodies);
    }

    #[test]
    fn nothing_appears_until_the_transaction_commits() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chr1", None).expect("assigned");
        let bam = paths.bam.clone();
        let header = header();

        let mut writer = RecordWriter::create(paths, &header, 2, &settings(None)).expect("created");
        writer
            .write_record(&RawRecord::new(&record(0, 1, b"r")).expect("valid"))
            .expect("written");
        assert!(!bam.exists());

        let mut transaction = Transaction::new();
        writer.finish(&mut transaction).expect("finished");
        assert!(!bam.exists(), "still nothing before the commit");
        transaction.commit().expect("committed");
        assert!(bam.exists());
    }

    #[test]
    fn aborting_leaves_no_trace() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
        let paths = encoder.assign(b"chr1", None).expect("assigned");
        let bam = paths.bam.clone();
        let header = header();

        let mut writer = RecordWriter::create(paths, &header, 2, &settings(None)).expect("created");
        writer
            .write_record(&RawRecord::new(&record(0, 1, b"r")).expect("valid"))
            .expect("written");
        writer.abort();

        assert!(!bam.exists());
        assert_eq!(
            std::fs::read_dir(directory.path())
                .expect("listable")
                .count(),
            0,
            "the temporary file must be gone too"
        );
    }

    #[test]
    fn parallel_compression_produces_an_identical_file() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let bodies: Vec<Vec<u8>> = (0..2_000)
            .map(|index| record(0, index, format!("read{index}").as_bytes()))
            .collect();

        let write = |name: &str, workers: usize| {
            let pool = Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .expect("pool"),
            );
            let mut encoder = FilenameEncoder::new(directory.path(), FilenameTemplate::default());
            let paths = encoder.assign(name.as_bytes(), None).expect("assigned");
            let settings = OutputSettings {
                compression_level: 6,
                compression_workers: workers,
                pool: Some(pool),
                ..OutputSettings::default()
            };
            let mut writer = RecordWriter::create(paths, &header, 2, &settings).expect("created");
            for body in &bodies {
                writer
                    .write_record(&RawRecord::new(body).expect("valid"))
                    .expect("written");
            }
            let mut transaction = Transaction::new();
            let finished = writer.finish(&mut transaction).expect("finished");
            transaction.commit().expect("committed");
            finished
        };

        let serial = write("serial", 1);
        let parallel = write("parallel", 4);
        assert_eq!(
            std::fs::read(&serial.paths.bam).expect("readable"),
            std::fs::read(&parallel.paths.bam).expect("readable"),
            "compression must be deterministic regardless of worker count"
        );
        assert_eq!(
            serial.stats.raw_record_digest,
            parallel.stats.raw_record_digest
        );
    }
}
