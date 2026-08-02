// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! A bounded, deterministic set of output BAMs.
//!
//! # What is bounded
//!
//! * **File descriptors.** At most `max_open_files` outputs hold an open handle.
//!   Beyond that the least recently used one is parked: its current BGZF block
//!   is flushed and its handle closed, and it resumes transparently on the next
//!   write. Parking is cheap — one block flush — and it is what lets `bamsplit
//!   tag` survive a 1 024-descriptor `ulimit` with 20 000 barcodes.
//! * **Memory.** A parked output keeps only its path, counters, and a partially
//!   built index. The index is the one part that can grow, which is why the
//!   spool engine exists for genuinely high-cardinality splits.
//!
//! # What is deterministic
//!
//! Outputs are stored in an [`indexmap::IndexMap`], so iteration order is first
//! creation order, not hash order. Every manifest, every `finish_all` sequence,
//! and every commit order is therefore reproducible run to run.
//!
//! # Rollback
//!
//! Each output commits as soon as it is finished, which keeps temporary-file
//! usage proportional to the number of *live* outputs rather than to the whole
//! run. The committed paths are remembered, so a later failure still removes
//! everything the run produced: the caller either gets a complete output
//! directory or an empty one.

use std::collections::VecDeque;
use std::path::PathBuf;

use indexmap::IndexMap;

use crate::bam::header::BamHeader;
use crate::bam::raw_record::RawRecord;
use crate::error::OutputError;
use crate::output::filename::{FilenameEncoder, FilenameTemplate};
use crate::output::transaction::{Transaction, prepare_directory};
use crate::output::writer::{FinishedOutput, OutputSettings, RecordWriter};

/// How an output set is bounded and named.
#[derive(Debug, Clone)]
pub struct OutputManagerOptions {
    /// The most outputs that may hold an open file descriptor at once.
    pub max_open_files: usize,
    /// Per-output writer settings.
    pub settings: OutputSettings,
    /// The filename template.
    pub template: FilenameTemplate,
}

impl Default for OutputManagerOptions {
    fn default() -> Self {
        Self {
            max_open_files: 64,
            settings: OutputSettings::default(),
            template: FilenameTemplate::default(),
        }
    }
}

/// Owns every output of a run.
pub struct OutputManager<'h> {
    header: &'h BamHeader,
    reference_count: usize,
    encoder: FilenameEncoder,
    options: OutputManagerOptions,
    live: IndexMap<Vec<u8>, RecordWriter>,
    /// Keys with an open handle, least recently used first.
    open_order: VecDeque<Vec<u8>>,
    finished: Vec<FinishedOutput>,
    committed: Vec<PathBuf>,
}

impl std::fmt::Debug for OutputManager<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OutputManager")
            .field("live", &self.live.len())
            .field("open", &self.open_order.len())
            .field("finished", &self.finished.len())
            .finish_non_exhaustive()
    }
}

impl<'h> OutputManager<'h> {
    /// Creates a manager writing into `directory`.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::BadOutputDirectory`] if the directory cannot be
    /// created or is not a directory.
    pub fn new(
        directory: impl Into<PathBuf>,
        header: &'h BamHeader,
        options: OutputManagerOptions,
    ) -> Result<Self, OutputError> {
        let directory = directory.into();
        prepare_directory(&directory)?;
        Ok(Self {
            header,
            reference_count: header.reference_count(),
            encoder: FilenameEncoder::new(directory, options.template.clone()),
            options,
            live: IndexMap::new(),
            open_order: VecDeque::new(),
            finished: Vec::new(),
            committed: Vec::new(),
        })
    }

    /// How many distinct output keys have been seen, live or finished.
    #[must_use]
    pub fn output_count(&self) -> usize {
        self.live.len() + self.finished.len()
    }

    /// Whether `key` has an output, live or finished.
    #[must_use]
    pub fn contains(&self, key: &[u8]) -> bool {
        self.live.contains_key(key)
            || self
                .finished
                .iter()
                .any(|output| output.paths.logical_key == key)
    }

    /// The outputs finished so far, in completion order.
    #[must_use]
    pub fn finished(&self) -> &[FinishedOutput] {
        &self.finished
    }

    /// Writes one record to `key`, creating the output if needed.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError`] for a naming collision, an existing output, or
    /// an I/O failure.
    pub fn write(&mut self, key: &[u8], record: &RawRecord<'_>) -> Result<(), OutputError> {
        self.ensure_live(key)?;
        self.touch(key)?;
        let writer = self
            .live
            .get_mut(key)
            .ok_or_else(|| internal("the output vanished between creation and use"))?;
        writer.write_record(record)
    }

    /// Creates a header-only output, as `--emit-empty` asks for.
    ///
    /// Does nothing when the key already has an output.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError`] for a naming collision or an I/O failure.
    pub fn create_empty(&mut self, key: &[u8]) -> Result<(), OutputError> {
        if self.contains(key) {
            return Ok(());
        }
        self.ensure_live(key)?;
        Ok(())
    }

    /// Finishes and commits one output.
    ///
    /// The stream engine calls this the moment a routing key changes, which is
    /// what keeps its resource usage flat across a whole genome.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError`] if the output cannot be finished or committed.
    pub fn finish_one(&mut self, key: &[u8]) -> Result<Option<&FinishedOutput>, OutputError> {
        let Some(writer) = self.live.shift_remove(key) else {
            return Ok(None);
        };
        self.open_order.retain(|open| open != key);

        let mut transaction = Transaction::new();
        let finished = writer.finish(&mut transaction)?;
        let paths = transaction.commit()?;
        self.committed.extend(paths);
        self.finished.push(finished);
        Ok(self.finished.last())
    }

    /// Finishes and commits every live output, in creation order.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError`] on the first failure, after rolling back
    /// everything the run has committed.
    pub fn finish_all(mut self) -> Result<Vec<FinishedOutput>, OutputError> {
        let keys: Vec<Vec<u8>> = self.live.keys().cloned().collect();
        for key in keys {
            if let Err(error) = self.finish_one(&key) {
                self.rollback();
                return Err(error);
            }
        }
        Ok(std::mem::take(&mut self.finished))
    }

    /// Discards every output, committed or not.
    ///
    /// Used on interruption and on any hard failure, so a failed run leaves the
    /// output directory as it found it.
    pub fn abort(mut self) {
        self.rollback();
    }

    fn rollback(&mut self) {
        for (_, writer) in std::mem::take(&mut self.live) {
            writer.abort();
        }
        self.open_order.clear();
        for path in std::mem::take(&mut self.committed) {
            let _ = std::fs::remove_file(path);
        }
        self.finished.clear();
    }

    /// Creates the output for `key` if it does not exist yet.
    fn ensure_live(&mut self, key: &[u8]) -> Result<(), OutputError> {
        if self.live.contains_key(key) {
            return Ok(());
        }
        if let Some(existing) = self
            .finished
            .iter()
            .find(|output| output.paths.logical_key == key)
        {
            // Re-opening a finished output would mean the input was not grouped
            // by routing key; the stream engine must catch that first.
            return Err(internal(&format!(
                "output {:?} was already finished and cannot be reopened",
                existing.paths.stem
            )));
        }

        self.make_room()?;
        let extension = self
            .options
            .settings
            .index_kind
            .map(crate::index::IndexKind::extension);
        let paths = self.encoder.assign(key, extension)?;
        let writer = RecordWriter::create(
            paths,
            self.header,
            self.reference_count,
            &self.options.settings,
        )?;
        self.live.insert(key.to_vec(), writer);
        self.open_order.push_back(key.to_vec());
        Ok(())
    }

    /// Parks outputs until another handle can be opened.
    fn make_room(&mut self) -> Result<(), OutputError> {
        if self.options.max_open_files == 0 {
            return Err(OutputError::TooManyOpenFiles { limit: 0 });
        }
        while self.open_order.len() >= self.options.max_open_files {
            let Some(victim) = self.open_order.pop_front() else {
                break;
            };
            if let Some(writer) = self.live.get_mut(&victim) {
                writer.park()?;
            }
        }
        Ok(())
    }

    /// Marks `key` as most recently used, reopening and parking as needed.
    fn touch(&mut self, key: &[u8]) -> Result<(), OutputError> {
        if let Some(position) = self.open_order.iter().position(|open| open == key) {
            if position + 1 != self.open_order.len() {
                let entry = self.open_order.remove(position);
                if let Some(entry) = entry {
                    self.open_order.push_back(entry);
                }
            }
            return Ok(());
        }
        self.make_room()?;
        self.open_order.push_back(key.to_vec());
        Ok(())
    }
}

fn internal(reason: &str) -> OutputError {
    OutputError::Validation {
        path: PathBuf::new(),
        reason: format!("internal error: {reason}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::header::ReferenceSequence;
    use crate::index::IndexKind;

    fn header() -> BamHeader {
        BamHeader::from_parts(
            b"@HD\tVN:1.6\tSO:coordinate\n".to_vec(),
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

    fn record(reference_id: i32, position: i32) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&reference_id.to_le_bytes());
        body.extend_from_slice(&position.to_le_bytes());
        body.extend_from_slice(&[2, 60]);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(b"r\0");
        body
    }

    fn options(max_open_files: usize) -> OutputManagerOptions {
        OutputManagerOptions {
            max_open_files,
            settings: OutputSettings {
                compression_level: 1,
                ..OutputSettings::default()
            },
            ..OutputManagerOptions::default()
        }
    }

    #[test]
    fn writes_and_commits_outputs_in_creation_order() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(8)).expect("created");

        for (key, reference_id) in [(&b"chr1"[..], 0), (b"chr2", 1)] {
            let body = record(reference_id, 10);
            manager
                .write(key, &RawRecord::new(&body).expect("valid"))
                .expect("written");
        }
        let finished = manager.finish_all().expect("finished");

        assert_eq!(finished.len(), 2);
        assert_eq!(finished[0].paths.stem, "chr1");
        assert_eq!(finished[1].paths.stem, "chr2");
        for output in &finished {
            assert!(output.paths.bam.exists());
            assert_eq!(output.stats.record_count, 1);
        }
    }

    #[test]
    fn descriptors_stay_within_the_limit() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(2)).expect("created");

        for index in 0..10u32 {
            let key = format!("key{index}").into_bytes();
            let body = record(0, index as i32);
            manager
                .write(&key, &RawRecord::new(&body).expect("valid"))
                .expect("written");
            assert!(
                manager.open_order.len() <= 2,
                "{} handles open",
                manager.open_order.len()
            );
        }
        // Interleave back to an early key to force a reopen.
        let body = record(0, 999);
        manager
            .write(b"key0", &RawRecord::new(&body).expect("valid"))
            .expect("written");

        let finished = manager.finish_all().expect("finished");
        assert_eq!(finished.len(), 10);
        let first = finished
            .iter()
            .find(|output| output.paths.stem == "key0")
            .expect("present");
        assert_eq!(first.stats.record_count, 2, "the reopened output kept both");
    }

    #[test]
    fn a_zero_descriptor_budget_is_refused() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(0)).expect("created");
        let body = record(0, 1);
        let error = manager
            .write(b"chr1", &RawRecord::new(&body).expect("valid"))
            .expect_err("must refuse");
        assert!(
            matches!(error, OutputError::TooManyOpenFiles { limit: 0 }),
            "{error}"
        );
    }

    #[test]
    fn emit_empty_creates_a_header_only_output() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager = OutputManager::new(
            directory.path(),
            &header,
            OutputManagerOptions {
                settings: OutputSettings {
                    compression_level: 1,
                    index_kind: Some(IndexKind::Bai),
                    ..OutputSettings::default()
                },
                ..options(8)
            },
        )
        .expect("created");

        manager.create_empty(b"chrEmpty").expect("created");
        let finished = manager.finish_all().expect("finished");
        assert_eq!(finished.len(), 1);
        assert_eq!(finished[0].stats.record_count, 0);
        assert!(finished[0].paths.bam.exists());
        assert!(
            finished[0]
                .index_path
                .as_deref()
                .is_some_and(std::path::Path::exists)
        );
    }

    #[test]
    fn finish_one_commits_immediately() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(8)).expect("created");

        let body = record(0, 1);
        manager
            .write(b"chr1", &RawRecord::new(&body).expect("valid"))
            .expect("written");
        let path = manager
            .finish_one(b"chr1")
            .expect("finished")
            .expect("some")
            .paths
            .bam
            .clone();
        assert!(path.exists());
        assert_eq!(manager.output_count(), 1);
        assert!(manager.finish_one(b"chr1").expect("no error").is_none());
    }

    #[test]
    fn reopening_a_finished_output_is_an_internal_error() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(8)).expect("created");

        let body = record(0, 1);
        manager
            .write(b"chr1", &RawRecord::new(&body).expect("valid"))
            .expect("written");
        manager.finish_one(b"chr1").expect("finished");
        let error = manager
            .write(b"chr1", &RawRecord::new(&body).expect("valid"))
            .expect_err("must refuse");
        assert!(matches!(error, OutputError::Validation { .. }), "{error}");
    }

    #[test]
    fn aborting_removes_committed_and_pending_outputs_alike() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(8)).expect("created");

        let body = record(0, 1);
        manager
            .write(b"chr1", &RawRecord::new(&body).expect("valid"))
            .expect("written");
        manager.finish_one(b"chr1").expect("finished");
        manager
            .write(b"chr2", &RawRecord::new(&body).expect("valid"))
            .expect("written");
        manager.abort();

        assert_eq!(
            std::fs::read_dir(directory.path())
                .expect("listable")
                .count(),
            0,
            "a failed run must leave the directory empty"
        );
    }

    #[test]
    fn an_existing_output_is_refused_without_force() {
        let directory = tempfile::tempdir().expect("temp dir");
        std::fs::write(directory.path().join("chr1.bam"), b"old").expect("seeded");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(8)).expect("created");
        let body = record(0, 1);
        let error = manager
            .write(b"chr1", &RawRecord::new(&body).expect("valid"))
            .expect_err("must refuse");
        assert!(
            matches!(error, OutputError::AlreadyExists { .. }),
            "{error}"
        );
    }

    #[test]
    fn unsafe_keys_cannot_escape_the_output_directory() {
        let directory = tempfile::tempdir().expect("temp dir");
        let header = header();
        let mut manager =
            OutputManager::new(directory.path(), &header, options(8)).expect("created");

        let body = record(0, 1);
        for key in [&b"../escape"[..], b"..", b"/etc/passwd", b""] {
            manager
                .write(key, &RawRecord::new(&body).expect("valid"))
                .expect("written");
        }
        let finished = manager.finish_all().expect("finished");
        for output in &finished {
            assert_eq!(
                output.paths.bam.parent(),
                Some(directory.path()),
                "{:?} escaped",
                output.paths.bam
            );
        }
    }
}
