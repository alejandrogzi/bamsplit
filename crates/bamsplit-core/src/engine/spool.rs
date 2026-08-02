// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The spool-and-finalize engine: for inputs whose keys are interleaved.
//!
//! # Why a second pass beats a thousand open BAMs
//!
//! An unsorted BAM, a query-name-sorted BAM, or any `--tag`/`--shard` split
//! visits its output keys in arbitrary order. Writing those directly would need
//! one live BGZF writer per key — each with a compression buffer and a partially
//! built index — which is where a naive splitter runs out of descriptors and
//! then out of memory.
//!
//! So routing and encoding are separated:
//!
//! ```text
//! phase 1   BAM ─▶ raw record ─▶ routing key ─▶ length-prefixed spool file
//! phase 2   spool ─▶ header + raw records ─▶ BGZF ─▶ index ─▶ atomic rename
//! ```
//!
//! Phase 1 writes plain, uncompressed, append-only files through an LRU handle
//! cache bounded by `--max-open-files`. Phase 2 replays each spool exactly once
//! and produces a finished BAM.
//!
//! Phase 2 fans out: every spool is independent, so finalization runs one task
//! per output across the thread budget. Paths are assigned first, on one thread,
//! because the filename encoder has to see every key in a deterministic order to
//! detect collisions; only the replay-and-write is parallel, and results are
//! collected in job order so the manifest does not depend on scheduling. Each
//! task compresses on its own thread rather than through a shared pool, which
//! keeps the total at `tasks` rather than `tasks * compression_workers`.
//!
//! A spool holds only record bodies — never a repeated BAM header — so
//! temporary usage is close to the uncompressed record volume, and each output's
//! index is built and released one at a time.
//!
//! # Spool format
//!
//! ```text
//! magic    "BSPL"                     4 bytes
//! version  u16 LE                     2 bytes
//! reserved u16 LE                     2 bytes  (zero; keeps the header 8-aligned)
//! repeated:
//!   length u32 LE                     4 bytes
//!   body   `length` bytes             the BAM record body, verbatim
//! footer:
//!   sentinel u32 LE = 0xffff_ffff     4 bytes  (an impossible record length)
//!   records  u64 LE                   8 bytes
//!   checksum u64 LE                   8 bytes  (xxh3 over every body)
//! ```
//!
//! The sentinel distinguishes a complete spool from a truncated one, and the
//! record count plus checksum catch a spool that was corrupted between the
//! phases — which is the failure mode that would otherwise silently drop
//! records.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use indexmap::IndexMap;

use crate::bam::raw_record::RawRecord;
use crate::engine::{
    EngineKind, ExecutionContext, ExecutionReport, INTERRUPT_CHECK_INTERVAL, SplitEngine,
    check_interrupt, open,
};
use crate::error::{EngineError, SpoolError};
use crate::output::transaction::{Transaction, prepare_directory};
use crate::routing::{Route, Router, RoutingKey};
use crate::stats::RunStats;

/// The spool magic number.
pub const SPOOL_MAGIC: [u8; 4] = *b"BSPL";

/// The spool format version this build writes and reads.
pub const SPOOL_VERSION: u16 = 1;

/// The end-of-records sentinel, chosen because no BAM record can be this long.
pub const SPOOL_SENTINEL: u32 = u32::MAX;

const SPOOL_HEADER_LEN: usize = 8;
const SPOOL_FOOTER_LEN: usize = 20;

/// The largest record body a spool will accept on replay.
const MAX_SPOOL_RECORD: u64 = 1 << 31;

/// One key's temporary record store.
struct Spool {
    path: PathBuf,
    handle: Option<BufWriter<File>>,
    records: u64,
    bytes: u64,
    hasher: xxhash_rust::xxh3::Xxh3,
    finished: bool,
}

impl std::fmt::Debug for Spool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Spool")
            .field("path", &self.path)
            .field("records", &self.records)
            .field("bytes", &self.bytes)
            .field("open", &self.handle.is_some())
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl Spool {
    fn create(directory: &Path, ordinal: usize) -> Result<Self, SpoolError> {
        // The file name is the *ordinal*, not the key: a key can be arbitrary
        // bytes, and encoding it here would duplicate the filename logic and
        // could collide. The mapping back to keys lives in memory.
        let path = directory.join(format!("spool-{ordinal:08}.bspl"));
        let mut handle = BufWriter::with_capacity(
            1 << 18,
            File::create(&path).map_err(|source| SpoolError::Io {
                path: path.clone(),
                source,
            })?,
        );
        handle
            .write_all(&SPOOL_MAGIC)
            .and_then(|()| handle.write_all(&SPOOL_VERSION.to_le_bytes()))
            .and_then(|()| handle.write_all(&0u16.to_le_bytes()))
            .map_err(|source| SpoolError::Io {
                path: path.clone(),
                source,
            })?;
        crate::output::transaction::registry::register(&path);
        Ok(Self {
            path,
            handle: Some(handle),
            records: 0,
            bytes: u64::try_from(SPOOL_HEADER_LEN).unwrap_or(0),
            hasher: xxhash_rust::xxh3::Xxh3::new(),
            finished: false,
        })
    }

    fn reopen(&mut self) -> Result<&mut BufWriter<File>, SpoolError> {
        if self.handle.is_none() {
            let file = std::fs::OpenOptions::new()
                .append(true)
                .open(&self.path)
                .map_err(|source| SpoolError::Io {
                    path: self.path.clone(),
                    source,
                })?;
            self.handle = Some(BufWriter::with_capacity(1 << 18, file));
        }
        self.handle.as_mut().ok_or_else(|| SpoolError::Io {
            path: self.path.clone(),
            source: std::io::Error::other("the spool handle vanished"),
        })
    }

    fn append(&mut self, body: &[u8]) -> Result<(), SpoolError> {
        let length = u32::try_from(body.len()).map_err(|_| SpoolError::RecordTooLarge {
            path: self.path.clone(),
            length: body.len() as u64,
            offset: self.bytes,
            limit: MAX_SPOOL_RECORD,
        })?;
        let path = self.path.clone();
        let handle = self.reopen()?;
        handle
            .write_all(&length.to_le_bytes())
            .and_then(|()| handle.write_all(body))
            .map_err(|source| SpoolError::Io { path, source })?;
        self.hasher.update(body);
        self.records += 1;
        self.bytes += 4 + u64::from(length);
        Ok(())
    }

    fn park(&mut self) -> Result<(), SpoolError> {
        if let Some(mut handle) = self.handle.take() {
            handle.flush().map_err(|source| SpoolError::Io {
                path: self.path.clone(),
                source,
            })?;
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<(), SpoolError> {
        if self.finished {
            return Ok(());
        }
        let checksum = self.hasher.digest();
        let records = self.records;
        let path = self.path.clone();
        let handle = self.reopen()?;
        handle
            .write_all(&SPOOL_SENTINEL.to_le_bytes())
            .and_then(|()| handle.write_all(&records.to_le_bytes()))
            .and_then(|()| handle.write_all(&checksum.to_le_bytes()))
            .and_then(|()| handle.flush())
            .map_err(|source| SpoolError::Io { path, source })?;
        self.handle = None;
        self.bytes += u64::try_from(SPOOL_FOOTER_LEN).unwrap_or(0);
        self.finished = true;
        Ok(())
    }
}

/// Replays a spool file, validating its framing, count, and checksum.
///
/// The callback receives each record body in the order it was spooled, so
/// output order matches input order.
///
/// # Errors
///
/// Returns [`SpoolError`] for bad magic, an unsupported version, truncation, an
/// implausible length, a record-count mismatch, or a checksum mismatch.
pub fn replay<F>(path: &Path, mut visit: F) -> Result<u64, SpoolError>
where
    F: FnMut(&[u8]) -> Result<(), SpoolError>,
{
    let file = File::open(path).map_err(|source| SpoolError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut reader = BufReader::with_capacity(1 << 18, file);

    let mut header = [0u8; SPOOL_HEADER_LEN];
    read_exact(&mut reader, &mut header, path, "the spool header")?;
    let magic: [u8; 4] = [header[0], header[1], header[2], header[3]];
    if magic != SPOOL_MAGIC {
        return Err(SpoolError::BadMagic {
            path: path.to_path_buf(),
            expected: SPOOL_MAGIC,
            found: magic,
        });
    }
    let version = u16::from_le_bytes([header[4], header[5]]);
    if version != SPOOL_VERSION {
        return Err(SpoolError::UnsupportedVersion {
            path: path.to_path_buf(),
            version,
            supported: SPOOL_VERSION,
        });
    }

    let mut hasher = xxhash_rust::xxh3::Xxh3::new();
    let mut body = Vec::new();
    let mut records = 0u64;
    let mut offset = SPOOL_HEADER_LEN as u64;

    loop {
        let mut prefix = [0u8; 4];
        read_exact(&mut reader, &mut prefix, path, "a record length prefix")?;
        let length = u32::from_le_bytes(prefix);
        offset += 4;
        if length == SPOOL_SENTINEL {
            break;
        }
        if u64::from(length) > MAX_SPOOL_RECORD {
            return Err(SpoolError::RecordTooLarge {
                path: path.to_path_buf(),
                length: u64::from(length),
                offset: offset - 4,
                limit: MAX_SPOOL_RECORD,
            });
        }
        body.clear();
        body.resize(length as usize, 0);
        read_exact(&mut reader, &mut body, path, "a record body")?;
        offset += u64::from(length);
        hasher.update(&body);
        records += 1;
        visit(&body)?;
    }

    let mut footer = [0u8; 16];
    read_exact(&mut reader, &mut footer, path, "the spool footer")?;
    let declared_records = u64::from_le_bytes(footer[..8].try_into().unwrap_or([0; 8]));
    let declared_checksum = u64::from_le_bytes(footer[8..].try_into().unwrap_or([0; 8]));

    if declared_records != records {
        return Err(SpoolError::RecordCountMismatch {
            path: path.to_path_buf(),
            expected: declared_records,
            actual: records,
        });
    }
    let checksum = hasher.digest();
    if declared_checksum != checksum {
        return Err(SpoolError::ChecksumMismatch {
            path: path.to_path_buf(),
            expected: declared_checksum,
            actual: checksum,
        });
    }
    Ok(records)
}

fn read_exact<R: Read>(
    reader: &mut R,
    buf: &mut [u8],
    path: &Path,
    what: &str,
) -> Result<(), SpoolError> {
    match reader.read_exact(buf) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => {
            Err(SpoolError::Truncated {
                path: path.to_path_buf(),
                reason: format!("the file ends inside {what}"),
            })
        }
        Err(source) => Err(SpoolError::Io {
            path: path.to_path_buf(),
            source,
        }),
    }
}

/// The spool-and-finalize engine.
#[derive(Debug, Clone, Copy, Default)]
pub struct SpoolEngine;

impl<R: Router> SplitEngine<R> for SpoolEngine {
    fn execute(
        &self,
        context: &ExecutionContext<'_>,
        router: &R,
    ) -> Result<ExecutionReport, EngineError> {
        let scratch = SpoolDirectory::new(context)?;
        let (input_header, mut input, effective_io) = open(
            &context.input,
            context.io,
            context.threads.decompression,
            context.max_record_size,
        )?;

        // Phase 1: route every record into a spool, in input order.
        let mut spools: IndexMap<Vec<u8>, Spool> = IndexMap::new();
        let mut open_order: VecDeque<Vec<u8>> = VecDeque::new();
        let mut stats = RunStats::default();
        let mut records = 0u64;
        let max_open = context.manager_options.max_open_files.max(1);

        let phase_one = (|| -> Result<(), EngineError> {
            while let Some(record) = input.next_record()? {
                records += 1;
                if records % INTERRUPT_CHECK_INTERVAL == 0 {
                    check_interrupt(context.interrupt, records)?;
                }
                let route = router
                    .route(&input_header, &record)
                    .map_err(|error| Box::new(error.relocate_routing(records)))?;

                match &route {
                    Route::Drop(reason) => {
                        if reason.is_unmatched() {
                            stats.observe_emissions(0);
                        } else {
                            stats.observe_drop();
                        }
                        continue;
                    }
                    Route::One(key) => spool_record(
                        &mut spools,
                        &mut open_order,
                        max_open,
                        scratch.path(),
                        key.logical(),
                        record.raw_bytes(),
                    )?,
                    Route::Many(keys) => {
                        for key in keys {
                            spool_record(
                                &mut spools,
                                &mut open_order,
                                max_open,
                                scratch.path(),
                                key.logical(),
                                record.raw_bytes(),
                            )?;
                        }
                    }
                }
                stats.observe_emissions(route.emission_count() as u64);
            }
            for spool in spools.values_mut() {
                spool.finish().map_err(Box::new)?;
            }
            Ok(())
        })();
        phase_one?;

        stats.temporary_bytes = spools.values().map(|spool| spool.bytes).sum();

        // Phase 2: replay each spool into a finished BAM.
        //
        // Every spool is independent — its own records, its own output, its own
        // index — so finalization fans out across the task budget. Paths are
        // assigned up front on one thread, because the filename encoder detects
        // collisions and must see every key in a deterministic order; only the
        // replay-and-write is parallel.
        //
        // Each task compresses on its own thread rather than through a shared
        // pool: `tasks * 1` stays inside the budget, whereas `tasks *
        // compression_workers` would multiply past it.
        // Phase 2 writes outputs directly rather than through the output
        // manager, so the directory it would have created has to be created
        // here.
        prepare_directory(&context.output_directory).map_err(Box::new)?;
        let mut encoder = crate::output::FilenameEncoder::new(
            &context.output_directory,
            context.manager_options.template.clone(),
        );
        let extension = context
            .manager_options
            .settings
            .index_kind
            .map(crate::index::IndexKind::extension);

        let mut jobs: Vec<FinalizeJob> = Vec::with_capacity(spools.len());
        for (key, spool) in &spools {
            let paths = encoder.assign(key, extension).map_err(Box::new)?;
            jobs.push(FinalizeJob {
                paths,
                spool: spool.path.clone(),
                expected: spool.records,
            });
        }

        let concurrency = context.threads.tasks.min(jobs.len()).max(1);
        let mut settings = context.manager_options.settings.clone();
        if concurrency > 1 {
            settings.compression_workers = 1;
            settings.pool = None;
        }

        let finalize = |job: &FinalizeJob| -> Result<FinalizedSpool, EngineError> {
            check_interrupt(context.interrupt, records)?;
            let mut writer = crate::output::RecordWriter::create(
                job.paths.clone(),
                context.header,
                input_header.reference_count(),
                &settings,
            )
            .map_err(Box::new)?;

            let count = replay(&job.spool, |body| {
                let record = RawRecord::from_validated_bytes(body);
                writer
                    .write_record(&record)
                    .map_err(|error| SpoolError::Io {
                        path: job.spool.clone(),
                        source: std::io::Error::other(error.to_string()),
                    })
            })
            .map_err(Box::new)?;

            if count != job.expected {
                return Err(EngineError::from(SpoolError::RecordCountMismatch {
                    path: job.spool.clone(),
                    expected: job.expected,
                    actual: count,
                }));
            }

            let mut transaction = Transaction::new();
            let finished = writer.finish(&mut transaction).map_err(Box::new)?;
            let committed = transaction.commit().map_err(Box::new)?;
            Ok(FinalizedSpool {
                finished,
                committed,
                records: count,
            })
        };

        let results: Vec<Result<FinalizedSpool, EngineError>> = if concurrency > 1 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(concurrency)
                .thread_name(|index| format!("bamsplit-f{index}"))
                .build()
                .map_err(|source| EngineError::WorkerFailure {
                    reason: source.to_string(),
                })?;
            pool.install(|| {
                use rayon::prelude::*;
                jobs.par_iter().map(&finalize).collect()
            })
        } else {
            jobs.iter().map(&finalize).collect()
        };

        // Outputs are collected in job order, not completion order, so the
        // manifest is identical however the pool happened to schedule them.
        let mut outputs = Vec::with_capacity(results.len());
        let mut committed = Vec::new();
        let mut replayed = 0u64;
        let mut failure = None;
        for result in results {
            match result {
                Ok(spool) => {
                    replayed += spool.records;
                    committed.extend(spool.committed);
                    outputs.push(spool.finished);
                }
                Err(error) => failure = failure.or(Some(error)),
            }
        }
        if let Some(error) = failure {
            for path in committed {
                let _ = std::fs::remove_file(path);
            }
            return Err(error);
        }

        let mut skipped_keys = Vec::new();
        for key in router.declared_keys(&input_header) {
            if spools.contains_key(key.logical()) {
                continue;
            }
            if context.emit_empty {
                match empty_output(
                    context,
                    &input_header,
                    &mut encoder,
                    extension,
                    key.logical(),
                ) {
                    Ok((finished, paths)) => {
                        committed.extend(paths);
                        outputs.push(finished);
                    }
                    Err(error) => {
                        for path in committed {
                            let _ = std::fs::remove_file(path);
                        }
                        return Err(error);
                    }
                }
            } else {
                skipped_keys.push(key.logical().to_vec());
            }
        }

        stats.total_output_emissions = outputs.iter().map(|output| output.stats.record_count).sum();

        let notes = vec![
            format!(
                "spooled {records} records into {} temporary files ({} bytes), then finalized each once",
                spools.len(),
                stats.temporary_bytes
            ),
            format!(
                "replayed {replayed} records with checksum verification across {concurrency} \
                 finalization task(s)"
            ),
        ];

        Ok(ExecutionReport {
            outputs,
            stats,
            skipped_keys,
            engine: EngineKind::Spool,
            io: effective_io,
            notes,
        })
    }
}

/// One spool waiting to be turned into a BAM.
#[derive(Debug)]
struct FinalizeJob {
    paths: crate::output::OutputPaths,
    spool: PathBuf,
    expected: u64,
}

/// One finished output, plus what it committed.
#[derive(Debug)]
struct FinalizedSpool {
    finished: crate::output::FinishedOutput,
    committed: Vec<PathBuf>,
    records: u64,
}

/// Creates a header-only output for a declared key that received no record.
fn empty_output(
    context: &ExecutionContext<'_>,
    input_header: &crate::bam::header::BamHeader,
    encoder: &mut crate::output::FilenameEncoder,
    extension: Option<&str>,
    key: &[u8],
) -> Result<(crate::output::FinishedOutput, Vec<PathBuf>), EngineError> {
    let paths = encoder.assign(key, extension).map_err(Box::new)?;
    let writer = crate::output::RecordWriter::create(
        paths,
        context.header,
        input_header.reference_count(),
        &context.manager_options.settings,
    )
    .map_err(Box::new)?;
    let mut transaction = Transaction::new();
    let finished = writer.finish(&mut transaction).map_err(Box::new)?;
    let committed = transaction.commit().map_err(Box::new)?;
    Ok((finished, committed))
}

/// Appends one record to a key's spool, parking handles to stay within the
/// descriptor budget.
fn spool_record(
    spools: &mut IndexMap<Vec<u8>, Spool>,
    open_order: &mut VecDeque<Vec<u8>>,
    max_open: usize,
    directory: &Path,
    key: &[u8],
    body: &[u8],
) -> Result<(), EngineError> {
    if !spools.contains_key(key) {
        let ordinal = spools.len();
        let spool = Spool::create(directory, ordinal).map_err(Box::new)?;
        spools.insert(key.to_vec(), spool);
    }

    if let Some(position) = open_order.iter().position(|open| open == key) {
        if position + 1 != open_order.len()
            && let Some(entry) = open_order.remove(position)
        {
            open_order.push_back(entry);
        }
    } else {
        while open_order.len() >= max_open {
            let Some(victim) = open_order.pop_front() else {
                break;
            };
            if let Some(spool) = spools.get_mut(&victim) {
                spool.park().map_err(Box::new)?;
            }
        }
        open_order.push_back(key.to_vec());
    }

    let spool = spools
        .get_mut(key)
        .ok_or_else(|| EngineError::WorkerFailure {
            reason: "internal error: a spool vanished after creation".to_string(),
        })?;
    spool.append(body).map_err(Box::new)?;
    Ok(())
}

/// The scratch directory for one run's spools.
///
/// Removed on drop unless `--keep-temp` was given, so an interrupted run does
/// not leave gigabytes behind.
#[derive(Debug)]
struct SpoolDirectory {
    path: PathBuf,
    keep: bool,
    owned: bool,
}

impl SpoolDirectory {
    fn new(context: &ExecutionContext<'_>) -> Result<Self, EngineError> {
        let base = context.temp_dir.clone().unwrap_or_else(std::env::temp_dir);
        prepare_directory(&base).map_err(Box::new)?;
        let path = base.join(format!(
            "bamsplit-spool-{}-{:x}",
            std::process::id(),
            xxhash_rust::xxh3::xxh3_64(context.output_directory.as_os_str().as_encoded_bytes())
        ));
        std::fs::create_dir_all(&path).map_err(|source| {
            EngineError::from(SpoolError::BadTempDirectory {
                path: path.clone(),
                reason: source.to_string(),
            })
        })?;
        Ok(Self {
            path,
            keep: context.keep_temp,
            owned: true,
        })
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for SpoolDirectory {
    fn drop(&mut self) {
        if !self.owned || self.keep {
            return;
        }
        if let Ok(entries) = std::fs::read_dir(&self.path) {
            for entry in entries.flatten() {
                crate::output::transaction::registry::unregister(&entry.path());
            }
        }
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spool_with(bodies: &[&[u8]]) -> (tempfile::TempDir, PathBuf) {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut spool = Spool::create(directory.path(), 0).expect("created");
        for body in bodies {
            spool.append(body).expect("appended");
        }
        spool.finish().expect("finished");
        let path = spool.path.clone();
        (directory, path)
    }

    #[test]
    fn a_spool_round_trips_bodies_in_order() {
        let bodies: Vec<Vec<u8>> = (0..200u32)
            .map(|index| vec![(index % 251) as u8; 1 + (index as usize % 97)])
            .collect();
        let refs: Vec<&[u8]> = bodies.iter().map(Vec::as_slice).collect();
        let (_directory, path) = spool_with(&refs);

        let mut seen = Vec::new();
        let count = replay(&path, |body| {
            seen.push(body.to_vec());
            Ok(())
        })
        .expect("replayed");
        assert_eq!(count, bodies.len() as u64);
        assert_eq!(seen, bodies);
    }

    #[test]
    fn an_empty_spool_is_valid() {
        let (_directory, path) = spool_with(&[]);
        let count = replay(&path, |_| Ok(())).expect("replayed");
        assert_eq!(count, 0);
    }

    #[test]
    fn a_spool_never_repeats_a_bam_header() {
        let (_directory, path) = spool_with(&[b"body-one", b"body-two"]);
        let bytes = std::fs::read(&path).expect("readable");
        assert_eq!(&bytes[..4], &SPOOL_MAGIC);
        // header(8) + 2 * (4 + 8) + footer(20): no room for a repeated header.
        assert_eq!(bytes.len(), 8 + 2 * (4 + 8) + 20);
    }

    #[test]
    fn bad_magic_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("bad.bspl");
        std::fs::write(&path, b"XXXX\x01\x00\x00\x00").expect("written");
        let error = replay(&path, |_| Ok(())).expect_err("must reject");
        assert!(matches!(error, SpoolError::BadMagic { .. }), "{error}");
    }

    #[test]
    fn an_unsupported_version_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("future.bspl");
        let mut bytes = SPOOL_MAGIC.to_vec();
        bytes.extend_from_slice(&99u16.to_le_bytes());
        bytes.extend_from_slice(&0u16.to_le_bytes());
        std::fs::write(&path, bytes).expect("written");
        let error = replay(&path, |_| Ok(())).expect_err("must reject");
        assert!(
            matches!(error, SpoolError::UnsupportedVersion { version: 99, .. }),
            "{error}"
        );
    }

    #[test]
    fn a_truncated_spool_is_rejected() {
        let (_directory, path) = spool_with(&[b"aaaa", b"bbbb"]);
        let mut bytes = std::fs::read(&path).expect("readable");
        bytes.truncate(bytes.len() - 6);
        std::fs::write(&path, bytes).expect("written");
        let error = replay(&path, |_| Ok(())).expect_err("must reject");
        assert!(matches!(error, SpoolError::Truncated { .. }), "{error}");
    }

    #[test]
    fn a_corrupted_body_fails_the_checksum() {
        let (_directory, path) = spool_with(&[b"aaaaaaaa"]);
        let mut bytes = std::fs::read(&path).expect("readable");
        // Flip a byte inside the record body, leaving the framing intact.
        let position = SPOOL_HEADER_LEN + 4;
        bytes[position] ^= 0xff;
        std::fs::write(&path, bytes).expect("written");
        let error = replay(&path, |_| Ok(())).expect_err("must reject");
        assert!(
            matches!(error, SpoolError::ChecksumMismatch { .. }),
            "{error}"
        );
    }

    #[test]
    fn a_tampered_record_count_is_rejected() {
        let (_directory, path) = spool_with(&[b"aaaa"]);
        let mut bytes = std::fs::read(&path).expect("readable");
        let count_at = bytes.len() - 16;
        bytes[count_at..count_at + 8].copy_from_slice(&99u64.to_le_bytes());
        std::fs::write(&path, bytes).expect("written");
        let error = replay(&path, |_| Ok(())).expect_err("must reject");
        assert!(
            matches!(
                error,
                SpoolError::RecordCountMismatch {
                    expected: 99,
                    actual: 1,
                    ..
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn an_implausible_length_prefix_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("huge.bspl");
        let mut bytes = SPOOL_MAGIC.to_vec();
        bytes.extend_from_slice(&SPOOL_VERSION.to_le_bytes());
        bytes.extend_from_slice(&0u16.to_le_bytes());
        // One below the sentinel, so it is read as a length rather than an end.
        bytes.extend_from_slice(&(u32::MAX - 1).to_le_bytes());
        std::fs::write(&path, bytes).expect("written");
        let error = replay(&path, |_| Ok(())).expect_err("must reject");
        assert!(
            matches!(error, SpoolError::RecordTooLarge { .. }),
            "{error}"
        );
    }

    #[test]
    fn parking_and_resuming_a_spool_preserves_records() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut spool = Spool::create(directory.path(), 0).expect("created");
        spool.append(b"first").expect("appended");
        spool.park().expect("parked");
        spool.append(b"second").expect("appended");
        spool.park().expect("parked");
        spool.finish().expect("finished");

        let mut seen = Vec::new();
        replay(&spool.path, |body| {
            seen.push(body.to_vec());
            Ok(())
        })
        .expect("replayed");
        assert_eq!(seen, vec![b"first".to_vec(), b"second".to_vec()]);
    }

    #[test]
    fn finishing_twice_is_harmless() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut spool = Spool::create(directory.path(), 0).expect("created");
        spool.append(b"x").expect("appended");
        spool.finish().expect("finished");
        spool.finish().expect("still fine");
        assert_eq!(replay(&spool.path, |_| Ok(())).expect("replayed"), 1);
    }
}
