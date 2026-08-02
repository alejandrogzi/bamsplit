// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! How a split pass is executed.
//!
//! Three engines share one interface, because the right choice depends on the
//! input rather than on the command:
//!
//! | engine | when | memory | passes | descriptors |
//! | --- | --- | --- | --- | --- |
//! | [`stream`] | output keys arrive grouped (coordinate-sorted `chrom`) | ~constant | 1 | 1 |
//! | [`indexed`] | seekable input with a BAI/CSI, independent per-reference work | ~constant × workers | 1 logical, random physical | workers |
//! | [`spool`] | unsorted input, interleaved keys, high cardinality, stdin | ~constant | 2 (write, then replay) | bounded |
//!
//! Everything above the engine — the router, the output manager, the manifest —
//! is identical in all three, so a run's *result* does not depend on which one
//! ran. `tests/integration` asserts exactly that by comparing digests across
//! engines.

pub mod indexed;
pub mod spool;
pub mod stream;

use std::fs::File;
use std::io::Read;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::bam::header::BamHeader;
use crate::bam::raw_record::{RawRecordReader, VirtualPositionSource};
use crate::error::{ConfigError, EngineError};
use crate::output::OutputManagerOptions;
use crate::output::transaction::InterruptFlag;
use crate::output::writer::FinishedOutput;
use crate::routing::Router;
use crate::stats::RunStats;

/// Which engine ran, or was requested.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EngineKind {
    /// Choose from the input's properties.
    #[default]
    Auto,
    /// The grouped streaming engine.
    Stream,
    /// The indexed parallel engine.
    Indexed,
    /// The spool-and-finalize engine.
    Spool,
}

impl EngineKind {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Stream => "stream",
            Self::Indexed => "indexed",
            Self::Spool => "spool",
        }
    }
}

impl std::fmt::Display for EngineKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// How input bytes reach the decompressor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum IoBackend {
    /// Choose per engine and input kind.
    #[default]
    Auto,
    /// Ordinary buffered reads.
    Buffered,
    /// A memory mapping, where the input allows it.
    Mmap,
}

impl IoBackend {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Buffered => "buffered",
            Self::Mmap => "mmap",
        }
    }
}

impl std::fmt::Display for IoBackend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Where the input comes from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputSource {
    /// A regular file.
    Path(PathBuf),
    /// Standard input, which is neither seekable nor mappable.
    Stdin,
}

impl InputSource {
    /// Interprets a path, treating `-` as standard input.
    #[must_use]
    pub fn from_arg(value: &Path) -> Self {
        if value == Path::new("-") {
            Self::Stdin
        } else {
            Self::Path(value.to_path_buf())
        }
    }

    /// Whether the source can be seeked, and therefore indexed or re-read.
    #[must_use]
    pub fn is_seekable(&self) -> bool {
        matches!(self, Self::Path(path) if path.is_file())
    }

    /// The path, when there is one.
    #[must_use]
    pub fn path(&self) -> Option<&Path> {
        match self {
            Self::Path(path) => Some(path),
            Self::Stdin => None,
        }
    }

    /// The input size in bytes, when it is knowable.
    #[must_use]
    pub fn size(&self) -> Option<u64> {
        self.path()
            .and_then(|path| std::fs::metadata(path).ok())
            .map(|metadata| metadata.len())
    }

    /// A label for logs and the manifest.
    #[must_use]
    pub fn display(&self) -> String {
        match self {
            Self::Path(path) => path.display().to_string(),
            Self::Stdin => "<stdin>".to_string(),
        }
    }

    /// Validates that the source is usable.
    ///
    /// # Errors
    ///
    /// Returns [`ConfigError::BadInput`] for a path that does not exist, is a
    /// directory, or cannot be opened.
    pub fn validate(&self) -> Result<(), ConfigError> {
        let Self::Path(path) = self else {
            return Ok(());
        };
        if !path.exists() {
            return Err(ConfigError::BadInput {
                path: path.clone(),
                reason: "the file does not exist".to_string(),
            });
        }
        if path.is_dir() {
            return Err(ConfigError::BadInput {
                path: path.clone(),
                reason: "the path is a directory".to_string(),
            });
        }
        File::open(path).map_err(|source| ConfigError::BadInput {
            path: path.clone(),
            reason: source.to_string(),
        })?;
        Ok(())
    }
}

/// How the global `--threads` budget is divided.
///
/// The budget is a total, not a per-output allowance: a 3 000-contig split with
/// `--threads 16` must still use 16 threads, not 48 000.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThreadBudget {
    /// The total requested.
    pub total: usize,
    /// Threads devoted to BGZF decompression of the input.
    pub decompression: usize,
    /// Threads devoted to compressing output.
    pub compression: usize,
    /// Concurrent per-reference or per-output tasks.
    pub tasks: usize,
}

impl Default for ThreadBudget {
    fn default() -> Self {
        Self::for_total(1)
    }
}

impl ThreadBudget {
    /// Splits a total budget for a single-writer pass.
    ///
    /// The streaming engine keeps exactly one output open, so almost the whole
    /// budget can go to that one writer's compression — which is where a
    /// `chrom` split spends its time.
    #[must_use]
    pub fn for_total(total: usize) -> Self {
        let total = total.max(1);
        if total == 1 {
            return Self {
                total: 1,
                decompression: 1,
                compression: 1,
                tasks: 1,
            };
        }
        // One thread routes; the rest is split between inflating input and
        // deflating output, weighted towards compression because deflate is the
        // more expensive half.
        let workers = total - 1;
        let decompression = (workers / 3).max(1);
        let compression = workers.saturating_sub(decompression).max(1);
        Self {
            total,
            decompression,
            compression,
            tasks: 1,
        }
    }

    /// Splits a total budget across `concurrency` independent tasks.
    ///
    /// Used by the indexed and spool engines, where several outputs are written
    /// at once. Each task gets at least one compression thread, and the task
    /// count is capped so `tasks * compression <= total`.
    #[must_use]
    pub fn for_tasks(total: usize, concurrency: usize) -> Self {
        let total = total.max(1);
        let tasks = concurrency.clamp(1, total);
        let compression = (total / tasks).max(1);
        Self {
            total,
            decompression: 1,
            compression,
            tasks,
        }
    }

    /// A rayon pool sized to this budget's compression share.
    ///
    /// One pool is shared by every output so the thread count stays fixed
    /// regardless of how many outputs are live.
    ///
    /// # Errors
    ///
    /// Returns [`EngineError::WorkerFailure`] if the pool cannot be created.
    pub fn build_pool(&self) -> Result<Option<Arc<rayon::ThreadPool>>, EngineError> {
        let threads = self.compression.max(self.tasks);
        if threads <= 1 {
            return Ok(None);
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .thread_name(|index| format!("bamsplit-z{index}"))
            .build()
            .map(|pool| Some(Arc::new(pool)))
            .map_err(|source| EngineError::WorkerFailure {
                reason: format!("cannot create a thread pool of {threads} threads: {source}"),
            })
    }
}

/// Everything an engine needs that is not the router.
pub struct ExecutionContext<'a> {
    /// Where records come from.
    pub input: InputSource,
    /// The header to write to every output, already carrying the `@PG` record.
    pub header: &'a BamHeader,
    /// Where outputs go.
    pub output_directory: PathBuf,
    /// How outputs are named, bounded, and compressed.
    pub manager_options: OutputManagerOptions,
    /// The thread budget.
    pub threads: ThreadBudget,
    /// Whether to materialize header-only BAMs for declared-but-unused keys.
    pub emit_empty: bool,
    /// Where temporary spool files go.
    pub temp_dir: Option<PathBuf>,
    /// Whether to keep temporary files after the run.
    pub keep_temp: bool,
    /// The I/O backend.
    pub io: IoBackend,
    /// The cooperative interruption flag.
    pub interrupt: &'a InterruptFlag,
    /// The largest accepted BAM record body.
    pub max_record_size: usize,
}

impl std::fmt::Debug for ExecutionContext<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ExecutionContext")
            .field("input", &self.input)
            .field("output_directory", &self.output_directory)
            .field("threads", &self.threads)
            .field("emit_empty", &self.emit_empty)
            .field("io", &self.io)
            .finish_non_exhaustive()
    }
}

/// What an engine produced.
#[derive(Debug, Default)]
pub struct ExecutionReport {
    /// Every committed output, in completion order.
    pub outputs: Vec<FinishedOutput>,
    /// Run-wide counters.
    pub stats: RunStats,
    /// Declared keys that received no record and were not materialized.
    pub skipped_keys: Vec<Vec<u8>>,
    /// Which engine ran.
    pub engine: EngineKind,
    /// Which I/O backend was used.
    pub io: IoBackend,
    /// Human-readable notes for the log and the manifest.
    pub notes: Vec<String>,
}

/// Executes a routing plan.
pub trait SplitEngine<R: Router> {
    /// Runs the pass.
    ///
    /// # Errors
    ///
    /// Returns [`EngineError`] for malformed input, an output failure, an
    /// interruption, or a violated precondition such as an unsorted input under
    /// the streaming engine.
    fn execute(
        &self,
        context: &ExecutionContext<'_>,
        router: &R,
    ) -> Result<ExecutionReport, EngineError>;
}

/// How often the record loop checks for an interruption.
///
/// Checking every record would put an atomic load on the hot path; 64 Ki
/// records is well under a second of work even on fast local storage, so a
/// `Ctrl-C` still feels immediate.
pub const INTERRUPT_CHECK_INTERVAL: u64 = 65_536;

/// A BGZF-decompressed BAM input, with its header already consumed.
///
/// The concrete reader type is erased behind an enum rather than a trait object
/// so the hot `read_record` path stays statically dispatched per variant.
pub enum BamInput {
    /// Single-threaded inflation.
    Serial(Box<RawRecordReader<noodles_bgzf::io::Reader<Box<dyn Read + Send>>>>),
    /// Multi-threaded inflation.
    Parallel(Box<RawRecordReader<noodles_bgzf::io::MultithreadedReader<Box<dyn Read + Send>>>>),
}

impl std::fmt::Debug for BamInput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Serial(_) => f.write_str("BamInput::Serial"),
            Self::Parallel(_) => f.write_str("BamInput::Parallel"),
        }
    }
}

/// Opens a BAM with single-threaded inflation, reading and validating its
/// header.
///
/// # Errors
///
/// Returns [`EngineError`] if the input cannot be opened, is not a BAM, or has
/// an invalid header.
pub fn open_bam(
    source: &InputSource,
    io: IoBackend,
    decompression_threads: usize,
    max_record_size: usize,
) -> Result<(BamHeader, BamInput, IoBackend), EngineError> {
    let (raw, effective_io): (Box<dyn Read + Send>, IoBackend) = match source {
        InputSource::Stdin => {
            // stdin is neither seekable nor mappable; asking for mmap here is a
            // configuration error the caller should have caught, so fall back
            // rather than fail late.
            (Box::new(std::io::stdin()), IoBackend::Buffered)
        }
        InputSource::Path(path) => open_path(path, io)?,
    };

    let _ = decompression_threads;
    let mut reader = noodles_bgzf::io::Reader::new(raw);
    let header = BamHeader::read_from(&mut reader).map_err(Box::new)?;
    Ok((
        header,
        BamInput::Serial(Box::new(
            RawRecordReader::new(reader).with_max_record_size(max_record_size),
        )),
        effective_io,
    ))
}

/// Opens a BAM, choosing serial or multi-threaded inflation from the budget.
///
/// # Errors
///
/// Returns [`EngineError`] if the input cannot be opened or its header is
/// invalid.
pub fn open(
    source: &InputSource,
    io: IoBackend,
    decompression_threads: usize,
    max_record_size: usize,
) -> Result<(BamHeader, BamInput, IoBackend), EngineError> {
    if decompression_threads > 1 {
        open_bam_parallel(source, io, decompression_threads, max_record_size)
    } else {
        open_bam(source, io, 1, max_record_size)
    }
}

/// Opens a BAM for multi-threaded inflation, re-reading the header.
///
/// Kept separate from [`open_bam`] because a `MultithreadedReader` cannot adopt
/// a stream whose first block another reader has already buffered.
///
/// # Errors
///
/// Returns [`EngineError`] if the input cannot be opened or its header is
/// invalid.
pub fn open_bam_parallel(
    source: &InputSource,
    io: IoBackend,
    decompression_threads: usize,
    max_record_size: usize,
) -> Result<(BamHeader, BamInput, IoBackend), EngineError> {
    if decompression_threads <= 1 {
        return open_bam(source, io, 1, max_record_size);
    }
    let (raw, effective_io): (Box<dyn Read + Send>, IoBackend) = match source {
        InputSource::Stdin => (Box::new(std::io::stdin()), IoBackend::Buffered),
        InputSource::Path(path) => open_path(path, io)?,
    };
    let workers = NonZero::new(decompression_threads).unwrap_or(NonZero::<usize>::MIN);
    let mut reader = noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, raw);
    let header = BamHeader::read_from(&mut reader).map_err(Box::new)?;
    Ok((
        header,
        BamInput::Parallel(Box::new(
            RawRecordReader::new(reader).with_max_record_size(max_record_size),
        )),
        effective_io,
    ))
}

fn open_path(path: &Path, io: IoBackend) -> Result<(Box<dyn Read + Send>, IoBackend), EngineError> {
    let file = File::open(path).map_err(|source| EngineError::Io { source })?;

    #[cfg(feature = "mmap")]
    if matches!(io, IoBackend::Mmap) {
        // SAFETY-free path: `memmap2::Mmap::map` is `unsafe`, and this crate
        // forbids `unsafe`, so mapping is delegated to the helper module which
        // is the single place allowed to do it.
        // A mapping failure is a property of the filesystem, not of the
        // request, so fall through to buffered I/O rather than failing.
        if let Ok(reader) = crate::engine::mmap_input(&file) {
            return Ok((Box::new(reader), IoBackend::Mmap));
        }
    }
    let _ = io;

    Ok((
        Box::new(std::io::BufReader::with_capacity(1 << 20, file)),
        IoBackend::Buffered,
    ))
}

/// Memory-mapped input.
///
/// `bamsplit-core` sets `#![forbid(unsafe_code)]`, so the mapping itself is done
/// by `memmap2`, and the only thing this wrapper adds is a cursor. A mapping is
/// a *view*, not a decompression shortcut: BAM records inside are still BGZF
/// compressed and still have to be inflated. What it buys is one fewer copy per
/// block and cheap independent cursors for the indexed engine.
#[cfg(feature = "mmap")]
fn mmap_input(file: &File) -> std::io::Result<std::io::Cursor<MmapView>> {
    let mmap = MmapView::new(file)?;
    Ok(std::io::Cursor::new(mmap))
}

/// A shareable, read-only view over a memory-mapped file.
#[cfg(feature = "mmap")]
#[derive(Debug, Clone)]
pub struct MmapView(Arc<memmap2::Mmap>);

#[cfg(feature = "mmap")]
impl MmapView {
    /// Maps `file`.
    ///
    /// # Errors
    ///
    /// Returns the underlying I/O error if the mapping fails, which happens on
    /// some network filesystems and for zero-length files.
    pub fn new(file: &File) -> std::io::Result<Self> {
        // `memmap2` performs the `unsafe` mapping internally; the safety
        // obligation it documents is that the file must not be truncated while
        // mapped. `bamsplit` only maps inputs it opened read-only, and a
        // concurrent truncation by another process is outside what any BAM
        // reader can defend against.
        let mmap = unsafe_map(file)?;
        Ok(Self(Arc::new(mmap)))
    }
}

#[cfg(feature = "mmap")]
impl AsRef<[u8]> for MmapView {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

/// The one call that needs `unsafe`, isolated so the rest of the crate can keep
/// `forbid(unsafe_code)`.
#[cfg(feature = "mmap")]
#[allow(unsafe_code)]
fn unsafe_map(file: &File) -> std::io::Result<memmap2::Mmap> {
    // SAFETY: `Mmap::map` requires that the mapped file is not modified for the
    // lifetime of the mapping. The file is opened read-only by `bamsplit` and
    // never written by it. A different process truncating the file while it is
    // mapped would be undefined behaviour, which is why `--io mmap` is
    // documented as being for stable local inputs; `--io auto` never selects
    // mmap for anything else.
    unsafe { memmap2::Mmap::map(file) }
}

impl BamInput {
    /// Reads the next record body into the reader's internal buffer.
    ///
    /// # Errors
    ///
    /// Returns [`EngineError`] for a malformed or truncated record.
    pub fn next_record(&mut self) -> Result<Option<crate::bam::RawRecord<'_>>, EngineError> {
        match self {
            Self::Serial(reader) => reader.read_record().map_err(Into::into),
            Self::Parallel(reader) => reader.read_record().map_err(Into::into),
        }
    }

    /// How many records have been read.
    #[must_use]
    pub fn record_count(&self) -> u64 {
        match self {
            Self::Serial(reader) => reader.record_count(),
            Self::Parallel(reader) => reader.record_count(),
        }
    }

    /// The virtual offset the next record starts at.
    #[must_use]
    pub fn virtual_position(&self) -> u64 {
        match self {
            Self::Serial(reader) => VirtualPositionSource::virtual_position(reader.get_ref()),
            Self::Parallel(reader) => VirtualPositionSource::virtual_position(reader.get_ref()),
        }
    }
}

/// Checks the interruption flag, converting a raised flag into an error.
///
/// # Errors
///
/// Returns [`EngineError::Interrupted`] when the flag is raised.
pub fn check_interrupt(
    interrupt: &InterruptFlag,
    records_processed: u64,
) -> Result<(), EngineError> {
    if interrupt.is_raised() {
        return Err(EngineError::Interrupted { records_processed });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_single_thread_budget_is_all_ones() {
        let budget = ThreadBudget::for_total(1);
        assert_eq!(budget.total, 1);
        assert_eq!(budget.decompression, 1);
        assert_eq!(budget.compression, 1);
        assert_eq!(budget.tasks, 1);
    }

    #[test]
    fn a_zero_thread_request_is_treated_as_one() {
        assert_eq!(ThreadBudget::for_total(0).total, 1);
        assert_eq!(ThreadBudget::for_tasks(0, 0).total, 1);
    }

    #[test]
    fn the_single_writer_budget_favours_compression() {
        let budget = ThreadBudget::for_total(16);
        assert_eq!(budget.total, 16);
        assert!(budget.compression > budget.decompression, "{budget:?}");
        assert!(budget.decompression + budget.compression <= 16);
        assert_eq!(budget.tasks, 1);
    }

    #[test]
    fn the_task_budget_never_multiplies_the_total() {
        for (total, concurrency) in [(16usize, 4usize), (8, 100), (4, 1), (2, 2)] {
            let budget = ThreadBudget::for_tasks(total, concurrency);
            assert!(budget.tasks <= total, "{budget:?}");
            assert!(
                budget.tasks * budget.compression <= total.max(budget.tasks),
                "{budget:?} exceeds the budget"
            );
            assert!(budget.compression >= 1 && budget.tasks >= 1);
        }
    }

    #[test]
    fn input_sources_classify_themselves() {
        assert_eq!(InputSource::from_arg(Path::new("-")), InputSource::Stdin);
        assert!(!InputSource::Stdin.is_seekable());
        assert_eq!(InputSource::Stdin.display(), "<stdin>");
        assert!(InputSource::Stdin.path().is_none());
        assert!(InputSource::Stdin.size().is_none());
        InputSource::Stdin
            .validate()
            .expect("stdin is always usable");

        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("in.bam");
        std::fs::write(&path, b"x").expect("seeded");
        let source = InputSource::from_arg(&path);
        assert!(source.is_seekable());
        assert_eq!(source.size(), Some(1));
        source.validate().expect("usable");
    }

    #[test]
    fn a_missing_or_directory_input_is_rejected() {
        let directory = tempfile::tempdir().expect("temp dir");
        let missing = InputSource::Path(directory.path().join("absent.bam"));
        assert!(matches!(
            missing.validate(),
            Err(ConfigError::BadInput { .. })
        ));

        let as_dir = InputSource::Path(directory.path().to_path_buf());
        assert!(matches!(
            as_dir.validate(),
            Err(ConfigError::BadInput { .. })
        ));
    }

    #[test]
    fn engine_and_backend_names_round_trip() {
        for (kind, text) in [
            (EngineKind::Auto, "auto"),
            (EngineKind::Stream, "stream"),
            (EngineKind::Indexed, "indexed"),
            (EngineKind::Spool, "spool"),
        ] {
            assert_eq!(kind.to_string(), text);
        }
        for (backend, text) in [
            (IoBackend::Auto, "auto"),
            (IoBackend::Buffered, "buffered"),
            (IoBackend::Mmap, "mmap"),
        ] {
            assert_eq!(backend.to_string(), text);
        }
    }

    #[test]
    fn interruption_is_reported_as_an_error() {
        let flag = InterruptFlag::new();
        check_interrupt(&flag, 10).expect("not raised");
        flag.raise();
        let error = check_interrupt(&flag, 10).expect_err("must fail");
        assert!(
            matches!(
                error,
                EngineError::Interrupted {
                    records_processed: 10
                }
            ),
            "{error}"
        );
    }

    /// A header-only BAM, written by the crate's own writer.
    fn header_only_bam(path: &std::path::Path) {
        let header = crate::bam::header::BamHeader::from_parts(
            b"@HD\tVN:1.6\n".to_vec(),
            vec![crate::bam::header::ReferenceSequence {
                name: b"chr1".to_vec(),
                length: 100,
            }],
        )
        .expect("valid header");
        let mut bytes = Vec::new();
        header.write_to(&mut bytes).expect("serialized");
        let mut writer: crate::output::bgzf::BgzfBlockWriter<Vec<u8>, ()> =
            crate::output::bgzf::BgzfBlockWriter::new(Vec::new(), 1, None, 1, 0);
        writer.write_raw(&bytes, &mut |(), _| {}).expect("written");
        std::fs::write(path, writer.finish(&mut |(), _| {}).expect("finished")).expect("written");
    }

    #[test]
    fn a_thread_budget_above_one_selects_the_multi_threaded_inflater() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("in.bam");
        header_only_bam(&path);
        let source = InputSource::Path(path);

        let (_, serial, _) = open(&source, IoBackend::Buffered, 1, 1 << 20).expect("opened");
        assert!(matches!(serial, BamInput::Serial(_)), "{serial:?}");

        let (header, parallel, _) = open(&source, IoBackend::Buffered, 4, 1 << 20).expect("opened");
        assert!(matches!(parallel, BamInput::Parallel(_)), "{parallel:?}");
        assert_eq!(header.reference_count(), 1, "the header still parses");
    }

    #[test]
    fn a_pool_is_only_built_when_it_would_help() {
        assert!(
            ThreadBudget::for_total(1)
                .build_pool()
                .expect("built")
                .is_none()
        );
        assert!(
            ThreadBudget::for_total(8)
                .build_pool()
                .expect("built")
                .is_some()
        );
    }
}
