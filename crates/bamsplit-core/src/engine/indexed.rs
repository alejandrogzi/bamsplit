// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The indexed parallel engine: independent per-reference queries.
//!
//! # Two query plans
//!
//! * **Per reference.** Each reference's records occupy a known set of BGZF
//!   chunks, so `chr1` and `chr2` can be extracted at the same time by two
//!   threads that never touch each other's bytes. This is what `bamsplit chrom`
//!   uses.
//! * **Per region.** For `bamsplit region` the interesting intervals are
//!   sub-ranges of a reference, so the plan queries *those* rather than the
//!   whole reference. A sparse annotation — a few hundred genes across a
//!   3 Gbp genome — then reads a few per-cent of the file instead of all of it,
//!   which is a different and much larger win than parallelism alone.
//!
//! Both plans partition by reference, which is what makes them safe: a region
//! belongs to exactly one reference, so every output is written by exactly one
//! task and no two tasks ever contend for a file.
//!
//! # When it wins, and when it does not
//!
//! It is **not** unconditionally faster, and `bamsplit` does not pretend
//! otherwise:
//!
//! * Chunk boundaries do not align with references, so a block holding the last
//!   `chr1` record and the first `chr2` record is inflated twice — once by each
//!   task. On an input with thousands of small contigs that duplicated work
//!   dominates.
//! * Random access defeats read-ahead. On a network filesystem the streaming
//!   engine's single sequential pass usually wins outright.
//! * The unplaced records still need a sequential tail scan.
//!
//! [`crate::planning`] weighs those factors; `--engine indexed` forces the issue
//! and is validated rather than silently ignored.
//!
//! # Correctness
//!
//! Index chunks are *conservative*: a chunk may contain records of a
//! neighbouring reference, and two chunks may overlap. So every record is
//! filtered by its actual `refID` after being read, and chunks are merged and
//! optimized before use so no record is read — and therefore emitted — twice.
//! `tests/integration` asserts the indexed engine's per-output digests equal the
//! streaming engine's.
//!
//! Querying a region's envelope is sufficient for every assignment mode: a
//! record that `start`, `midpoint`, `contained`, `best-overlap`, or `overlap`
//! could assign to that region must intersect its envelope, and a binning index
//! returns every record intersecting the queried interval.
//!
//! # What the record counts mean here
//!
//! Under the per-region plan the engine deliberately never reads records
//! outside any region, so the manifest's `input_records` counts *records
//! examined* rather than records in the file, and `unmatched_records` counts
//! only those that were read and matched nothing. Conservation still holds
//! among the records it saw. The report says so in its notes, because the same
//! split under `--engine stream` reports the file's true totals and the two
//! numbers would otherwise look like a discrepancy.

use std::fs::File;
use std::path::Path;

use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::bam::header::BamHeader;
use crate::bam::raw_record::RawRecordReader;
use crate::engine::{
    EngineKind, ExecutionContext, ExecutionReport, INTERRUPT_CHECK_INTERVAL, IoBackend,
    SplitEngine, check_interrupt,
};
use crate::error::{ConfigError, EngineError, IndexError};
use crate::index::{bai, csi};
use crate::interval::Interval;
use crate::output::OutputManager;
use crate::output::transaction::Transaction;
use crate::output::writer::{FinishedOutput, RecordWriter};
use crate::routing::{Route, Router, RoutingKey};
use crate::stats::RunStats;

/// An input index, in whichever format was found.
pub enum InputIndex {
    /// A `.bai` alongside the BAM.
    Bai(Box<bai::Index>),
    /// A `.csi` alongside the BAM.
    Csi(Box<csi::Index>),
}

impl std::fmt::Debug for InputIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bai(_) => f.write_str("InputIndex::Bai"),
            Self::Csi(_) => f.write_str("InputIndex::Csi"),
        }
    }
}

impl InputIndex {
    /// The format name, for the manifest.
    #[must_use]
    pub const fn kind(&self) -> &'static str {
        match self {
            Self::Bai(_) => "bai",
            Self::Csi(_) => "csi",
        }
    }

    /// The underlying binning index.
    #[must_use]
    pub fn as_binning_index(&self) -> &dyn BinningIndex {
        match self {
            Self::Bai(index) => index.as_ref(),
            Self::Csi(index) => index.as_ref(),
        }
    }

    /// How many references the index describes.
    #[must_use]
    pub fn reference_count(&self) -> usize {
        BinningIndex::reference_sequences(self.as_binning_index()).count()
    }
}

/// Looks for an index next to `bam`.
///
/// Both `sample.bam.bai` and `sample.bai` are checked, in that order, matching
/// what `samtools index` can produce.
#[must_use]
pub fn find_index(bam: &Path) -> Option<std::path::PathBuf> {
    let mut suffixed = bam.as_os_str().to_os_string();
    suffixed.push(".bai");
    let candidates = [
        std::path::PathBuf::from(&suffixed),
        {
            let mut csi = bam.as_os_str().to_os_string();
            csi.push(".csi");
            std::path::PathBuf::from(csi)
        },
        bam.with_extension("bai"),
        bam.with_extension("csi"),
    ];
    candidates.into_iter().find(|path| path.is_file())
}

/// Loads an index and checks it against the header.
///
/// # Errors
///
/// Returns [`IndexError`] if the file cannot be parsed, or if it describes a
/// different number of references than the header declares — which means the
/// index is stale and using it would silently lose records.
pub fn load_index(path: &Path, header: &BamHeader) -> Result<InputIndex, IndexError> {
    let index = match path.extension().and_then(std::ffi::OsStr::to_str) {
        Some("csi") => InputIndex::Csi(Box::new(csi::read(path)?)),
        _ => InputIndex::Bai(Box::new(bai::read(path)?)),
    };
    let declared = index.reference_count();
    if declared != header.reference_count() {
        return Err(IndexError::IndexHeaderMismatch {
            path: path.to_path_buf(),
            index_references: declared,
            header_references: header.reference_count(),
        });
    }
    Ok(index)
}

/// The intervals a per-region plan will query, grouped by reference.
///
/// Intervals are 0-based half-open, sorted, and merged, so the chunk sets they
/// produce are as small as the index allows.
#[derive(Debug, Clone, Default)]
pub struct QueryTargets {
    by_reference: std::collections::BTreeMap<usize, Vec<Interval>>,
}

impl QueryTargets {
    /// Builds targets from `(reference id, interval)` pairs.
    ///
    /// Overlapping and abutting intervals on one reference are merged: querying
    /// them separately would return overlapping chunk sets, and merging here is
    /// cheaper than deduplicating chunks later.
    #[must_use]
    pub fn from_intervals(pairs: impl IntoIterator<Item = (usize, Interval)>) -> Self {
        let mut by_reference: std::collections::BTreeMap<usize, Vec<Interval>> =
            std::collections::BTreeMap::new();
        for (reference_id, interval) in pairs {
            if !interval.is_empty() {
                by_reference.entry(reference_id).or_default().push(interval);
            }
        }
        for intervals in by_reference.values_mut() {
            crate::interval::merge_in_place(intervals);
        }
        Self { by_reference }
    }

    /// The references that have at least one target.
    pub fn references(&self) -> impl Iterator<Item = usize> + '_ {
        self.by_reference.keys().copied()
    }

    /// The intervals on one reference.
    #[must_use]
    pub fn intervals(&self, reference_id: usize) -> Option<&[Interval]> {
        self.by_reference
            .get(&reference_id)
            .map(std::vec::Vec::as_slice)
    }

    /// How many intervals will be queried in total.
    #[must_use]
    pub fn len(&self) -> usize {
        self.by_reference.values().map(Vec::len).sum()
    }

    /// Whether there is nothing to query.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.by_reference.is_empty()
    }
}

/// The indexed parallel engine.
#[derive(Debug, Default)]
pub struct IndexedEngine {
    /// The index to use, loaded by the planner.
    pub index: Option<InputIndex>,
    /// The intervals to restrict queries to, for the per-region plan.
    ///
    /// [`None`] means the per-reference plan: read every reference in full.
    pub targets: Option<QueryTargets>,
}

impl IndexedEngine {
    /// An engine that extracts each reference in full.
    ///
    /// Requires a router whose key is a function of `refID`, so that every
    /// output belongs to exactly one task.
    #[must_use]
    pub fn per_reference(index: InputIndex) -> Self {
        Self {
            index: Some(index),
            targets: None,
        }
    }

    /// An engine that reads only the given intervals.
    ///
    /// Safe for any router whose keys are confined to one reference each, which
    /// region routing satisfies: a region names one reference.
    #[must_use]
    pub fn per_region(index: InputIndex, targets: QueryTargets) -> Self {
        Self {
            index: Some(index),
            targets: Some(targets),
        }
    }

    /// Whether this engine is reading only selected intervals.
    #[must_use]
    pub const fn is_per_region(&self) -> bool {
        self.targets.is_some()
    }
}

impl<R: Router> SplitEngine<R> for IndexedEngine {
    fn execute(
        &self,
        context: &ExecutionContext<'_>,
        router: &R,
    ) -> Result<ExecutionReport, EngineError> {
        let Some(index) = self.index.as_ref() else {
            return Err(EngineError::from(ConfigError::IncompatibleEngine {
                engine: "indexed",
                reason: "no BAI or CSI index was found for the input".to_string(),
            }));
        };
        let Some(bam_path) = context.input.path().map(Path::to_path_buf) else {
            return Err(EngineError::from(ConfigError::IncompatibleEngine {
                engine: "indexed",
                reason: "standard input cannot be seeked, so an index is unusable".to_string(),
            }));
        };
        // The per-reference plan needs every output to belong to one reference,
        // which only a `refID`-derived key guarantees. The per-region plan
        // carries that guarantee in its targets instead.
        if self.targets.is_none() && !router.is_grouped_by_coordinate() {
            return Err(EngineError::from(ConfigError::IncompatibleEngine {
                engine: "indexed",
                reason: format!(
                    "`{}` routing cannot be decomposed into independent reference queries",
                    router.mode()
                ),
            }));
        }

        let (input_header, _, _) = crate::engine::open_bam(
            &context.input,
            IoBackend::Buffered,
            1,
            context.max_record_size,
        )?;
        crate::output::transaction::prepare_directory(&context.output_directory)
            .map_err(Box::new)?;

        // One task per reference, bounded by the thread budget. Tasks are
        // completely independent: separate file handles, separate outputs, no
        // shared mutable state, so no work-stealing coordination is needed
        // beyond the pool itself.
        let candidates: Vec<usize> = match &self.targets {
            Some(targets) => targets.references().collect(),
            None => (0..input_header.reference_count()).collect(),
        };
        let mut tasks: Vec<ReferenceTask> = Vec::with_capacity(candidates.len());
        for reference_id in candidates {
            let Some(length) = input_header.reference_length(reference_id) else {
                continue;
            };
            let intervals = self
                .targets
                .as_ref()
                .and_then(|targets| targets.intervals(reference_id));
            let chunks = merged_chunks(index, reference_id, u64::from(length), intervals)?;
            if !chunks.is_empty() {
                tasks.push(ReferenceTask {
                    reference_id,
                    chunks,
                });
            }
        }

        let mut outputs = Vec::new();
        let mut stats = RunStats::default();
        let mut notes = match &self.targets {
            Some(targets) => vec![
                format!(
                    "queried {} interval(s) across {} reference sequence(s) from the {} index",
                    targets.len(),
                    tasks.len(),
                    index.kind()
                ),
                "record counts cover only the records these queries examined, not the whole \
                 input; the same split under `--engine stream` reports the file's totals"
                    .to_string(),
            ],
            None => vec![format!(
                "queried {} reference sequences from the {} index",
                tasks.len(),
                index.kind()
            )],
        };

        let results: Vec<Result<TaskOutcome<R::Key>, EngineError>> = if context.threads.tasks > 1 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(context.threads.tasks)
                .thread_name(|index| format!("bamsplit-q{index}"))
                .build()
                .map_err(|source| EngineError::WorkerFailure {
                    reason: source.to_string(),
                })?;
            pool.install(|| {
                use rayon::prelude::*;
                tasks
                    .par_iter()
                    .map(|task| run_task(context, router, &input_header, &bam_path, task))
                    .collect()
            })
        } else {
            tasks
                .iter()
                .map(|task| run_task(context, router, &input_header, &bam_path, task))
                .collect()
        };

        let mut pending: Vec<TaskOutcome<R::Key>> = Vec::with_capacity(results.len());
        for result in results {
            match result {
                Ok(outcome) => pending.push(outcome),
                Err(error) => {
                    for outcome in pending {
                        for path in outcome.committed {
                            let _ = std::fs::remove_file(path);
                        }
                    }
                    return Err(error);
                }
            }
        }

        // Ordering by reference id rather than completion order keeps the
        // manifest reproducible regardless of how the pool scheduled the work.
        pending.sort_by_key(|outcome| outcome.reference_id);
        for outcome in pending {
            stats.merge(&outcome.stats);
            outputs.extend(outcome.outputs);
        }

        // The unplaced records are not addressable by a binning index, so they
        // need one sequential tail scan. The index tells us whether there are
        // any at all, which lets the common case skip the scan entirely — and
        // the per-region plan skips it always, because an unplaced record has no
        // coordinate and so can never fall inside a region.
        let unplaced_declared = if self.targets.is_some() {
            0
        } else {
            index
                .as_binning_index()
                .unplaced_unmapped_record_count()
                .unwrap_or(0)
        };
        if unplaced_declared > 0 {
            notes.push(format!(
                "scanned the tail of the input for {unplaced_declared} unplaced records"
            ));
            let outcome = run_unplaced(context, router, &input_header, &bam_path, index)?;
            stats.merge(&outcome.stats);
            outputs.extend(outcome.outputs);
        }

        stats.total_output_emissions = outputs.iter().map(|output| output.stats.record_count).sum();

        let mut skipped_keys = Vec::new();
        for key in router.declared_keys(&input_header) {
            if outputs
                .iter()
                .any(|output| output.paths.logical_key == key.logical())
            {
                continue;
            }
            if context.emit_empty {
                let mut encoder = crate::output::FilenameEncoder::new(
                    &context.output_directory,
                    context.manager_options.template.clone(),
                );
                let paths = encoder
                    .assign(
                        key.logical(),
                        context
                            .manager_options
                            .settings
                            .index_kind
                            .map(crate::index::IndexKind::extension),
                    )
                    .map_err(Box::new)?;
                let writer = RecordWriter::create(
                    paths,
                    context.header,
                    input_header.reference_count(),
                    &context.manager_options.settings,
                )
                .map_err(Box::new)?;
                let mut transaction = Transaction::new();
                let finished = writer.finish(&mut transaction).map_err(Box::new)?;
                transaction.commit().map_err(Box::new)?;
                outputs.push(finished);
            } else {
                skipped_keys.push(key.logical().to_vec());
            }
        }

        Ok(ExecutionReport {
            outputs,
            stats,
            skipped_keys,
            engine: EngineKind::Indexed,
            io: IoBackend::Buffered,
            notes,
        })
    }
}

#[derive(Debug)]
struct ReferenceTask {
    reference_id: usize,
    chunks: Vec<Chunk>,
}

/// The paths one task committed, so another task's failure can roll them back.
fn committed_paths(outputs: &[FinishedOutput]) -> Vec<std::path::PathBuf> {
    outputs
        .iter()
        .flat_map(|output| {
            std::iter::once(output.paths.bam.clone()).chain(output.index_path.clone())
        })
        .collect()
}

struct TaskOutcome<K> {
    reference_id: usize,
    outputs: Vec<FinishedOutput>,
    stats: RunStats,
    committed: Vec<std::path::PathBuf>,
    _key: std::marker::PhantomData<K>,
}

/// The chunks covering a reference, or just the given intervals of it.
///
/// Every query bound is clamped to what the index can address. `noodles` rejects
/// a bound beyond `2^(min_shift + 3*depth) - 1` outright, so asking for
/// `Position::MAX` would make every query fail — and, because a failed query is
/// indistinguishable from an empty one at this level, would silently produce no
/// output at all. Clamping is what keeps the whole reference covered.
///
/// Chunks from several intervals are merged into one non-overlapping set. That
/// matters for correctness, not just speed: two neighbouring regions routinely
/// share a BGZF block, and reading it twice would emit its records twice.
fn merged_chunks(
    index: &InputIndex,
    reference_id: usize,
    reference_length: u64,
    targets: Option<&[Interval]>,
) -> Result<Vec<Chunk>, IndexError> {
    use noodles_core::Position;
    use noodles_core::region::Interval as CoreInterval;

    let binning = index.as_binning_index();
    let addressable = csi::addressable(binning.min_shift(), binning.depth());
    let limit = reference_length.max(1).min(addressable);

    // 1-based inclusive query bounds, clamped to the addressable space.
    let bounds: Vec<(u64, u64)> = match targets {
        None => vec![(1, limit)],
        Some(intervals) => intervals
            .iter()
            .filter(|interval| !interval.is_empty())
            .map(|interval| {
                let (start, end) = interval.one_based_inclusive();
                let start = (start.max(1) as u64).min(limit);
                let end = (end.max(1) as u64).min(limit);
                (start, end.max(start))
            })
            .collect(),
    };

    let mut raw = Vec::new();
    for (start, end) in bounds {
        let (Some(start), Some(end)) = (
            usize::try_from(start).ok().and_then(Position::new),
            usize::try_from(end).ok().and_then(Position::new),
        ) else {
            continue;
        };
        raw.extend(
            binning
                .query(reference_id, CoreInterval::from(start..=end))
                .map_err(|source| IndexError::NotIndexable {
                    path: std::path::PathBuf::new(),
                    reason: format!("querying reference {reference_id} failed: {source}"),
                })?,
        );
    }

    // `merge_chunks` sorts internally, but sorting here too makes the input
    // order irrelevant to the result and therefore reproducible.
    raw.sort_unstable_by_key(|chunk| (u64::from(chunk.start()), u64::from(chunk.end())));
    Ok(noodles_csi::binning_index::merge_chunks(&raw))
}

/// Extracts one reference — or the requested intervals of it — into outputs.
///
/// Each task owns its own [`OutputManager`], so descriptors stay inside
/// `--max-open-files` *per task* rather than growing with the number of regions
/// on a reference. Two tasks never share an output, because a region and a
/// reference key both belong to exactly one reference.
fn run_task<R: Router>(
    context: &ExecutionContext<'_>,
    router: &R,
    header: &BamHeader,
    bam_path: &Path,
    task: &ReferenceTask,
) -> Result<TaskOutcome<R::Key>, EngineError> {
    let mut stats = RunStats::default();
    let mut manager = OutputManager::new(
        &context.output_directory,
        context.header,
        context.manager_options.clone(),
    )
    .map_err(Box::new)?;

    let file = File::open(bam_path).map_err(|source| EngineError::Io { source })?;
    let mut reader = RawRecordReader::new(noodles_bgzf::io::Reader::new(file))
        .with_max_record_size(context.max_record_size);
    let mut records = 0u64;

    let outcome = (|| -> Result<(), EngineError> {
        for chunk in &task.chunks {
            reader
                .get_mut()
                .seek(chunk.start())
                .map_err(|source| EngineError::Io { source })?;
            loop {
                if u64::from(reader.get_mut().virtual_position()) >= u64::from(chunk.end()) {
                    break;
                }
                let Some(record) = reader.read_record().map_err(Box::new)? else {
                    break;
                };
                records += 1;
                if records.is_multiple_of(INTERRUPT_CHECK_INTERVAL) {
                    check_interrupt(context.interrupt, records)?;
                }
                // Chunks are conservative, so a record from a neighbouring
                // reference can appear here; skip it without counting, because
                // the task that owns it will see it too.
                if record.reference_sequence_id().map_err(Box::new)?
                    != Some(task.reference_id as i32)
                {
                    continue;
                }
                let route = router
                    .route(header, &record)
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
                    Route::One(_) | Route::Many(_) => {
                        for key in route.keys() {
                            manager.write(key.logical(), &record).map_err(Box::new)?;
                        }
                    }
                }
                stats.observe_emissions(route.emission_count() as u64);
            }
        }
        Ok(())
    })();

    if let Err(error) = outcome {
        manager.abort();
        return Err(error);
    }

    let outputs = manager.finish_all().map_err(Box::new)?;
    let committed = committed_paths(&outputs);
    Ok(TaskOutcome {
        reference_id: task.reference_id,
        outputs,
        stats,
        committed,
        _key: std::marker::PhantomData,
    })
}

/// Scans the tail of the input for unplaced records.
fn run_unplaced<R: Router>(
    context: &ExecutionContext<'_>,
    router: &R,
    header: &BamHeader,
    bam_path: &Path,
    index: &InputIndex,
) -> Result<TaskOutcome<R::Key>, EngineError> {
    let mut stats = RunStats::default();
    let mut manager = OutputManager::new(
        &context.output_directory,
        context.header,
        context.manager_options.clone(),
    )
    .map_err(Box::new)?;

    let file = File::open(bam_path).map_err(|source| EngineError::Io { source })?;
    let mut reader = RawRecordReader::new(noodles_bgzf::io::Reader::new(file))
        .with_max_record_size(context.max_record_size);

    // Start from the last reference's first record: everything unplaced sorts
    // after every placed record in a coordinate-sorted BAM, so this skips the
    // bulk of the file without needing a full scan.
    let start = if let Some(position) = index.as_binning_index().last_first_record_start_position()
    {
        position
    } else {
        // No placed records at all, so start just after the header.
        let mut probe = noodles_bgzf::io::Reader::new(
            File::open(bam_path).map_err(|source| EngineError::Io { source })?,
        );
        BamHeader::read_from(&mut probe).map_err(Box::new)?;
        probe.virtual_position()
    };
    reader
        .get_mut()
        .seek(start)
        .map_err(|source| EngineError::Io { source })?;

    let mut records = 0u64;
    let outcome = (|| -> Result<(), EngineError> {
        while let Some(record) = reader.read_record().map_err(Box::new)? {
            records += 1;
            if records.is_multiple_of(INTERRUPT_CHECK_INTERVAL) {
                check_interrupt(context.interrupt, records)?;
            }
            if record.reference_sequence_id().map_err(Box::new)?.is_some() {
                continue;
            }
            let route = router
                .route(header, &record)
                .map_err(|error| Box::new(error.relocate_routing(records)))?;
            if let Route::Drop(reason) = &route {
                if reason.is_unmatched() {
                    stats.observe_emissions(0);
                } else {
                    stats.observe_drop();
                }
                continue;
            }
            for key in route.keys() {
                manager.write(key.logical(), &record).map_err(Box::new)?;
            }
            stats.observe_emissions(route.emission_count() as u64);
        }
        Ok(())
    })();

    if let Err(error) = outcome {
        manager.abort();
        return Err(error);
    }

    let outputs = manager.finish_all().map_err(Box::new)?;
    let committed = committed_paths(&outputs);
    Ok(TaskOutcome {
        reference_id: usize::MAX,
        outputs,
        stats,
        committed,
        _key: std::marker::PhantomData,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_discovery_prefers_the_suffixed_name() {
        let directory = tempfile::tempdir().expect("temp dir");
        let bam = directory.path().join("sample.bam");
        std::fs::write(&bam, b"x").expect("seeded");
        assert!(find_index(&bam).is_none());

        let suffixed = directory.path().join("sample.bam.bai");
        std::fs::write(&suffixed, b"x").expect("seeded");
        assert_eq!(find_index(&bam), Some(suffixed.clone()));

        std::fs::remove_file(&suffixed).expect("removed");
        let csi = directory.path().join("sample.bam.csi");
        std::fs::write(&csi, b"x").expect("seeded");
        assert_eq!(find_index(&bam), Some(csi));
    }

    #[test]
    fn a_stale_index_is_rejected() {
        use crate::bam::header::ReferenceSequence;
        use crate::index::{IndexBuilder, IndexKind};

        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("sample.bam.bai");
        // An index over two references.
        IndexBuilder::new(IndexKind::Bai, 2)
            .finish(&path)
            .expect("written");

        // A header declaring three.
        let header = BamHeader::from_parts(
            b"@HD\tVN:1.6\n".to_vec(),
            (0..3)
                .map(|index| ReferenceSequence {
                    name: format!("chr{index}").into_bytes(),
                    length: 1000,
                })
                .collect(),
        )
        .expect("valid header");

        let error = load_index(&path, &header).expect_err("must reject");
        assert!(
            matches!(
                error,
                IndexError::IndexHeaderMismatch {
                    index_references: 2,
                    header_references: 3,
                    ..
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn an_index_matching_the_header_is_accepted() {
        use crate::bam::header::ReferenceSequence;
        use crate::index::{IndexBuilder, IndexKind};

        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("sample.bam.bai");
        IndexBuilder::new(IndexKind::Bai, 2)
            .finish(&path)
            .expect("written");

        let header = BamHeader::from_parts(
            b"@HD\tVN:1.6\n".to_vec(),
            (0..2)
                .map(|index| ReferenceSequence {
                    name: format!("chr{index}").into_bytes(),
                    length: 1000,
                })
                .collect(),
        )
        .expect("valid header");

        let index = load_index(&path, &header).expect("accepted");
        assert_eq!(index.kind(), "bai");
        assert_eq!(index.reference_count(), 2);
    }
}
