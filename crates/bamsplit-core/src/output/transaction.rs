// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Transactional output: nothing appears at its final path until it is
//! complete.
//!
//! # The guarantee
//!
//! A `bamsplit` run must never leave behind a file that *looks* finished but is
//! not. A truncated `chr1.bam` in a directory of good ones is worse than no
//! output at all, because the next pipeline stage will happily consume it.
//!
//! So every output is written to
//!
//! ```text
//! <final-name>.part.<pid>.<counter><nonce>
//! ```
//!
//! and renamed into place only after the BAM header, all records, the BGZF
//! end-of-file marker, the index, and the statistics are all done and
//! validated. `rename(2)` within one directory is atomic on every filesystem
//! `bamsplit` targets, so a reader sees either the old file or the complete new
//! one.
//!
//! # Cleanup
//!
//! Three layers, deliberately overlapping:
//!
//! 1. [`PendingFile`] removes its temporary file on `Drop` unless it was
//!    committed, so an early `return` or a `?` cannot leak.
//! 2. A process-wide [`registry`] tracks every live temporary path, so a signal
//!    handler can sweep paths whose owning value never got to run `Drop`.
//! 3. [`InterruptFlag`] lets the record loop notice a signal and unwind
//!    normally, which is what actually triggers layer 1.
//!
//! Layer 2 exists because a signal handler cannot safely unwind the stack;
//! layer 1 exists because most failures are not signals.

use std::collections::BTreeSet;
use std::fs::{File, OpenOptions};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};

use crate::error::OutputError;

/// A process-wide set of temporary paths that are currently live.
///
/// Only ever holds paths this process created, and only for as long as they
/// exist, so sweeping it can never delete a user's file.
pub mod registry {
    use super::{BTreeSet, Mutex, OnceLock, PathBuf};

    fn store() -> &'static Mutex<BTreeSet<PathBuf>> {
        static STORE: OnceLock<Mutex<BTreeSet<PathBuf>>> = OnceLock::new();
        STORE.get_or_init(|| Mutex::new(BTreeSet::new()))
    }

    /// Records a temporary path.
    pub fn register(path: &std::path::Path) {
        if let Ok(mut paths) = store().lock() {
            paths.insert(path.to_path_buf());
        }
    }

    /// Forgets a temporary path, because it was committed or already removed.
    pub fn unregister(path: &std::path::Path) {
        if let Ok(mut paths) = store().lock() {
            paths.remove(path);
        }
    }

    /// Removes every registered temporary path.
    ///
    /// Returns how many were removed. Errors are swallowed: this runs during
    /// teardown, when there is nowhere useful to report them, and a path that
    /// cannot be removed is not a reason to stop removing the others.
    ///
    /// This is process-wide by design — a signal handler has no way to reach
    /// the values that own those paths — so it must only ever be called while
    /// tearing the process down.
    #[must_use]
    pub fn sweep() -> usize {
        let Ok(mut paths) = store().lock() else {
            return 0;
        };
        remove_all(std::mem::take(&mut *paths))
    }

    /// Removes a specific set of paths and forgets them.
    ///
    /// Exists so the sweep logic can be tested without a global sweep, which
    /// would delete the temporary files of every other test running in the same
    /// process.
    pub fn remove_all<I>(paths: I) -> usize
    where
        I: IntoIterator<Item = PathBuf>,
    {
        let mut removed = 0;
        let mut store_guard = store().lock().ok();
        for path in paths {
            if std::fs::remove_file(&path).is_ok() {
                removed += 1;
            }
            if let Some(paths) = store_guard.as_mut() {
                paths.remove(&path);
            }
        }
        removed
    }

    /// A snapshot of the live temporary paths, for tests and diagnostics.
    #[must_use]
    pub fn snapshot() -> Vec<PathBuf> {
        store()
            .lock()
            .map(|paths| paths.iter().cloned().collect())
            .unwrap_or_default()
    }
}

/// A cooperative interruption flag.
///
/// A signal handler can only do async-signal-safe work, and unwinding is not
/// that. So the handler sets this flag and the record loop, which checks it
/// every few thousand records, returns an error — at which point ordinary
/// `Drop`-based cleanup runs.
#[derive(Debug, Default)]
pub struct InterruptFlag {
    raised: AtomicBool,
}

impl InterruptFlag {
    /// Creates a clear flag.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            raised: AtomicBool::new(false),
        }
    }

    /// Raises the flag. Safe to call from a signal handler.
    pub fn raise(&self) {
        self.raised.store(true, Ordering::SeqCst);
    }

    /// Whether the flag has been raised.
    #[must_use]
    pub fn is_raised(&self) -> bool {
        self.raised.load(Ordering::Relaxed)
    }

    /// Clears the flag, for tests.
    pub fn clear(&self) {
        self.raised.store(false, Ordering::SeqCst);
    }
}

/// The process-wide interruption flag the CLI's signal handler drives.
#[must_use]
pub fn interrupt_flag() -> &'static InterruptFlag {
    static FLAG: InterruptFlag = InterruptFlag::new();
    &FLAG
}

/// Builds a temporary path next to `final_path`.
///
/// The name is unique per process and per call, and it stays in the same
/// directory so the eventual rename is a same-filesystem operation.
#[must_use]
pub fn temporary_path(final_path: &Path) -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);

    let counter = COUNTER.fetch_add(1, Ordering::Relaxed);
    // A per-process nonce keeps two concurrent runs, which share a pid space
    // only across containers, from colliding on the same output directory.
    let nonce = nonce();
    let pid = std::process::id();

    let name = final_path.file_name().map_or_else(
        || std::ffi::OsString::from("output"),
        std::ffi::OsStr::to_os_string,
    );
    let mut temporary = name;
    temporary.push(format!(".part.{pid}.{counter:x}{nonce:08x}"));
    final_path.with_file_name(temporary)
}

fn nonce() -> u32 {
    static NONCE: OnceLock<u32> = OnceLock::new();
    *NONCE.get_or_init(|| {
        let seed = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |elapsed| elapsed.as_nanos() as u64);
        (xxhash_rust::xxh3::xxh3_64(&seed.to_le_bytes()) >> 32) as u32
    })
}

/// A file being written to a temporary path, to be renamed on success.
///
/// Dropping without [`commit`](Self::commit) removes the temporary file.
#[derive(Debug)]
pub struct PendingFile {
    temporary: PathBuf,
    final_path: PathBuf,
    force: bool,
    committed: bool,
}

impl PendingFile {
    /// Creates the temporary file and returns it alongside the handle.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::AlreadyExists`] when `final_path` exists and
    /// `force` is not set, and [`OutputError::Io`] if the temporary file
    /// cannot be created.
    pub fn create(final_path: &Path, force: bool) -> Result<(File, Self), OutputError> {
        if !force && final_path.exists() {
            return Err(OutputError::AlreadyExists {
                path: final_path.to_path_buf(),
            });
        }

        let temporary = temporary_path(final_path);
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|source| OutputError::Io {
                path: temporary.clone(),
                source,
            })?;
        registry::register(&temporary);

        Ok((
            file,
            Self {
                temporary,
                final_path: final_path.to_path_buf(),
                force,
                committed: false,
            },
        ))
    }

    /// Reopens the temporary file for appending.
    ///
    /// Used when an output is parked to stay inside `--max-open-files` and then
    /// resumed.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::Io`] if the file cannot be reopened.
    pub fn reopen_append(&self) -> Result<File, OutputError> {
        OpenOptions::new()
            .append(true)
            .open(&self.temporary)
            .map_err(|source| OutputError::Io {
                path: self.temporary.clone(),
                source,
            })
    }

    /// The temporary path.
    #[must_use]
    pub fn temporary_path(&self) -> &Path {
        &self.temporary
    }

    /// The path the file will be renamed to.
    #[must_use]
    pub fn final_path(&self) -> &Path {
        &self.final_path
    }

    /// The size of the temporary file, in bytes.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::Io`] if the file cannot be inspected.
    pub fn size(&self) -> Result<u64, OutputError> {
        std::fs::metadata(&self.temporary)
            .map(|metadata| metadata.len())
            .map_err(|source| OutputError::Io {
                path: self.temporary.clone(),
                source,
            })
    }

    /// Renames the temporary file into place.
    ///
    /// Re-checks for an existing destination, because another process may have
    /// created it since [`create`](Self::create) looked. That check is
    /// inherently racy on POSIX — `rename(2)` has no "fail if exists" mode
    /// portable enough to rely on — but it turns the common case of two
    /// concurrent runs into a clear error instead of a silent overwrite.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::AlreadyExists`] or [`OutputError::Commit`].
    pub fn commit(mut self) -> Result<PathBuf, OutputError> {
        if !self.force && self.final_path.exists() {
            return Err(OutputError::AlreadyExists {
                path: self.final_path.clone(),
            });
        }
        std::fs::rename(&self.temporary, &self.final_path).map_err(|source| {
            OutputError::Commit {
                temporary: self.temporary.clone(),
                final_path: self.final_path.clone(),
                source,
            }
        })?;
        registry::unregister(&self.temporary);
        self.committed = true;
        Ok(self.final_path.clone())
    }

    /// Removes the temporary file now, rather than at drop.
    pub fn abort(mut self) {
        self.remove();
        self.committed = true;
    }

    fn remove(&mut self) {
        let _ = std::fs::remove_file(&self.temporary);
        registry::unregister(&self.temporary);
    }
}

impl Drop for PendingFile {
    fn drop(&mut self) {
        if !self.committed {
            self.remove();
        }
    }
}

/// A group of files that must appear together or not at all.
///
/// A BAM and its index are one unit: an index without its BAM is useless, and a
/// BAM whose index is stale is dangerous. `commit` renames them in a fixed
/// order — index first, then BAM — so that a crash between the two renames
/// leaves an orphan index rather than an unindexed BAM. An orphan `.bai` is
/// harmless and obvious; a BAM that appears finished but has no index is the
/// failure mode that costs a pipeline an hour.
#[derive(Debug, Default)]
pub struct Transaction {
    files: Vec<PendingFile>,
}

impl Transaction {
    /// Creates an empty transaction.
    #[must_use]
    pub const fn new() -> Self {
        Self { files: Vec::new() }
    }

    /// Adds a pending file. Files commit in the order they are added.
    pub fn push(&mut self, file: PendingFile) {
        self.files.push(file);
    }

    /// How many files are pending.
    #[must_use]
    pub fn len(&self) -> usize {
        self.files.len()
    }

    /// Whether the transaction holds no files.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.files.is_empty()
    }

    /// Commits every file.
    ///
    /// On failure, the files that were already renamed are rolled back by
    /// removing them, and the rest are dropped, which removes their temporary
    /// files. The output directory is then exactly as it was.
    ///
    /// # Errors
    ///
    /// Returns the first [`OutputError`] encountered.
    pub fn commit(self) -> Result<Vec<PathBuf>, OutputError> {
        let mut committed = Vec::with_capacity(self.files.len());
        for file in self.files {
            match file.commit() {
                Ok(path) => committed.push(path),
                Err(error) => {
                    for path in &committed {
                        let _ = std::fs::remove_file(path);
                    }
                    return Err(error);
                }
            }
        }
        Ok(committed)
    }

    /// Discards every file.
    pub fn abort(self) {
        for file in self.files {
            file.abort();
        }
    }
}

/// Prepares an output directory.
///
/// Creates it, including parents, and rejects a path that exists but is not a
/// directory.
///
/// # Errors
///
/// Returns [`OutputError::BadOutputDirectory`].
pub fn prepare_directory(path: &Path) -> Result<(), OutputError> {
    if path.exists() {
        if !path.is_dir() {
            return Err(OutputError::BadOutputDirectory {
                path: path.to_path_buf(),
                reason: "the path exists but is not a directory".to_string(),
            });
        }
        return Ok(());
    }
    std::fs::create_dir_all(path).map_err(|source| OutputError::BadOutputDirectory {
        path: path.to_path_buf(),
        reason: source.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    #[test]
    fn temporary_paths_are_unique_and_local() {
        let final_path = Path::new("/out/chr1.bam");
        let first = temporary_path(final_path);
        let second = temporary_path(final_path);
        assert_ne!(first, second);
        assert_eq!(first.parent(), final_path.parent());
        let name = first.file_name().and_then(|n| n.to_str()).expect("utf-8");
        assert!(name.starts_with("chr1.bam.part."), "{name}");
    }

    #[test]
    fn a_committed_file_appears_at_its_final_path() {
        let directory = tempfile::tempdir().expect("temp dir");
        let final_path = directory.path().join("chr1.bam");

        let (mut file, pending) = PendingFile::create(&final_path, false).expect("created");
        let temporary = pending.temporary_path().to_path_buf();
        file.write_all(b"payload").expect("written");
        drop(file);

        assert!(temporary.exists());
        assert!(!final_path.exists(), "nothing appears before the commit");
        assert!(registry::snapshot().contains(&temporary));

        pending.commit().expect("committed");
        assert!(final_path.exists());
        assert!(!temporary.exists());
        assert!(!registry::snapshot().contains(&temporary));
        assert_eq!(std::fs::read(&final_path).expect("readable"), b"payload");
    }

    #[test]
    fn dropping_without_committing_removes_the_temporary_file() {
        let directory = tempfile::tempdir().expect("temp dir");
        let final_path = directory.path().join("chr1.bam");

        let temporary = {
            let (mut file, pending) = PendingFile::create(&final_path, false).expect("created");
            file.write_all(b"partial").expect("written");
            let temporary = pending.temporary_path().to_path_buf();
            assert!(temporary.exists());
            temporary
        };

        assert!(!temporary.exists(), "the temporary file must be removed");
        assert!(!final_path.exists());
        assert!(!registry::snapshot().contains(&temporary));
    }

    #[test]
    fn an_existing_output_is_refused_without_force() {
        let directory = tempfile::tempdir().expect("temp dir");
        let final_path = directory.path().join("chr1.bam");
        std::fs::write(&final_path, b"old").expect("seeded");

        let error = PendingFile::create(&final_path, false).expect_err("must refuse");
        assert!(
            matches!(error, OutputError::AlreadyExists { .. }),
            "{error}"
        );
        assert_eq!(std::fs::read(&final_path).expect("readable"), b"old");
    }

    #[test]
    fn force_replaces_transactionally() {
        let directory = tempfile::tempdir().expect("temp dir");
        let final_path = directory.path().join("chr1.bam");
        std::fs::write(&final_path, b"old").expect("seeded");

        let (mut file, pending) = PendingFile::create(&final_path, true).expect("created");
        file.write_all(b"new").expect("written");
        drop(file);
        // The old file is still intact right up to the rename.
        assert_eq!(std::fs::read(&final_path).expect("readable"), b"old");
        pending.commit().expect("committed");
        assert_eq!(std::fs::read(&final_path).expect("readable"), b"new");
    }

    #[test]
    fn a_transaction_commits_every_file() {
        let directory = tempfile::tempdir().expect("temp dir");
        let bam = directory.path().join("chr1.bam");
        let index = directory.path().join("chr1.bam.bai");

        let mut transaction = Transaction::new();
        for path in [&index, &bam] {
            let (mut file, pending) = PendingFile::create(path, false).expect("created");
            file.write_all(b"x").expect("written");
            transaction.push(pending);
        }
        assert_eq!(transaction.len(), 2);

        let committed = transaction.commit().expect("committed");
        assert_eq!(committed, vec![index.clone(), bam.clone()]);
        assert!(bam.exists() && index.exists());
    }

    #[test]
    fn a_failed_transaction_leaves_nothing_behind() {
        let directory = tempfile::tempdir().expect("temp dir");
        let good = directory.path().join("good.bam");
        let blocked = directory.path().join("blocked.bam");
        std::fs::write(&blocked, b"existing").expect("seeded");

        let mut transaction = Transaction::new();
        let (mut file, pending) = PendingFile::create(&good, false).expect("created");
        file.write_all(b"x").expect("written");
        transaction.push(pending);

        // Simulate a competing writer creating the second output after the
        // pre-flight check but before the commit.
        let (mut file, pending) = PendingFile::create(&blocked, true).expect("created");
        file.write_all(b"y").expect("written");
        let mut pending = pending;
        pending.force = false;
        transaction.push(pending);

        let error = transaction.commit().expect_err("must fail");
        assert!(
            matches!(error, OutputError::AlreadyExists { .. }),
            "{error}"
        );
        assert!(!good.exists(), "the first commit must be rolled back");
        assert_eq!(std::fs::read(&blocked).expect("readable"), b"existing");
    }

    #[test]
    fn aborting_a_transaction_removes_every_temporary_file() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut transaction = Transaction::new();
        let mut temporaries = Vec::new();
        for name in ["a.bam", "b.bam"] {
            let (_, pending) =
                PendingFile::create(&directory.path().join(name), false).expect("created");
            temporaries.push(pending.temporary_path().to_path_buf());
            transaction.push(pending);
        }
        transaction.abort();
        for temporary in temporaries {
            assert!(!temporary.exists(), "{temporary:?}");
        }
    }

    #[test]
    fn a_parked_file_can_be_reopened_for_appending() {
        let directory = tempfile::tempdir().expect("temp dir");
        let final_path = directory.path().join("chr1.bam");
        let (mut file, pending) = PendingFile::create(&final_path, false).expect("created");
        file.write_all(b"first").expect("written");
        drop(file);

        let mut reopened = pending.reopen_append().expect("reopened");
        reopened.write_all(b"second").expect("written");
        drop(reopened);

        assert_eq!(pending.size().expect("sized"), 11);
        pending.commit().expect("committed");
        assert_eq!(
            std::fs::read(&final_path).expect("readable"),
            b"firstsecond"
        );
    }

    #[test]
    fn the_registry_sweeps_orphans() {
        let directory = tempfile::tempdir().expect("temp dir");
        let final_path = directory.path().join("orphan.bam");
        let (file, pending) = PendingFile::create(&final_path, false).expect("created");
        drop(file);
        let temporary = pending.temporary_path().to_path_buf();
        // Forget the handle without running `Drop`, as a signal would.
        std::mem::forget(pending);

        assert!(temporary.exists());
        assert!(registry::snapshot().contains(&temporary));
        assert_eq!(registry::remove_all([temporary.clone()]), 1);
        assert!(!temporary.exists());
        assert!(!registry::snapshot().contains(&temporary));
    }

    #[test]
    fn the_interrupt_flag_is_observable() {
        let flag = InterruptFlag::new();
        assert!(!flag.is_raised());
        flag.raise();
        assert!(flag.is_raised());
        flag.clear();
        assert!(!flag.is_raised());
    }

    #[test]
    fn directories_are_created_and_validated() {
        let directory = tempfile::tempdir().expect("temp dir");
        let nested = directory.path().join("a/b/c");
        prepare_directory(&nested).expect("created");
        assert!(nested.is_dir());
        // Idempotent.
        prepare_directory(&nested).expect("still fine");

        let file = directory.path().join("file");
        std::fs::write(&file, b"x").expect("seeded");
        let error = prepare_directory(&file).expect_err("must reject");
        assert!(
            matches!(error, OutputError::BadOutputDirectory { .. }),
            "{error}"
        );
    }
}
