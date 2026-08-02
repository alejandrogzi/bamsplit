// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Lossless, high-throughput partitioning of BAM files.
//!
//! `bamsplit-core` is the library behind the `bamsplit` command-line tool. It
//! is deliberately free of argument parsing and terminal presentation, so it is
//! equally usable from a Rust program.
//!
//! # The shape of a run
//!
//! ```text
//! input BAM ──▶ header ──▶ router ──▶ engine ──▶ output manager ──▶ manifest
//!                            │           │            │
//!                    Route::One/Many   stream      transactional
//!                    Route::Drop       indexed     .part → rename
//!                                      spool       streaming BAI/CSI
//! ```
//!
//! * [`bam`] turns bytes into borrowed, bounds-checked views.
//! * [`routing`] decides *where* a record goes, and never touches a file.
//! * [`engine`] decides *how* the pass is executed.
//! * [`output`] owns filenames, temporary paths, and atomic commits.
//! * [`manifest`] records what happened and proves no record was lost.
//!
//! # Guarantees
//!
//! * **Lossless.** Record bodies are written back byte-for-byte.
//! * **Bounded.** Memory and open file descriptors stay within configured
//!   limits regardless of output cardinality.
//! * **Deterministic.** The same input and options produce the same outputs,
//!   the same filenames, and the same manifest.
//! * **Transactional.** An output appears at its final path only after it is
//!   complete and validated.
//! * **Non-panicking.** Malformed input, bad configuration, integer overflow,
//!   and filesystem failures all become typed errors. The only `panic!`s are on
//!   internal invariants proven by construction, and each is documented at its
//!   site.

// `deny` rather than `forbid`: memory-mapped input needs exactly one `unsafe`
// call, which `engine::unsafe_map` isolates behind an `#[allow]` and a `SAFETY`
// comment. Everything else in the crate is safe Rust.
#![deny(unsafe_code)]
#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(
    // The crate reports counts and coordinates as `u64`/`i64` and converts at
    // well-defined boundaries; blanket-denying these lints would bury the few
    // conversions that genuinely need attention.
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    // Long, explicit error enums and option structs are the point.
    clippy::module_name_repetitions,
    clippy::struct_excessive_bools,
    clippy::too_many_lines,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc
)]

#[cfg(feature = "annotation")]
pub mod annotation;
pub mod bam;
pub mod engine;
pub mod error;
pub mod index;
pub mod interval;
pub mod manifest;
pub mod output;
pub mod planning;
pub mod routing;
pub mod stats;

mod api;

pub use api::Stats;
#[cfg(feature = "annotation")]
pub use api::{AnnotationReport, RegionOptions, split_by_region};
pub use api::{
    ChromOptions, InspectOptions, InspectReport, RunOptions, ScanReport, ShardOptions, SplitReport,
    TagOptions, inspect, split_by_chromosome, split_by_shard, split_by_tag,
};
pub use bam::{BamHeader, RawRecord, RawRecordReader};
pub use error::{Classify, Error, ExitCode, Result};
pub use interval::Interval;
pub use manifest::{Manifest, ManifestFormat};

/// The crate version, as declared in `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The program name written into `@PG PN` and the manifest.
pub const PROGRAM_NAME: &str = "bamsplit";
