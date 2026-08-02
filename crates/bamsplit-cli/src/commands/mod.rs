// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! One module per subcommand, each translating parsed arguments into a
//! [`bamsplit_core`] call and rendering the result.
//!
//! Rendering is deliberately separate from execution: the core returns a
//! manifest, and these modules decide how to show it. That is why `--quiet` and
//! `--manifest none` cannot change what was written.

pub mod chrom;
pub mod inspect;
pub mod region;
pub mod shard;
pub mod tag;

use bamsplit_core::SplitReport;

/// Prints the one-line summary a user sees after a successful split.
pub fn report_summary(report: &SplitReport) {
    let manifest = &report.manifest;
    tracing::info!(
        "wrote {} outputs from {} records in {:.2}s using the {} engine",
        report.created_outputs(),
        manifest.input_records,
        manifest.elapsed_seconds,
        manifest.selected_engine
    );
    for note in &manifest.notes {
        tracing::debug!("{note}");
    }
    for path in &report.manifest_paths {
        tracing::info!("manifest: {}", path.display());
    }
}
