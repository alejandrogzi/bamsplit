// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `bamsplit chrom`.

use bamsplit_core::{ChromOptions, Result};

use crate::cli::{ChromArgs, GlobalOptions};

/// Runs a per-reference split.
///
/// # Errors
///
/// Propagates every [`bamsplit_core::Error`]; `main` maps it to an exit code.
pub fn run(globals: &GlobalOptions, args: ChromArgs) -> Result<()> {
    let run = globals.to_run_options(
        crate::cli::command_line(),
        args.index,
        args.filename_template.clone(),
    );
    let options = ChromOptions {
        out_dir: args.out_dir,
        placed_unmapped: args.placed_unmapped.into(),
        unplaced: args.unplaced.into(),
        unmapped_name: args.unmapped_name,
        emit_empty: args.emit_empty,
        include: args.include,
        exclude: args.exclude,
        reference_list: args.reference_list,
        ignore_missing_references: args.ignore_missing_references,
    };
    let report = bamsplit_core::split_by_chromosome(&args.input, &options, &run)?;
    super::report_summary(&report);
    Ok(())
}
