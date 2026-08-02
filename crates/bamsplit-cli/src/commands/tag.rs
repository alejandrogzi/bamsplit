// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `bamsplit tag`.

use bamsplit_core::{Result, TagOptions};

use crate::cli::{GlobalOptions, TagArgs};

/// Runs a per-tag-value split.
///
/// # Errors
///
/// Propagates every [`bamsplit_core::Error`], including the cardinality guard.
pub fn run(globals: &GlobalOptions, args: TagArgs) -> Result<()> {
    let run = globals.to_run_options(crate::cli::command_line(), args.index, None);
    let options = TagOptions {
        out_dir: args.out_dir,
        tag: args.tag,
        field: args.field.map(Into::into),
        missing: args.missing.into(),
        missing_name: args.missing_name,
        unknown_read_group: args.unknown_read_group.into(),
        max_outputs: args.max_outputs,
        allow_high_cardinality: args.allow_high_cardinality,
        allow_array_tags: args.allow_array_tags,
    };
    let report = bamsplit_core::split_by_tag(&args.input, &options, &run)?;
    super::report_summary(&report);
    Ok(())
}
