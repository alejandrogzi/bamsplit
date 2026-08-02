// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The `bamsplit` command-line front-end.
//!
//! Exposed as a library so the integration tests can drive argument parsing and
//! dispatch directly, without spawning a process for every case.

#![warn(missing_docs)]

pub mod cli;
pub mod commands;

use bamsplit_core::error::{Classify, ExitCode, render_chain};

/// Parses arguments, runs the requested command, and returns an exit code.
///
/// Never panics on user input and never returns [`Result`]: every failure is
/// already classified into an [`ExitCode`] by the time it gets here.
#[must_use]
pub fn run() -> ExitCode {
    use clap::Parser as _;

    let cli = match cli::Cli::try_parse() {
        Ok(cli) => cli,
        Err(error) => {
            // `clap` writes its own message, including `--help` and `--version`,
            // which exit successfully rather than as an argument error.
            let _ = error.print();
            return if error.use_stderr() {
                ExitCode::InvalidArguments
            } else {
                ExitCode::Success
            };
        }
    };

    init_logging(&cli.globals);
    install_interrupt_handler();

    let globals = cli.globals.clone();
    let outcome = match cli.command {
        cli::Command::Chrom(args) => commands::chrom::run(&globals, args),
        cli::Command::Shard(args) => commands::shard::run(&globals, args),
        cli::Command::Tag(args) => commands::tag::run(&globals, args),
        cli::Command::Region(args) => commands::region::run(&globals, args),
        cli::Command::Inspect(args) => commands::inspect::run(&globals, args),
    };

    match outcome {
        Ok(()) => ExitCode::Success,
        Err(error) => {
            let code = error.exit_code();
            tracing::error!("{}", render_chain(&error));
            // A signal leaves temporary files whose owners never ran `Drop`,
            // because unwinding cannot be triggered from a handler. Sweeping
            // here is the last line of defence.
            if code == ExitCode::Interrupted {
                let removed = bamsplit_core::output::registry::sweep();
                if removed > 0 {
                    tracing::warn!("removed {removed} partial output(s) after interruption");
                }
            }
            code
        }
    }
}

fn init_logging(globals: &cli::GlobalOptions) {
    use tracing_subscriber::EnvFilter;

    if globals.quiet {
        return;
    }
    let filter = std::env::var("BAMSPLIT_LOG").ok().map_or_else(
        || EnvFilter::new(globals.log_level.as_filter()),
        EnvFilter::new,
    );
    // Logs go to stderr so stdout stays free for machine-readable output.
    let _ = tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_target(false)
        .without_time()
        .try_init();
}

fn install_interrupt_handler() {
    // The handler does the minimum a signal handler safely can: set a flag. The
    // record loop polls it and returns an error, which unwinds and triggers the
    // ordinary `Drop`-based cleanup.
    let result = ctrlc::set_handler(|| {
        bamsplit_core::output::interrupt_flag().raise();
    });
    if let Err(error) = result {
        tracing::warn!("cannot install a signal handler, so Ctrl-C will not clean up: {error}");
    }
}
