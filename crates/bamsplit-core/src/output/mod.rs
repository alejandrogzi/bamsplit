// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Naming, writing, and atomically committing output BAMs.
//!
//! ```text
//! logical key ──filename── stem ──writer── .part.<pid>.<nonce> ──commit── final path
//!                  │                 │                              │
//!            reversible,        BGZF blocks +                  rename(2),
//!         traversal-safe      streaming index                after validation
//! ```
//!
//! Routers produce logical keys; nothing below this module knows or cares that
//! a key came from a reference name, a tag value, or a transcript identifier.

pub mod bgzf;
pub mod filename;
pub mod manager;
pub mod transaction;
pub mod writer;

pub use bgzf::{BGZF_EOF, BgzfBlockWriter, MAX_BLOCK_UNCOMPRESSED, VirtualRange};
pub use filename::{FilenameEncoder, FilenameTemplate, OutputPaths};
pub use manager::{OutputManager, OutputManagerOptions};
pub use transaction::{
    InterruptFlag, PendingFile, Transaction, interrupt_flag, prepare_directory, registry,
};
pub use writer::{FinishedOutput, OutputSettings, RecordWriter};
