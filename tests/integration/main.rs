// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! End-to-end tests: real BAMs in, real BAMs out, checked for losslessness.
//!
//! These exercise the whole stack — reader, router, engine, writer, index,
//! manifest — against the synthetic fixtures in the `bamsplit` crate. Nothing
//! here needs `samtools`; the differential suite handles cross-checking against
//! an external oracle.

mod support;

mod cli;
mod conservation;
mod engines;
mod properties;
mod regions;
mod robustness;
