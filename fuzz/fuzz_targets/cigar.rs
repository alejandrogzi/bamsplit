// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! CIGAR decoding and block derivation must never panic or overflow.

#![no_main]

use bamsplit_core::bam::cigar::{CigarOps, DeletionPolicy};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let ops = CigarOps::new(data);
    let _ = ops.reference_span();
    let _ = ops.read_length();
    let _ = ops.to_sam_string();
    for start in [0i64, 1, i64::MAX - 1, i64::MIN + 1] {
        let _ = ops.reference_interval(start);
        for policy in [
            DeletionPolicy::ExcludeFromBlocks,
            DeletionPolicy::IncludeInBlocks,
        ] {
            if let Ok(blocks) = ops.aligned_blocks(start, policy) {
                // Blocks must be sorted, non-empty, and non-overlapping.
                for window in blocks.windows(2) {
                    assert!(window[0].end <= window[1].start);
                }
                assert!(blocks.iter().all(|block| !block.is_empty()));
            }
        }
    }
});
