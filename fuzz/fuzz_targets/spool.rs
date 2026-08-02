// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Spool replay must reject corruption rather than dropping records.

#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Replay reads from a path, so the input has to reach the filesystem. A
    // single reused file keeps the fuzzer fast.
    let Ok(directory) = std::env::var("BAMSPLIT_FUZZ_DIR") else {
        return;
    };
    let path = std::path::Path::new(&directory).join("fuzz.bspl");
    if std::fs::write(&path, data).is_err() {
        return;
    }

    let mut records = 0u64;
    let outcome = bamsplit_core::engine::spool::replay(&path, |body| {
        records += 1;
        // A replayed body must be usable as a record or rejected as one; either
        // way, no panic.
        let _ = bamsplit_core::bam::RawRecord::new(body);
        Ok(())
    });

    if let Ok(count) = outcome {
        // A successful replay must agree with what the callback saw.
        assert_eq!(count, records);
    }
});
