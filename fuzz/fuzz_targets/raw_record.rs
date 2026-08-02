// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Record framing and field access must never panic.

#![no_main]

use bamsplit_core::bam::{RawRecord, RawRecordReader};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // The framed path: `block_size` prefixes come from the fuzzer too.
    let mut reader = RawRecordReader::new(data).with_max_record_size(1 << 20);
    let mut budget = 0;
    while budget < 10_000 {
        budget += 1;
        match reader.read_record() {
            Ok(Some(record)) => {
                let _ = record.reference_sequence_id();
                let _ = record.alignment_start();
                let _ = record.alignment_end();
                let _ = record.qname();
                let _ = record.flags();
                let _ = record.mapping_quality();
                let _ = record.tag(*b"RG");
                let _ = record.validate_deeply(true);
            }
            Ok(None) | Err(_) => break,
        }
    }

    // The unframed path: the whole input treated as one record body.
    if let Ok(record) = RawRecord::new(data) {
        let _ = record.reference_span();
        let _ = record.data();
    }
});
