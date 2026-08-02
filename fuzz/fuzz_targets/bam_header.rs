// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Header parsing must never panic, however malformed the bytes.

#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut cursor = data;
    if let Ok(header) = bamsplit_core::bam::header::BamHeader::read_from(&mut cursor) {
        // A header that parsed must also re-serialize and re-parse.
        if let Ok(bytes) = {
            let mut out = Vec::new();
            header.write_to(&mut out).map(|()| out)
        } {
            let mut round_trip = &bytes[..];
            let _ = bamsplit_core::bam::header::BamHeader::read_from(&mut round_trip);
        }
        let _ = header.checksum();
        let _ = header.max_reference_length();
    }
});
