// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Auxiliary-tag decoding must never panic or over-allocate.

#![no_main]

use bamsplit_core::bam::tags::{TagReader, find_tag, validate_section};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut reader = TagReader::new(data);
    while let Some(field) = reader.next() {
        match field {
            Ok((_, value)) => {
                let _ = value.to_routing_bytes();
                if let Some(array) = value.as_array() {
                    // Iteration must respect the validated payload length, not
                    // the declared element count.
                    let observed = array.iter_f64().count();
                    assert!(observed <= array.len() as usize);
                }
            }
            Err(_) => break,
        }
    }
    let _ = find_tag(data, *b"RG");
    let _ = validate_section(data, true);
});
