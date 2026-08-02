// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Filename encoding must be containment-safe for every input, and reversible
//! for every input short enough to escape the length cap.

#![no_main]

use bamsplit_core::output::filename::{
    FilenameTemplate, MAX_STEM_LEN, decode, encode, is_safe_component,
};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let encoded = encode(data);

    // Containment holds unconditionally.
    assert!(is_safe_component(&encoded), "unsafe stem: {encoded:?}");
    assert!(encoded.len() <= MAX_STEM_LEN);
    let joined = std::path::Path::new("/output").join(&encoded);
    assert!(joined.starts_with("/output"));
    assert_eq!(joined.components().count(), 3);

    // Reversibility holds below the length cap.
    if encoded.len() < MAX_STEM_LEN {
        assert_eq!(decode(&encoded).as_deref(), Some(data));
    }

    // Templates must reject anything they cannot render safely.
    if let Ok(text) = std::str::from_utf8(data)
        && let Ok(template) = FilenameTemplate::parse(text)
        && let Ok(rendered) = template.render(&encoded, 0)
    {
        assert!(is_safe_component(&rendered), "unsafe render: {rendered:?}");
    }
});
