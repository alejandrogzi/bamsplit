// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The `bamsplit` binary.

fn main() -> std::process::ExitCode {
    let code = bamsplit_cli::run();
    // `ExitCode::from` takes a `u8`; every `bamsplit` code is 0..=6.
    std::process::ExitCode::from(u8::try_from(code.code()).unwrap_or(1))
}
