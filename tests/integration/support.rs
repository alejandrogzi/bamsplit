// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Shared helpers for the integration suite.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use bamsplit_core::bam::RawRecordReader;
use bamsplit_core::bam::header::BamHeader;
use bamsplit_core::{Manifest, RunOptions};
use bamsplit_fixtures::Fixture;

/// A temporary working directory holding one fixture.
pub struct Workspace {
    pub directory: tempfile::TempDir,
    pub input: PathBuf,
}

impl Workspace {
    /// Writes `fixture` into a fresh temporary directory.
    pub fn with(fixture: &Fixture) -> Self {
        let directory = tempfile::tempdir().expect("temp dir");
        let input = directory.path().join("input.bam");
        fixture.write(&input).expect("fixture written");
        Self { directory, input }
    }

    /// A path inside the workspace.
    pub fn path(&self, name: &str) -> PathBuf {
        self.directory.path().join(name)
    }
}

/// Options that keep tests fast: level-1 compression and a single thread unless
/// a test says otherwise.
pub fn run_options() -> RunOptions {
    RunOptions {
        compression_level: 1,
        command_line: "bamsplit test".to_string(),
        ..RunOptions::default()
    }
}

/// Reads every record body from a BAM, plus its header.
pub fn read_bam(path: &Path) -> (BamHeader, Vec<Vec<u8>>) {
    let file = std::fs::File::open(path).unwrap_or_else(|error| panic!("{path:?}: {error}"));
    let mut reader = noodles_bgzf::io::Reader::new(file);
    let header = BamHeader::read_from(&mut reader).expect("valid header");
    let mut records = RawRecordReader::new(reader);
    let mut bodies = Vec::new();
    while let Some(record) = records.read_record().expect("valid record") {
        bodies.push(record.raw_bytes().to_vec());
    }
    (header, bodies)
}

/// Maps each output stem to its record bodies.
pub fn collect_outputs(directory: &Path) -> BTreeMap<String, Vec<Vec<u8>>> {
    let mut outputs = BTreeMap::new();
    let entries = std::fs::read_dir(directory).expect("listable");
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|extension| extension.to_str()) != Some("bam") {
            continue;
        }
        let stem = path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or_default()
            .to_string();
        let (_, bodies) = read_bam(&path);
        outputs.insert(stem, bodies);
    }
    outputs
}

/// Loads the JSON manifest a run wrote.
pub fn read_manifest(directory: &Path) -> Manifest {
    let path = directory.join("bamsplit.manifest.json");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{path:?}: {error}"));
    serde_json::from_str(&text).expect("valid manifest JSON")
}

/// Every record body in the input fixture, in input order.
pub fn fixture_bodies(fixture: &Fixture) -> Vec<Vec<u8>> {
    fixture
        .records
        .iter()
        .map(bamsplit_fixtures::RecordSpec::encode)
        .collect()
}

/// Asserts that the union of every output's records equals the input's, as a
/// multiset. This is the losslessness property, stated once.
pub fn assert_lossless(fixture: &Fixture, directory: &Path) {
    let mut expected = fixture_bodies(fixture);
    let mut actual: Vec<Vec<u8>> = collect_outputs(directory).into_values().flatten().collect();
    expected.sort();
    actual.sort();
    assert_eq!(
        actual.len(),
        expected.len(),
        "record count changed: {} in, {} out",
        expected.len(),
        actual.len()
    );
    assert_eq!(actual, expected, "record bodies were altered");
}

/// The path to the built `bamsplit` binary.
pub fn binary() -> PathBuf {
    // `CARGO_BIN_EXE_` is only set for binaries of the same package, so the
    // path is derived from the test executable's own location instead.
    let mut path = std::env::current_exe().expect("test executable path");
    path.pop();
    if path.ends_with("deps") {
        path.pop();
    }
    path.join(if cfg!(windows) {
        "bamsplit.exe"
    } else {
        "bamsplit"
    })
}
