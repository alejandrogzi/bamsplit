// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Whole-run benchmarks: engine against engine, thread count against thread
//! count.
//!
//! ```console
//! cargo bench --bench end_to_end
//! ```
//!
//! These write real files into a temporary directory, so they measure the thing
//! a user experiences — including index construction and the atomic commit — and
//! not just the in-memory pipeline. `scripts/benchmark.sh` compares the same
//! workload against `samtools view` per chromosome.

use bamsplit_core::engine::EngineKind;
use bamsplit_core::index::IndexMode;
use bamsplit_core::{ChromOptions, RunOptions, ShardOptions};
use bamsplit_fixtures::{Fixture, RecordSpec};
use criterion::{Criterion, criterion_group, criterion_main};

/// A coordinate-sorted BAM with `references` references and `per_reference`
/// records each.
fn fixture(references: usize, per_reference: i32) -> Fixture {
    let names: Vec<String> = (0..references).map(|index| format!("chr{index}")).collect();
    let entries: Vec<(&str, u32)> = names
        .iter()
        .map(|name| (name.as_str(), 250_000_000u32))
        .collect();

    // Seeded payloads, not constant ones: an all-zero sequence with "quality
    // unavailable" deflates to nothing, and a whole-run benchmark whose output
    // costs nothing to compress measures the wrong thing entirely.
    let mut records = Vec::with_capacity(references * per_reference as usize);
    for reference in 0..references {
        for index in 0..per_reference {
            records.push(RecordSpec {
                sequence_length: 150,
                ..RecordSpec::mapped(
                    &format!("read{reference}_{index:07}"),
                    i32::try_from(reference).unwrap_or(0),
                    index * 100,
                )
                .with_cigar(vec![(150, 0)])
                .with_string_tag(b"RG", "rg0")
                .with_int_tag(b"NM", 1)
                .with_payload_seed(
                    0x9e37_79b9_7f4a_7c15 ^ ((reference as u64) << 32) ^ index as u64,
                )
            });
        }
    }
    Fixture::new(&entries, "coordinate").with_records(records)
}

struct Input {
    _directory: tempfile::TempDir,
    path: std::path::PathBuf,
    records: u64,
}

fn write_input(fixture: &Fixture) -> Input {
    let directory = tempfile::tempdir().expect("temp dir");
    let path = directory.path().join("input.bam");
    fixture.write(&path).expect("fixture written");
    Input {
        _directory: directory,
        path,
        records: fixture.records.len() as u64,
    }
}

fn run_chrom(input: &std::path::Path, run: &RunOptions) {
    let directory = tempfile::tempdir().expect("temp dir");
    bamsplit_core::split_by_chromosome(
        input,
        &ChromOptions::new(directory.path().join("out")),
        run,
    )
    .expect("split succeeds");
}

fn threads(criterion: &mut Criterion) {
    let input = write_input(&fixture(8, 4_000));

    let mut group = criterion.benchmark_group("chrom_threads");
    group.sample_size(10);
    group.throughput(criterion::Throughput::Elements(input.records));
    for count in [1usize, 2, 4, 8, 16] {
        group.bench_function(format!("{count}"), |bencher| {
            bencher.iter(|| {
                run_chrom(
                    &input.path,
                    &RunOptions {
                        threads: count,
                        engine: EngineKind::Stream,
                        compression_level: 6,
                        manifest: bamsplit_core::ManifestFormat::None,
                        ..RunOptions::default()
                    },
                );
            });
        });
    }
    group.finish();
}

fn engines(criterion: &mut Criterion) {
    let input = write_input(&fixture(8, 4_000));

    let mut group = criterion.benchmark_group("chrom_engines");
    group.sample_size(10);
    group.throughput(criterion::Throughput::Elements(input.records));
    for (name, engine) in [("stream", EngineKind::Stream), ("spool", EngineKind::Spool)] {
        group.bench_function(name, |bencher| {
            bencher.iter(|| {
                run_chrom(
                    &input.path,
                    &RunOptions {
                        threads: 4,
                        engine,
                        compression_level: 1,
                        manifest: bamsplit_core::ManifestFormat::None,
                        ..RunOptions::default()
                    },
                );
            });
        });
    }
    group.finish();
}

fn indexing(criterion: &mut Criterion) {
    let input = write_input(&fixture(4, 5_000));

    let mut group = criterion.benchmark_group("chrom_indexing");
    group.sample_size(10);
    group.throughput(criterion::Throughput::Elements(input.records));
    for (name, index) in [("bai", IndexMode::Bai), ("none", IndexMode::None)] {
        group.bench_function(name, |bencher| {
            bencher.iter(|| {
                run_chrom(
                    &input.path,
                    &RunOptions {
                        threads: 4,
                        engine: EngineKind::Stream,
                        compression_level: 1,
                        index,
                        manifest: bamsplit_core::ManifestFormat::None,
                        ..RunOptions::default()
                    },
                );
            });
        });
    }
    group.finish();
}

fn contig_count(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("chrom_contigs");
    group.sample_size(10);
    for references in [4usize, 64, 512] {
        let input = write_input(&fixture(references, 40_000 / references as i32));
        group.throughput(criterion::Throughput::Elements(input.records));
        group.bench_function(format!("{references}"), |bencher| {
            bencher.iter(|| {
                run_chrom(
                    &input.path,
                    &RunOptions {
                        threads: 4,
                        compression_level: 1,
                        manifest: bamsplit_core::ManifestFormat::None,
                        ..RunOptions::default()
                    },
                );
            });
        });
    }
    group.finish();
}

fn sharding(criterion: &mut Criterion) {
    let input = write_input(&fixture(4, 5_000));

    let mut group = criterion.benchmark_group("shard");
    group.sample_size(10);
    group.throughput(criterion::Throughput::Elements(input.records));
    for shards in [8u32, 64] {
        group.bench_function(format!("{shards}"), |bencher| {
            bencher.iter(|| {
                let directory = tempfile::tempdir().expect("temp dir");
                bamsplit_core::split_by_shard(
                    &input.path,
                    &ShardOptions::new(directory.path().join("out"), shards),
                    &RunOptions {
                        threads: 4,
                        compression_level: 1,
                        manifest: bamsplit_core::ManifestFormat::None,
                        ..RunOptions::default()
                    },
                )
                .expect("split succeeds");
            });
        });
    }
    group.finish();
}

criterion_group!(benches, threads, engines, indexing, contig_count, sharding);
criterion_main!(benches);
