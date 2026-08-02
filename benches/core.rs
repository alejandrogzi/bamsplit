// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Micro-benchmarks for the hot path.
//!
//! The question these answer is where a `bamsplit chrom` run actually spends its
//! time, so the answer decides what is worth optimizing. Run with:
//!
//! ```console
//! cargo bench --bench core
//! ```
//!
//! The `noodles` decode/encode benchmark is the comparison that justifies the
//! raw-transfer path: if it were not materially slower, the raw path would be
//! complexity for nothing.

use bamsplit_core::bam::cigar::{CigarOps, DeletionPolicy};
use bamsplit_core::bam::tags::find_tag;
use bamsplit_core::bam::{RawRecord, RawRecordReader};
use bamsplit_core::output::bgzf::BgzfBlockWriter;
use bamsplit_core::output::filename::encode;
use bamsplit_core::routing::shard::{ShardRouter, ShardRouterOptions};
use bamsplit_fixtures::{Fixture, RecordSpec};
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

/// A record with a realistic auxiliary section and a spliced CIGAR.
fn sample_record() -> Vec<u8> {
    RecordSpec {
        sequence_length: 150,
        ..RecordSpec::mapped("HWI-ST745:123:C0FVAACXX:5:1101:1234:5678", 0, 1_000_000)
            .with_cigar(vec![(75, 0), (5_000, 3), (75, 0)])
            .with_string_tag(b"RG", "rg0")
            .with_int_tag(b"NM", 2)
            .with_string_tag(b"MD", "75A74")
            .with_string_tag(b"CB", "ACGTACGTACGTACGT-1")
    }
    .encode()
}

fn parsing(criterion: &mut Criterion) {
    let body = sample_record();
    let record = RawRecord::new(&body).expect("valid");

    let mut group = criterion.benchmark_group("record");
    group.bench_function("validate_and_wrap", |bencher| {
        bencher.iter(|| RawRecord::new(black_box(&body)).expect("valid"));
    });
    group.bench_function("fixed_core", |bencher| {
        bencher.iter(|| {
            (
                black_box(&record).reference_sequence_id().expect("valid"),
                record.alignment_start().expect("valid"),
                record.flags().expect("valid"),
                record.mapping_quality().expect("valid"),
            )
        });
    });
    group.bench_function("qname", |bencher| {
        bencher.iter(|| black_box(&record).qname().expect("valid"));
    });
    group.bench_function("tag_lookup_first", |bencher| {
        bencher.iter(|| black_box(&record).tag(*b"RG").expect("valid"));
    });
    group.bench_function("tag_lookup_last", |bencher| {
        bencher.iter(|| black_box(&record).tag(*b"CB").expect("valid"));
    });
    group.bench_function("tag_lookup_absent", |bencher| {
        bencher.iter(|| black_box(&record).tag(*b"ZZ").expect("valid"));
    });
    group.bench_function("alignment_end", |bencher| {
        bencher.iter(|| black_box(&record).alignment_end().expect("valid"));
    });
    group.finish();
}

fn cigar(criterion: &mut Criterion) {
    let mut packed = Vec::new();
    for (length, kind) in [(75u32, 0u32), (5_000, 3), (10, 2), (65, 0), (10, 4)] {
        packed.extend_from_slice(&((length << 4) | kind).to_le_bytes());
    }
    let ops = CigarOps::new(&packed);
    let mut blocks = Vec::with_capacity(8);

    let mut group = criterion.benchmark_group("cigar");
    group.bench_function("reference_span", |bencher| {
        bencher.iter(|| black_box(ops).reference_span().expect("valid"));
    });
    group.bench_function("aligned_blocks_collect", |bencher| {
        bencher.iter(|| {
            black_box(ops)
                .aligned_blocks(1_000, DeletionPolicy::ExcludeFromBlocks)
                .expect("valid")
        });
    });
    group.bench_function("aligned_blocks_reuse_buffer", |bencher| {
        bencher.iter(|| {
            blocks.clear();
            black_box(ops)
                .for_each_aligned_block(1_000, DeletionPolicy::ExcludeFromBlocks, |block| {
                    blocks.push(block);
                })
                .expect("valid");
        });
    });
    group.finish();
}

fn tags(criterion: &mut Criterion) {
    let body = sample_record();
    let record = RawRecord::new(&body).expect("valid");
    let data = record.data().expect("valid");

    let mut group = criterion.benchmark_group("tags");
    group.bench_function("find_tag", |bencher| {
        bencher.iter(|| find_tag(black_box(data), *b"MD").expect("valid"));
    });
    group.bench_function("normalize_scalar", |bencher| {
        let value = find_tag(data, *b"NM").expect("valid").expect("present");
        bencher.iter(|| black_box(value).to_routing_bytes());
    });
    group.bench_function("validate_section", |bencher| {
        bencher.iter(|| {
            bamsplit_core::bam::tags::validate_section(black_box(data), false).expect("valid")
        });
    });
    group.finish();
}

fn routing(criterion: &mut Criterion) {
    let router = ShardRouter::new(ShardRouterOptions {
        shards: 64,
        ..ShardRouterOptions::default()
    })
    .expect("built");
    let name = b"HWI-ST745:123:C0FVAACXX:5:1101:1234:5678";

    let mut group = criterion.benchmark_group("routing");
    group.bench_function("xxh3_shard", |bencher| {
        bencher.iter(|| black_box(&router).shard_of(black_box(name)));
    });
    group.bench_function("filename_encode_plain", |bencher| {
        bencher.iter(|| encode(black_box(b"chr1")));
    });
    group.bench_function("filename_encode_escaped", |bencher| {
        bencher.iter(|| encode(black_box(b"chr1/alternate:1-2")));
    });
    group.finish();
}

/// The whole read-route-write loop, in memory.
///
/// This is the number that matters: it includes BGZF inflation, record framing,
/// routing, and BGZF deflation, so it is directly comparable across changes.
fn round_trip(criterion: &mut Criterion) {
    // The payload is seeded rather than constant. An all-zero sequence with
    // "quality unavailable" deflates to almost nothing, which would make the
    // level comparison below measure a run of zeroes instead of real data.
    let fixture = {
        let base = Fixture::new(&[("chr1", 250_000_000)], "coordinate");
        base.with_records((0..20_000i32).map(|index| {
            RecordSpec {
                sequence_length: 150,
                ..RecordSpec::mapped(&format!("read{index:08}"), 0, index * 50)
                    .with_cigar(vec![(150, 0)])
                    .with_string_tag(b"RG", "rg0")
                    .with_int_tag(b"NM", 1)
                    .with_payload_seed(0x9e37_79b9_7f4a_7c15 ^ index as u64)
            }
        }))
    };
    let bytes = fixture.to_bytes().expect("encodes");

    let mut group = criterion.benchmark_group("round_trip");
    group.throughput(criterion::Throughput::Elements(fixture.records.len() as u64));

    group.bench_function("read_only", |bencher| {
        bencher.iter(|| {
            let mut reader = noodles_bgzf::io::Reader::new(&bytes[..]);
            bamsplit_core::bam::header::BamHeader::read_from(&mut reader).expect("valid");
            let mut records = RawRecordReader::new(reader);
            let mut total = 0u64;
            while let Some(record) = records.read_record().expect("valid") {
                total += u64::from(record.reference_sequence_id().expect("valid").is_some());
            }
            black_box(total)
        });
    });

    for level in [0u32, 1, 6] {
        group.bench_function(format!("read_and_write_level_{level}"), |bencher| {
            bencher.iter(|| {
                let mut reader = noodles_bgzf::io::Reader::new(&bytes[..]);
                bamsplit_core::bam::header::BamHeader::read_from(&mut reader).expect("valid");
                let mut records = RawRecordReader::new(reader);
                let mut writer: BgzfBlockWriter<Vec<u8>, ()> =
                    BgzfBlockWriter::new(Vec::new(), level, None, 1, 0);
                while let Some(record) = records.read_record().expect("valid") {
                    writer
                        .write_record(record.raw_bytes(), (), &mut |(), _| {})
                        .expect("written");
                }
                black_box(writer.finish(&mut |(), _| {}).expect("finished").len())
            });
        });
    }
    group.finish();
}

/// Raw transfer against a full decode/re-encode through `noodles`.
///
/// The raw path exists only if this gap is real; the benchmark is the evidence.
fn raw_versus_decoded(criterion: &mut Criterion) {
    let bodies: Vec<Vec<u8>> = (0..5_000i32)
        .map(|index| {
            RecordSpec {
                sequence_length: 150,
                ..RecordSpec::mapped(&format!("read{index:08}"), 0, index * 50)
                    .with_cigar(vec![(150, 0)])
                    .with_string_tag(b"RG", "rg0")
                    .with_int_tag(b"NM", 1)
            }
            .encode()
        })
        .collect();

    let mut group = criterion.benchmark_group("transfer");
    group.throughput(criterion::Throughput::Elements(bodies.len() as u64));

    group.bench_function("raw_bytes", |bencher| {
        bencher.iter(|| {
            let mut total = 0usize;
            for body in &bodies {
                let record = RawRecord::new(black_box(body)).expect("valid");
                total += record.raw_bytes().len();
            }
            black_box(total)
        });
    });

    group.bench_function("parse_every_field", |bencher| {
        // The closest honest comparison without a second BAM library: parse
        // every field the raw path deliberately skips.
        bencher.iter(|| {
            let mut total = 0u64;
            for body in &bodies {
                let record = RawRecord::new(black_box(body)).expect("valid");
                total += u64::from(record.flags().expect("valid"));
                total += record.qname().expect("valid").len() as u64;
                total += u64::from(record.reference_span().expect("valid").total());
                let mut tags = record.tags().expect("valid");
                while let Some(field) = tags.next() {
                    let (_, value) = field.expect("valid");
                    total += value.to_routing_bytes().len() as u64;
                }
            }
            black_box(total)
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    parsing,
    cigar,
    tags,
    routing,
    round_trip,
    raw_versus_decoded
);
criterion_main!(benches);
