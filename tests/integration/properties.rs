// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Property tests over generated inputs.
//!
//! Example-based tests pin the cases we thought of. These pin the ones we did
//! not: `proptest` generates record sets, reference names, and byte strings, and
//! each property below must hold for all of them.

use bamsplit_core::output::filename::{decode, encode, is_safe_component};
use bamsplit_core::routing::shard::{ShardRouter, ShardRouterOptions};
use bamsplit_core::{ChromOptions, ShardOptions};
use bamsplit_fixtures::{Fixture, RecordSpec};
use proptest::prelude::*;

use crate::support::*;

/// Reference names that are legal SAM `SN` values.
fn reference_name() -> impl Strategy<Value = String> {
    proptest::string::string_regex("[0-9A-Za-z!#$%&+./:;?@^_|~-]{1,12}").expect("a valid regex")
}

/// A record set over `reference_count` references, in coordinate order.
fn sorted_records(reference_count: usize) -> impl Strategy<Value = Vec<RecordSpec>> {
    proptest::collection::vec((0..reference_count, 0..10_000u32, 0..4u16), 0..40usize).prop_map(
        move |mut triples| {
            triples.sort();
            triples
                .into_iter()
                .enumerate()
                .map(|(index, (reference_id, position, flag_pick))| {
                    let flags = match flag_pick {
                        0 => 0,
                        1 => bamsplit_fixtures::flags::SECONDARY,
                        2 => bamsplit_fixtures::flags::SUPPLEMENTARY,
                        _ => bamsplit_fixtures::flags::UNMAPPED,
                    };
                    RecordSpec::mapped(
                        &format!("read{index:04}"),
                        i32::try_from(reference_id).unwrap_or(0),
                        i32::try_from(position).unwrap_or(0),
                    )
                    .with_flags(flags)
                })
                .collect()
        },
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    /// Splitting by chromosome conserves every raw record body, for any set of
    /// coordinate-ordered records.
    #[test]
    fn chromosome_splitting_conserves_record_bodies(records in sorted_records(3)) {
        let fixture = Fixture::new(
            &[("chrA", 20_000), ("chrB", 20_000), ("chrC", 20_000)],
            "coordinate",
        )
        .with_records(records);

        let workspace = Workspace::with(&fixture);
        let out = workspace.path("out");
        bamsplit_core::split_by_chromosome(
            &workspace.input,
            &ChromOptions::new(&out),
            &run_options(),
        )
        .expect("split succeeds");

        assert_lossless(&fixture, &out);
        read_manifest(&out).validate().expect("conservation holds");
    }

    /// Every output BAM keeps the complete input reference dictionary.
    #[test]
    fn outputs_keep_the_full_reference_dictionary(records in sorted_records(3)) {
        let fixture = Fixture::new(
            &[("chrA", 20_000), ("chrB", 20_000), ("chrC", 20_000)],
            "coordinate",
        )
        .with_records(records);
        let workspace = Workspace::with(&fixture);
        let out = workspace.path("out");
        bamsplit_core::split_by_chromosome(
            &workspace.input,
            &ChromOptions::new(&out),
            &run_options(),
        )
        .expect("split succeeds");

        for entry in std::fs::read_dir(&out).expect("listable").flatten() {
            let path = entry.path();
            if path.extension().and_then(|e| e.to_str()) != Some("bam") {
                continue;
            }
            let (header, _) = read_bam(&path);
            prop_assert_eq!(header.reference_count(), 3);
            prop_assert_eq!(header.reference_name(0), Some(&b"chrA"[..]));
        }
    }

    /// A record is emitted at most once by every non-duplicating mode.
    #[test]
    fn non_duplicating_modes_emit_each_record_once(records in sorted_records(2)) {
        let fixture = Fixture::new(&[("chrA", 20_000), ("chrB", 20_000)], "coordinate")
            .with_records(records);
        let workspace = Workspace::with(&fixture);

        for (name, shards) in [("chrom", None), ("shard", Some(4u32))] {
            let out = workspace.path(name);
            if let Some(shards) = shards {
                bamsplit_core::split_by_shard(
                    &workspace.input,
                    &ShardOptions::new(&out, shards),
                    &run_options(),
                )
                .expect("split succeeds");
            } else {
                bamsplit_core::split_by_chromosome(
                    &workspace.input,
                    &ChromOptions::new(&out),
                    &run_options(),
                )
                .expect("split succeeds");
            }
            let manifest = read_manifest(&out);
            prop_assert_eq!(manifest.duplicate_emissions, 0);
            prop_assert_eq!(
                manifest.total_output_emissions,
                manifest.unique_emitted_records
            );
        }
    }

    /// All records sharing a QNAME land in the same shard, for any shard count.
    #[test]
    fn qname_sharding_keeps_templates_together(
        names in proptest::collection::vec("[a-z]{1,10}", 1..30),
        shards in 1..64u32,
    ) {
        let router = ShardRouter::new(ShardRouterOptions {
            shards,
            ..ShardRouterOptions::default()
        })
        .expect("built");

        for name in &names {
            let first = router.shard_of(name.as_bytes());
            // Recomputing must agree, and must be inside the shard range.
            prop_assert_eq!(first, router.shard_of(name.as_bytes()));
            prop_assert!(first < shards);
        }
    }

    /// Filename encoding is reversible for any byte string short enough to
    /// escape the length cap.
    #[test]
    fn filename_encoding_round_trips(bytes in proptest::collection::vec(any::<u8>(), 0..60)) {
        let encoded = encode(&bytes);
        prop_assert!(is_safe_component(&encoded), "{}", encoded);
        let decoded = decode(&encoded);
        prop_assert_eq!(decoded.as_deref(), Some(&bytes[..]));
    }

    /// An encoded name can never escape its directory, whatever the input.
    #[test]
    fn encoded_names_cannot_escape(bytes in proptest::collection::vec(any::<u8>(), 0..300)) {
        let encoded = encode(&bytes);
        let joined = std::path::Path::new("/output").join(&encoded);
        prop_assert!(joined.starts_with("/output"));
        prop_assert_eq!(joined.components().count(), 3);
        prop_assert!(!encoded.contains('/'));
        prop_assert!(encoded.len() <= bamsplit_core::output::filename::MAX_STEM_LEN);
    }

    /// Distinct reference names always produce distinct output files.
    #[test]
    fn distinct_reference_names_produce_distinct_files(
        names in proptest::collection::hash_set(reference_name(), 1..6),
    ) {
        let names: Vec<String> = names.into_iter().collect();
        let references: Vec<(&str, u32)> =
            names.iter().map(|name| (name.as_str(), 1_000u32)).collect();
        let fixture = Fixture::new(&references, "coordinate").with_records(
            (0..names.len()).map(|index| {
                RecordSpec::mapped(
                    &format!("r{index}"),
                    i32::try_from(index).unwrap_or(0),
                    10,
                )
            }),
        );

        let workspace = Workspace::with(&fixture);
        let out = workspace.path("out");
        bamsplit_core::split_by_chromosome(
            &workspace.input,
            &ChromOptions::new(&out),
            &run_options(),
        )
        .expect("split succeeds");

        // One BAM per reference plus the unplaced output that was declared but
        // skipped; every one distinct.
        let outputs = collect_outputs(&out);
        prop_assert_eq!(outputs.len(), names.len());
        assert_lossless(&fixture, &out);
    }

    /// Output ordering in the manifest is deterministic.
    #[test]
    fn manifest_output_order_is_deterministic(records in sorted_records(3)) {
        let fixture = Fixture::new(
            &[("chrA", 20_000), ("chrB", 20_000), ("chrC", 20_000)],
            "coordinate",
        )
        .with_records(records);
        let workspace = Workspace::with(&fixture);

        let order = |name: &str| {
            let out = workspace.path(name);
            bamsplit_core::split_by_chromosome(
                &workspace.input,
                &ChromOptions::new(&out),
                &run_options(),
            )
            .expect("split succeeds");
            read_manifest(&out)
                .outputs
                .into_iter()
                .map(|output| output.logical_key)
                .collect::<Vec<_>>()
        };
        prop_assert_eq!(order("first"), order("second"));
    }
}
