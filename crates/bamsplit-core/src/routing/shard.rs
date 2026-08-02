// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Deterministic computational sharding.
//!
//! # The hash
//!
//! ```text
//! shard = xxh3_64_with_seed(key_bytes, seed) % shard_count
//! ```
//!
//! with [`DEFAULT_SEED`] unless `--seed` overrides it. XXH3 is used because it
//! is fast enough to disappear next to BGZF decompression, and because its
//! output is fully specified — the same bytes and seed give the same shard on
//! every platform, architecture, and `bamsplit` version.
//!
//! # Template coherence
//!
//! With `--key qname` every record sharing a QNAME lands in the same shard:
//! both mates, every secondary alignment, and every supplementary alignment,
//! because the shard depends only on the name. No QNAME table is kept, so
//! memory does not grow with the number of templates.
//!
//! `--key record` deliberately gives up that property; see
//! [`ShardKeySource::Record`].

use crate::bam::header::{BamHeader, ReadGroupField};
use crate::bam::raw_record::{RawRecord, RecordLocation};
use crate::error::{ConfigError, RoutingError};
use crate::routing::{ByteKey, DropReason, MissingPolicy, Route, Router};

/// The default hash seed.
///
/// A fixed, documented constant so two runs of the same version — and two
/// different tools implementing the same rule — agree. Changing it would
/// silently reshuffle every existing shard layout, so it is part of the public
/// interface.
pub const DEFAULT_SEED: u64 = 0x6261_6d73_706c_6974; // "bamsplit" in ASCII

/// What is hashed to pick a shard.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShardKeySource {
    /// The complete QNAME. Preserves template, mate, and alignment coherence.
    QName,
    /// The record's own bytes.
    ///
    /// Spreads records evenly regardless of template structure, and therefore
    /// **does not** keep mates, secondary alignments, or supplementary
    /// alignments together. Use it only when downstream work is per-record.
    Record,
    /// An auxiliary tag's value.
    Tag([u8; 2]),
    /// A `@RG`-derived field.
    ReadGroupField(ReadGroupField),
}

impl ShardKeySource {
    /// The option value a user writes.
    #[must_use]
    pub fn as_string(&self) -> String {
        match self {
            Self::QName => "qname".to_string(),
            Self::Record => "record".to_string(),
            Self::Tag(tag) => format!("tag:{}", crate::bam::tags::render_tag(*tag)),
            Self::ReadGroupField(field) => field.as_str().to_string(),
        }
    }

    /// Whether this source keeps every record of a template together.
    #[must_use]
    pub const fn preserves_template_coherence(&self) -> bool {
        !matches!(self, Self::Record)
    }
}

/// How a shard split is configured.
#[derive(Debug, Clone)]
pub struct ShardRouterOptions {
    /// How many shards to produce.
    pub shards: u32,
    /// What to hash.
    pub key: ShardKeySource,
    /// The hash seed.
    pub seed: u64,
    /// What to do when the key is absent.
    pub missing: MissingPolicy,
    /// The stem of the output for records with no key.
    pub missing_name: Vec<u8>,
    /// Zero-padding width; derived from `shards` when [`None`].
    pub shard_width: Option<usize>,
}

impl Default for ShardRouterOptions {
    fn default() -> Self {
        Self {
            shards: 1,
            key: ShardKeySource::QName,
            seed: DEFAULT_SEED,
            missing: MissingPolicy::File,
            missing_name: b"no-key".to_vec(),
            shard_width: None,
        }
    }
}

/// Routes records into a fixed number of deterministic shards.
#[derive(Debug)]
pub struct ShardRouter {
    options: ShardRouterOptions,
    width: usize,
    /// Pre-rendered shard names, so routing never formats a string.
    names: Vec<ByteKey>,
    missing_key: ByteKey,
}

impl ShardRouter {
    /// Builds a router.
    ///
    /// # Errors
    ///
    /// Returns [`ConfigError::OutOfRange`] if `shards` is zero or the requested
    /// padding width cannot hold the largest shard number.
    pub fn new(options: ShardRouterOptions) -> Result<Self, ConfigError> {
        if options.shards == 0 {
            return Err(ConfigError::OutOfRange {
                option: "--shards",
                constraint: "at least 1",
                value: options.shards.to_string(),
            });
        }
        let natural = (options.shards - 1).to_string().len();
        let width = options.shard_width.unwrap_or(natural.max(4));
        if width < natural {
            return Err(ConfigError::OutOfRange {
                option: "--shard-width",
                constraint: "wide enough for the largest shard number",
                value: width.to_string(),
            });
        }

        let names = (0..options.shards)
            .map(|index| ByteKey::new(format!("shard-{index:0width$}").into_bytes()))
            .collect();

        Ok(Self {
            missing_key: ByteKey::new(options.missing_name.clone()),
            options,
            width,
            names,
        })
    }

    /// The zero-padding width in use.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    /// The shard a byte string maps to.
    ///
    /// Exposed so callers and tests can verify coherence without constructing a
    /// record.
    #[must_use]
    pub fn shard_of(&self, bytes: &[u8]) -> u32 {
        let hash = xxhash_rust::xxh3::xxh3_64_with_seed(bytes, self.options.seed);
        // `shards >= 1` was checked at construction.
        (hash % u64::from(self.options.shards)) as u32
    }

    /// The bytes this router hashes for a record, or [`None`] when absent.
    fn key_bytes<'a>(
        &self,
        header: &'a BamHeader,
        record: &RawRecord<'a>,
    ) -> Result<Option<std::borrow::Cow<'a, [u8]>>, RoutingError> {
        use std::borrow::Cow;

        Ok(match &self.options.key {
            ShardKeySource::QName => Some(Cow::Borrowed(record.qname().map_err(Box::new)?)),
            ShardKeySource::Record => Some(Cow::Borrowed(record.raw_bytes())),
            ShardKeySource::Tag(tag) => record
                .tag(*tag)
                .map_err(Box::new)?
                .map(|value| Cow::Owned(value.to_routing_bytes())),
            ShardKeySource::ReadGroupField(field) => {
                let Some(group) = read_group(header, record)? else {
                    return Ok(None);
                };
                group.field(*field).map(Cow::Borrowed)
            }
        })
    }
}

/// Resolves a record's `RG` tag against the header.
///
/// Returns [`None`] when the record has no `RG` or the header does not declare
/// it; the caller applies the configured policy.
fn read_group<'a>(
    header: &'a BamHeader,
    record: &RawRecord<'_>,
) -> Result<Option<&'a crate::bam::header::ReadGroupInfo>, RoutingError> {
    let Some(value) = record.tag(*b"RG").map_err(Box::new)? else {
        return Ok(None);
    };
    let Some(id) = value.as_bytes() else {
        // `RG` must be a `Z` string; anything else is a malformed header/record
        // pairing rather than a missing value.
        return Err(RoutingError::MissingReadGroup {
            location: RecordLocation::UNKNOWN,
        });
    };
    Ok(header.read_group(id))
}

impl Router for ShardRouter {
    type Key = ByteKey;

    fn route(
        &self,
        header: &BamHeader,
        record: &RawRecord<'_>,
    ) -> Result<Route<Self::Key>, RoutingError> {
        let Some(bytes) = self.key_bytes(header, record)? else {
            return Ok(match self.options.missing {
                MissingPolicy::File => Route::One(self.missing_key.clone()),
                MissingPolicy::Drop => Route::Drop(missing_reason(&self.options.key)),
                MissingPolicy::Error => return Err(missing_error(&self.options.key)),
            });
        };
        let shard = self.shard_of(&bytes);
        // `shard < shards == names.len()`.
        self.names.get(shard as usize).cloned().map_or_else(
            || {
                Err(RoutingError::InvalidReferenceId {
                    id: shard as i32,
                    reference_count: self.names.len(),
                    location: RecordLocation::UNKNOWN,
                })
            },
            |key| Ok(Route::One(key)),
        )
    }

    fn mode(&self) -> &'static str {
        "shard"
    }

    fn declared_keys(&self, _header: &BamHeader) -> Vec<Self::Key> {
        let mut keys = self.names.clone();
        if self.options.missing == MissingPolicy::File {
            keys.push(self.missing_key.clone());
        }
        keys
    }
}

fn missing_reason(key: &ShardKeySource) -> DropReason {
    match key {
        ShardKeySource::ReadGroupField(_) => DropReason::UnknownReadGroup,
        // QNAME and record bytes always exist for a validated record, so those
        // arms are unreachable in practice; `MissingTag` is the honest label for
        // "the value this router hashes was not there".
        ShardKeySource::Tag(_) | ShardKeySource::QName | ShardKeySource::Record => {
            DropReason::MissingTag
        }
    }
}

fn missing_error(key: &ShardKeySource) -> RoutingError {
    let location = RecordLocation::UNKNOWN;
    match key {
        ShardKeySource::Tag(tag) => RoutingError::MissingTag {
            tag: crate::bam::tags::render_tag(*tag),
            location,
        },
        ShardKeySource::ReadGroupField(field) => RoutingError::MissingReadGroupField {
            id: "<unresolved>".to_string(),
            field: field.sam_tag(),
        },
        ShardKeySource::QName | ShardKeySource::Record => {
            RoutingError::MissingReadGroup { location }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::header::ReferenceSequence;
    use crate::routing::RoutingKey as _;

    fn header_with_read_groups() -> BamHeader {
        BamHeader::from_parts(
            b"@HD\tVN:1.6\n@RG\tID:rg0\tSM:sample0\tLB:lib0\tPU:unit0\n\
              @RG\tID:rg1\tSM:sample1\tLB:lib1\n"
                .to_vec(),
            vec![ReferenceSequence {
                name: b"chr1".to_vec(),
                length: 1000,
            }],
        )
        .expect("valid header")
    }

    fn record(name: &[u8], data: &[u8]) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.push(u8::try_from(name.len() + 1).expect("short"));
        body.extend_from_slice(&[60]);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(name);
        body.push(0);
        body.extend_from_slice(data);
        body
    }

    fn shard_label(router: &ShardRouter, header: &BamHeader, body: &[u8]) -> String {
        router
            .route(header, &RawRecord::new(body).expect("valid"))
            .expect("routable")
            .keys()[0]
            .label()
            .into_owned()
    }

    #[test]
    fn shard_names_are_zero_padded_from_the_shard_count() {
        for (shards, expected) in [(4u32, "shard-0000"), (10_000, "shard-0000")] {
            let router = ShardRouter::new(ShardRouterOptions {
                shards,
                ..ShardRouterOptions::default()
            })
            .expect("built");
            assert_eq!(
                router.declared_keys(&header_with_read_groups())[0]
                    .label()
                    .into_owned(),
                expected
            );
        }
        let wide = ShardRouter::new(ShardRouterOptions {
            shards: 100_000,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        assert_eq!(wide.width(), 5);
    }

    #[test]
    fn an_explicit_width_is_honoured_and_validated() {
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 8,
            shard_width: Some(2),
            ..ShardRouterOptions::default()
        })
        .expect("built");
        assert_eq!(
            router.declared_keys(&header_with_read_groups())[0]
                .label()
                .into_owned(),
            "shard-00"
        );

        let error = ShardRouter::new(ShardRouterOptions {
            shards: 1_000,
            shard_width: Some(2),
            ..ShardRouterOptions::default()
        })
        .expect_err("must reject");
        assert!(matches!(error, ConfigError::OutOfRange { .. }), "{error}");
    }

    #[test]
    fn zero_shards_is_rejected() {
        let error = ShardRouter::new(ShardRouterOptions {
            shards: 0,
            ..ShardRouterOptions::default()
        })
        .expect_err("must reject");
        assert!(matches!(error, ConfigError::OutOfRange { .. }), "{error}");
    }

    #[test]
    fn every_record_of_a_template_lands_in_one_shard() {
        let header = header_with_read_groups();
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 32,
            ..ShardRouterOptions::default()
        })
        .expect("built");

        // Mate, secondary, and supplementary records of one template differ in
        // flags and auxiliary data but share a QNAME.
        let name = b"template-000123";
        let variants = [
            record(name, b""),
            record(name, b"NMi\x01\x00\x00\x00"),
            record(name, b"RGZrg0\x00"),
        ];
        let expected = shard_label(&router, &header, &variants[0]);
        for variant in &variants[1..] {
            assert_eq!(shard_label(&router, &header, variant), expected);
        }
    }

    #[test]
    fn hashing_is_deterministic_across_routers() {
        let first = ShardRouter::new(ShardRouterOptions {
            shards: 64,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        let second = ShardRouter::new(ShardRouterOptions {
            shards: 64,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        for name in ["a", "read-1", "x".repeat(300).as_str()] {
            assert_eq!(
                first.shard_of(name.as_bytes()),
                second.shard_of(name.as_bytes())
            );
        }
    }

    #[test]
    fn the_seed_changes_the_assignment() {
        let default = ShardRouter::new(ShardRouterOptions {
            shards: 64,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        let seeded = ShardRouter::new(ShardRouterOptions {
            shards: 64,
            seed: 12345,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        let differences = (0..200u32)
            .filter(|index| {
                let name = format!("read{index}");
                default.shard_of(name.as_bytes()) != seeded.shard_of(name.as_bytes())
            })
            .count();
        assert!(differences > 150, "only {differences} of 200 differed");
    }

    #[test]
    fn shards_are_reasonably_balanced() {
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 16,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        let mut counts = [0u32; 16];
        for index in 0..16_000u32 {
            let name = format!("read{index}");
            counts[router.shard_of(name.as_bytes()) as usize] += 1;
        }
        for count in counts {
            assert!((700..1300).contains(&count), "shard sizes {counts:?}");
        }
    }

    #[test]
    fn a_single_shard_collects_everything() {
        let header = header_with_read_groups();
        let router = ShardRouter::new(ShardRouterOptions::default()).expect("built");
        for index in 0..20u32 {
            let body = record(format!("r{index}").as_bytes(), b"");
            assert_eq!(shard_label(&router, &header, &body), "shard-0000");
        }
    }

    #[test]
    fn record_based_sharding_does_not_preserve_coherence() {
        let header = header_with_read_groups();
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 64,
            key: ShardKeySource::Record,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        assert!(!ShardKeySource::Record.preserves_template_coherence());

        let name = b"template-1";
        let first = shard_label(&router, &header, &record(name, b""));
        let second = shard_label(&router, &header, &record(name, b"NMi\x01\x00\x00\x00"));
        assert_ne!(first, second, "distinct bodies should usually diverge");
    }

    #[test]
    fn tag_based_sharding_uses_the_canonical_value() {
        let header = header_with_read_groups();
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 16,
            key: ShardKeySource::Tag(*b"CB"),
            ..ShardRouterOptions::default()
        })
        .expect("built");
        let first = shard_label(&router, &header, &record(b"a", b"CBZbarcode1\x00"));
        let second = shard_label(&router, &header, &record(b"b", b"CBZbarcode1\x00"));
        assert_eq!(first, second, "the same barcode must share a shard");
        assert_eq!(first, format!("shard-{:04}", router.shard_of(b"barcode1")),);
    }

    #[test]
    fn header_derived_sharding_resolves_the_read_group() {
        let header = header_with_read_groups();
        for (field, value) in [
            (ReadGroupField::Sample, &b"sample0"[..]),
            (ReadGroupField::Library, b"lib0"),
            (ReadGroupField::PlatformUnit, b"unit0"),
            (ReadGroupField::ReadGroup, b"rg0"),
        ] {
            let router = ShardRouter::new(ShardRouterOptions {
                shards: 16,
                key: ShardKeySource::ReadGroupField(field),
                ..ShardRouterOptions::default()
            })
            .expect("built");
            let body = record(b"a", b"RGZrg0\x00");
            assert_eq!(
                shard_label(&router, &header, &body),
                format!("shard-{:04}", router.shard_of(value)),
                "{field}"
            );
        }
    }

    #[test]
    fn a_missing_key_follows_the_policy() {
        let header = header_with_read_groups();
        let options = |missing| ShardRouterOptions {
            shards: 4,
            key: ShardKeySource::Tag(*b"CB"),
            missing,
            ..ShardRouterOptions::default()
        };
        let body = record(b"a", b"");

        let filing = ShardRouter::new(options(MissingPolicy::File)).expect("built");
        assert_eq!(shard_label(&filing, &header, &body), "no-key");

        let dropping = ShardRouter::new(options(MissingPolicy::Drop)).expect("built");
        assert_eq!(
            dropping
                .route(&header, &RawRecord::new(&body).expect("valid"))
                .expect("routable"),
            Route::Drop(DropReason::MissingTag)
        );

        let strict = ShardRouter::new(options(MissingPolicy::Error)).expect("built");
        let error = strict
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        assert!(matches!(error, RoutingError::MissingTag { .. }), "{error}");
    }

    #[test]
    fn an_unknown_read_group_counts_as_missing() {
        let header = header_with_read_groups();
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 4,
            key: ShardKeySource::ReadGroupField(ReadGroupField::Sample),
            ..ShardRouterOptions::default()
        })
        .expect("built");
        let body = record(b"a", b"RGZabsent\x00");
        assert_eq!(shard_label(&router, &header, &body), "no-key");
    }

    #[test]
    fn a_read_group_missing_the_requested_field_counts_as_missing() {
        let header = header_with_read_groups();
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 4,
            key: ShardKeySource::ReadGroupField(ReadGroupField::PlatformUnit),
            ..ShardRouterOptions::default()
        })
        .expect("built");
        // rg1 has no PU.
        let body = record(b"a", b"RGZrg1\x00");
        assert_eq!(shard_label(&router, &header, &body), "no-key");
    }

    #[test]
    fn declared_keys_cover_every_shard_plus_the_missing_output() {
        let header = header_with_read_groups();
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 3,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        let labels: Vec<String> = router
            .declared_keys(&header)
            .iter()
            .map(|key| key.label().into_owned())
            .collect();
        assert_eq!(labels, ["shard-0000", "shard-0001", "shard-0002", "no-key"]);
    }

    #[test]
    fn key_sources_render_their_option_values() {
        assert_eq!(ShardKeySource::QName.as_string(), "qname");
        assert_eq!(ShardKeySource::Record.as_string(), "record");
        assert_eq!(ShardKeySource::Tag(*b"CB").as_string(), "tag:CB");
        assert_eq!(
            ShardKeySource::ReadGroupField(ReadGroupField::PlatformUnit).as_string(),
            "platform-unit"
        );
    }

    #[test]
    fn the_default_seed_is_pinned() {
        // Changing this constant reshuffles every existing shard layout, so it
        // is asserted rather than merely documented.
        assert_eq!(DEFAULT_SEED, 0x6261_6d73_706c_6974);
        let router = ShardRouter::new(ShardRouterOptions {
            shards: 1_000,
            ..ShardRouterOptions::default()
        })
        .expect("built");
        assert_eq!(
            router.shard_of(b"template-000123"),
            (xxhash_rust::xxh3::xxh3_64_with_seed(b"template-000123", DEFAULT_SEED) % 1_000) as u32
        );
    }
}
