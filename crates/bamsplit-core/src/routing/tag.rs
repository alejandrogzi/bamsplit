// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Routing by an auxiliary tag or a `@RG`-derived field.
//!
//! # Cardinality is the hazard
//!
//! `--tag RG` produces a handful of outputs. `--tag CB` on a single-cell BAM
//! produces hundreds of thousands, one file per cell barcode, and a naive
//! implementation would exhaust file descriptors, inodes, and patience.
//!
//! So this router counts distinct keys as it discovers them and refuses to
//! exceed `--max-outputs`. The failure is deliberate and early: the run stops,
//! every partial file is removed, and the error names the alternative —
//! `bamsplit shard --key tag:CB`, which gives the same grouping guarantee in a
//! bounded number of files.
//!
//! The counter lives behind a lock because the indexed engine routes from
//! several threads. It is touched once per *new* key, not once per record, so
//! contention is negligible.
//!
//! # Array-valued tags
//!
//! A `B` array is rejected by default. `CB:B:i,1,2,3` has no obvious filename
//! and no obvious equality: is it the same key as `CB:B:i,3,2,1`? Rather than
//! pick an answer silently, `bamsplit` asks for `--allow-array-tags`, and then
//! uses the canonical `<subtype>,<v0>,<v1>,…` rendering, which is
//! order-sensitive.
//!
//! # Floating-point values
//!
//! Never formatted through anything locale-sensitive. See
//! [`crate::bam::tags::RawTagValue::to_routing_bytes`].

use std::collections::HashSet;
use std::sync::Mutex;

use crate::bam::header::{BamHeader, ReadGroupField};
use crate::bam::raw_record::{RawRecord, RecordLocation};
use crate::error::{ConfigError, RoutingError};
use crate::routing::{ByteKey, DropReason, MissingPolicy, Route, Router};

/// What the router reads from each record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagSource {
    /// An auxiliary tag, used verbatim.
    Tag([u8; 2]),
    /// A field resolved from the record's read group.
    Field(ReadGroupField),
}

impl TagSource {
    /// The option value a user writes.
    #[must_use]
    pub fn as_string(&self) -> String {
        match self {
            Self::Tag(tag) => crate::bam::tags::render_tag(*tag),
            Self::Field(field) => field.as_str().to_string(),
        }
    }
}

/// How a tag split is configured.
#[derive(Debug, Clone)]
pub struct TagRouterOptions {
    /// What to route on.
    pub source: TagSource,
    /// What to do when the value is absent.
    pub missing: MissingPolicy,
    /// The stem of the output for records with no value.
    pub missing_name: Vec<u8>,
    /// What to do when a record's read group is not in the header.
    pub unknown_read_group: MissingPolicy,
    /// The most distinct outputs allowed.
    pub max_outputs: usize,
    /// Whether to lift the cardinality limit entirely.
    pub allow_high_cardinality: bool,
    /// Whether `B` arrays may be used as keys.
    pub allow_array_tags: bool,
}

impl Default for TagRouterOptions {
    fn default() -> Self {
        Self {
            source: TagSource::Field(ReadGroupField::ReadGroup),
            missing: MissingPolicy::File,
            missing_name: b"no-tag".to_vec(),
            unknown_read_group: MissingPolicy::File,
            max_outputs: 1_000,
            allow_high_cardinality: false,
            allow_array_tags: false,
        }
    }
}

/// Routes records to one output per distinct tag value.
#[derive(Debug)]
pub struct TagRouter {
    options: TagRouterOptions,
    missing_key: ByteKey,
    /// Distinct keys seen so far, for the cardinality check.
    seen: Mutex<HashSet<Vec<u8>>>,
}

impl TagRouter {
    /// Builds a router.
    ///
    /// # Errors
    ///
    /// Returns [`ConfigError::OutOfRange`] if `max_outputs` is zero while the
    /// limit is in force.
    pub fn new(options: TagRouterOptions) -> Result<Self, ConfigError> {
        if options.max_outputs == 0 && !options.allow_high_cardinality {
            return Err(ConfigError::OutOfRange {
                option: "--max-outputs",
                constraint: "at least 1, or use `--allow-high-cardinality`",
                value: "0".to_string(),
            });
        }
        Ok(Self {
            missing_key: ByteKey::new(options.missing_name.clone()),
            options,
            seen: Mutex::new(HashSet::new()),
        })
    }

    /// How many distinct keys have been discovered.
    #[must_use]
    pub fn discovered_keys(&self) -> usize {
        self.seen.lock().map_or(0, |seen| seen.len())
    }

    /// Records a key and enforces the cardinality limit.
    fn admit(&self, key: &[u8]) -> Result<(), RoutingError> {
        let Ok(mut seen) = self.seen.lock() else {
            // A poisoned lock means another thread panicked mid-routing. Rather
            // than mask that as a routing failure, stop counting and let the
            // run continue to its own error; the manifest's conservation check
            // is the backstop.
            return Ok(());
        };
        if seen.contains(key) {
            return Ok(());
        }
        if !self.options.allow_high_cardinality && seen.len() >= self.options.max_outputs {
            return Err(RoutingError::CardinalityExceeded {
                discovered: seen.len() + 1,
                limit: self.options.max_outputs,
                hint: self.options.source.as_string(),
            });
        }
        seen.insert(key.to_vec());
        Ok(())
    }

    fn resolve(
        &self,
        header: &BamHeader,
        record: &RawRecord<'_>,
    ) -> Result<Resolution, RoutingError> {
        let location = RecordLocation::UNKNOWN;
        match &self.options.source {
            TagSource::Tag(tag) => {
                let Some(value) = record.tag(*tag).map_err(Box::new)? else {
                    return Ok(Resolution::Missing);
                };
                if value.is_array() && !self.options.allow_array_tags {
                    return Err(RoutingError::ArrayValuedTag {
                        tag: crate::bam::tags::render_tag(*tag),
                        location,
                    });
                }
                Ok(Resolution::Found(value.to_routing_bytes()))
            }
            TagSource::Field(field) => {
                let Some(value) = record.tag(*b"RG").map_err(Box::new)? else {
                    return Ok(Resolution::Missing);
                };
                let Some(id) = value.as_bytes() else {
                    return Err(RoutingError::MissingReadGroup { location });
                };
                let Some(group) = header.read_group(id) else {
                    return Ok(Resolution::UnknownReadGroup(
                        String::from_utf8_lossy(id).into_owned(),
                    ));
                };
                match group.field(*field) {
                    Some(value) => Ok(Resolution::Found(value.to_vec())),
                    None => Ok(Resolution::MissingField {
                        id: String::from_utf8_lossy(id).into_owned(),
                        field: field.sam_tag(),
                    }),
                }
            }
        }
    }
}

enum Resolution {
    Found(Vec<u8>),
    Missing,
    UnknownReadGroup(String),
    MissingField { id: String, field: &'static str },
}

impl Router for TagRouter {
    type Key = ByteKey;

    fn route(
        &self,
        header: &BamHeader,
        record: &RawRecord<'_>,
    ) -> Result<Route<Self::Key>, RoutingError> {
        let location = RecordLocation::UNKNOWN;
        match self.resolve(header, record)? {
            Resolution::Found(bytes) => {
                self.admit(&bytes)?;
                Ok(Route::One(ByteKey::new(bytes)))
            }
            Resolution::Missing => Ok(match self.options.missing {
                MissingPolicy::File => Route::One(self.missing_key.clone()),
                MissingPolicy::Drop => Route::Drop(DropReason::MissingTag),
                MissingPolicy::Error => {
                    return Err(RoutingError::MissingTag {
                        tag: self.options.source.as_string(),
                        location,
                    });
                }
            }),
            Resolution::UnknownReadGroup(id) => Ok(match self.options.unknown_read_group {
                MissingPolicy::File => Route::One(self.missing_key.clone()),
                MissingPolicy::Drop => Route::Drop(DropReason::UnknownReadGroup),
                MissingPolicy::Error => {
                    return Err(RoutingError::UnknownReadGroup { id, location });
                }
            }),
            Resolution::MissingField { id, field } => Ok(match self.options.missing {
                MissingPolicy::File => Route::One(self.missing_key.clone()),
                MissingPolicy::Drop => Route::Drop(DropReason::MissingReadGroupField),
                MissingPolicy::Error => {
                    return Err(RoutingError::MissingReadGroupField { id, field });
                }
            }),
        }
    }

    fn mode(&self) -> &'static str {
        "tag"
    }

    fn declared_keys(&self, header: &BamHeader) -> Vec<Self::Key> {
        // Only header-derived routing can know its keys in advance; an
        // auxiliary tag's values are only discoverable by reading the records.
        let mut keys: Vec<ByteKey> = match &self.options.source {
            TagSource::Field(field) => {
                let mut values: Vec<Vec<u8>> = header
                    .read_groups()
                    .iter()
                    .filter_map(|group| group.field(*field).map(<[u8]>::to_vec))
                    .collect();
                values.sort_unstable();
                values.dedup();
                values.into_iter().map(ByteKey::new).collect()
            }
            TagSource::Tag(_) => Vec::new(),
        };
        if self.options.missing == MissingPolicy::File
            || self.options.unknown_read_group == MissingPolicy::File
        {
            keys.push(self.missing_key.clone());
        }
        keys
    }
}

/// Validates that exactly one of `--tag` and `--field` was given.
///
/// # Errors
///
/// Returns [`ConfigError::ExactlyOneRequired`] when both or neither is present.
pub fn require_exactly_one_source(
    tag: Option<&str>,
    field: Option<&str>,
) -> Result<(), ConfigError> {
    if tag.is_some() == field.is_some() {
        return Err(ConfigError::ExactlyOneRequired {
            options: "`--tag`, `--field`",
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::header::ReferenceSequence;
    use crate::routing::RoutingKey as _;

    fn header() -> BamHeader {
        BamHeader::from_parts(
            b"@HD\tVN:1.6\n@RG\tID:rg0\tSM:sample0\tLB:lib0\tPU:unit0\tPL:ILLUMINA\tCN:centre\n\
              @RG\tID:rg1\tSM:sample1\n"
                .to_vec(),
            vec![ReferenceSequence {
                name: b"chr1".to_vec(),
                length: 1000,
            }],
        )
        .expect("valid header")
    }

    fn record(data: &[u8]) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&[2, 60]);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(b"r\0");
        body.extend_from_slice(data);
        body
    }

    fn label(router: &TagRouter, header: &BamHeader, body: &[u8]) -> String {
        router
            .route(header, &RawRecord::new(body).expect("valid"))
            .expect("routable")
            .keys()[0]
            .label()
            .into_owned()
    }

    #[test]
    fn routes_on_a_string_tag() {
        let header = header();
        let router = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"RG"),
            ..TagRouterOptions::default()
        })
        .expect("built");
        assert_eq!(label(&router, &header, &record(b"RGZrg0\x00")), "rg0");
        assert_eq!(router.discovered_keys(), 1);
    }

    #[test]
    fn routes_on_every_scalar_tag_type() {
        let header = header();
        let router = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"XX"),
            ..TagRouterOptions::default()
        })
        .expect("built");

        let mut character = b"XXA".to_vec();
        character.push(b'+');
        let mut int32 = b"XXi".to_vec();
        int32.extend_from_slice(&(-42i32).to_le_bytes());
        let mut float = b"XXf".to_vec();
        float.extend_from_slice(&0.5f32.to_le_bytes());

        assert_eq!(label(&router, &header, &record(&character)), "+");
        assert_eq!(label(&router, &header, &record(&int32)), "-42");
        assert_eq!(label(&router, &header, &record(&float)), "0.5");
        assert_eq!(label(&router, &header, &record(b"XXZtext\x00")), "text");
        assert_eq!(label(&router, &header, &record(b"XXHDEAD\x00")), "DEAD");
    }

    #[test]
    fn array_tags_are_rejected_by_default_and_allowed_on_request() {
        let header = header();
        let mut data = b"XXB".to_vec();
        data.push(b'i');
        data.extend_from_slice(&2u32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes());
        data.extend_from_slice(&2i32.to_le_bytes());
        let body = record(&data);

        let strict = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"XX"),
            ..TagRouterOptions::default()
        })
        .expect("built");
        let error = strict
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        assert!(
            matches!(error, RoutingError::ArrayValuedTag { .. }),
            "{error}"
        );

        let permissive = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"XX"),
            allow_array_tags: true,
            ..TagRouterOptions::default()
        })
        .expect("built");
        assert_eq!(label(&permissive, &header, &body), "i,1,2");
    }

    #[test]
    fn header_derived_fields_resolve_through_the_read_group() {
        let header = header();
        for (field, expected) in [
            (ReadGroupField::Sample, "sample0"),
            (ReadGroupField::Library, "lib0"),
            (ReadGroupField::PlatformUnit, "unit0"),
            (ReadGroupField::Platform, "ILLUMINA"),
            (ReadGroupField::SequencingCenter, "centre"),
            (ReadGroupField::ReadGroup, "rg0"),
        ] {
            let router = TagRouter::new(TagRouterOptions {
                source: TagSource::Field(field),
                ..TagRouterOptions::default()
            })
            .expect("built");
            assert_eq!(label(&router, &header, &record(b"RGZrg0\x00")), expected);
        }
    }

    #[test]
    fn a_missing_tag_follows_the_policy() {
        let header = header();
        let build = |missing| {
            TagRouter::new(TagRouterOptions {
                source: TagSource::Tag(*b"CB"),
                missing,
                ..TagRouterOptions::default()
            })
            .expect("built")
        };
        let body = record(b"");

        assert_eq!(label(&build(MissingPolicy::File), &header, &body), "no-tag");
        assert_eq!(
            build(MissingPolicy::Drop)
                .route(&header, &RawRecord::new(&body).expect("valid"))
                .expect("routable"),
            Route::Drop(DropReason::MissingTag)
        );
        let error = build(MissingPolicy::Error)
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        assert!(matches!(error, RoutingError::MissingTag { .. }), "{error}");
    }

    #[test]
    fn an_unknown_read_group_has_its_own_policy() {
        let header = header();
        let build = |policy| {
            TagRouter::new(TagRouterOptions {
                source: TagSource::Field(ReadGroupField::Sample),
                unknown_read_group: policy,
                ..TagRouterOptions::default()
            })
            .expect("built")
        };
        let body = record(b"RGZabsent\x00");

        assert_eq!(label(&build(MissingPolicy::File), &header, &body), "no-tag");
        assert_eq!(
            build(MissingPolicy::Drop)
                .route(&header, &RawRecord::new(&body).expect("valid"))
                .expect("routable"),
            Route::Drop(DropReason::UnknownReadGroup)
        );
        let error = build(MissingPolicy::Error)
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        assert!(
            matches!(error, RoutingError::UnknownReadGroup { .. }),
            "{error}"
        );
    }

    #[test]
    fn a_read_group_lacking_the_requested_field_follows_the_missing_policy() {
        let header = header();
        let router = TagRouter::new(TagRouterOptions {
            source: TagSource::Field(ReadGroupField::PlatformUnit),
            missing: MissingPolicy::Error,
            ..TagRouterOptions::default()
        })
        .expect("built");
        // rg1 declares SM but not PU.
        let body = record(b"RGZrg1\x00");
        let error = router
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        assert!(
            matches!(
                error,
                RoutingError::MissingReadGroupField { field: "PU", .. }
            ),
            "{error}"
        );
    }

    #[test]
    fn the_cardinality_limit_fails_with_an_actionable_message() {
        let header = header();
        let router = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"CB"),
            max_outputs: 3,
            ..TagRouterOptions::default()
        })
        .expect("built");

        for index in 0..3u32 {
            let mut data = b"CBZ".to_vec();
            data.extend_from_slice(format!("barcode{index}").as_bytes());
            data.push(0);
            label(&router, &header, &record(&data));
        }
        assert_eq!(router.discovered_keys(), 3);

        let body = record(b"CBZone-too-many\x00");
        let error = router
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        let rendered = error.to_string();
        assert!(
            matches!(error, RoutingError::CardinalityExceeded { limit: 3, .. }),
            "{rendered}"
        );
        assert!(
            rendered.contains("bamsplit shard --key tag:CB"),
            "{rendered}"
        );
    }

    #[test]
    fn repeated_keys_do_not_consume_cardinality() {
        let header = header();
        let router = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"CB"),
            max_outputs: 1,
            ..TagRouterOptions::default()
        })
        .expect("built");
        for _ in 0..100 {
            label(&router, &header, &record(b"CBZsame\x00"));
        }
        assert_eq!(router.discovered_keys(), 1);
    }

    #[test]
    fn the_limit_can_be_lifted() {
        let header = header();
        let router = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"CB"),
            max_outputs: 1,
            allow_high_cardinality: true,
            ..TagRouterOptions::default()
        })
        .expect("built");
        for index in 0..50u32 {
            let mut data = b"CBZ".to_vec();
            data.extend_from_slice(format!("bc{index}").as_bytes());
            data.push(0);
            label(&router, &header, &record(&data));
        }
        assert_eq!(router.discovered_keys(), 50);
    }

    #[test]
    fn a_zero_limit_is_rejected_unless_lifted() {
        let error = TagRouter::new(TagRouterOptions {
            max_outputs: 0,
            ..TagRouterOptions::default()
        })
        .expect_err("must reject");
        assert!(matches!(error, ConfigError::OutOfRange { .. }), "{error}");

        TagRouter::new(TagRouterOptions {
            max_outputs: 0,
            allow_high_cardinality: true,
            ..TagRouterOptions::default()
        })
        .expect("allowed");
    }

    #[test]
    fn header_derived_keys_are_declared_in_advance() {
        let header = header();
        let router = TagRouter::new(TagRouterOptions {
            source: TagSource::Field(ReadGroupField::Sample),
            ..TagRouterOptions::default()
        })
        .expect("built");
        let labels: Vec<String> = router
            .declared_keys(&header)
            .iter()
            .map(|key| key.label().into_owned())
            .collect();
        assert_eq!(labels, ["sample0", "sample1", "no-tag"]);

        // A raw auxiliary tag cannot be enumerated from the header.
        let opaque = TagRouter::new(TagRouterOptions {
            source: TagSource::Tag(*b"CB"),
            ..TagRouterOptions::default()
        })
        .expect("built");
        assert_eq!(opaque.declared_keys(&header).len(), 1);
    }

    #[test]
    fn exactly_one_source_is_required() {
        assert!(require_exactly_one_source(Some("RG"), None).is_ok());
        assert!(require_exactly_one_source(None, Some("sample")).is_ok());
        assert!(matches!(
            require_exactly_one_source(Some("RG"), Some("sample")),
            Err(ConfigError::ExactlyOneRequired { .. })
        ));
        assert!(matches!(
            require_exactly_one_source(None, None),
            Err(ConfigError::ExactlyOneRequired { .. })
        ));
    }
}
