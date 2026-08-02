// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Deciding *where* a record goes, and nothing else.
//!
//! A [`Router`] borrows a record and returns a [`Route`]. It never opens a
//! file, never learns an output path, and never counts anything. That
//! separation is what lets the same router run under three different engines —
//! streaming, indexed, and spooled — without change, and it is why the routing
//! tests can be pure functions over synthetic records.
//!
//! # Keys
//!
//! A [`RoutingKey`] is the identity of an output. It carries raw bytes, because
//! reference names and tag values are byte strings and are not required to be
//! UTF-8, and it derives everything else — sort order, filename, display label,
//! manifest fields — from those bytes plus whatever the router knows.
//!
//! # Fan-out
//!
//! [`Route::Many`] exists for exactly one documented case: `bamsplit region
//! --assignment overlap`, where one alignment can legitimately belong to
//! several annotated features. Every other router returns [`Route::One`] or
//! [`Route::Drop`], and the manifest's conservation check is what enforces that.

pub mod chrom;
#[cfg(feature = "annotation")]
pub mod region;
pub mod shard;
pub mod tag;

use std::borrow::Cow;
use std::hash::Hash;

use smallvec::SmallVec;

use crate::bam::header::BamHeader;
use crate::bam::raw_record::RawRecord;
use crate::error::RoutingError;

/// The identity of one output.
pub trait RoutingKey: Clone + Eq + Ord + Hash + Send + Sync + std::fmt::Debug {
    /// The logical value, as raw bytes. Not necessarily UTF-8.
    fn logical(&self) -> &[u8];

    /// A human-readable label, for logs and progress.
    ///
    /// Lossy by design: a label is for people, and the manifest keeps the exact
    /// bytes.
    fn label(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(self.logical())
    }

    /// The filesystem-safe encoded name.
    fn encoded(&self) -> String {
        crate::output::filename::encode(self.logical())
    }

    /// Extra manifest fields specific to this key kind.
    fn manifest_fields(&self) -> Vec<(&'static str, serde_json::Value)> {
        Vec::new()
    }
}

/// A routing key that is just a byte string.
///
/// Used by `chrom`, `shard`, and `tag`; ordering is lexicographic over the
/// bytes, which makes output ordering reproducible without depending on hash
/// iteration order.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteKey(Vec<u8>);

impl ByteKey {
    /// Creates a key.
    #[must_use]
    pub fn new(bytes: impl Into<Vec<u8>>) -> Self {
        Self(bytes.into())
    }

    /// Consumes the key, returning its bytes.
    #[must_use]
    pub fn into_bytes(self) -> Vec<u8> {
        self.0
    }
}

impl RoutingKey for ByteKey {
    fn logical(&self) -> &[u8] {
        &self.0
    }
}

impl From<&[u8]> for ByteKey {
    fn from(bytes: &[u8]) -> Self {
        Self(bytes.to_vec())
    }
}

/// Why a record was not written anywhere.
///
/// Every variant is counted separately in the manifest, because "dropped by
/// policy" and "matched nothing" are very different things to a user checking
/// that a split was lossless.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DropReason {
    /// `refID == -1` and `--unplaced drop`.
    Unplaced,
    /// The reference was excluded by `--include`/`--exclude`/`--reference-list`.
    ExcludedReference,
    /// The routing tag is absent and the policy is `drop`.
    MissingTag,
    /// The record's read group is not in the header and the policy is `drop`.
    UnknownReadGroup,
    /// The requested `@RG` sub-field is absent and the policy is `drop`.
    MissingReadGroupField,
    /// No annotation region matched.
    Unmatched,
}

impl DropReason {
    /// A short label for the manifest.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Unplaced => "unplaced",
            Self::ExcludedReference => "excluded-reference",
            Self::MissingTag => "missing-tag",
            Self::UnknownReadGroup => "unknown-read-group",
            Self::MissingReadGroupField => "missing-read-group-field",
            Self::Unmatched => "unmatched",
        }
    }

    /// Whether this drop should be counted as "unmatched" rather than
    /// "dropped".
    ///
    /// The conservation equation treats the two separately, and only region
    /// routing can produce an unmatched record.
    #[must_use]
    pub const fn is_unmatched(self) -> bool {
        matches!(self, Self::Unmatched)
    }
}

impl std::fmt::Display for DropReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Where one record goes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Route<K> {
    /// Exactly one output.
    One(K),
    /// Several outputs. Only region routing in `overlap` mode produces this.
    Many(SmallVec<[K; 2]>),
    /// Nowhere, for the given reason.
    Drop(DropReason),
}

impl<K> Route<K> {
    /// How many outputs this record will be written to.
    #[must_use]
    pub fn emission_count(&self) -> usize {
        match self {
            Self::One(_) => 1,
            Self::Many(keys) => keys.len(),
            Self::Drop(_) => 0,
        }
    }

    /// The keys, as a slice.
    #[must_use]
    pub fn keys(&self) -> &[K] {
        match self {
            Self::One(key) => std::slice::from_ref(key),
            Self::Many(keys) => keys,
            Self::Drop(_) => &[],
        }
    }
}

/// Decides where records go.
///
/// `Send + Sync` because the indexed engine routes from several threads at once
/// over a shared router; every implementation here is effectively immutable
/// after construction, and the one that is not ([`tag::TagRouter`]'s cardinality
/// counter) uses a lock.
pub trait Router: Send + Sync {
    /// The key type this router produces.
    type Key: RoutingKey;

    /// Routes one record.
    ///
    /// # Errors
    ///
    /// Returns [`RoutingError`] when the record is malformed in a way that
    /// prevents routing, or when a policy is set to `error` and its condition
    /// is met.
    fn route(
        &self,
        header: &BamHeader,
        record: &RawRecord<'_>,
    ) -> Result<Route<Self::Key>, RoutingError>;

    /// A short name for logs and the manifest, e.g. `"chrom"`.
    fn mode(&self) -> &'static str;

    /// Whether this router can emit a record to more than one output.
    ///
    /// The manifest picks its conservation equation from this.
    fn may_duplicate(&self) -> bool {
        false
    }

    /// Keys that must exist even if no record routes to them.
    ///
    /// `--emit-empty` materializes these as header-only BAMs; without it they
    /// are listed in the manifest as skipped.
    fn declared_keys(&self, _header: &BamHeader) -> Vec<Self::Key> {
        Vec::new()
    }

    /// Whether the routing keys of a coordinate-sorted input arrive grouped.
    ///
    /// Only a router whose key is a function of the reference id can promise
    /// this, and only that promise makes the streaming engine usable.
    fn is_grouped_by_coordinate(&self) -> bool {
        false
    }
}

/// What to do with a record whose routing value is absent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MissingPolicy {
    /// Send it to a dedicated output.
    #[default]
    File,
    /// Discard it.
    Drop,
    /// Fail the run.
    Error,
}

impl MissingPolicy {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::File => "file",
            Self::Drop => "drop",
            Self::Error => "error",
        }
    }
}

impl std::fmt::Display for MissingPolicy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Parses and validates a two-character BAM tag.
///
/// # Errors
///
/// Returns [`crate::error::ConfigError::InvalidTag`] unless the value is
/// exactly `[A-Za-z][A-Za-z0-9]`, which is what the SAM specification allows.
pub fn parse_tag(value: &str) -> Result<[u8; 2], crate::error::ConfigError> {
    let bytes = value.as_bytes();
    let invalid = || crate::error::ConfigError::InvalidTag {
        tag: value.to_string(),
    };
    if bytes.len() != 2 {
        return Err(invalid());
    }
    if !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
        return Err(invalid());
    }
    Ok([bytes[0], bytes[1]])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn byte_keys_derive_everything_from_their_bytes() {
        let key = ByteKey::new(b"chr1/alt".to_vec());
        assert_eq!(key.logical(), b"chr1/alt");
        assert_eq!(key.label(), "chr1/alt");
        assert_eq!(key.encoded(), "chr1%2Falt");
        assert!(key.manifest_fields().is_empty());
    }

    #[test]
    fn byte_keys_sort_lexicographically() {
        let mut keys = [
            ByteKey::new(b"chr2".to_vec()),
            ByteKey::new(b"chr10".to_vec()),
            ByteKey::new(b"chr1".to_vec()),
        ];
        keys.sort();
        let labels: Vec<_> = keys.iter().map(|key| key.label().into_owned()).collect();
        assert_eq!(labels, ["chr1", "chr10", "chr2"]);
    }

    #[test]
    fn routes_report_their_emission_count() {
        let one: Route<ByteKey> = Route::One(ByteKey::new(b"a".to_vec()));
        assert_eq!(one.emission_count(), 1);
        assert_eq!(one.keys().len(), 1);

        let many: Route<ByteKey> = Route::Many(SmallVec::from_vec(vec![
            ByteKey::new(b"a".to_vec()),
            ByteKey::new(b"b".to_vec()),
        ]));
        assert_eq!(many.emission_count(), 2);

        let dropped: Route<ByteKey> = Route::Drop(DropReason::Unplaced);
        assert_eq!(dropped.emission_count(), 0);
        assert!(dropped.keys().is_empty());
    }

    #[test]
    fn drop_reasons_distinguish_unmatched_from_dropped() {
        assert!(DropReason::Unmatched.is_unmatched());
        for reason in [
            DropReason::Unplaced,
            DropReason::ExcludedReference,
            DropReason::MissingTag,
            DropReason::UnknownReadGroup,
            DropReason::MissingReadGroupField,
        ] {
            assert!(!reason.is_unmatched(), "{reason}");
        }
    }

    #[test]
    fn tags_are_validated_against_the_sam_grammar() {
        assert_eq!(parse_tag("RG").expect("valid"), *b"RG");
        assert_eq!(parse_tag("X1").expect("valid"), *b"X1");
        for bad in ["", "R", "RGB", "1R", "R-", "  "] {
            assert!(parse_tag(bad).is_err(), "{bad:?} must be rejected");
        }
    }
}
