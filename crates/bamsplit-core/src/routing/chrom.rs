// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Routing by reference sequence — the flagship split.
//!
//! # The rules, exactly
//!
//! | record | goes to |
//! | --- | --- |
//! | valid `refID` | that reference |
//! | `refID == -1` | the unplaced output (`--unplaced keep`), nowhere (`drop`), or an error |
//! | `UNMAPPED` flag **and** a valid `refID` | that reference (`--placed-unmapped by-reference`) or the unplaced output |
//! | `refID` outside the dictionary | an input error |
//! | secondary / supplementary | their own `refID`, never their primary's |
//!
//! Mate fields are **never** rewritten. A read whose mate is on another
//! chromosome ends up in a different file from its mate, and its
//! `next_ref_id`/`next_pos` still resolve, because every output keeps the whole
//! reference dictionary.
//!
//! # Why placed-unmapped defaults to `by-reference`
//!
//! A coordinate-sorted BAM parks an unmapped read immediately after its mapped
//! mate, sharing the mate's `refID` and position. Sending those to
//! `unmapped.bam` would split pairs across files for no benefit and would break
//! the grouping the streaming engine depends on. Keeping them with their
//! reference preserves both.

use crate::bam::header::BamHeader;
use crate::bam::raw_record::RawRecord;
use crate::bam::validation::check_reference_ids;
use crate::error::{ConfigError, RoutingError};
use crate::routing::{ByteKey, DropReason, Route, Router};

/// Where a placed-but-unmapped record goes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PlacedUnmapped {
    /// Keep it with its reference. The default.
    #[default]
    ByReference,
    /// Send it to the unplaced output.
    Unmapped,
}

impl PlacedUnmapped {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ByReference => "by-reference",
            Self::Unmapped => "unmapped",
        }
    }
}

/// What to do with a record that has no reference at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UnplacedPolicy {
    /// Write it to the unplaced output. The default.
    #[default]
    Keep,
    /// Discard it.
    Drop,
    /// Fail the run.
    Error,
}

impl UnplacedPolicy {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Keep => "keep",
            Self::Drop => "drop",
            Self::Error => "error",
        }
    }
}

/// Which references a run selects, and what to call the unplaced output.
#[derive(Debug, Clone)]
pub struct ChromRouterOptions {
    /// Where placed-unmapped records go.
    pub placed_unmapped: PlacedUnmapped,
    /// What to do with unplaced records.
    pub unplaced: UnplacedPolicy,
    /// The stem of the unplaced output.
    pub unmapped_name: Vec<u8>,
    /// Only these references, if non-empty.
    pub include: Vec<Vec<u8>>,
    /// Never these references.
    pub exclude: Vec<Vec<u8>>,
    /// Whether an unknown name in `include`/`exclude` is tolerated.
    pub ignore_missing_references: bool,
}

impl Default for ChromRouterOptions {
    fn default() -> Self {
        Self {
            placed_unmapped: PlacedUnmapped::ByReference,
            unplaced: UnplacedPolicy::Keep,
            unmapped_name: b"unmapped".to_vec(),
            include: Vec::new(),
            exclude: Vec::new(),
            ignore_missing_references: false,
        }
    }
}

/// Routes records to one output per reference sequence.
#[derive(Debug)]
pub struct ChromRouter {
    options: ChromRouterOptions,
    /// Indexed by `refID`: the key for that reference, or [`None`] when the
    /// reference is not selected.
    ///
    /// A flat `Vec` lookup keeps routing to one bounds-checked index per record,
    /// which matters when this runs a billion times.
    selected: Vec<Option<ByteKey>>,
    unmapped_key: ByteKey,
    /// Names that were requested but do not exist, kept for reporting.
    missing_references: Vec<Vec<u8>>,
}

impl ChromRouter {
    /// Builds a router for `header`.
    ///
    /// `--include`, `--exclude`, and `--reference-list` are resolved here
    /// against the **original logical reference names**, never against encoded
    /// filenames, so a user writes `chr1/alt` and not `chr1%2Falt`.
    ///
    /// # Errors
    ///
    /// Returns [`ConfigError::UnknownReferenceNames`] for a requested name that
    /// is not in the dictionary, unless `ignore_missing_references` is set, and
    /// [`ConfigError::EmptySelection`] if nothing is left to write.
    pub fn new(header: &BamHeader, options: ChromRouterOptions) -> Result<Self, ConfigError> {
        let mut missing_references = Vec::new();
        let mut check = |names: &[Vec<u8>], option: &'static str| -> Result<(), ConfigError> {
            let unknown: Vec<Vec<u8>> = names
                .iter()
                .filter(|name| header.reference_id(name).is_none())
                .cloned()
                .collect();
            if unknown.is_empty() {
                return Ok(());
            }
            if !options.ignore_missing_references {
                return Err(ConfigError::UnknownReferenceNames {
                    plural: if unknown.len() == 1 { "" } else { "s" },
                    names: render_names(&unknown),
                    option,
                });
            }
            missing_references.extend(unknown);
            Ok(())
        };
        check(&options.include, "--include")?;
        check(&options.exclude, "--exclude")?;

        let selected: Vec<Option<ByteKey>> = header
            .references()
            .iter()
            .map(|reference| {
                let included =
                    options.include.is_empty() || options.include.contains(&reference.name);
                let excluded = options.exclude.contains(&reference.name);
                (included && !excluded).then(|| ByteKey::new(reference.name.clone()))
            })
            .collect();

        if selected.iter().all(Option::is_none) && !header.references().is_empty() {
            return Err(ConfigError::EmptySelection {
                reason: if options.include.is_empty() {
                    "`--exclude` removed every reference sequence".to_string()
                } else {
                    "`--include` and `--exclude` together select no reference sequence".to_string()
                },
            });
        }

        Ok(Self {
            unmapped_key: ByteKey::new(options.unmapped_name.clone()),
            options,
            selected,
            missing_references,
        })
    }

    /// Names requested via `--include`/`--exclude` that the header does not
    /// declare, tolerated because `--ignore-missing-references` was set.
    #[must_use]
    pub fn missing_references(&self) -> &[Vec<u8>] {
        &self.missing_references
    }

    /// Whether the unplaced output can receive records at all.
    #[must_use]
    pub const fn writes_unplaced(&self) -> bool {
        matches!(self.options.unplaced, UnplacedPolicy::Keep)
            || matches!(self.options.placed_unmapped, PlacedUnmapped::Unmapped)
    }

    /// The key of the unplaced output.
    #[must_use]
    pub const fn unmapped_key(&self) -> &ByteKey {
        &self.unmapped_key
    }
}

fn render_names(names: &[Vec<u8>]) -> String {
    names
        .iter()
        .map(|name| format!("{:?}", String::from_utf8_lossy(name)))
        .collect::<Vec<_>>()
        .join(", ")
}

impl Router for ChromRouter {
    type Key = ByteKey;

    fn route(
        &self,
        header: &BamHeader,
        record: &RawRecord<'_>,
    ) -> Result<Route<Self::Key>, RoutingError> {
        let location = crate::bam::raw_record::RecordLocation::UNKNOWN;
        let reference_id = record.reference_sequence_id().map_err(Box::new)?;

        let Some(reference_id) = reference_id else {
            return Ok(match self.options.unplaced {
                UnplacedPolicy::Keep => Route::One(self.unmapped_key.clone()),
                UnplacedPolicy::Drop => Route::Drop(DropReason::Unplaced),
                UnplacedPolicy::Error => {
                    return Err(RoutingError::UnplacedRecord { location });
                }
            });
        };

        // A `refID` the header does not declare is an input error, not something
        // to route somewhere arbitrary. The mate id is checked too, because a
        // dangling mate reference means the BAM and its header disagree.
        check_reference_ids(header, record, location).map_err(Box::new)?;

        if record.is_unmapped().map_err(Box::new)?
            && self.options.placed_unmapped == PlacedUnmapped::Unmapped
        {
            return Ok(Route::One(self.unmapped_key.clone()));
        }

        let index =
            usize::try_from(reference_id).map_err(|_| RoutingError::InvalidReferenceId {
                id: reference_id,
                reference_count: header.reference_count(),
                location,
            })?;
        match self.selected.get(index) {
            Some(Some(key)) => Ok(Route::One(key.clone())),
            Some(None) => Ok(Route::Drop(DropReason::ExcludedReference)),
            None => Err(RoutingError::InvalidReferenceId {
                id: reference_id,
                reference_count: header.reference_count(),
                location,
            }),
        }
    }

    fn mode(&self) -> &'static str {
        "chrom"
    }

    fn declared_keys(&self, _header: &BamHeader) -> Vec<Self::Key> {
        let mut keys: Vec<ByteKey> = self.selected.iter().flatten().cloned().collect();
        if self.writes_unplaced() {
            keys.push(self.unmapped_key.clone());
        }
        keys
    }

    fn is_grouped_by_coordinate(&self) -> bool {
        // The key is a pure function of `refID`, and a coordinate-sorted BAM
        // orders by `refID` first, so keys arrive in contiguous runs. The one
        // subtlety is the unplaced output, which such a BAM places in a single
        // run at the end — also contiguous.
        true
    }
}

/// Reads a `--reference-list` file: one logical reference name per line.
///
/// Blank lines and `#` comments are ignored. Names are taken verbatim, so a
/// name containing a space is written as-is.
///
/// # Errors
///
/// Returns [`ConfigError::ReferenceList`] if the file cannot be read.
pub fn read_reference_list(path: &std::path::Path) -> Result<Vec<Vec<u8>>, ConfigError> {
    let contents = std::fs::read(path).map_err(|source| ConfigError::ReferenceList {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(contents
        .split(|byte| *byte == b'\n')
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
        .filter(|line| !line.is_empty() && !line.starts_with(b"#"))
        .map(<[u8]>::to_vec)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::header::ReferenceSequence;
    use crate::routing::RoutingKey as _;

    fn header() -> BamHeader {
        BamHeader::from_parts(
            b"@HD\tVN:1.6\tSO:coordinate\n".to_vec(),
            ["chr1", "chr2", "chrX"]
                .iter()
                .map(|name| ReferenceSequence {
                    name: name.as_bytes().to_vec(),
                    length: 1_000_000,
                })
                .collect(),
        )
        .expect("valid header")
    }

    fn record(reference_id: i32, flags: u16) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&reference_id.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&[2, 60]);
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&0u16.to_le_bytes());
        body.extend_from_slice(&flags.to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&(-1i32).to_le_bytes());
        body.extend_from_slice(&0i32.to_le_bytes());
        body.extend_from_slice(b"r\0");
        body
    }

    fn route(router: &ChromRouter, header: &BamHeader, body: &[u8]) -> Route<ByteKey> {
        router
            .route(header, &RawRecord::new(body).expect("valid"))
            .expect("routable")
    }

    fn label(route: &Route<ByteKey>) -> String {
        route.keys()[0].label().into_owned()
    }

    #[test]
    fn a_placed_record_goes_to_its_reference() {
        let header = header();
        let router = ChromRouter::new(&header, ChromRouterOptions::default()).expect("built");
        assert_eq!(label(&route(&router, &header, &record(0, 0))), "chr1");
        assert_eq!(label(&route(&router, &header, &record(2, 0))), "chrX");
    }

    #[test]
    fn secondary_and_supplementary_records_follow_their_own_reference() {
        let header = header();
        let router = ChromRouter::new(&header, ChromRouterOptions::default()).expect("built");
        assert_eq!(label(&route(&router, &header, &record(1, 0x100))), "chr2");
        assert_eq!(label(&route(&router, &header, &record(1, 0x800))), "chr2");
    }

    #[test]
    fn placed_unmapped_defaults_to_staying_with_its_reference() {
        let header = header();
        let router = ChromRouter::new(&header, ChromRouterOptions::default()).expect("built");
        assert_eq!(label(&route(&router, &header, &record(0, 0x4))), "chr1");
    }

    #[test]
    fn placed_unmapped_can_be_diverted() {
        let header = header();
        let router = ChromRouter::new(
            &header,
            ChromRouterOptions {
                placed_unmapped: PlacedUnmapped::Unmapped,
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        assert_eq!(label(&route(&router, &header, &record(0, 0x4))), "unmapped");
    }

    #[test]
    fn unplaced_records_default_to_the_unmapped_output() {
        let header = header();
        let router = ChromRouter::new(&header, ChromRouterOptions::default()).expect("built");
        assert_eq!(
            label(&route(&router, &header, &record(-1, 0x4))),
            "unmapped"
        );
    }

    #[test]
    fn unplaced_records_can_be_dropped_or_rejected() {
        let header = header();
        let dropping = ChromRouter::new(
            &header,
            ChromRouterOptions {
                unplaced: UnplacedPolicy::Drop,
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        assert_eq!(
            route(&dropping, &header, &record(-1, 0x4)),
            Route::Drop(DropReason::Unplaced)
        );

        let failing = ChromRouter::new(
            &header,
            ChromRouterOptions {
                unplaced: UnplacedPolicy::Error,
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        let body = record(-1, 0x4);
        let error = failing
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        assert!(
            matches!(error, RoutingError::UnplacedRecord { .. }),
            "{error}"
        );
    }

    #[test]
    fn the_unmapped_output_can_be_renamed() {
        let header = header();
        let router = ChromRouter::new(
            &header,
            ChromRouterOptions {
                unmapped_name: b"no-coordinate".to_vec(),
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        assert_eq!(
            label(&route(&router, &header, &record(-1, 0x4))),
            "no-coordinate"
        );
    }

    #[test]
    fn an_invalid_reference_id_is_an_input_error() {
        let header = header();
        let router = ChromRouter::new(&header, ChromRouterOptions::default()).expect("built");
        let body = record(9, 0);
        let error = router
            .route(&header, &RawRecord::new(&body).expect("valid"))
            .expect_err("must reject");
        assert!(
            matches!(
                error,
                RoutingError::Record(_) | RoutingError::InvalidReferenceId { .. }
            ),
            "{error}"
        );
    }

    #[test]
    fn include_and_exclude_operate_on_logical_names() {
        let header = header();
        let router = ChromRouter::new(
            &header,
            ChromRouterOptions {
                include: vec![b"chr1".to_vec(), b"chrX".to_vec()],
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        assert_eq!(label(&route(&router, &header, &record(0, 0))), "chr1");
        assert_eq!(
            route(&router, &header, &record(1, 0)),
            Route::Drop(DropReason::ExcludedReference)
        );
        assert_eq!(label(&route(&router, &header, &record(2, 0))), "chrX");
    }

    #[test]
    fn exclude_wins_over_include() {
        let header = header();
        let router = ChromRouter::new(
            &header,
            ChromRouterOptions {
                include: vec![b"chr1".to_vec(), b"chr2".to_vec()],
                exclude: vec![b"chr2".to_vec()],
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        assert_eq!(
            route(&router, &header, &record(1, 0)),
            Route::Drop(DropReason::ExcludedReference)
        );
    }

    #[test]
    fn an_unknown_requested_reference_is_an_error_by_default() {
        let header = header();
        let error = ChromRouter::new(
            &header,
            ChromRouterOptions {
                include: vec![b"chr99".to_vec()],
                ..ChromRouterOptions::default()
            },
        )
        .expect_err("must reject");
        assert!(
            matches!(error, ConfigError::UnknownReferenceNames { .. }),
            "{error}"
        );
    }

    #[test]
    fn unknown_references_can_be_ignored_explicitly() {
        let header = header();
        let router = ChromRouter::new(
            &header,
            ChromRouterOptions {
                include: vec![b"chr1".to_vec(), b"chr99".to_vec()],
                ignore_missing_references: true,
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        assert_eq!(router.missing_references(), [b"chr99".to_vec()]);
        assert_eq!(label(&route(&router, &header, &record(0, 0))), "chr1");
    }

    #[test]
    fn an_empty_selection_is_rejected() {
        let header = header();
        let error = ChromRouter::new(
            &header,
            ChromRouterOptions {
                exclude: vec![b"chr1".to_vec(), b"chr2".to_vec(), b"chrX".to_vec()],
                ..ChromRouterOptions::default()
            },
        )
        .expect_err("must reject");
        assert!(
            matches!(error, ConfigError::EmptySelection { .. }),
            "{error}"
        );
    }

    #[test]
    fn declared_keys_cover_every_selected_reference_plus_unplaced() {
        let header = header();
        let router = ChromRouter::new(&header, ChromRouterOptions::default()).expect("built");
        let labels: Vec<String> = router
            .declared_keys(&header)
            .iter()
            .map(|key| key.label().into_owned())
            .collect();
        assert_eq!(labels, ["chr1", "chr2", "chrX", "unmapped"]);

        let no_unplaced = ChromRouter::new(
            &header,
            ChromRouterOptions {
                unplaced: UnplacedPolicy::Drop,
                ..ChromRouterOptions::default()
            },
        )
        .expect("built");
        assert_eq!(no_unplaced.declared_keys(&header).len(), 3);
    }

    #[test]
    fn the_router_promises_coordinate_grouping() {
        let header = header();
        let router = ChromRouter::new(&header, ChromRouterOptions::default()).expect("built");
        assert!(router.is_grouped_by_coordinate());
        assert!(!router.may_duplicate());
        assert_eq!(router.mode(), "chrom");
    }

    #[test]
    fn reference_lists_ignore_blanks_and_comments() {
        let directory = tempfile::tempdir().expect("temp dir");
        let path = directory.path().join("refs.txt");
        std::fs::write(&path, "# comment\nchr1\n\nchr2\r\n").expect("written");
        assert_eq!(
            read_reference_list(&path).expect("read"),
            vec![b"chr1".to_vec(), b"chr2".to_vec()]
        );
    }

    #[test]
    fn a_missing_reference_list_is_a_configuration_error() {
        let error = read_reference_list(std::path::Path::new("/nonexistent/refs.txt"))
            .expect_err("must fail");
        assert!(
            matches!(error, ConfigError::ReferenceList { .. }),
            "{error}"
        );
    }
}
