// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Choosing an engine and an I/O backend, and saying why.
//!
//! # The rules
//!
//! `--engine auto` is not a guess. Every decision below is derived from
//! observable input properties, and the reason string it produces goes into both
//! the log and the manifest, so a surprising choice can be argued with.
//!
//! | command | input | choice |
//! | --- | --- | --- |
//! | `chrom` | coordinate-sorted | **stream** — one pass, one output, no index needed |
//! | `chrom` | coordinate-sorted, indexed, seekable, few large references, threads > 1 | **indexed** |
//! | `chrom` | not coordinate-sorted | **spool** |
//! | `shard` | any | **spool** — hashed keys are interleaved by construction |
//! | `tag` | any | **spool** — tag values are interleaved by construction |
//! | `region` | coordinate-sorted, dense annotation, few regions | **stream** |
//! | `region` | indexed local file, sparse annotation | **indexed** |
//! | `region` | otherwise | **spool** |
//!
//! # Why `indexed` is not the default even when it could run
//!
//! Random access duplicates work at chunk boundaries and defeats read-ahead. It
//! pays off when references are few and large, so the per-reference setup and
//! the duplicated boundary blocks are amortized. With thousands of small
//! contigs, or on network storage, the streaming engine wins — so the planner
//! requires *both* a small reference count and a local file before it will pick
//! `indexed`.
//!
//! An explicit `--engine` is validated, never silently overridden: asking for
//! something impossible is an error, because quietly doing something else is how
//! a benchmark ends up measuring the wrong thing.

use std::path::Path;

use crate::bam::header::{BamHeader, SortOrder};
use crate::engine::{EngineKind, InputSource, IoBackend};
use crate::error::ConfigError;

/// The reference count above which per-reference parallelism stops paying.
///
/// Chosen so a human genome with its primary assembly plus alternate contigs
/// (~200–600 entries) stays on the streaming engine, while a 24-chromosome
/// reference can use the indexed one.
pub const MAX_REFERENCES_FOR_INDEXED: usize = 64;

/// The input size below which parallel decompression is not worth its setup.
pub const MIN_SIZE_FOR_PARALLEL_IO: u64 = 8 * 1024 * 1024;

/// The most regions the streaming engine will hold open for region routing.
pub const MAX_REGIONS_FOR_STREAM: usize = 256;

/// The most regions the indexed engine will plan queries for.
///
/// Each per-reference task bounds its own descriptors through the output
/// manager, so the real limit is the number of index queries and the memory
/// their partially built indexes occupy. A whole-transcriptome annotation is
/// past that, and belongs on the spool engine.
pub const MAX_REGIONS_FOR_INDEXED: usize = 8_192;

/// What the planner observed about the input.
#[derive(Debug, Clone)]
pub struct InputFacts {
    /// The declared sort order.
    pub sort_order: SortOrder,
    /// How many references the header declares.
    pub reference_count: usize,
    /// The input size in bytes, when known.
    pub size: Option<u64>,
    /// Whether the input can be seeked.
    pub seekable: bool,
    /// Whether an index was found next to the input.
    pub indexed: bool,
    /// Whether the input is a local regular file.
    pub local_file: bool,
}

impl InputFacts {
    /// Gathers facts about an input.
    #[must_use]
    pub fn gather(source: &InputSource, header: &BamHeader) -> Self {
        let path = source.path();
        Self {
            sort_order: header.sort_order(),
            reference_count: header.reference_count(),
            size: source.size(),
            seekable: source.is_seekable(),
            indexed: path.is_some_and(|path| crate::engine::indexed::find_index(path).is_some()),
            local_file: path.is_some_and(Path::is_file),
        }
    }
}

/// What the planner decided.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Plan {
    /// The engine to run.
    pub engine: EngineKind,
    /// Why.
    pub engine_reason: String,
    /// The I/O backend to use.
    pub io: IoBackend,
    /// Why.
    pub io_reason: String,
}

/// Which routing mode is being planned for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoutingMode {
    /// One output per reference sequence.
    Chrom,
    /// Hashed shards.
    Shard,
    /// One output per tag value.
    Tag,
    /// One output per annotated region.
    Region {
        /// Whether the annotation covers most of the genome, in which case a
        /// single sequential pass beats many index queries.
        dense: bool,
        /// How many regions the annotation produced.
        ///
        /// The streaming engine holds every live output open at once for region
        /// routing — regions overlap, so keys are not contiguous — and each live
        /// output carries a partially built index. Past a few hundred that stops
        /// being bounded, and the spool engine is the right answer.
        regions: usize,
    },
}

impl RoutingMode {
    /// The manifest label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Chrom => "chrom",
            Self::Shard => "shard",
            Self::Tag => "tag",
            Self::Region { .. } => "region",
        }
    }
}

/// Chooses an engine and backend.
///
/// # Errors
///
/// Returns [`ConfigError::IncompatibleEngine`] or
/// [`ConfigError::IncompatibleIoBackend`] when an explicit request cannot be
/// honoured.
pub fn plan(
    requested_engine: EngineKind,
    requested_io: IoBackend,
    mode: RoutingMode,
    facts: &InputFacts,
    threads: usize,
) -> Result<Plan, ConfigError> {
    let (engine, engine_reason) = match requested_engine {
        EngineKind::Auto => choose_engine(mode, facts, threads),
        explicit => {
            validate_engine(explicit, mode, facts)?;
            (
                explicit,
                format!("`--engine {explicit}` was requested explicitly"),
            )
        }
    };

    let (io, io_reason) = match requested_io {
        IoBackend::Auto => choose_io(engine, facts),
        IoBackend::Buffered => (
            IoBackend::Buffered,
            "`--io buffered` was requested explicitly".to_string(),
        ),
        IoBackend::Mmap => {
            if !facts.local_file {
                return Err(ConfigError::IncompatibleIoBackend {
                    backend: "mmap",
                    reason: "the input is not a local regular file, so it cannot be mapped"
                        .to_string(),
                });
            }
            (
                IoBackend::Mmap,
                "`--io mmap` was requested explicitly".to_string(),
            )
        }
    };

    Ok(Plan {
        engine,
        engine_reason,
        io,
        io_reason,
    })
}

fn choose_engine(mode: RoutingMode, facts: &InputFacts, threads: usize) -> (EngineKind, String) {
    match mode {
        RoutingMode::Chrom => {
            if !facts.sort_order.is_coordinate() {
                return (
                    EngineKind::Spool,
                    format!(
                        "the input declares `SO:{}`, so output keys are interleaved and must be \
                         spooled",
                        facts.sort_order
                    ),
                );
            }
            if facts.indexed
                && facts.seekable
                && facts.local_file
                && threads > 1
                && facts.reference_count <= MAX_REFERENCES_FOR_INDEXED
                && facts
                    .size
                    .is_some_and(|size| size >= MIN_SIZE_FOR_PARALLEL_IO)
            {
                return (
                    EngineKind::Indexed,
                    format!(
                        "the input is a local indexed file with {} reference sequences and \
                         {threads} threads are available, so per-reference queries can run in \
                         parallel",
                        facts.reference_count
                    ),
                );
            }
            (
                EngineKind::Stream,
                format!(
                    "the input is coordinate-sorted, so one sequential pass suffices ({} \
                     reference sequences, index {})",
                    facts.reference_count,
                    if facts.indexed { "present" } else { "absent" }
                ),
            )
        }
        RoutingMode::Shard => (
            EngineKind::Spool,
            "hashed shard keys are interleaved throughout the input by construction".to_string(),
        ),
        RoutingMode::Tag => (
            EngineKind::Spool,
            "tag values are interleaved throughout the input, and their cardinality is not \
             knowable in advance"
                .to_string(),
        ),
        RoutingMode::Region { dense, regions } => {
            if dense && regions <= MAX_REGIONS_FOR_STREAM && facts.sort_order.is_coordinate() {
                (
                    EngineKind::Stream,
                    format!(
                        "the annotation covers most of the genome in only {regions} region(s), so \
                         one sequential pass is cheaper than spooling"
                    ),
                )
            } else if !dense
                && facts.indexed
                && facts.seekable
                && facts.local_file
                && regions <= MAX_REGIONS_FOR_INDEXED
            {
                (
                    EngineKind::Indexed,
                    format!(
                        "the annotation is sparse, so the {regions} region(s) can be queried \
                         directly and the rest of the input never read"
                    ),
                )
            } else {
                (
                    EngineKind::Spool,
                    format!(
                        "{regions} region(s) can be visited in any order and may overlap, so \
                         records are spooled per region"
                    ),
                )
            }
        }
    }
}

fn validate_engine(
    engine: EngineKind,
    mode: RoutingMode,
    facts: &InputFacts,
) -> Result<(), ConfigError> {
    match engine {
        EngineKind::Auto | EngineKind::Spool => Ok(()),
        EngineKind::Stream => {
            // The streaming engine only requires that keys arrive grouped. It
            // detects a violation at runtime and fails loudly rather than
            // guessing here, so the only up-front rejection is a mode whose keys
            // are interleaved by construction.
            match mode {
                RoutingMode::Shard | RoutingMode::Tag => Err(ConfigError::IncompatibleEngine {
                    engine: "stream",
                    reason: format!(
                        "`{}` keys are interleaved throughout the input, so a single grouped pass \
                         cannot finalize any output; use `--engine spool`",
                        mode.as_str()
                    ),
                }),
                RoutingMode::Chrom | RoutingMode::Region { .. } => Ok(()),
            }
        }
        EngineKind::Indexed => {
            if !facts.seekable {
                return Err(ConfigError::IncompatibleEngine {
                    engine: "indexed",
                    reason: "the input is not seekable, so an index cannot be used".to_string(),
                });
            }
            if !facts.indexed {
                return Err(ConfigError::IncompatibleEngine {
                    engine: "indexed",
                    reason: "no `.bai` or `.csi` index was found next to the input; create one \
                             with `samtools index` or use `--engine stream`"
                        .to_string(),
                });
            }
            match mode {
                // `chrom` decomposes by reference; `region` decomposes by the
                // intervals its annotation names, which still belong to one
                // reference each.
                RoutingMode::Chrom | RoutingMode::Region { .. } => Ok(()),
                RoutingMode::Shard | RoutingMode::Tag => Err(ConfigError::IncompatibleEngine {
                    engine: "indexed",
                    reason: format!(
                        "`{}` keys are a function of the record, not of its position, so they \
                         cannot be decomposed into index queries",
                        mode.as_str()
                    ),
                }),
            }
        }
    }
}

fn choose_io(engine: EngineKind, facts: &InputFacts) -> (IoBackend, String) {
    // A mapping is only worth it where cursors are created repeatedly and
    // randomly, which is the indexed engine's access pattern. For a sequential
    // pass, buffered reads with read-ahead are at least as fast and do not risk
    // the mapping failing on an unusual filesystem.
    if engine == EngineKind::Indexed && facts.local_file {
        return (
            IoBackend::Mmap,
            "the indexed engine creates many independent cursors over a local file, which a \
             mapping serves more cheaply than repeated opens"
                .to_string(),
        );
    }
    (
        IoBackend::Buffered,
        "the pass is sequential, so buffered reads with read-ahead are used".to_string(),
    )
}

/// Whether an annotation is dense enough to prefer a single pass.
///
/// "Dense" is deliberately crude — a fraction of the addressable genome — because
/// the decision only needs to separate "a few hundred genes" from "every exon in
/// the genome".
#[must_use]
pub fn annotation_is_dense(covered_bases: u64, header: &BamHeader) -> bool {
    let genome: u64 = header
        .references()
        .iter()
        .map(|reference| u64::from(reference.length))
        .sum();
    genome > 0 && covered_bases * 4 >= genome
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::header::ReferenceSequence;

    fn facts(sort_order: SortOrder, reference_count: usize, indexed: bool) -> InputFacts {
        InputFacts {
            sort_order,
            reference_count,
            size: Some(64 * 1024 * 1024),
            seekable: true,
            indexed,
            local_file: true,
        }
    }

    #[test]
    fn a_coordinate_sorted_chrom_split_streams() {
        let plan = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &facts(SortOrder::Coordinate, 3_000, false),
            8,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Stream);
        assert_eq!(plan.io, IoBackend::Buffered);
        assert!(plan.engine_reason.contains("coordinate-sorted"), "{plan:?}");
    }

    #[test]
    fn an_unsorted_chrom_split_spools() {
        for sort_order in [
            SortOrder::Unsorted,
            SortOrder::QueryName,
            SortOrder::Unknown,
        ] {
            let plan = plan(
                EngineKind::Auto,
                IoBackend::Auto,
                RoutingMode::Chrom,
                &facts(sort_order, 24, false),
                8,
            )
            .expect("planned");
            assert_eq!(plan.engine, EngineKind::Spool, "{sort_order}");
        }
    }

    #[test]
    fn an_indexed_local_input_with_few_references_uses_the_indexed_engine() {
        let plan = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &facts(SortOrder::Coordinate, 24, true),
            16,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Indexed);
        assert_eq!(plan.io, IoBackend::Mmap, "{plan:?}");
    }

    #[test]
    fn many_contigs_keep_the_streaming_engine_even_when_indexed() {
        let plan = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &facts(SortOrder::Coordinate, MAX_REFERENCES_FOR_INDEXED + 1, true),
            16,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Stream);
    }

    #[test]
    fn a_single_thread_keeps_the_streaming_engine() {
        let plan = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &facts(SortOrder::Coordinate, 24, true),
            1,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Stream);
    }

    #[test]
    fn a_small_input_keeps_the_streaming_engine() {
        let mut small = facts(SortOrder::Coordinate, 24, true);
        small.size = Some(1024);
        let plan = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &small,
            16,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Stream);
    }

    #[test]
    fn shard_and_tag_always_spool() {
        for mode in [RoutingMode::Shard, RoutingMode::Tag] {
            let plan = plan(
                EngineKind::Auto,
                IoBackend::Auto,
                mode,
                &facts(SortOrder::Coordinate, 24, true),
                16,
            )
            .expect("planned");
            assert_eq!(plan.engine, EngineKind::Spool, "{}", mode.as_str());
        }
    }

    #[test]
    fn region_planning_follows_density_and_region_count() {
        let dense = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Region {
                dense: true,
                regions: 24,
            },
            &facts(SortOrder::Coordinate, 24, true),
            8,
        )
        .expect("planned");
        assert_eq!(dense.engine, EngineKind::Stream);

        // Dense but with far too many regions to hold open at once.
        let many = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Region {
                dense: true,
                regions: MAX_REGIONS_FOR_STREAM + 1,
            },
            &facts(SortOrder::Coordinate, 24, true),
            8,
        )
        .expect("planned");
        assert_eq!(many.engine, EngineKind::Spool);

        // Sparse and indexed: query the regions directly.
        let sparse = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Region {
                dense: false,
                regions: 24,
            },
            &facts(SortOrder::Coordinate, 24, true),
            8,
        )
        .expect("planned");
        assert_eq!(sparse.engine, EngineKind::Indexed);

        // Sparse but with no index to query.
        let unindexed = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Region {
                dense: false,
                regions: 24,
            },
            &facts(SortOrder::Coordinate, 24, false),
            8,
        )
        .expect("planned");
        assert_eq!(unindexed.engine, EngineKind::Spool);
    }

    #[test]
    fn the_indexed_engine_serves_region_routing() {
        let plan = plan(
            EngineKind::Indexed,
            IoBackend::Auto,
            RoutingMode::Region {
                dense: false,
                regions: 10,
            },
            &facts(SortOrder::Coordinate, 24, true),
            8,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Indexed);
    }

    #[test]
    fn too_many_regions_fall_back_to_the_spool_engine() {
        let plan = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Region {
                dense: false,
                regions: MAX_REGIONS_FOR_INDEXED + 1,
            },
            &facts(SortOrder::Coordinate, 24, true),
            8,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Spool);
    }

    #[test]
    fn an_explicit_engine_is_validated_not_overridden() {
        // Requesting `stream` on an unsorted input is allowed: the engine
        // detects the violation itself and reports it precisely.
        let plan = plan(
            EngineKind::Stream,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &facts(SortOrder::Unsorted, 24, false),
            8,
        )
        .expect("planned");
        assert_eq!(plan.engine, EngineKind::Stream);
        assert!(plan.engine_reason.contains("explicitly"), "{plan:?}");
    }

    #[test]
    fn an_impossible_engine_request_is_an_error() {
        let no_index = plan(
            EngineKind::Indexed,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &facts(SortOrder::Coordinate, 24, false),
            8,
        )
        .expect_err("must reject");
        assert!(
            matches!(no_index, ConfigError::IncompatibleEngine { .. }),
            "{no_index}"
        );

        let not_seekable = plan(
            EngineKind::Indexed,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &InputFacts {
                seekable: false,
                local_file: false,
                ..facts(SortOrder::Coordinate, 24, true)
            },
            8,
        )
        .expect_err("must reject");
        assert!(matches!(
            not_seekable,
            ConfigError::IncompatibleEngine { .. }
        ));

        let wrong_mode = plan(
            EngineKind::Stream,
            IoBackend::Auto,
            RoutingMode::Tag,
            &facts(SortOrder::Coordinate, 24, true),
            8,
        )
        .expect_err("must reject");
        assert!(matches!(wrong_mode, ConfigError::IncompatibleEngine { .. }));

        let indexed_tag = plan(
            EngineKind::Indexed,
            IoBackend::Auto,
            RoutingMode::Tag,
            &facts(SortOrder::Coordinate, 24, true),
            8,
        )
        .expect_err("must reject");
        assert!(matches!(
            indexed_tag,
            ConfigError::IncompatibleEngine { .. }
        ));
    }

    #[test]
    fn mmap_is_refused_for_a_non_file_input() {
        let error = plan(
            EngineKind::Auto,
            IoBackend::Mmap,
            RoutingMode::Chrom,
            &InputFacts {
                local_file: false,
                ..facts(SortOrder::Coordinate, 24, false)
            },
            8,
        )
        .expect_err("must reject");
        assert!(
            matches!(error, ConfigError::IncompatibleIoBackend { .. }),
            "{error}"
        );
    }

    #[test]
    fn density_is_measured_against_the_whole_dictionary() {
        let header = BamHeader::from_parts(
            b"@HD\tVN:1.6\n".to_vec(),
            vec![ReferenceSequence {
                name: b"chr1".to_vec(),
                length: 1_000_000,
            }],
        )
        .expect("valid header");
        assert!(annotation_is_dense(500_000, &header));
        assert!(annotation_is_dense(250_000, &header));
        assert!(!annotation_is_dense(100_000, &header));

        let empty =
            BamHeader::from_parts(b"@HD\tVN:1.6\n".to_vec(), Vec::new()).expect("valid header");
        assert!(!annotation_is_dense(0, &empty));
    }

    #[test]
    fn every_reason_string_is_populated() {
        let plan = plan(
            EngineKind::Auto,
            IoBackend::Auto,
            RoutingMode::Chrom,
            &facts(SortOrder::Coordinate, 24, false),
            1,
        )
        .expect("planned");
        assert!(!plan.engine_reason.is_empty());
        assert!(!plan.io_reason.is_empty());
    }
}
