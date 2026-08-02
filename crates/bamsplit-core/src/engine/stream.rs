// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The grouped streaming engine — the flagship `bamsplit chrom` path.
//!
//! # The pipeline
//!
//! ```text
//! BGZF input ─▶ raw record ─▶ route ─▶ key changed? ─▶ finalize previous
//!                                          │                    │
//!                                          ▼                    ▼
//!                                    raw record write ◀── open next output
//!                                          │
//!                                          ▼
//!                              streaming index + statistics
//! ```
//!
//! One sequential pass. One active output. No input index. No re-decompression.
//! Memory is a handful of buffers plus the index of the *one* chromosome being
//! written, so a 3 000-contig genome costs the same as a two-contig one.
//!
//! # Why grouping is checked, not assumed
//!
//! `@HD SO:coordinate` is a claim, not a guarantee — plenty of BAMs in the wild
//! carry it after a filtering step that reordered records. If this engine
//! trusted the header and the input turned out to be interleaved, it would
//! finalize `chr1.bam`, later meet another `chr1` record, and either overwrite
//! the finished file or write a second partial one. Both are silent data loss.
//!
//! So the engine tracks the keys it has finalized, and a reappearance is
//! [`EngineError::UngroupedInput`]. For `chrom` it additionally verifies
//! monotonic `(refID, pos)` and reports [`EngineError::SortOrderViolation`].
//! Either way every partial output is removed, and the caller — the planner
//! under `--engine auto` — retries with the spool engine when the input is
//! seekable.

use crate::bam::validation::CoordinateOrderTracker;
use crate::engine::{
    EngineKind, ExecutionContext, ExecutionReport, INTERRUPT_CHECK_INTERVAL, SplitEngine,
    check_interrupt, open,
};
use crate::error::EngineError;
use crate::output::OutputManager;
use crate::routing::{Route, Router, RoutingKey};
use crate::stats::RunStats;

/// The grouped streaming engine.
#[derive(Debug, Clone, Copy, Default)]
pub struct StreamEngine {
    /// Whether to verify monotonic `(refID, position)` as well as key grouping.
    ///
    /// Enabled for `chrom`, where the header's sort-order claim is load-bearing.
    /// Disabled for routers whose keys are not coordinate-derived, where
    /// out-of-order coordinates are expected and harmless.
    pub verify_coordinate_order: bool,
    /// Whether to finalize an output the moment its key stops appearing.
    ///
    /// Only a router whose key is a function of `refID` can promise that a
    /// coordinate-sorted input visits keys in contiguous runs. Region routing
    /// cannot — regions overlap and nest — so for it the engine writes through
    /// the output manager's bounded LRU and finalizes everything at the end,
    /// trading the one-open-output property for correctness.
    pub finalize_on_key_change: bool,
}

impl StreamEngine {
    /// An engine that verifies coordinate order and finalizes eagerly.
    #[must_use]
    pub const fn with_order_verification() -> Self {
        Self {
            verify_coordinate_order: true,
            finalize_on_key_change: true,
        }
    }
}

impl<R: Router> SplitEngine<R> for StreamEngine {
    fn execute(
        &self,
        context: &ExecutionContext<'_>,
        router: &R,
    ) -> Result<ExecutionReport, EngineError> {
        let (input_header, mut input, effective_io) = open(
            &context.input,
            context.io,
            context.threads.decompression,
            context.max_record_size,
        )?;

        let mut manager = OutputManager::new(
            &context.output_directory,
            context.header,
            context.manager_options.clone(),
        )
        .map_err(Box::new)?;

        let mut stats = RunStats::default();
        let mut notes = Vec::new();
        let mut order = CoordinateOrderTracker::new();
        // Keys already finalized. A `Vec` is right here: the streaming engine
        // finalizes at most one key per transition, so this is small, and a
        // linear scan over a handful of `Vec<u8>` beats hashing every record.
        let mut finalized: Vec<Vec<u8>> = Vec::new();
        let mut active: Option<R::Key> = None;
        let mut records = 0u64;

        let outcome = (|| -> Result<(), EngineError> {
            while let Some(record) = input.next_record()? {
                records += 1;
                if records.is_multiple_of(INTERRUPT_CHECK_INTERVAL) {
                    check_interrupt(context.interrupt, records)?;
                }

                if self.verify_coordinate_order
                    && !order.observe(&record, records).map_err(Box::new)?
                    && let Some(violation) = order.first_violation()
                {
                    return Err(EngineError::SortOrderViolation {
                        record_number: violation.record_number,
                        reference_id: violation.reference_id,
                        position: violation.position,
                        previous_reference_id: violation.previous_reference_id,
                        previous_position: violation.previous_position,
                    });
                }

                let route = router
                    .route(&input_header, &record)
                    .map_err(|error| Box::new(error.relocate_routing(records)))?;

                match &route {
                    Route::Drop(reason) => {
                        if reason.is_unmatched() {
                            stats.observe_emissions(0);
                        } else {
                            stats.observe_drop();
                        }
                        continue;
                    }
                    Route::One(key) => {
                        transition(&mut manager, &mut active, &mut finalized, key, records)?;
                        manager.write(key.logical(), &record).map_err(Box::new)?;
                    }
                    Route::Many(keys) => {
                        // Fan-out defeats the one-active-output invariant, so
                        // the manager's parking is what bounds descriptors here.
                        for key in keys {
                            if self.finalize_on_key_change
                                && finalized.iter().any(|seen| seen == key.logical())
                            {
                                return Err(EngineError::UngroupedInput {
                                    key: key.label().into_owned(),
                                    record_number: records,
                                });
                            }
                            manager.write(key.logical(), &record).map_err(Box::new)?;
                        }
                        active = None;
                    }
                }
                stats.observe_emissions(route.emission_count() as u64);
            }
            Ok(())
        })();

        if let Err(error) = outcome {
            manager.abort();
            return Err(error);
        }

        // Materialize or record the keys that never saw a record.
        let mut skipped_keys = Vec::new();
        for key in router.declared_keys(&input_header) {
            if manager.contains(key.logical()) {
                continue;
            }
            if context.emit_empty {
                if let Err(error) = manager.create_empty(key.logical()).map_err(Box::new) {
                    manager.abort();
                    return Err(error.into());
                }
            } else {
                skipped_keys.push(key.logical().to_vec());
            }
        }

        let outputs = manager.finish_all().map_err(Box::new)?;
        stats.total_output_emissions = outputs.iter().map(|output| output.stats.record_count).sum();

        if self.verify_coordinate_order && order.is_sorted() {
            notes.push("verified monotonic (reference, position) ordering".to_string());
        }
        notes.push(if self.finalize_on_key_change {
            format!("one sequential pass over {records} records, one active output at a time")
        } else {
            format!(
                "one sequential pass over {records} records, outputs bounded by the descriptor \
                 budget rather than finalized on key change"
            )
        });

        Ok(ExecutionReport {
            outputs,
            stats,
            skipped_keys,
            engine: EngineKind::Stream,
            io: effective_io,
            notes,
        })
    }
}

/// Finalizes the previous output when the routing key changes.
fn transition<K: RoutingKey>(
    manager: &mut OutputManager<'_>,
    active: &mut Option<K>,
    finalized: &mut Vec<Vec<u8>>,
    key: &K,
    record_number: u64,
) -> Result<(), EngineError> {
    if active.as_ref().is_some_and(|current| current == key) {
        return Ok(());
    }
    if finalized.iter().any(|seen| seen == key.logical()) {
        return Err(EngineError::UngroupedInput {
            key: key.label().into_owned(),
            record_number,
        });
    }
    if let Some(previous) = active.take() {
        manager.finish_one(previous.logical()).map_err(Box::new)?;
        finalized.push(previous.logical().to_vec());
    }
    *active = Some(key.clone());
    Ok(())
}

impl crate::error::RoutingError {
    /// Attaches a record ordinal to a routing failure raised without one.
    ///
    /// Routers see only a borrowed record and cannot know where it came from;
    /// the engine can, so it fills the gap here rather than threading the
    /// location through every router signature.
    #[must_use]
    pub fn relocate_routing(self, record_number: u64) -> Self {
        use crate::bam::raw_record::RecordLocation;
        use crate::error::RoutingError as E;

        let located = RecordLocation::new(record_number, 0);
        match self {
            E::Record(inner) => E::Record(Box::new(inner.relocate(located))),
            E::InvalidReferenceId {
                id,
                reference_count,
                location,
            } if location.is_unknown() => E::InvalidReferenceId {
                id,
                reference_count,
                location: located,
            },
            E::MissingReadGroup { location } if location.is_unknown() => {
                E::MissingReadGroup { location: located }
            }
            E::UnknownReadGroup { id, location } if location.is_unknown() => E::UnknownReadGroup {
                id,
                location: located,
            },
            E::MissingTag { tag, location } if location.is_unknown() => E::MissingTag {
                tag,
                location: located,
            },
            E::ArrayValuedTag { tag, location } if location.is_unknown() => E::ArrayValuedTag {
                tag,
                location: located,
            },
            E::UnplacedRecord { location } if location.is_unknown() => {
                E::UnplacedRecord { location: located }
            }
            E::UnmatchedRecord { location } if location.is_unknown() => {
                E::UnmatchedRecord { location: located }
            }
            E::Geometry { source, location } if location.is_unknown() => E::Geometry {
                source,
                location: located,
            },
            other => other,
        }
    }
}
