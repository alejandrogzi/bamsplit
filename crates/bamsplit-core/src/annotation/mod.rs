// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Annotation-driven region routing.
//!
//! ```text
//! annotation file
//!       │
//!       ▼  detect      format (BED/GTF/GFF), BED width, compression
//!       │              — content beats extension
//!       ▼  load        typed genepred readers, Bed3..Bed12 / Gtf / Gff
//!       │              — no parser of our own
//!       ▼  extract     span · exon · intron · cds · utr · 5'utr · 3'utr
//!       │              — one record, one output key
//!       ▼  index       per-reference sorted envelopes + prefix maxima
//!       │
//!       ▼  route       start · midpoint · contained · best-overlap · overlap
//! ```
//!
//! Every stage is independently testable, and none of them touches a BAM: the
//! output of this module is a [`interval_index::IntervalIndex`], which
//! [`crate::routing::region`] then queries per record.

pub mod detect;
pub mod features;
pub mod genepred_reader;
pub mod interval_index;

pub use detect::{AnnotationFormat, BedType, Compression, Detection, detect};
pub use features::{
    ExtractionStats, FeatureOptions, FeatureType, LogicalRegion, MissingFeaturePolicy, RegionKey,
    RegionMetadata, Strand, extract_regions, generate_windows, parse_window_size,
};
pub use genepred_reader::{LoadOptions, LoadedAnnotations, load};
pub use interval_index::IntervalIndex;
