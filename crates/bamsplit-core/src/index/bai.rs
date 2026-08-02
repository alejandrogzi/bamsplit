// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! BAI specifics: the fixed binning parameters and the coordinate ceiling.
//!
//! BAI is not parameterized. Its R-tree is hard-wired to `min_shift = 14` and
//! `depth = 5`, which means the deepest bin covers 2<sup>14</sup> = 16 384
//! bases and the whole index addresses
//!
//! ```text
//! 2^(14 + 3*5) - 1  =  2^29 - 1  =  536 870 911
//! ```
//!
//! 1-based positions. Any reference longer than that — a few plant and
//! amphibian genomes, and every "one contig per genome" assembly — needs CSI.
//! `bamsplit` decides between the two from the header, before a single record
//! is written, so the choice is deterministic and never has to be revisited
//! mid-stream.

use std::path::Path;

use crate::error::IndexError;

/// The BAI minimum shift: the deepest bin spans 2^14 bases.
pub const MIN_SHIFT: u8 = 14;

/// The BAI R-tree depth.
pub const DEPTH: u8 = 5;

/// The largest 1-based position a BAI can address, `2^29 - 1`.
pub const MAX_POSITION: u64 = (1 << (MIN_SHIFT as u32 + 3 * DEPTH as u32)) - 1;

/// The filename extension, without the dot.
pub const EXTENSION: &str = "bai";

/// The concrete `noodles` index type BAI builds.
pub type Index = noodles_bam::bai::Index;

/// The `noodles` indexer that produces a [`Index`].
pub type Indexer = noodles_csi::binning_index::Indexer<
    noodles_csi::binning_index::index::reference_sequence::index::LinearIndex,
>;

/// Whether every reference in the dictionary fits inside BAI's coordinate
/// space.
#[must_use]
pub fn can_represent(max_reference_length: u64) -> bool {
    max_reference_length <= MAX_POSITION
}

/// Writes an index to `path`.
///
/// # Errors
///
/// Returns [`IndexError::Io`] if the file cannot be created or written.
pub fn write(path: &Path, index: &Index) -> Result<(), IndexError> {
    noodles_bam::bai::fs::write(path, index).map_err(|source| IndexError::Io {
        path: path.to_path_buf(),
        source,
    })
}

/// Reads an index from `path`.
///
/// # Errors
///
/// Returns [`IndexError::ReadInput`] if the file cannot be read or parsed.
pub fn read(path: &Path) -> Result<Index, IndexError> {
    noodles_bam::bai::fs::read(path).map_err(|source| IndexError::ReadInput {
        path: path.to_path_buf(),
        source,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_coordinate_ceiling_matches_the_specification() {
        assert_eq!(MAX_POSITION, 536_870_911);
        assert_eq!(MAX_POSITION, (1u64 << 29) - 1);
    }

    #[test]
    fn representability_is_decided_by_the_longest_reference() {
        assert!(can_represent(1));
        assert!(can_represent(248_956_422)); // human chr1
        assert!(can_represent(MAX_POSITION));
        assert!(!can_represent(MAX_POSITION + 1));
        assert!(!can_represent(1_000_000_000)); // e.g. Ambystoma chromosomes
    }
}
