// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! CSI specifics: choosing a depth that covers the reference dictionary.
//!
//! Unlike BAI, CSI carries its own `min_shift` and `depth`, so it can address
//! any reference length. `bamsplit` keeps `min_shift = 14` — the same 16 KiB
//! leaf bin BAI and `htslib` use, so bin granularity is familiar — and raises
//! only the depth, to the smallest value that covers the longest reference:
//!
//! ```text
//! addressable(depth) = 2^(14 + 3*depth) - 1
//!
//! depth 5 -> 536 870 911          (BAI parity)
//! depth 6 -> 4 294 967 295        (covers any u32 reference length)
//! ```
//!
//! Depth 6 is therefore always sufficient for a BAM, whose `l_ref` is a signed
//! 32-bit field. The general computation is kept anyway because the ceiling is
//! a property of the format, not of the current caller, and because a
//! needlessly deep index costs bins.

use std::path::Path;

use crate::error::IndexError;

/// The CSI minimum shift `bamsplit` writes.
pub const MIN_SHIFT: u8 = 14;

/// The smallest depth `bamsplit` will use, matching BAI.
pub const MIN_DEPTH: u8 = 5;

/// The largest depth the CSI bin arithmetic supports.
///
/// `noodles` asserts `depth <= 10` when computing the bin limit.
pub const MAX_DEPTH: u8 = 10;

/// The filename extension, without the dot.
pub const EXTENSION: &str = "csi";

/// The concrete `noodles` index type CSI builds.
pub type Index = noodles_csi::Index;

/// The `noodles` indexer that produces a [`Index`].
pub type Indexer = noodles_csi::binning_index::Indexer<
    noodles_csi::binning_index::index::reference_sequence::index::BinnedIndex,
>;

/// The largest 1-based position addressable at `min_shift`/`depth`.
#[must_use]
pub fn addressable(min_shift: u8, depth: u8) -> u64 {
    let bits = u32::from(min_shift) + 3 * u32::from(depth);
    if bits >= 64 {
        return u64::MAX;
    }
    (1u64 << bits) - 1
}

/// The smallest depth that covers `max_reference_length`.
///
/// Returns [`MAX_DEPTH`] when even that is not enough, which cannot happen for
/// a BAM: `l_ref` is 32 bits and depth 6 already covers `2^32 - 1`.
#[must_use]
pub fn depth_for(max_reference_length: u64) -> u8 {
    (MIN_DEPTH..=MAX_DEPTH)
        .find(|depth| addressable(MIN_SHIFT, *depth) >= max_reference_length)
        .unwrap_or(MAX_DEPTH)
}

/// Writes an index to `path`.
///
/// # Errors
///
/// Returns [`IndexError::Io`] if the file cannot be created or written.
pub fn write(path: &Path, index: &Index) -> Result<(), IndexError> {
    noodles_csi::fs::write(path, index).map_err(|source| IndexError::Io {
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
    noodles_csi::fs::read(path).map_err(|source| IndexError::ReadInput {
        path: path.to_path_buf(),
        source,
    })
}

#[cfg(test)]
#[allow(clippy::items_after_statements)]
mod tests {
    use super::*;

    #[test]
    fn addressable_space_grows_by_three_bits_per_level() {
        assert_eq!(addressable(14, 5), 536_870_911);
        assert_eq!(addressable(14, 6), 4_294_967_295);
        assert_eq!(addressable(14, 7), 34_359_738_367);
    }

    #[test]
    fn addressable_saturates_instead_of_overflowing() {
        assert_eq!(addressable(14, 20), u64::MAX);
        assert_eq!(addressable(64, 0), u64::MAX);
    }

    #[test]
    fn depth_is_the_smallest_that_fits() {
        assert_eq!(depth_for(1), MIN_DEPTH);
        assert_eq!(depth_for(248_956_422), 5);
        assert_eq!(depth_for(536_870_911), 5);
        assert_eq!(depth_for(536_870_912), 6);
        assert_eq!(depth_for(u64::from(u32::MAX)), 6);
    }

    #[test]
    fn depth_is_capped_at_the_supported_maximum() {
        assert_eq!(depth_for(u64::MAX), MAX_DEPTH);
        const _: () = assert!(MAX_DEPTH <= 10, "noodles asserts depth <= 10");
    }
}
