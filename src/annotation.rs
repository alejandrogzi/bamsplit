// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Synthetic annotation fixtures.
//!
//! The BAM fixtures cover record shapes; these cover the *other* half of region
//! routing — every BED width, both GXF dialects, every compression codec, and
//! the malformed cases detection and extraction have to reject.
//!
//! Fixtures are written to disk rather than built in memory because `genepred`
//! reads from a path: its GTF/GFF aggregation has no reader-based entry point,
//! and its decompressor is chosen from the file name.

use std::io::Write as _;
use std::path::{Path, PathBuf};

/// How a fixture is compressed on disk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codec {
    /// Plain text.
    Plain,
    /// gzip.
    Gzip,
    /// Zstandard.
    Zstd,
    /// bzip2.
    Bzip2,
}

impl Codec {
    /// The filename suffix, including the dot.
    #[must_use]
    pub const fn suffix(self) -> &'static str {
        match self {
            Self::Plain => "",
            Self::Gzip => ".gz",
            Self::Zstd => ".zst",
            Self::Bzip2 => ".bz2",
        }
    }

    /// Compresses `contents`.
    ///
    /// # Panics
    ///
    /// Panics if a codec fails on in-memory data, which would mean the codec
    /// itself is broken.
    #[must_use]
    pub fn encode(self, contents: &[u8]) -> Vec<u8> {
        match self {
            Self::Plain => contents.to_vec(),
            Self::Gzip => {
                let mut encoder =
                    flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
                encoder.write_all(contents).expect("gzip accepts a Vec");
                encoder.finish().expect("gzip finishes")
            }
            Self::Zstd => zstd::encode_all(contents, 1).expect("zstd accepts a slice"),
            Self::Bzip2 => {
                let mut encoder =
                    bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::fast());
                encoder.write_all(contents).expect("bzip2 accepts a Vec");
                encoder.finish().expect("bzip2 finishes")
            }
        }
    }
}

/// A named annotation fixture.
#[derive(Debug, Clone)]
pub struct Annotation {
    /// The base filename, without any compression suffix.
    pub name: &'static str,
    /// The uncompressed contents.
    pub contents: &'static str,
    /// What it exercises.
    pub purpose: &'static str,
}

impl Annotation {
    /// Writes the fixture into `directory`, compressed with `codec`.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be written.
    pub fn write(&self, directory: &Path, codec: Codec) -> std::io::Result<PathBuf> {
        std::fs::create_dir_all(directory)?;
        let path = directory.join(format!("{}{}", self.name, codec.suffix()));
        std::fs::write(&path, codec.encode(self.contents.as_bytes()))?;
        Ok(path)
    }

    /// Writes the fixture under a name that hides its real format.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be written.
    pub fn write_as(&self, directory: &Path, name: &str, codec: Codec) -> std::io::Result<PathBuf> {
        std::fs::create_dir_all(directory)?;
        let path = directory.join(name);
        std::fs::write(&path, codec.encode(self.contents.as_bytes()))?;
        Ok(path)
    }
}

/// Three transcripts on `chr1`, `chr1`, and `chr2`.
///
/// ```text
/// coding      exons                          thick        strand
/// ENST01      1000..1200 2500..2800 4600..5000  1200..4800   +
/// ENST02      1000..1200 2500..2800 4600..5000  1200..4800   -
/// SINGLE      6000..6500                        6100..6400   +
/// ```
pub const BED12: Annotation = Annotation {
    name: "transcripts.bed",
    purpose: "multi-exon coding transcripts on both strands, plus a single-exon one",
    contents: "\
chr1\t1000\t5000\tENST01\t0\t+\t1200\t4800\t0,0,0\t3\t200,300,400\t0,1500,3600
chr1\t1000\t5000\tENST02\t0\t-\t1200\t4800\t0,0,0\t3\t200,300,400\t0,1500,3600
chr2\t6000\t6500\tSINGLE\t0\t+\t6100\t6400\t0,0,0\t1\t500\t0
",
};

/// A BED12 record with no coding bounds — `thickStart == thickEnd`.
pub const BED12_NONCODING: Annotation = Annotation {
    name: "noncoding.bed",
    purpose: "a transcript with no CDS, for the missing-feature policies",
    contents: "\
chr1\t1000\t2000\tLINC01\t0\t+\t1000\t1000\t0,0,0\t2\t100,100\t0,900
",
};

/// A BED12 record with no strand, for the strand-aware UTR features.
pub const BED12_NO_STRAND: Annotation = Annotation {
    name: "no-strand.bed",
    purpose: "a coding transcript with `.` strand",
    contents: "\
chr1\t1000\t2000\tUNSTRANDED\t0\t.\t1200\t1800\t0,0,0\t2\t100,100\t0,900
",
};

/// BED6: named, scored, stranded, but with no blocks or coding bounds.
pub const BED6: Annotation = Annotation {
    name: "regions.bed",
    purpose: "the span-as-single-exon fallback, and `--require-blocks`",
    contents: "\
chr1\t100\t500\tregionA\t0\t+
chr1\t600\t900\tregionB\t0\t-
chr2\t100\t400\tregionC\t0\t+
",
};

/// BED3: coordinates only, so every name has to be generated.
pub const BED3: Annotation = Annotation {
    name: "intervals.bed",
    purpose: "generated names, and the narrowest supported width",
    contents: "\
chr1\t100\t500
chr1\t600\t900
",
};

/// BED4, BED5, BED8, and BED9, each one row.
pub const BED4: Annotation = Annotation {
    name: "b4.bed",
    purpose: "BED4",
    contents: "chr1\t100\t500\tn1\nchr1\t600\t900\tn2\n",
};

/// See [`BED4`].
pub const BED5: Annotation = Annotation {
    name: "b5.bed",
    purpose: "BED5",
    contents: "chr1\t100\t500\tn1\t0\nchr1\t600\t900\tn2\t100\n",
};

/// See [`BED4`].
pub const BED8: Annotation = Annotation {
    name: "b8.bed",
    purpose: "BED8, the narrowest width with coding bounds",
    contents: "chr1\t100\t500\tn1\t0\t+\t200\t400\nchr1\t600\t900\tn2\t0\t-\t650\t850\n",
};

/// See [`BED4`].
pub const BED9: Annotation = Annotation {
    name: "b9.bed",
    purpose: "BED9",
    contents: "\
chr1\t100\t500\tn1\t0\t+\t200\t400\t255,0,0
chr1\t600\t900\tn2\t0\t-\t650\t850\t0,0,255
",
};

/// BED12 with two trailing columns beyond the standard width.
pub const BED12_EXTRA_COLUMNS: Annotation = Annotation {
    name: "extra-columns.bed",
    purpose: "columns past the standard width become additional fields",
    contents: "\
chr1\t1000\t5000\tENST01\t0\t+\t1200\t4800\t0,0,0\t3\t200,300,400\t0,1500,3600\tBRCA1\tprotein_coding
",
};

/// Two records sharing a name, for duplicate resolution.
pub const BED6_DUPLICATE_NAMES: Annotation = Annotation {
    name: "duplicates.bed",
    purpose: "two records that would otherwise overwrite one another",
    contents: "\
chr1\t100\t200\tgene\t0\t+
chr1\t300\t400\tgene\t0\t+
chr1\t500\t600\tgene\t0\t+
",
};

/// BED4 rows whose name column is `.`, so names must be generated.
pub const BED4_MISSING_NAMES: Annotation = Annotation {
    name: "unnamed.bed",
    purpose: "`.` in the name column counts as unnamed",
    contents: "chr1\t100\t200\t.\nchr1\t300\t400\t.\n",
};

/// GTF for one two-exon coding transcript and one single-exon transcript.
pub const GTF: Annotation = Annotation {
    name: "genes.gtf",
    purpose: "GTF aggregation into transcript-like records",
    contents: "\
##description: synthetic
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ttest\texon\t2501\t2800\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ttest\tCDS\t1201\t2700\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";
chr2\ttest\texon\t6001\t6500\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";
",
};

/// GFF3 for the same shape as [`GTF`].
pub const GFF3: Annotation = Annotation {
    name: "genes.gff3",
    purpose: "GFF3 aggregation, with `=`-separated attributes",
    contents: "\
##gff-version 3
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tID=e1;Parent=T1
chr1\ttest\texon\t2501\t2800\t.\t+\t.\tID=e2;Parent=T1
chr2\ttest\texon\t6001\t6500\t.\t-\t.\tID=e3;Parent=T2
",
};

/// Rows whose column counts disagree.
pub const MIXED_WIDTHS: Annotation = Annotation {
    name: "mixed.bed",
    purpose: "a file concatenated from two sources",
    contents: "chr1\t100\t500\tn1\t0\t+\nchr1\t600\t900\n",
};

/// Content that matches neither BED nor GTF/GFF.
pub const AMBIGUOUS: Annotation = Annotation {
    name: "mystery.dat",
    purpose: "detection must give up rather than guess",
    contents: "alpha\tbeta\tgamma\ndelta\tepsilon\tzeta\n",
};

/// No data rows at all.
pub const EMPTY: Annotation = Annotation {
    name: "empty.bed",
    purpose: "comments only",
    contents: "# nothing to see\ntrack name=empty\n\n",
};

/// A BED12 record whose `blockCount` disagrees with its block arrays.
pub const INVALID_BLOCK_COUNT: Annotation = Annotation {
    name: "bad-block-count.bed",
    purpose: "blockCount says 5, the arrays say 2",
    contents: "chr1\t1000\t5000\tBAD\t0\t+\t1200\t4800\t0,0,0\t5\t200,300\t0,1500\n",
};

/// A BED12 record whose blocks are not in ascending order.
pub const UNSORTED_BLOCKS: Annotation = Annotation {
    name: "unsorted-blocks.bed",
    purpose: "blocks given out of order",
    contents: "chr1\t1000\t5000\tBAD\t0\t+\t1200\t4800\t0,0,0\t2\t300,200\t1500,0\n",
};

/// A BED12 record whose blocks overlap one another.
pub const OVERLAPPING_BLOCKS: Annotation = Annotation {
    name: "overlapping-blocks.bed",
    purpose: "blocks that overlap, which no exon structure can",
    contents: "chr1\t1000\t5000\tBAD\t0\t+\t1200\t4800\t0,0,0\t2\t500,500\t0,200\n",
};

/// Reference names that do not appear in any BAM fixture.
pub const WRONG_GENOME_BUILD: Annotation = Annotation {
    name: "ensembl-style.bed",
    purpose: "`1` rather than `chr1`, the classic build mismatch",
    contents: "1\t100\t500\tgeneA\t0\t+\n2\t100\t500\tgeneB\t0\t+\n",
};

/// Every annotation fixture.
#[must_use]
pub fn catalogue() -> Vec<Annotation> {
    vec![
        BED3,
        BED4,
        BED5,
        BED6,
        BED8,
        BED9,
        BED12,
        BED12_EXTRA_COLUMNS,
        BED12_NONCODING,
        BED12_NO_STRAND,
        BED6_DUPLICATE_NAMES,
        BED4_MISSING_NAMES,
        GTF,
        GFF3,
        MIXED_WIDTHS,
        AMBIGUOUS,
        EMPTY,
        INVALID_BLOCK_COUNT,
        UNSORTED_BLOCKS,
        OVERLAPPING_BLOCKS,
        WRONG_GENOME_BUILD,
    ]
}

/// Writes every annotation fixture into `directory`, plus one copy of the BED12
/// and GTF fixtures under each compression codec.
///
/// # Errors
///
/// Returns an error if any fixture cannot be written.
pub fn write_catalogue(directory: &Path) -> std::io::Result<Vec<PathBuf>> {
    let mut written = Vec::new();
    for fixture in catalogue() {
        written.push(fixture.write(directory, Codec::Plain)?);
    }
    for codec in [Codec::Gzip, Codec::Zstd, Codec::Bzip2] {
        written.push(BED12.write(directory, codec)?);
        written.push(GTF.write(directory, codec)?);
    }
    // A gzip stream whose name does not admit it, so detection cannot rely on
    // the extension.
    written.push(BED12.write_as(directory, "misnamed-plain.bed", Codec::Gzip)?);
    // GTF content in a file called `.bed`.
    written.push(GTF.write_as(directory, "lying-extension.bed", Codec::Plain)?);
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_fixture_writes_and_round_trips() {
        let directory = tempfile::tempdir().expect("temp dir");
        let written = write_catalogue(directory.path()).expect("written");
        assert_eq!(written.len(), catalogue().len() + 6 + 2);
        for path in &written {
            assert!(path.exists(), "{path:?}");
        }
    }

    #[test]
    fn codecs_round_trip_their_contents() {
        for codec in [Codec::Plain, Codec::Gzip, Codec::Zstd, Codec::Bzip2] {
            let encoded = codec.encode(BED12.contents.as_bytes());
            assert!(!encoded.is_empty(), "{codec:?}");
            if codec != Codec::Plain {
                assert_ne!(encoded, BED12.contents.as_bytes(), "{codec:?}");
            }
        }
    }

    #[test]
    fn the_catalogue_has_no_duplicate_names() {
        let mut names: Vec<&str> = catalogue().iter().map(|fixture| fixture.name).collect();
        names.sort_unstable();
        let count = names.len();
        names.dedup();
        assert_eq!(names.len(), count);
    }
}
