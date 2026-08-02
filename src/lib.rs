// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Synthetic BAM fixtures for the `bamsplit` test-suite and benchmarks.
//!
//! # Why hand-built fixtures
//!
//! The behaviour `bamsplit` has to get right is mostly *unusual* input:
//! placed-unmapped records, cross-chromosome mates, references named `chr1/alt`,
//! a truncated BGZF stream, an invalid `refID`, a `CG` long-CIGAR record. None of
//! that appears in a normal BAM, and `samtools` will not write most of it because
//! it validates on the way out.
//!
//! So the fixtures are assembled byte by byte, using `bamsplit`'s own BGZF writer
//! and header serializer. That has a second benefit: a fixture that
//! `samtools quickcheck` accepts is independent evidence that the writer is
//! correct, and the differential suite leans on exactly that.
//!
//! # Example
//!
//! ```
//! use bamsplit_fixtures::Fixture;
//!
//! let directory = tempfile::tempdir()?;
//! let path = directory.path().join("sorted.bam");
//! Fixture::coordinate_sorted().write(&path)?;
//! assert!(path.exists());
//! # Ok::<_, Box<dyn std::error::Error>>(())
//! ```

#![warn(missing_docs)]

pub mod annotation;

pub use bamsplit_core as core;

use std::path::Path;

use bamsplit_core::bam::header::{BamHeader, ReferenceSequence};
use bamsplit_core::output::bgzf::BgzfBlockWriter;

/// SAM flag shorthands, so a fixture reads as prose.
pub mod flags {
    pub use bamsplit_core::bam::flags::*;
}

/// One record to synthesize.
///
/// Everything has a sensible default, so a fixture only states what it is
/// actually testing.
#[derive(Debug, Clone)]
pub struct RecordSpec {
    /// The read name.
    pub name: Vec<u8>,
    /// The reference id, or `-1` for unplaced.
    pub reference_id: i32,
    /// The 0-based position, or `-1` for absent.
    pub position: i32,
    /// The SAM flags.
    pub flags: u16,
    /// The mapping quality.
    pub mapping_quality: u8,
    /// CIGAR operations as `(length, op_code)` pairs, with `M = 0`.
    pub cigar: Vec<(u32, u32)>,
    /// How many bases the read has.
    pub sequence_length: usize,
    /// The mate's reference id.
    pub mate_reference_id: i32,
    /// The mate's 0-based position.
    pub mate_position: i32,
    /// The template length.
    pub template_length: i32,
    /// The auxiliary section, already encoded.
    pub data: Vec<u8>,
    /// When set, fill the sequence and quality with a pseudo-random stream
    /// derived from this seed rather than with constants.
    ///
    /// The default fixtures use `=` bases and "quality unavailable", which
    /// compress to almost nothing. That is fine for correctness but useless for
    /// a benchmark: it would measure DEFLATE on a run of zeroes. Seeding the
    /// payload gives compression ratios in the range real sequencing data
    /// produces.
    pub payload_seed: Option<u64>,
}

impl Default for RecordSpec {
    fn default() -> Self {
        Self {
            name: b"read".to_vec(),
            reference_id: 0,
            position: 0,
            flags: 0,
            mapping_quality: 60,
            cigar: vec![(10, 0)],
            sequence_length: 10,
            mate_reference_id: -1,
            mate_position: -1,
            template_length: 0,
            data: Vec::new(),
            payload_seed: None,
        }
    }
}

impl RecordSpec {
    /// A mapped record at `(reference_id, position)`.
    #[must_use]
    pub fn mapped(name: &str, reference_id: i32, position: i32) -> Self {
        Self {
            name: name.as_bytes().to_vec(),
            reference_id,
            position,
            ..Self::default()
        }
    }

    /// An unplaced, unmapped record.
    #[must_use]
    pub fn unplaced(name: &str) -> Self {
        Self {
            name: name.as_bytes().to_vec(),
            reference_id: -1,
            position: -1,
            flags: flags::UNMAPPED,
            mapping_quality: 0,
            cigar: Vec::new(),
            ..Self::default()
        }
    }

    /// A record that carries a reference and position but is flagged unmapped —
    /// the case a coordinate-sorted BAM produces for the unmapped mate of a
    /// mapped read.
    #[must_use]
    pub fn placed_unmapped(name: &str, reference_id: i32, position: i32) -> Self {
        Self {
            flags: flags::UNMAPPED,
            mapping_quality: 0,
            cigar: Vec::new(),
            ..Self::mapped(name, reference_id, position)
        }
    }

    /// Adds flag bits.
    #[must_use]
    pub const fn with_flags(mut self, bits: u16) -> Self {
        self.flags |= bits;
        self
    }

    /// Overrides the CIGAR.
    #[must_use]
    pub fn with_cigar(mut self, cigar: Vec<(u32, u32)>) -> Self {
        self.cigar = cigar;
        self
    }

    /// Points the record at a mate.
    #[must_use]
    pub const fn with_mate(mut self, reference_id: i32, position: i32) -> Self {
        self.mate_reference_id = reference_id;
        self.mate_position = position;
        self
    }

    /// Appends a `Z`-typed auxiliary field.
    #[must_use]
    pub fn with_string_tag(mut self, tag: &[u8; 2], value: &str) -> Self {
        self.data.extend_from_slice(tag);
        self.data.push(b'Z');
        self.data.extend_from_slice(value.as_bytes());
        self.data.push(0);
        self
    }

    /// Appends an `i`-typed auxiliary field.
    #[must_use]
    pub fn with_int_tag(mut self, tag: &[u8; 2], value: i32) -> Self {
        self.data.extend_from_slice(tag);
        self.data.push(b'i');
        self.data.extend_from_slice(&value.to_le_bytes());
        self
    }

    /// Appends a raw, pre-encoded auxiliary field.
    #[must_use]
    pub fn with_raw_tag(mut self, bytes: &[u8]) -> Self {
        self.data.extend_from_slice(bytes);
        self
    }

    /// Fills the sequence and quality with a pseudo-random stream.
    #[must_use]
    pub const fn with_payload_seed(mut self, seed: u64) -> Self {
        self.payload_seed = Some(seed);
        self
    }

    /// Encodes the record body, exactly as BAM stores it.
    #[must_use]
    pub fn encode(&self) -> Vec<u8> {
        let mut body = Vec::with_capacity(64 + self.data.len());
        body.extend_from_slice(&self.reference_id.to_le_bytes());
        body.extend_from_slice(&self.position.to_le_bytes());
        body.push(u8::try_from(self.name.len() + 1).unwrap_or(u8::MAX));
        body.push(self.mapping_quality);
        // `bin` is advisory and recomputed by every reader that cares; a fixture
        // writes the unmapped sentinel so a mismatch can never mask a bug in
        // code that (wrongly) trusted it.
        body.extend_from_slice(&4680u16.to_le_bytes());
        body.extend_from_slice(&u16::try_from(self.cigar.len()).unwrap_or(0).to_le_bytes());
        body.extend_from_slice(&self.flags.to_le_bytes());
        body.extend_from_slice(
            &i32::try_from(self.sequence_length)
                .unwrap_or(0)
                .to_le_bytes(),
        );
        body.extend_from_slice(&self.mate_reference_id.to_le_bytes());
        body.extend_from_slice(&self.mate_position.to_le_bytes());
        body.extend_from_slice(&self.template_length.to_le_bytes());
        body.extend_from_slice(&self.name);
        body.push(0);
        for (length, kind) in &self.cigar {
            body.extend_from_slice(&((length << 4) | kind).to_le_bytes());
        }
        match self.payload_seed {
            None => {
                // `=` in the 4-bit alphabet is 0, so an all-zero sequence is
                // `====…`, which every reader accepts.
                body.extend(std::iter::repeat_n(0u8, self.sequence_length.div_ceil(2)));
                // 0xff means "quality unavailable".
                body.extend(std::iter::repeat_n(0xffu8, self.sequence_length));
            }
            Some(seed) => {
                let mut state = seed | 1;
                let mut next = move || {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    state
                };
                // Four random bases per byte pair, from the 1..=8 range of the
                // 4-bit alphabet (A, C, G, T and friends) rather than `=`.
                for _ in 0..self.sequence_length.div_ceil(2) {
                    let value = next();
                    let high = 1u8 << (value % 4) as u8;
                    let low = 1u8 << ((value >> 8) % 4) as u8;
                    body.push((high << 4) | low);
                }
                // Phred scores clustered high, as a real run produces.
                for _ in 0..self.sequence_length {
                    body.push(25 + (next() % 15) as u8);
                }
            }
        }
        body.extend_from_slice(&self.data);
        body
    }
}

/// How a fixture's BGZF stream should be damaged, for the robustness tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Damage {
    /// Write a well-formed file.
    #[default]
    None,
    /// Omit the 28-byte end-of-file marker.
    MissingEof,
    /// Cut the last 64 bytes, leaving a partial BGZF block.
    TruncatedBlock,
    /// Cut the final record in half.
    TruncatedRecord,
    /// Write a negative `block_size` for the last record.
    NegativeBlockSize,
}

/// A synthetic BAM.
#[derive(Debug, Clone)]
pub struct Fixture {
    /// The reference dictionary.
    pub references: Vec<ReferenceSequence>,
    /// The SAM header text.
    pub header_text: String,
    /// The records, in the order they will be written.
    pub records: Vec<RecordSpec>,
    /// How to damage the output.
    pub damage: Damage,
    /// The DEFLATE level; a low level keeps fixtures quick to build.
    pub compression_level: u32,
}

impl Fixture {
    /// An empty fixture with the given references.
    #[must_use]
    pub fn new(references: &[(&str, u32)], sort_order: &str) -> Self {
        let mut header_text = format!("@HD\tVN:1.6\tSO:{sort_order}\n");
        for (name, length) in references {
            header_text.push_str(&format!("@SQ\tSN:{name}\tLN:{length}\n"));
        }
        Self {
            references: references
                .iter()
                .map(|(name, length)| ReferenceSequence {
                    name: name.as_bytes().to_vec(),
                    length: *length,
                })
                .collect(),
            header_text,
            records: Vec::new(),
            damage: Damage::None,
            compression_level: 1,
        }
    }

    /// Appends raw header lines, such as `@RG` or `@PG`.
    ///
    /// The lines are emitted in SAM's canonical group order when the header is
    /// serialized, so they can be given here in any order.
    #[must_use]
    pub fn with_header_lines(mut self, lines: &[&str]) -> Self {
        for line in lines {
            self.header_text.push_str(line);
            if !line.ends_with('\n') {
                self.header_text.push('\n');
            }
        }
        self
    }

    /// Appends a record.
    #[must_use]
    pub fn with_record(mut self, record: RecordSpec) -> Self {
        self.records.push(record);
        self
    }

    /// Appends several records.
    #[must_use]
    pub fn with_records(mut self, records: impl IntoIterator<Item = RecordSpec>) -> Self {
        self.records.extend(records);
        self
    }

    /// Damages the output.
    #[must_use]
    pub const fn with_damage(mut self, damage: Damage) -> Self {
        self.damage = damage;
        self
    }

    /// The parsed header this fixture will write.
    ///
    /// # Errors
    ///
    /// Returns an error if the fixture's header text is not valid SAM.
    pub fn header(&self) -> Result<BamHeader, Box<dyn std::error::Error>> {
        Ok(BamHeader::from_parts(
            self.header_text.clone().into_bytes(),
            self.references.clone(),
        )?)
    }

    /// Encodes the fixture to bytes.
    ///
    /// # Errors
    ///
    /// Returns an error if the header is invalid or the BGZF stream cannot be
    /// built.
    pub fn to_bytes(&self) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
        let header = self.header()?;
        let mut header_bytes = Vec::new();
        header.write_to(&mut header_bytes)?;

        let mut writer: BgzfBlockWriter<Vec<u8>, ()> =
            BgzfBlockWriter::new(Vec::new(), self.compression_level, None, 1, 0);
        writer.write_raw(&header_bytes, &mut |(), _| {})?;

        let last = self.records.len().saturating_sub(1);
        for (index, record) in self.records.iter().enumerate() {
            let body = record.encode();
            if index == last && self.damage == Damage::NegativeBlockSize {
                // Bypass `write_record`, which computes a correct prefix.
                let mut framed = (-1i32).to_le_bytes().to_vec();
                framed.extend_from_slice(&body);
                writer.write_raw(&framed, &mut |(), _| {})?;
            } else {
                writer.write_record(&body, (), &mut |(), _| {})?;
            }
        }

        let mut bytes = match self.damage {
            Damage::MissingEof | Damage::TruncatedBlock => {
                let (sink, _) = writer.park(&mut |(), _| {})?;
                sink
            }
            Damage::None | Damage::NegativeBlockSize | Damage::TruncatedRecord => {
                writer.finish(&mut |(), _| {})?
            }
        };

        match self.damage {
            Damage::TruncatedBlock => {
                let cut = bytes.len().saturating_sub(64);
                bytes.truncate(cut.max(1));
            }
            Damage::TruncatedRecord => {
                // Cut inside the final BGZF block so the stream inflates but the
                // last record's body runs out mid-way.
                let cut = bytes
                    .len()
                    .saturating_sub(bamsplit_core::output::BGZF_EOF.len() + 12);
                bytes.truncate(cut.max(1));
            }
            Damage::None | Damage::MissingEof | Damage::NegativeBlockSize => {}
        }
        Ok(bytes)
    }

    /// Writes the fixture to `path`.
    ///
    /// # Errors
    ///
    /// Returns an error if encoding or the write fails.
    pub fn write(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::write(path, self.to_bytes()?)?;
        Ok(())
    }

    // ---------------------------------------------------------------------
    // The named fixtures the suites use.
    // ---------------------------------------------------------------------

    /// The workhorse: coordinate-sorted, three references, and one of every
    /// interesting record shape.
    #[must_use]
    pub fn coordinate_sorted() -> Self {
        Self::new(
            &[("chr1", 10_000), ("chr2", 20_000), ("chrX", 5_000)],
            "coordinate",
        )
        .with_header_lines(&[
            "@RG\tID:rg0\tSM:sample0\tLB:lib0\tPU:unit0\tPL:ILLUMINA\tCN:centre",
            "@RG\tID:rg1\tSM:sample1\tLB:lib1",
            // An existing program record, so the `@PG` chain has something to
            // link onto.
            "@PG\tID:aligner\tPN:aligner\tVN:1.0",
        ])
        .with_records([
            // A proper pair on chr1.
            RecordSpec::mapped("pair1", 0, 100)
                .with_flags(flags::SEGMENTED | flags::PROPERLY_SEGMENTED | flags::FIRST_SEGMENT)
                .with_mate(0, 200)
                .with_string_tag(b"RG", "rg0")
                .with_int_tag(b"NM", 0),
            RecordSpec::mapped("pair1", 0, 200)
                .with_flags(flags::SEGMENTED | flags::PROPERLY_SEGMENTED | flags::LAST_SEGMENT)
                .with_mate(0, 100)
                .with_string_tag(b"RG", "rg0"),
            // A spliced RNA-seq alignment: 5M 100N 5M.
            RecordSpec::mapped("spliced", 0, 300)
                .with_cigar(vec![(5, 0), (100, 3), (5, 0)])
                .with_string_tag(b"RG", "rg1"),
            // No read group at all.
            RecordSpec::mapped("nogroup", 0, 400),
            // Placed but unmapped, sharing its mate's coordinate.
            RecordSpec::placed_unmapped("halfmapped", 0, 400)
                .with_flags(flags::SEGMENTED)
                .with_string_tag(b"RG", "rg0"),
            // A cross-chromosome mate: chr2 pointing at chrX.
            RecordSpec::mapped("crossmate", 1, 100)
                .with_flags(flags::SEGMENTED | flags::FIRST_SEGMENT)
                .with_mate(2, 50)
                .with_string_tag(b"RG", "rg0"),
            // Secondary and supplementary alignments.
            RecordSpec::mapped("secondary", 1, 150)
                .with_flags(flags::SECONDARY)
                .with_string_tag(b"RG", "rg1"),
            RecordSpec::mapped("supplementary", 1, 160)
                .with_cigar(vec![(5, 0), (5, 4)])
                .with_flags(flags::SUPPLEMENTARY)
                .with_string_tag(b"RG", "rg0"),
            // Duplicate-flagged and QC-fail records.
            RecordSpec::mapped("duplicate", 2, 50)
                .with_flags(flags::DUPLICATE)
                .with_string_tag(b"RG", "rg1"),
            RecordSpec::mapped("qcfail", 2, 60)
                .with_flags(flags::QC_FAIL)
                .with_string_tag(b"RG", "unknown-group"),
            // Unplaced records, which a coordinate-sorted BAM puts last.
            RecordSpec::unplaced("unplaced1").with_string_tag(b"RG", "rg0"),
            RecordSpec::unplaced("unplaced2").with_flags(flags::SEGMENTED | flags::MATE_UNMAPPED),
        ])
    }

    /// The same records, deliberately out of coordinate order.
    ///
    /// Note that *reversing* a grouped stream leaves it grouped — just in
    /// descending reference order — so this fixture exercises the sort-order
    /// check but not the grouping check. Use [`Fixture::interleaved`] for that.
    #[must_use]
    pub fn unsorted() -> Self {
        let sorted = Self::coordinate_sorted();
        let mut records = sorted.records.clone();
        records.reverse();
        Self {
            header_text: sorted.header_text.replace("SO:coordinate", "SO:unsorted"),
            records,
            ..sorted
        }
    }

    /// References visited round-robin, so no routing key is ever contiguous.
    ///
    /// This is the input that defeats the streaming engine: it would have to
    /// finalize `chr1`, then meet another `chr1` record. The engine detects that
    /// and the planner retries with the spool engine.
    #[must_use]
    pub fn interleaved() -> Self {
        Self::new(
            &[("chr1", 10_000), ("chr2", 20_000), ("chrX", 5_000)],
            "unsorted",
        )
        .with_records((0..24i32).map(|index| {
            RecordSpec::mapped(&format!("read{index:03}"), index % 3, 100 + index * 10)
        }))
    }

    /// Records grouped by name rather than by coordinate.
    #[must_use]
    pub fn query_name_sorted() -> Self {
        let sorted = Self::coordinate_sorted();
        let mut records = sorted.records.clone();
        records.sort_by(|left, right| left.name.cmp(&right.name));
        Self {
            header_text: sorted.header_text.replace("SO:coordinate", "SO:queryname"),
            records,
            ..sorted
        }
    }

    /// A header and no records.
    #[must_use]
    pub fn header_only() -> Self {
        Self::new(&[("chr1", 1_000), ("chr2", 2_000)], "coordinate")
    }

    /// A header with no references at all.
    #[must_use]
    pub fn no_references() -> Self {
        Self::new(&[], "unsorted").with_record(RecordSpec::unplaced("only"))
    }

    /// Reference names that are legal SAM but need filename encoding.
    ///
    /// The SAM grammar for `SN` is
    /// `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`, so a space,
    /// a comma, and a bracket are *not* representable and are covered by the
    /// filename-encoder unit tests instead. Everything here is a name a real
    /// aligner can emit and that would be dangerous to concatenate into a path.
    #[must_use]
    pub fn hostile_reference_names() -> Self {
        Self::new(
            &[
                ("chr1/alternate", 1_000),
                ("HLA-A*01:01", 1_000),
                ("100%", 1_000),
                ("..", 1_000),
                ("CON", 1_000),
                ("pipe|name", 1_000),
            ],
            "coordinate",
        )
        .with_records((0..6).map(|index| RecordSpec::mapped(&format!("r{index}"), index, 10)))
    }

    /// A reference longer than BAI can address, forcing CSI.
    #[must_use]
    pub fn beyond_bai_limit() -> Self {
        Self::new(&[("chrHuge", 600_000_000)], "coordinate").with_record(RecordSpec::mapped(
            "far",
            0,
            550_000_000,
        ))
    }

    /// Every scalar auxiliary type, plus a `B` array.
    #[must_use]
    pub fn every_tag_type() -> Self {
        let mut data = Vec::new();
        data.extend_from_slice(b"AaA*");
        data.extend_from_slice(b"Ccc");
        data.push(0xff);
        data.extend_from_slice(b"CuC");
        data.push(0xff);
        data.extend_from_slice(b"Sis");
        data.extend_from_slice(&(-2i16).to_le_bytes());
        data.extend_from_slice(b"SuS");
        data.extend_from_slice(&65_535u16.to_le_bytes());
        data.extend_from_slice(b"Iii");
        data.extend_from_slice(&(-3i32).to_le_bytes());
        data.extend_from_slice(b"IuI");
        data.extend_from_slice(&4_000_000_000u32.to_le_bytes());
        data.extend_from_slice(b"Flf");
        data.extend_from_slice(&0.5f32.to_le_bytes());
        data.extend_from_slice(b"ZzZhello\x00");
        data.extend_from_slice(b"HxHDEADBEEF\x00");
        data.extend_from_slice(b"BaB");
        data.push(b'i');
        data.extend_from_slice(&2u32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes());
        data.extend_from_slice(&(-1i32).to_le_bytes());

        Self::new(&[("chr1", 1_000)], "coordinate")
            .with_record(RecordSpec::mapped("tagged", 0, 10).with_raw_tag(&data))
    }

    /// A record using the `CG` long-CIGAR convention.
    #[must_use]
    pub fn long_cigar() -> Self {
        // The real CIGAR is 70 000 operations, so it cannot fit `n_cigar_op`.
        let operations = 70_000u32;
        let mut data = b"CGB".to_vec();
        data.push(b'I');
        data.extend_from_slice(&operations.to_le_bytes());
        for _ in 0..operations {
            data.extend_from_slice(&(1u32 << 4).to_le_bytes()); // 1M
        }

        // The placeholder is `<l_seq>S <reference span>N`, and the real CIGAR is
        // `70000M`, so the read is 70 000 bases long and spans 70 000 reference
        // bases. `htslib` validates that the expanded CIGAR's read length equals
        // `l_seq`, so those three numbers have to agree or the fixture is
        // rejected before it can test anything.
        Self::new(&[("chr1", 200_000)], "coordinate").with_record(RecordSpec {
            sequence_length: operations as usize,
            ..RecordSpec::mapped("longcigar", 0, 100)
                .with_cigar(vec![(operations, 4), (operations, 3)])
                .with_raw_tag(&data)
        })
    }

    /// Many small contigs, to exercise the descriptor budget.
    #[must_use]
    pub fn many_contigs(count: usize) -> Self {
        let names: Vec<String> = (0..count)
            .map(|index| format!("contig{index:05}"))
            .collect();
        let references: Vec<(&str, u32)> =
            names.iter().map(|name| (name.as_str(), 1_000u32)).collect();
        Self::new(&references, "coordinate").with_records((0..count).map(|index| {
            RecordSpec::mapped(&format!("r{index}"), i32::try_from(index).unwrap_or(0), 10)
        }))
    }

    /// High tag cardinality, for the `--max-outputs` guard.
    #[must_use]
    pub fn high_tag_cardinality(count: usize) -> Self {
        Self::new(&[("chr1", 1_000_000)], "coordinate").with_records((0..count).map(|index| {
            RecordSpec::mapped(&format!("r{index}"), 0, i32::try_from(index).unwrap_or(0))
                .with_string_tag(b"CB", &format!("barcode-{index:06}"))
        }))
    }

    /// A record whose `refID` is not in the dictionary.
    #[must_use]
    pub fn invalid_reference_id() -> Self {
        Self::new(&[("chr1", 1_000)], "coordinate").with_record(RecordSpec::mapped("bad", 99, 10))
    }

    /// A record whose mate `refID` is not in the dictionary.
    #[must_use]
    pub fn invalid_mate_reference_id() -> Self {
        Self::new(&[("chr1", 1_000)], "coordinate")
            .with_record(RecordSpec::mapped("bad", 0, 10).with_mate(99, 10))
    }

    /// A record with a malformed auxiliary section.
    #[must_use]
    pub fn malformed_tag() -> Self {
        Self::new(&[("chr1", 1_000)], "coordinate")
            .with_record(RecordSpec::mapped("bad", 0, 10).with_raw_tag(b"XXq\x00"))
    }

    /// A large, realistically compressible BAM, for benchmarking.
    ///
    /// `references` chromosomes of `records_per_reference` 150-base reads each,
    /// coordinate-sorted, with a read group and two auxiliary tags. The payload
    /// is pseudo-random so DEFLATE has real work to do; the whole fixture is
    /// deterministic, so two runs benchmark the same bytes.
    #[must_use]
    pub fn benchmark(references: usize, records_per_reference: usize) -> Self {
        let names: Vec<String> = (0..references).map(|index| format!("chr{index}")).collect();
        let entries: Vec<(&str, u32)> = names
            .iter()
            .map(|name| (name.as_str(), 250_000_000u32))
            .collect();

        let mut records = Vec::with_capacity(references * records_per_reference);
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        for reference in 0..references {
            for index in 0..records_per_reference {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                records.push(RecordSpec {
                    sequence_length: 150,
                    ..RecordSpec::mapped(
                        &format!("SIM:{reference}:{index:09}"),
                        i32::try_from(reference).unwrap_or(0),
                        i32::try_from(index.saturating_mul(37) % 240_000_000).unwrap_or(0),
                    )
                    .with_cigar(vec![(150, 0)])
                    .with_string_tag(b"RG", "rg0")
                    .with_int_tag(b"NM", (seed % 4) as i32)
                    .with_payload_seed(seed)
                });
            }
            // Coordinate order within a reference.
            let start = reference * records_per_reference;
            records[start..].sort_by_key(|record| record.position);
        }

        Self::new(&entries, "coordinate")
            .with_header_lines(&["@RG\tID:rg0\tSM:sample0\tLB:lib0\tPL:ILLUMINA"])
            .with_records(records)
    }

    /// Every named fixture, for the generator binary.
    #[must_use]
    pub fn catalogue() -> Vec<(&'static str, Self)> {
        vec![
            ("coordinate-sorted", Self::coordinate_sorted()),
            ("unsorted", Self::unsorted()),
            ("interleaved", Self::interleaved()),
            ("query-name-sorted", Self::query_name_sorted()),
            ("header-only", Self::header_only()),
            ("no-references", Self::no_references()),
            ("hostile-reference-names", Self::hostile_reference_names()),
            ("beyond-bai-limit", Self::beyond_bai_limit()),
            ("every-tag-type", Self::every_tag_type()),
            ("long-cigar", Self::long_cigar()),
            ("many-contigs", Self::many_contigs(200)),
            ("high-tag-cardinality", Self::high_tag_cardinality(500)),
            ("invalid-reference-id", Self::invalid_reference_id()),
            (
                "invalid-mate-reference-id",
                Self::invalid_mate_reference_id(),
            ),
            ("malformed-tag", Self::malformed_tag()),
            (
                "missing-eof",
                Self::coordinate_sorted().with_damage(Damage::MissingEof),
            ),
            (
                "truncated-block",
                Self::coordinate_sorted().with_damage(Damage::TruncatedBlock),
            ),
            (
                "truncated-record",
                Self::coordinate_sorted().with_damage(Damage::TruncatedRecord),
            ),
            (
                "negative-block-size",
                Self::coordinate_sorted().with_damage(Damage::NegativeBlockSize),
            ),
        ]
    }
}

/// Writes every fixture in the catalogue into `directory`.
///
/// # Errors
///
/// Returns an error if any fixture cannot be written.
pub fn write_catalogue(
    directory: &Path,
) -> Result<Vec<std::path::PathBuf>, Box<dyn std::error::Error>> {
    std::fs::create_dir_all(directory)?;
    let mut written = Vec::new();
    for (name, fixture) in Fixture::catalogue() {
        let path = directory.join(format!("{name}.bam"));
        fixture.write(&path)?;
        written.push(path);
    }
    written.extend(annotation::write_catalogue(&directory.join("annotations"))?);
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_fixture_encodes() {
        for (name, fixture) in Fixture::catalogue() {
            let bytes = fixture
                .to_bytes()
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            assert!(!bytes.is_empty(), "{name} produced no bytes");
        }
    }

    #[test]
    fn a_well_formed_fixture_reads_back_completely() {
        let fixture = Fixture::coordinate_sorted();
        let bytes = fixture.to_bytes().expect("encodes");
        let mut reader = noodles_bgzf::io::Reader::new(&bytes[..]);
        let header = BamHeader::read_from(&mut reader).expect("valid header");
        assert_eq!(header.reference_count(), 3);
        assert_eq!(header.read_groups().len(), 2);

        let mut records = bamsplit_core::bam::RawRecordReader::new(reader);
        let mut count = 0;
        while records.read_record().expect("valid record").is_some() {
            count += 1;
        }
        assert_eq!(count, fixture.records.len());
    }

    #[test]
    fn record_bodies_round_trip_through_the_fixture_writer() {
        let fixture = Fixture::every_tag_type();
        let expected: Vec<Vec<u8>> = fixture.records.iter().map(RecordSpec::encode).collect();
        let bytes = fixture.to_bytes().expect("encodes");

        let mut reader = noodles_bgzf::io::Reader::new(&bytes[..]);
        BamHeader::read_from(&mut reader).expect("valid header");
        let mut records = bamsplit_core::bam::RawRecordReader::new(reader);
        let mut actual = Vec::new();
        while let Some(record) = records.read_record().expect("valid record") {
            actual.push(record.raw_bytes().to_vec());
        }
        assert_eq!(actual, expected);
    }

    #[test]
    fn the_missing_eof_fixture_lacks_the_marker() {
        let bytes = Fixture::coordinate_sorted()
            .with_damage(Damage::MissingEof)
            .to_bytes()
            .expect("encodes");
        let marker = bamsplit_core::output::BGZF_EOF;
        assert_ne!(&bytes[bytes.len() - marker.len()..], &marker[..]);
    }

    #[test]
    fn the_damaged_fixtures_are_rejected() {
        for damage in [
            Damage::TruncatedBlock,
            Damage::TruncatedRecord,
            Damage::NegativeBlockSize,
        ] {
            let bytes = Fixture::coordinate_sorted()
                .with_damage(damage)
                .to_bytes()
                .expect("encodes");
            let mut reader = noodles_bgzf::io::Reader::new(&bytes[..]);
            let Ok(_) = BamHeader::read_from(&mut reader) else {
                continue;
            };
            let mut records = bamsplit_core::bam::RawRecordReader::new(reader);
            let mut failed = false;
            loop {
                match records.read_record() {
                    Ok(None) => break,
                    Ok(Some(_)) => {}
                    Err(_) => {
                        failed = true;
                        break;
                    }
                }
            }
            assert!(failed, "{damage:?} should have been rejected");
        }
    }

    #[test]
    fn the_long_cigar_fixture_resolves_through_its_cg_tag() {
        let fixture = Fixture::long_cigar();
        let bytes = fixture.to_bytes().expect("encodes");
        let mut reader = noodles_bgzf::io::Reader::new(&bytes[..]);
        BamHeader::read_from(&mut reader).expect("valid header");
        let mut records = bamsplit_core::bam::RawRecordReader::new(reader);
        let record = records.read_record().expect("valid").expect("one record");
        assert_eq!(record.reference_span().expect("valid").total(), 70_000);
    }

    #[test]
    fn hostile_reference_names_all_encode_safely() {
        let fixture = Fixture::hostile_reference_names();
        let header = fixture.header().expect("valid header");
        for reference in header.references() {
            let encoded = bamsplit_core::output::filename::encode(&reference.name);
            assert!(
                bamsplit_core::output::filename::is_safe_component(&encoded),
                "{encoded}"
            );
        }
    }

    #[test]
    fn the_catalogue_writes_to_disk() {
        let directory = tempfile::tempdir().expect("temp dir");
        let written = write_catalogue(directory.path()).expect("written");
        assert_eq!(
            written.len(),
            Fixture::catalogue().len() + annotation::catalogue().len() + 8,
            "every BAM fixture, every annotation, and the compressed variants"
        );
        for path in written {
            assert!(path.exists(), "{path:?}");
        }
    }
}
