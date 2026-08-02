// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The run manifest: what was produced, and proof that nothing was lost.
//!
//! # Conservation
//!
//! The manifest is not just a listing. It carries the equations that make
//! losslessness checkable rather than merely claimed.
//!
//! For every non-duplicating mode:
//!
//! ```text
//! input_records = unique_emitted_records + dropped_records + unmatched_records
//! ```
//!
//! For `region --assignment overlap`, where a record may legitimately go to
//! several outputs:
//!
//! ```text
//! total_output_emissions = unique_emitted_records + duplicate_emissions
//! ```
//!
//! Both are checked before the manifest is written, and a failure is
//! [`ExitCode::ValidationFailure`](crate::ExitCode::ValidationFailure) — the run
//! is reported as failed even though every byte was written, because a splitter
//! that cannot account for its own records has no business claiming success.
//!
//! # Formats
//!
//! JSON is the default and carries everything. TSV carries one row per output
//! with the columns a workflow manager needs, so a Nextflow channel can be built
//! with `splitCsv(sep: '\t', header: true)` and no JSON parsing.

use std::io::Write;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::ManifestError;
use crate::stats::{OutputStats, RunStats};

/// Which manifests to write.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ManifestFormat {
    /// `bamsplit.manifest.json` only. The default.
    #[default]
    Json,
    /// `bamsplit.manifest.tsv` only.
    Tsv,
    /// Both.
    Both,
    /// Neither.
    None,
}

impl ManifestFormat {
    /// The option value a user writes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Json => "json",
            Self::Tsv => "tsv",
            Self::Both => "both",
            Self::None => "none",
        }
    }

    /// Whether a JSON manifest is wanted.
    #[must_use]
    pub const fn wants_json(self) -> bool {
        matches!(self, Self::Json | Self::Both)
    }

    /// Whether a TSV manifest is wanted.
    #[must_use]
    pub const fn wants_tsv(self) -> bool {
        matches!(self, Self::Tsv | Self::Both)
    }
}

impl std::fmt::Display for ManifestFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The base name of the manifest files.
pub const MANIFEST_STEM: &str = "bamsplit.manifest";

/// Region-specific manifest fields, present only for `bamsplit region`.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct RegionManifest {
    /// The annotation file.
    pub annotation_source: String,
    /// The detected or requested format.
    pub annotation_format: String,
    /// The detected or requested BED width, when the input was BED.
    pub bed_type: Option<String>,
    /// The requested feature type.
    pub feature_type: String,
    /// The requested assignment mode.
    pub assignment_mode: String,
    /// The alignment geometry used for overlap testing.
    pub alignment_geometry: String,
    /// How many annotation records were read.
    pub annotation_records: u64,
    /// How many genomic segments those records produced.
    pub derived_segment_count: u64,
    /// Records lacking the requested feature information.
    pub missing_feature_records: u64,
    /// Names that had to be generated.
    pub generated_names: u64,
    /// Duplicate logical names that were disambiguated.
    pub duplicate_name_resolutions: u64,
    /// Assignments that needed tie-breaking.
    pub ambiguous_assignments: u64,
    /// BED records whose whole span was treated as a single exon.
    pub span_as_exon_fallbacks: u64,
    /// Regions that produced no segments at all.
    pub empty_regions: u64,
}

/// One output's manifest entry.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OutputEntry {
    /// The logical key, rendered lossily for humans.
    pub logical_key: String,
    /// The logical key as hex, so a non-UTF-8 key survives round-tripping.
    pub logical_key_hex: String,
    /// The key after duplicate-name resolution.
    pub resolved_key: String,
    /// The filesystem-safe encoded key.
    pub encoded_key: String,
    /// The BAM path, relative to the output directory.
    pub bam_path: String,
    /// The index path, when one was written.
    pub index_path: Option<String>,
    /// The index format, when one was written.
    pub index_type: Option<String>,
    /// Why no index was written, when one was requested but skipped.
    pub index_note: Option<String>,
    /// Whether this output was listed but not created.
    pub skipped: bool,
    /// The counters and digest.
    #[serde(flatten)]
    pub stats: OutputStats,
    /// Region-specific fields.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub region: Option<OutputRegionEntry>,
}

/// Region-specific per-output fields.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct OutputRegionEntry {
    /// The annotation record's 0-based input ordinal.
    pub annotation_ordinal: u64,
    /// How many segments this region contributed.
    pub derived_segment_count: u64,
    /// Whether the whole span stood in for a missing exon structure.
    pub span_as_exon_fallback: bool,
    /// The original logical name, before duplicate resolution.
    pub original_name: String,
}

/// The complete manifest.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Manifest {
    /// Always `"bamsplit"`.
    pub program: String,
    /// The crate version.
    pub version: String,
    /// The full command line.
    pub command: String,
    /// RFC-3339 start time.
    pub start_time: String,
    /// RFC-3339 end time.
    pub end_time: String,
    /// Wall-clock duration.
    pub elapsed_seconds: f64,
    /// The input, as given.
    pub input_path: String,
    /// The input size in bytes, when known.
    pub input_size: Option<u64>,
    /// A digest of the input header.
    pub input_header_checksum: String,
    /// The input index format, when one was used.
    pub input_index_type: Option<String>,
    /// The input's declared sort order.
    pub input_sort_order: String,
    /// The engine that ran.
    pub selected_engine: String,
    /// Why that engine was chosen.
    pub engine_selection_reason: String,
    /// The I/O backend used.
    pub io_backend: String,
    /// The thread budget.
    pub threads: usize,
    /// The DEFLATE level.
    pub compression_level: u32,
    /// The routing mode, e.g. `"chrom"`.
    pub routing_mode: String,
    /// The requested index mode.
    pub index_mode: String,
    /// Why that index format was chosen.
    pub index_selection_reason: String,
    /// Temporary bytes written.
    pub temporary_bytes: u64,
    /// Records read from the input.
    pub input_records: u64,
    /// Records emitted at least once.
    pub unique_emitted_records: u64,
    /// The sum of every output's record count.
    pub total_output_emissions: u64,
    /// Records dropped by policy.
    pub dropped_records: u64,
    /// Records that matched nothing.
    pub unmatched_records: u64,
    /// Emissions beyond the first for a record.
    pub duplicate_emissions: u64,
    /// The largest number of outputs one record reached.
    pub max_emissions_for_one_record: u64,
    /// Whether the routing mode may duplicate records.
    pub may_duplicate: bool,
    /// Notes worth surfacing to a user.
    pub notes: Vec<String>,
    /// Region-specific fields.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub region: Option<RegionManifest>,
    /// One entry per output, including skipped ones.
    pub outputs: Vec<OutputEntry>,
}

impl Manifest {
    /// Validates the conservation equations.
    ///
    /// # Errors
    ///
    /// Returns [`ManifestError::ConservationFailed`],
    /// [`ManifestError::EmissionConservationFailed`], or
    /// [`ManifestError::PerOutputAccountingFailed`].
    pub fn validate(&self) -> Result<(), ManifestError> {
        let accounted = self.unique_emitted_records + self.dropped_records + self.unmatched_records;
        if accounted != self.input_records {
            return Err(ManifestError::ConservationFailed {
                input_records: self.input_records,
                unique_emitted: self.unique_emitted_records,
                dropped: self.dropped_records,
                unmatched: self.unmatched_records,
                difference: i128::from(self.input_records) - i128::from(accounted),
            });
        }

        let emitted = self.unique_emitted_records + self.duplicate_emissions;
        if emitted != self.total_output_emissions {
            return Err(ManifestError::EmissionConservationFailed {
                total_emissions: self.total_output_emissions,
                unique_emitted: self.unique_emitted_records,
                duplicates: self.duplicate_emissions,
                difference: i128::from(self.total_output_emissions) - i128::from(emitted),
            });
        }

        if !self.may_duplicate && self.duplicate_emissions != 0 {
            return Err(ManifestError::EmissionConservationFailed {
                total_emissions: self.total_output_emissions,
                unique_emitted: self.unique_emitted_records,
                duplicates: self.duplicate_emissions,
                difference: i128::from(self.duplicate_emissions),
            });
        }

        for output in &self.outputs {
            if !output.stats.is_consistent() {
                return Err(ManifestError::PerOutputAccountingFailed {
                    key: output.logical_key.clone(),
                    record_count: output.stats.record_count,
                    mapped: output.stats.mapped_count,
                    placed_unmapped: output.stats.placed_unmapped_count,
                    unplaced_unmapped: output.stats.unplaced_unmapped_count,
                });
            }
        }
        Ok(())
    }

    /// Renders the manifest as pretty-printed JSON.
    ///
    /// # Errors
    ///
    /// Returns [`ManifestError::Serialize`] if serialization fails.
    pub fn to_json(&self) -> Result<String, ManifestError> {
        serde_json::to_string_pretty(self).map_err(|source| ManifestError::Serialize {
            format: "json",
            source: Box::new(source),
        })
    }

    /// Renders the per-output table as TSV, with a header row.
    ///
    /// The columns are chosen so a workflow manager can build a
    /// `[meta, key, bam, index]` tuple straight from a row.
    ///
    /// # Errors
    ///
    /// Returns [`ManifestError::Serialize`] if a row cannot be written.
    pub fn to_tsv(&self) -> Result<String, ManifestError> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());

        let serialize = |source: csv::Error| ManifestError::Serialize {
            format: "tsv",
            source: Box::new(source),
        };

        writer
            .write_record([
                "logical_key",
                "resolved_key",
                "encoded_key",
                "bam",
                "index",
                "index_type",
                "records",
                "mapped",
                "placed_unmapped",
                "unplaced_unmapped",
                "primary",
                "secondary",
                "supplementary",
                "compressed_bytes",
                "coordinate_sorted",
                "raw_record_digest",
                "skipped",
            ])
            .map_err(serialize)?;

        for output in &self.outputs {
            writer
                .write_record([
                    output.logical_key.as_str(),
                    output.resolved_key.as_str(),
                    output.encoded_key.as_str(),
                    output.bam_path.as_str(),
                    output.index_path.as_deref().unwrap_or(""),
                    output.index_type.as_deref().unwrap_or(""),
                    &output.stats.record_count.to_string(),
                    &output.stats.mapped_count.to_string(),
                    &output.stats.placed_unmapped_count.to_string(),
                    &output.stats.unplaced_unmapped_count.to_string(),
                    &output.stats.primary_count.to_string(),
                    &output.stats.secondary_count.to_string(),
                    &output.stats.supplementary_count.to_string(),
                    &output.stats.compressed_bytes.to_string(),
                    &output.stats.coordinate_sorted.to_string(),
                    output.stats.raw_record_digest.as_str(),
                    &output.skipped.to_string(),
                ])
                .map_err(serialize)?;
        }

        let bytes = writer
            .into_inner()
            .map_err(|source| ManifestError::Serialize {
                format: "tsv",
                source: Box::new(source),
            })?;
        String::from_utf8(bytes).map_err(|source| ManifestError::Serialize {
            format: "tsv",
            source: Box::new(source),
        })
    }

    /// Writes the requested manifests into `directory`.
    ///
    /// Returns the paths written. Validation runs first, so a manifest that
    /// fails its conservation check is never written — a bad manifest on disk
    /// would be worse than none.
    ///
    /// # Errors
    ///
    /// Returns [`ManifestError`] if validation, serialization, or the write
    /// fails.
    pub fn write(
        &self,
        directory: &Path,
        format: ManifestFormat,
    ) -> Result<Vec<PathBuf>, ManifestError> {
        self.validate()?;
        let mut written = Vec::new();
        if format.wants_json() {
            let path = directory.join(format!("{MANIFEST_STEM}.json"));
            write_atomically(&path, self.to_json()?.as_bytes())?;
            written.push(path);
        }
        if format.wants_tsv() {
            let path = directory.join(format!("{MANIFEST_STEM}.tsv"));
            write_atomically(&path, self.to_tsv()?.as_bytes())?;
            written.push(path);
        }
        Ok(written)
    }

    /// Folds run-wide counters in.
    pub fn apply_run_stats(&mut self, stats: &RunStats) {
        self.input_records = stats.input_records;
        self.unique_emitted_records = stats.unique_emitted_records;
        self.total_output_emissions = stats.total_output_emissions;
        self.dropped_records = stats.dropped_records;
        self.unmatched_records = stats.unmatched_records;
        self.duplicate_emissions = stats.duplicate_emissions;
        self.max_emissions_for_one_record = stats.max_emissions_for_one_record;
        self.temporary_bytes = stats.temporary_bytes;
    }
}

/// The manifest itself is written transactionally, for the same reason the BAMs
/// are: a half-written manifest read by a workflow manager is worse than a
/// missing one.
fn write_atomically(path: &Path, bytes: &[u8]) -> Result<(), ManifestError> {
    let temporary = crate::output::transaction::temporary_path(path);
    let io = |source| ManifestError::Io {
        path: temporary.clone(),
        source,
    };
    {
        let mut file = std::fs::File::create(&temporary).map_err(io)?;
        file.write_all(bytes).map_err(io)?;
        file.sync_all().map_err(io)?;
    }
    std::fs::rename(&temporary, path).map_err(|source| ManifestError::Io {
        path: path.to_path_buf(),
        source,
    })
}

/// Formats a [`std::time::SystemTime`] as an RFC-3339 UTC timestamp.
///
/// Written by hand rather than pulling in a date-time crate: the manifest needs
/// exactly one format, and a 40-line civil-date conversion is cheaper than a
/// dependency with its own leap-second and time-zone database.
#[must_use]
pub fn format_timestamp(time: std::time::SystemTime) -> String {
    let Ok(elapsed) = time.duration_since(std::time::UNIX_EPOCH) else {
        return "1970-01-01T00:00:00Z".to_string();
    };
    let seconds = elapsed.as_secs();
    let (days, time_of_day) = (seconds / 86_400, seconds % 86_400);
    let (hour, minute, second) = (
        time_of_day / 3_600,
        (time_of_day % 3_600) / 60,
        time_of_day % 60,
    );
    let (year, month, day) = civil_from_days(days as i64);
    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z")
}

/// Converts days since the Unix epoch to a civil date.
///
/// This is Howard Hinnant's `civil_from_days`, which is exact for the whole
/// proleptic Gregorian calendar and needs no tables.
fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let day_of_era = z - era * 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let year = year_of_era + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let shifted_month = (5 * day_of_year + 2) / 153;
    let day = (day_of_year - (153 * shifted_month + 2) / 5 + 1) as u32;
    let month = if shifted_month < 10 {
        shifted_month + 3
    } else {
        shifted_month - 9
    } as u32;
    (if month <= 2 { year + 1 } else { year }, month, day)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn manifest() -> Manifest {
        Manifest {
            program: "bamsplit".to_string(),
            version: "0.1.0".to_string(),
            command: "bamsplit chrom in.bam".to_string(),
            start_time: "2025-01-01T00:00:00Z".to_string(),
            end_time: "2025-01-01T00:00:01Z".to_string(),
            elapsed_seconds: 1.0,
            input_path: "in.bam".to_string(),
            input_size: Some(1024),
            input_header_checksum: "0000000000000000".to_string(),
            input_index_type: None,
            input_sort_order: "coordinate".to_string(),
            selected_engine: "stream".to_string(),
            engine_selection_reason: "coordinate-sorted input".to_string(),
            io_backend: "buffered".to_string(),
            threads: 4,
            compression_level: 6,
            routing_mode: "chrom".to_string(),
            index_mode: "auto".to_string(),
            index_selection_reason: "within BAI limits".to_string(),
            temporary_bytes: 0,
            input_records: 10,
            unique_emitted_records: 8,
            total_output_emissions: 8,
            dropped_records: 1,
            unmatched_records: 1,
            duplicate_emissions: 0,
            max_emissions_for_one_record: 1,
            may_duplicate: false,
            notes: vec!["one pass".to_string()],
            region: None,
            outputs: Vec::new(),
        }
    }

    fn entry(key: &str, records: u64) -> OutputEntry {
        OutputEntry {
            logical_key: key.to_string(),
            logical_key_hex: hex(key.as_bytes()),
            resolved_key: key.to_string(),
            encoded_key: key.to_string(),
            bam_path: format!("{key}.bam"),
            index_path: Some(format!("{key}.bam.bai")),
            index_type: Some("bai".to_string()),
            index_note: None,
            skipped: false,
            stats: OutputStats {
                record_count: records,
                mapped_count: records,
                primary_count: records,
                coordinate_sorted: true,
                raw_record_digest: "0".repeat(16),
                ..OutputStats::new()
            },
            region: None,
        }
    }

    fn hex(bytes: &[u8]) -> String {
        use std::fmt::Write as _;
        bytes.iter().fold(String::new(), |mut out, byte| {
            let _ = write!(out, "{byte:02x}");
            out
        })
    }

    #[test]
    fn a_conserving_manifest_validates() {
        let mut m = manifest();
        m.outputs = vec![entry("chr1", 5), entry("chr2", 3)];
        m.validate().expect("conserves");
    }

    #[test]
    fn a_lost_record_fails_validation() {
        let mut m = manifest();
        m.unique_emitted_records = 7; // one record vanished
        m.total_output_emissions = 7;
        let error = m.validate().expect_err("must fail");
        assert!(
            matches!(error, ManifestError::ConservationFailed { .. }),
            "{error}"
        );
    }

    #[test]
    fn an_unexpected_duplicate_fails_validation() {
        let mut m = manifest();
        m.duplicate_emissions = 2;
        m.total_output_emissions = 10;
        let error = m.validate().expect_err("must fail");
        assert!(
            matches!(error, ManifestError::EmissionConservationFailed { .. }),
            "{error}"
        );
    }

    #[test]
    fn duplicates_are_accepted_in_a_duplicating_mode() {
        let mut m = manifest();
        m.may_duplicate = true;
        m.duplicate_emissions = 2;
        m.total_output_emissions = 10;
        m.validate().expect("overlap mode conserves");
    }

    #[test]
    fn emission_accounting_must_add_up() {
        let mut m = manifest();
        m.may_duplicate = true;
        m.total_output_emissions = 99;
        let error = m.validate().expect_err("must fail");
        assert!(
            matches!(error, ManifestError::EmissionConservationFailed { .. }),
            "{error}"
        );
    }

    #[test]
    fn an_inconsistent_output_fails_validation() {
        let mut m = manifest();
        let mut broken = entry("chr1", 8);
        broken.stats.mapped_count = 99;
        m.outputs = vec![broken];
        let error = m.validate().expect_err("must fail");
        assert!(
            matches!(error, ManifestError::PerOutputAccountingFailed { .. }),
            "{error}"
        );
    }

    #[test]
    fn json_round_trips() {
        let mut m = manifest();
        m.outputs = vec![entry("chr1", 5), entry("chr2", 3)];
        let json = m.to_json().expect("serializes");
        let parsed: Manifest = serde_json::from_str(&json).expect("deserializes");
        assert_eq!(parsed, m);
        assert!(json.contains("\"raw_record_digest\""));
        assert!(json.contains("\"unique_emitted_records\": 8"));
    }

    #[test]
    fn tsv_has_a_header_and_one_row_per_output() {
        let mut m = manifest();
        m.outputs = vec![entry("chr1", 5), entry("chr2", 3)];
        let tsv = m.to_tsv().expect("serializes");
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("logical_key\tresolved_key\tencoded_key\tbam\tindex"));
        assert!(lines[1].starts_with("chr1\tchr1\tchr1\tchr1.bam\tchr1.bam.bai\tbai\t5"));
    }

    #[test]
    fn writing_is_gated_on_validation() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut m = manifest();
        m.unique_emitted_records = 0;
        let error = m
            .write(directory.path(), ManifestFormat::Both)
            .expect_err("must fail");
        assert!(
            matches!(error, ManifestError::ConservationFailed { .. }),
            "{error}"
        );
        assert_eq!(
            std::fs::read_dir(directory.path())
                .expect("listable")
                .count(),
            0,
            "an invalid manifest must not be written"
        );
    }

    #[test]
    fn both_formats_are_written() {
        let directory = tempfile::tempdir().expect("temp dir");
        let mut m = manifest();
        m.outputs = vec![entry("chr1", 8)];
        let written = m
            .write(directory.path(), ManifestFormat::Both)
            .expect("written");
        assert_eq!(written.len(), 2);
        assert!(directory.path().join("bamsplit.manifest.json").exists());
        assert!(directory.path().join("bamsplit.manifest.tsv").exists());

        let none = m
            .write(directory.path(), ManifestFormat::None)
            .expect("written");
        assert!(none.is_empty());
    }

    #[test]
    fn run_stats_are_folded_in() {
        let mut m = manifest();
        m.apply_run_stats(&RunStats {
            input_records: 100,
            unique_emitted_records: 90,
            total_output_emissions: 95,
            dropped_records: 5,
            unmatched_records: 5,
            duplicate_emissions: 5,
            max_emissions_for_one_record: 3,
            temporary_bytes: 42,
            ambiguous_assignments: 1,
        });
        assert_eq!(m.input_records, 100);
        assert_eq!(m.temporary_bytes, 42);
        m.may_duplicate = true;
        m.validate().expect("conserves");
    }

    #[test]
    fn timestamps_are_rfc_3339() {
        assert_eq!(
            format_timestamp(std::time::UNIX_EPOCH),
            "1970-01-01T00:00:00Z"
        );
        assert_eq!(
            format_timestamp(std::time::UNIX_EPOCH + std::time::Duration::from_secs(1_700_000_000)),
            "2023-11-14T22:13:20Z"
        );
        // A leap day, to exercise the civil-date conversion.
        assert_eq!(
            format_timestamp(std::time::UNIX_EPOCH + std::time::Duration::from_secs(1_709_164_800)),
            "2024-02-29T00:00:00Z"
        );
    }

    #[test]
    fn manifest_formats_render_their_option_values() {
        for (format, text, json, tsv) in [
            (ManifestFormat::Json, "json", true, false),
            (ManifestFormat::Tsv, "tsv", false, true),
            (ManifestFormat::Both, "both", true, true),
            (ManifestFormat::None, "none", false, false),
        ] {
            assert_eq!(format.to_string(), text);
            assert_eq!(format.wants_json(), json);
            assert_eq!(format.wants_tsv(), tsv);
        }
    }
}
