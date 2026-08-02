// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The command-line surface.
//!
//! Everything here is presentation: parsing, validation that `clap` can express,
//! and translation into [`bamsplit_core`] option structs. No BAM logic lives in
//! this crate.
//!
//! # Streams
//!
//! Logs and progress go to **stderr**. Manifests go to files, and requested
//! tabular output goes to **stdout**. That split is what lets
//! `bamsplit inspect --full | jq` work while `--log-level debug` is on.

use std::path::PathBuf;

use bamsplit_core::bam::header::ReadGroupField;
use bamsplit_core::engine::{EngineKind, IoBackend};
use bamsplit_core::index::IndexMode;
use bamsplit_core::manifest::ManifestFormat;
use bamsplit_core::routing::chrom::{PlacedUnmapped, UnplacedPolicy};
use bamsplit_core::routing::MissingPolicy;
use clap::{Args, Parser, Subcommand, ValueEnum};

/// Lossless, high-performance partitioning of BAM files.
#[derive(Debug, Parser)]
#[command(
    name = "bamsplit",
    version,
    about = "Split a BAM by chromosome, shard, tag, or annotated region — losslessly",
    long_about = None,
    propagate_version = true,
    after_help = "\
EXIT CODES
  0  success
  1  general runtime failure
  2  invalid command-line arguments
  3  invalid or malformed input
  4  output conflict
  5  validation failure
  6  interrupted execution

EXAMPLES
  bamsplit chrom sample.bam --out-dir chromosomes --threads 16 --index auto
  bamsplit shard sample.bam --shards 32 --key qname --out-dir shards
  bamsplit tag sample.bam --field sample --out-dir by-sample
  bamsplit region sample.bam --regions genes.gtf.gz --feature cds \\
      --assignment best-overlap --out-dir by-cds
  bamsplit inspect sample.bam --full"
)]
pub struct Cli {
    /// Options that apply to every subcommand.
    #[command(flatten)]
    pub globals: GlobalOptions,

    /// The subcommand to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Options that apply to every subcommand.
#[derive(Debug, Args, Clone)]
pub struct GlobalOptions {
    /// Total thread budget, shared across decompression and compression.
    #[arg(
        long,
        short = 't',
        global = true,
        default_value_t = 1,
        value_name = "INT"
    )]
    pub threads: usize,

    /// DEFLATE level for output BGZF blocks.
    #[arg(long, global = true, default_value_t = 6, value_name = "INT")]
    pub compression_level: u32,

    /// Execution engine.
    #[arg(long, global = true, default_value = "auto", value_name = "ENGINE")]
    pub engine: EngineArg,

    /// How input bytes are read.
    #[arg(long, global = true, default_value = "auto", value_name = "MODE")]
    pub io: IoArg,

    /// Most outputs holding an open file descriptor at once.
    #[arg(long, global = true, default_value_t = 64, value_name = "INT")]
    pub max_open_files: usize,

    /// Where temporary spool files go.
    #[arg(long, global = true, value_name = "PATH")]
    pub temp_dir: Option<PathBuf>,

    /// Keep temporary spool files after the run, for debugging.
    #[arg(long, global = true)]
    pub keep_temp: bool,

    /// Replace existing outputs instead of refusing to.
    #[arg(long, global = true)]
    pub force: bool,

    /// Log verbosity.
    #[arg(long, global = true, default_value = "info", value_name = "LEVEL")]
    pub log_level: LogLevel,

    /// Suppress all logging.
    #[arg(long, short = 'q', global = true)]
    pub quiet: bool,

    /// Do not render a progress indicator.
    #[arg(long, global = true)]
    pub no_progress: bool,

    /// Do not append an `@PG` record to output headers.
    #[arg(long, global = true)]
    pub no_pg: bool,

    /// Which manifests to write.
    #[arg(long, global = true, default_value = "json", value_name = "FORMAT")]
    pub manifest: ManifestArg,

    /// Largest accepted BAM record body, in bytes.
    #[arg(long, global = true, value_name = "BYTES")]
    pub max_record_size: Option<usize>,
}

/// The subcommands.
#[derive(Debug, Subcommand)]
pub enum Command {
    /// One output per reference sequence. The flagship, most optimized split.
    Chrom(ChromArgs),
    /// A fixed number of deterministic, hash-assigned shards.
    Shard(ShardArgs),
    /// One output per auxiliary-tag value or `@RG`-derived field.
    Tag(TagArgs),
    /// One output per annotated region, or per generated window.
    Region(RegionArgs),
    /// Describe a BAM and what a split of it would look like. Writes no BAMs.
    Inspect(InspectArgs),
}

/// `bamsplit chrom`.
#[derive(Debug, Args)]
pub struct ChromArgs {
    /// The input BAM, or `-` for standard input.
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// Where outputs go.
    #[arg(long, default_value = "chromosomes", value_name = "PATH")]
    pub out_dir: PathBuf,

    /// Where a record that is flagged unmapped but has a reference goes.
    #[arg(long, default_value = "by-reference", value_name = "MODE")]
    pub placed_unmapped: PlacedUnmappedArg,

    /// What to do with a record that has no reference at all.
    #[arg(long, default_value = "keep", value_name = "MODE")]
    pub unplaced: UnplacedArg,

    /// Stem of the output holding unplaced records.
    #[arg(long, default_value = "unmapped", value_name = "STRING")]
    pub unmapped_name: String,

    /// Write a valid header-only BAM for references with no records.
    #[arg(long)]
    pub emit_empty: bool,

    /// Only these reference names, by their logical name.
    #[arg(long, value_delimiter = ',', value_name = "NAME,...")]
    pub include: Vec<String>,

    /// Never these reference names, by their logical name.
    #[arg(long, value_delimiter = ',', value_name = "NAME,...")]
    pub exclude: Vec<String>,

    /// A file of reference names, one per line, merged into `--include`.
    #[arg(long, value_name = "FILE")]
    pub reference_list: Option<PathBuf>,

    /// Tolerate requested reference names the header does not declare.
    #[arg(long)]
    pub ignore_missing_references: bool,

    /// Which output index to build.
    #[arg(long, default_value = "auto", value_name = "MODE")]
    pub index: IndexArg,

    /// Output filename template; `{key}` is required, `{index}` is available.
    #[arg(long, value_name = "TEMPLATE")]
    pub filename_template: Option<String>,
}

/// `bamsplit shard`.
#[derive(Debug, Args)]
pub struct ShardArgs {
    /// The input BAM, or `-` for standard input.
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// How many shards to produce.
    #[arg(long, default_value_t = 8, value_name = "INT")]
    pub shards: u32,

    /// What to hash: `qname`, `record`, `tag:XX`, `read-group`, `sample`,
    /// `library`, or `platform-unit`.
    #[arg(long, default_value = "qname", value_name = "KEY")]
    pub key: String,

    /// Hash seed. The default is fixed and documented, so shards are stable.
    #[arg(long, value_name = "U64")]
    pub seed: Option<u64>,

    /// What to do when the key is absent.
    #[arg(long, default_value = "file", value_name = "MODE")]
    pub missing: MissingArg,

    /// Stem of the output holding records with no key.
    #[arg(long, default_value = "no-key", value_name = "STRING")]
    pub missing_name: String,

    /// Zero-padding width for shard numbers.
    #[arg(long, value_name = "INT")]
    pub shard_width: Option<usize>,

    /// Where outputs go.
    #[arg(long, default_value = "shards", value_name = "PATH")]
    pub out_dir: PathBuf,

    /// Which output index to build.
    #[arg(long, default_value = "auto", value_name = "MODE")]
    pub index: IndexArg,
}

/// `bamsplit tag`.
#[derive(Debug, Args)]
pub struct TagArgs {
    /// The input BAM, or `-` for standard input.
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// Route on this auxiliary tag. Mutually exclusive with `--field`.
    #[arg(long, value_name = "XX", conflicts_with = "field")]
    pub tag: Option<String>,

    /// Route on this `@RG`-derived field. Mutually exclusive with `--tag`.
    #[arg(long, value_name = "FIELD")]
    pub field: Option<FieldArg>,

    /// What to do when the value is absent.
    #[arg(long, default_value = "file", value_name = "MODE")]
    pub missing: MissingArg,

    /// Stem of the output holding records with no value.
    #[arg(long, default_value = "no-tag", value_name = "STRING")]
    pub missing_name: String,

    /// What to do when a record's read group is not in the header.
    #[arg(long, default_value = "file", value_name = "MODE")]
    pub unknown_read_group: MissingArg,

    /// Refuse to produce more than this many outputs.
    #[arg(long, default_value_t = 1_000, value_name = "INT")]
    pub max_outputs: usize,

    /// Lift the `--max-outputs` limit entirely.
    #[arg(long)]
    pub allow_high_cardinality: bool,

    /// Allow `B`-array tag values as routing keys.
    #[arg(long)]
    pub allow_array_tags: bool,

    /// Where outputs go.
    #[arg(long, default_value = "by-tag", value_name = "PATH")]
    pub out_dir: PathBuf,

    /// Which output index to build.
    #[arg(long, default_value = "auto", value_name = "MODE")]
    pub index: IndexArg,
}

/// `bamsplit region`.
#[derive(Debug, Args)]
pub struct RegionArgs {
    /// The input BAM, or `-` for standard input.
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// A BED, GTF, or GFF annotation. Mutually exclusive with `--window-size`.
    #[arg(long, value_name = "PATH", conflicts_with = "window_size")]
    pub regions: Option<PathBuf>,

    /// Generate fixed windows instead of reading an annotation.
    #[arg(long, value_name = "SIZE")]
    pub window_size: Option<String>,

    /// BED width override. Applies to BED input only.
    #[arg(long = "type", default_value = "auto", value_name = "WIDTH")]
    pub bed_type: String,

    /// Annotation format override.
    #[arg(long, default_value = "auto", value_name = "FORMAT")]
    pub format: String,

    /// Which genomic segments of each annotation record participate.
    #[arg(long, default_value = "span", value_name = "FEATURE")]
    pub feature: String,

    /// How an alignment is assigned to a region.
    #[arg(long, default_value = "start", value_name = "MODE")]
    pub assignment: String,

    /// An annotation attribute to use as the logical name.
    #[arg(long, value_name = "STRING")]
    pub name_field: Option<String>,

    /// Prefix for generated names.
    #[arg(long, default_value = "region", value_name = "STRING")]
    pub unnamed_prefix: String,

    /// What to do with a record lacking the requested feature information.
    #[arg(long, default_value = "skip", value_name = "MODE")]
    pub missing_feature: String,

    /// Turn the BED span-as-single-exon fallback into an error.
    #[arg(long)]
    pub require_blocks: bool,

    /// Compare annotation segments against aligned blocks or the outer span.
    #[arg(long, default_value = "blocks", value_name = "MODE")]
    pub alignment_geometry: String,

    /// Write a valid header-only BAM for regions with no records.
    #[arg(long)]
    pub emit_empty: bool,

    /// Where outputs go.
    #[arg(long, default_value = "regions", value_name = "PATH")]
    pub out_dir: PathBuf,

    /// Which output index to build.
    #[arg(long, default_value = "auto", value_name = "MODE")]
    pub index: IndexArg,
}

/// `bamsplit inspect`.
#[derive(Debug, Args)]
pub struct InspectArgs {
    /// The input BAM, or `-` for standard input.
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// Scan every record, not just the header.
    #[arg(long)]
    pub full: bool,

    /// Predict the output count for this routing mode.
    #[arg(long, value_name = "MODE")]
    pub by: Option<ByArg>,

    /// Measure this auxiliary tag's cardinality during a full scan.
    #[arg(long, value_name = "XX")]
    pub tag: Option<String>,

    /// An annotation to inspect alongside the BAM.
    #[arg(long, value_name = "PATH")]
    pub regions: Option<PathBuf>,

    /// The feature type to check the annotation for.
    #[arg(long, default_value = "span", value_name = "FEATURE")]
    pub feature: String,

    /// Emit JSON instead of a human-readable report.
    #[arg(long)]
    pub json: bool,
}

macro_rules! value_enum {
    (
        $(#[$meta:meta])*
        $name:ident -> $target:ty {
            $($variant:ident = $value:literal => $mapped:expr),* $(,)?
        }
    ) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
        pub enum $name {
            $(
                /// See the option's documentation on the argument struct.
                #[value(name = $value)]
                $variant,
            )*
        }

        impl From<$name> for $target {
            fn from(value: $name) -> Self {
                match value {
                    $($name::$variant => $mapped,)*
                }
            }
        }
    };
}

value_enum! {
    /// `--engine`.
    EngineArg -> EngineKind {
        Auto = "auto" => EngineKind::Auto,
        Stream = "stream" => EngineKind::Stream,
        Indexed = "indexed" => EngineKind::Indexed,
        Spool = "spool" => EngineKind::Spool,
    }
}

value_enum! {
    /// `--io`.
    IoArg -> IoBackend {
        Auto = "auto" => IoBackend::Auto,
        Buffered = "buffered" => IoBackend::Buffered,
        Mmap = "mmap" => IoBackend::Mmap,
    }
}

value_enum! {
    /// `--index`.
    IndexArg -> IndexMode {
        Auto = "auto" => IndexMode::Auto,
        Bai = "bai" => IndexMode::Bai,
        Csi = "csi" => IndexMode::Csi,
        None = "none" => IndexMode::None,
    }
}

value_enum! {
    /// `--manifest`.
    ManifestArg -> ManifestFormat {
        Json = "json" => ManifestFormat::Json,
        Tsv = "tsv" => ManifestFormat::Tsv,
        Both = "both" => ManifestFormat::Both,
        None = "none" => ManifestFormat::None,
    }
}

value_enum! {
    /// `--placed-unmapped`.
    PlacedUnmappedArg -> PlacedUnmapped {
        ByReference = "by-reference" => PlacedUnmapped::ByReference,
        Unmapped = "unmapped" => PlacedUnmapped::Unmapped,
    }
}

value_enum! {
    /// `--unplaced`.
    UnplacedArg -> UnplacedPolicy {
        Keep = "keep" => UnplacedPolicy::Keep,
        Drop = "drop" => UnplacedPolicy::Drop,
        Error = "error" => UnplacedPolicy::Error,
    }
}

value_enum! {
    /// `--missing` and `--unknown-read-group`.
    MissingArg -> MissingPolicy {
        File = "file" => MissingPolicy::File,
        Drop = "drop" => MissingPolicy::Drop,
        Error = "error" => MissingPolicy::Error,
    }
}

value_enum! {
    /// `--field`.
    FieldArg -> ReadGroupField {
        ReadGroup = "read-group" => ReadGroupField::ReadGroup,
        Sample = "sample" => ReadGroupField::Sample,
        Library = "library" => ReadGroupField::Library,
        Platform = "platform" => ReadGroupField::Platform,
        PlatformUnit = "platform-unit" => ReadGroupField::PlatformUnit,
        SequencingCenter = "sequencing-center" => ReadGroupField::SequencingCenter,
    }
}

/// `--by`, for `inspect`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ByArg {
    /// Predict a per-reference split.
    #[value(name = "chrom")]
    Chrom,
    /// Predict a shard split.
    #[value(name = "shard")]
    Shard,
    /// Predict a per-tag split.
    #[value(name = "tag")]
    Tag,
    /// Predict a per-region split.
    #[value(name = "region")]
    Region,
}

/// `--log-level`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum LogLevel {
    /// Only errors.
    #[value(name = "error")]
    Error,
    /// Errors and warnings.
    #[value(name = "warn")]
    Warn,
    /// The default: progress and decisions.
    #[value(name = "info")]
    Info,
    /// Per-output detail.
    #[value(name = "debug")]
    Debug,
    /// Everything, including per-record decisions.
    #[value(name = "trace")]
    Trace,
}

impl LogLevel {
    /// The `tracing` filter string.
    #[must_use]
    pub const fn as_filter(self) -> &'static str {
        match self {
            Self::Error => "error",
            Self::Warn => "warn",
            Self::Info => "info",
            Self::Debug => "debug",
            Self::Trace => "trace",
        }
    }
}

impl GlobalOptions {
    /// Translates the global flags into core options.
    ///
    /// `command_line` is recorded verbatim in `@PG CL` and the manifest, so a
    /// finished BAM says exactly how it was produced.
    #[must_use]
    pub fn to_run_options(
        &self,
        command_line: String,
        index: IndexArg,
        filename_template: Option<String>,
    ) -> bamsplit_core::RunOptions {
        bamsplit_core::RunOptions {
            threads: self.threads,
            compression_level: self.compression_level,
            engine: self.engine.into(),
            io: self.io.into(),
            max_open_files: self.max_open_files,
            temp_dir: self.temp_dir.clone(),
            keep_temp: self.keep_temp,
            force: self.force,
            no_pg: self.no_pg,
            manifest: self.manifest.into(),
            index: index.into(),
            command_line,
            filename_template,
            max_record_size: self
                .max_record_size
                .unwrap_or(bamsplit_core::bam::raw_record::DEFAULT_MAX_RECORD_SIZE),
        }
    }
}

/// Renders the invocation as a single line, for `@PG CL`.
#[must_use]
pub fn command_line() -> String {
    std::env::args()
        .map(|argument| {
            if argument.contains(' ') || argument.is_empty() {
                format!("{argument:?}")
            } else {
                argument
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn the_command_tree_is_well_formed() {
        Cli::command().debug_assert();
    }

    #[test]
    fn chrom_parses_its_documented_invocation() {
        let cli = Cli::try_parse_from([
            "bamsplit",
            "chrom",
            "input.bam",
            "--out-dir",
            "chromosomes",
            "--threads",
            "16",
            "--index",
            "auto",
        ])
        .expect("parses");
        assert_eq!(cli.globals.threads, 16);
        let Command::Chrom(args) = cli.command else {
            panic!("expected the chrom subcommand");
        };
        assert_eq!(args.out_dir, PathBuf::from("chromosomes"));
        assert_eq!(args.index, IndexArg::Auto);
        assert_eq!(args.unmapped_name, "unmapped");
        assert_eq!(args.placed_unmapped, PlacedUnmappedArg::ByReference);
        assert_eq!(args.unplaced, UnplacedArg::Keep);
    }

    #[test]
    fn include_and_exclude_accept_comma_separated_names() {
        let cli = Cli::try_parse_from([
            "bamsplit",
            "chrom",
            "in.bam",
            "--include",
            "chr1,chr2",
            "--exclude",
            "chrM",
        ])
        .expect("parses");
        let Command::Chrom(args) = cli.command else {
            panic!("expected chrom");
        };
        assert_eq!(args.include, ["chr1", "chr2"]);
        assert_eq!(args.exclude, ["chrM"]);
    }

    #[test]
    fn shard_parses_its_documented_invocation() {
        let cli = Cli::try_parse_from([
            "bamsplit",
            "shard",
            "input.bam",
            "--shards",
            "32",
            "--key",
            "qname",
            "--out-dir",
            "shards",
        ])
        .expect("parses");
        let Command::Shard(args) = cli.command else {
            panic!("expected shard");
        };
        assert_eq!(args.shards, 32);
        assert_eq!(args.key, "qname");
    }

    #[test]
    fn tag_and_field_are_mutually_exclusive() {
        Cli::try_parse_from(["bamsplit", "tag", "in.bam", "--tag", "RG"]).expect("tag alone");
        Cli::try_parse_from(["bamsplit", "tag", "in.bam", "--field", "sample"])
            .expect("field alone");
        assert!(
            Cli::try_parse_from(["bamsplit", "tag", "in.bam", "--tag", "RG", "--field", "sample"])
                .is_err(),
            "both must be refused"
        );
    }

    #[test]
    fn regions_and_window_size_are_mutually_exclusive() {
        assert!(
            Cli::try_parse_from([
                "bamsplit",
                "region",
                "in.bam",
                "--regions",
                "a.bed",
                "--window-size",
                "50M",
            ])
            .is_err(),
            "both must be refused"
        );
    }

    #[test]
    fn every_global_option_is_accepted_after_the_subcommand() {
        let cli = Cli::try_parse_from([
            "bamsplit",
            "chrom",
            "in.bam",
            "--threads",
            "4",
            "--compression-level",
            "1",
            "--engine",
            "spool",
            "--io",
            "buffered",
            "--max-open-files",
            "8",
            "--temp-dir",
            "/tmp/x",
            "--keep-temp",
            "--force",
            "--log-level",
            "debug",
            "--quiet",
            "--no-progress",
            "--no-pg",
            "--manifest",
            "both",
        ])
        .expect("parses");
        let globals = cli.globals;
        assert_eq!(globals.threads, 4);
        assert_eq!(globals.compression_level, 1);
        assert_eq!(globals.engine, EngineArg::Spool);
        assert_eq!(globals.io, IoArg::Buffered);
        assert_eq!(globals.max_open_files, 8);
        assert_eq!(globals.temp_dir, Some(PathBuf::from("/tmp/x")));
        assert!(globals.keep_temp && globals.force && globals.quiet);
        assert!(globals.no_progress && globals.no_pg);
        assert_eq!(globals.log_level, LogLevel::Debug);
        assert_eq!(globals.manifest, ManifestArg::Both);
    }

    #[test]
    fn global_options_translate_to_core_options() {
        let cli =
            Cli::try_parse_from(["bamsplit", "chrom", "in.bam", "--threads", "8"]).expect("parses");
        let run =
            cli.globals
                .to_run_options("bamsplit chrom in.bam".to_string(), IndexArg::Csi, None);
        assert_eq!(run.threads, 8);
        assert_eq!(run.index, IndexMode::Csi);
        assert_eq!(run.command_line, "bamsplit chrom in.bam");
    }

    #[test]
    fn an_unknown_enum_value_is_rejected() {
        assert!(Cli::try_parse_from(["bamsplit", "chrom", "in.bam", "--index", "tbi"]).is_err());
        assert!(Cli::try_parse_from(["bamsplit", "chrom", "in.bam", "--engine", "magic"]).is_err());
        assert!(Cli::try_parse_from(["bamsplit", "nonesuch", "in.bam"]).is_err());
    }

    #[test]
    fn inspect_accepts_its_documented_invocations() {
        for arguments in [
            vec!["bamsplit", "inspect", "in.bam"],
            vec!["bamsplit", "inspect", "in.bam", "--full"],
            vec!["bamsplit", "inspect", "in.bam", "--by", "chrom"],
            vec!["bamsplit", "inspect", "in.bam", "--tag", "RG"],
            vec![
                "bamsplit",
                "inspect",
                "in.bam",
                "--regions",
                "a.gtf",
                "--feature",
                "exon",
            ],
        ] {
            Cli::try_parse_from(arguments.clone())
                .unwrap_or_else(|error| panic!("{arguments:?}: {error}"));
        }
    }

    #[test]
    fn log_levels_map_to_filters() {
        assert_eq!(LogLevel::Error.as_filter(), "error");
        assert_eq!(LogLevel::Trace.as_filter(), "trace");
    }
}
