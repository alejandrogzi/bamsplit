// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! The high-level entry points: one function per subcommand.
//!
//! Everything the CLI does, a Rust program can do too, with the same
//! guarantees and the same manifest. That is the point of keeping argument
//! parsing out of this crate.
//!
//! Each entry point runs the same seven steps:
//!
//! ```text
//! validate input ─▶ read header ─▶ plan engine + index ─▶ amend header (@PG)
//!                ─▶ build router ─▶ execute ─▶ validate and write the manifest
//! ```

use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime};

use crate::bam::header::{BamHeader, ReadGroupField, SortOrder};
use crate::engine::{
    EngineKind, ExecutionContext, ExecutionReport, InputSource, IoBackend, SplitEngine,
    ThreadBudget,
};
use crate::error::{ConfigError, Error, Result};
use crate::index::{IndexMode, IndexPlan};
use crate::manifest::{Manifest, ManifestFormat, OutputEntry, format_timestamp};
use crate::output::filename::FilenameTemplate;
use crate::output::{OutputManagerOptions, OutputSettings};
use crate::planning::{InputFacts, RoutingMode};
use crate::routing::chrom::{ChromRouter, ChromRouterOptions, PlacedUnmapped, UnplacedPolicy};
use crate::routing::shard::{ShardKeySource, ShardRouter, ShardRouterOptions};
use crate::routing::tag::{TagRouter, TagRouterOptions, TagSource};
use crate::routing::{MissingPolicy, Router};
use crate::stats::RunStats;

/// Options shared by every subcommand — the CLI's global flags.
#[derive(Debug, Clone)]
pub struct RunOptions {
    /// The total thread budget.
    pub threads: usize,
    /// The DEFLATE level, `0..=9`.
    pub compression_level: u32,
    /// Which engine to use.
    pub engine: EngineKind,
    /// Which I/O backend to use.
    pub io: IoBackend,
    /// The most outputs that may hold an open descriptor at once.
    pub max_open_files: usize,
    /// Where spool files go.
    pub temp_dir: Option<PathBuf>,
    /// Whether to keep spool files after the run.
    pub keep_temp: bool,
    /// Whether an existing output may be replaced.
    pub force: bool,
    /// Whether to suppress the `@PG` record.
    pub no_pg: bool,
    /// Which manifests to write.
    pub manifest: ManifestFormat,
    /// Which output index to build.
    pub index: IndexMode,
    /// The command line, recorded in `@PG CL` and the manifest.
    pub command_line: String,
    /// An optional output filename template.
    pub filename_template: Option<String>,
    /// The largest accepted BAM record body.
    pub max_record_size: usize,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            compression_level: 6,
            engine: EngineKind::Auto,
            io: IoBackend::Auto,
            max_open_files: 64,
            temp_dir: None,
            keep_temp: false,
            force: false,
            no_pg: false,
            manifest: ManifestFormat::Json,
            index: IndexMode::Auto,
            command_line: "bamsplit".to_string(),
            filename_template: None,
            max_record_size: crate::bam::raw_record::DEFAULT_MAX_RECORD_SIZE,
        }
    }
}

impl RunOptions {
    /// Sets the thread budget.
    #[must_use]
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Validates the numeric options.
    ///
    /// # Errors
    ///
    /// Returns [`ConfigError::OutOfRange`] for a compression level above 9 or a
    /// zero descriptor budget.
    pub fn validate(&self) -> std::result::Result<(), ConfigError> {
        if self.compression_level > 9 {
            return Err(ConfigError::OutOfRange {
                option: "--compression-level",
                constraint: "between 0 and 9",
                value: self.compression_level.to_string(),
            });
        }
        if self.max_open_files == 0 {
            return Err(ConfigError::OutOfRange {
                option: "--max-open-files",
                constraint: "at least 1",
                value: "0".to_string(),
            });
        }
        if self.max_record_size < crate::bam::FIXED_CORE_SIZE {
            return Err(ConfigError::OutOfRange {
                option: "--max-record-size",
                constraint: "at least 32, the size of the BAM fixed core",
                value: self.max_record_size.to_string(),
            });
        }
        Ok(())
    }
}

/// `bamsplit chrom` options.
#[derive(Debug, Clone)]
pub struct ChromOptions {
    /// Where outputs go.
    pub out_dir: PathBuf,
    /// Where placed-unmapped records go.
    pub placed_unmapped: PlacedUnmapped,
    /// What to do with unplaced records.
    pub unplaced: UnplacedPolicy,
    /// The stem of the unplaced output.
    pub unmapped_name: String,
    /// Whether to write header-only BAMs for references with no records.
    pub emit_empty: bool,
    /// Only these references.
    pub include: Vec<String>,
    /// Never these references.
    pub exclude: Vec<String>,
    /// A file of reference names, one per line, merged into `include`.
    pub reference_list: Option<PathBuf>,
    /// Whether an unknown requested reference is tolerated.
    pub ignore_missing_references: bool,
}

impl ChromOptions {
    /// Options writing into `out_dir`, with every default.
    pub fn new(out_dir: impl Into<PathBuf>) -> Self {
        Self {
            out_dir: out_dir.into(),
            placed_unmapped: PlacedUnmapped::ByReference,
            unplaced: UnplacedPolicy::Keep,
            unmapped_name: "unmapped".to_string(),
            emit_empty: false,
            include: Vec::new(),
            exclude: Vec::new(),
            reference_list: None,
            ignore_missing_references: false,
        }
    }
}

/// `bamsplit shard` options.
#[derive(Debug, Clone)]
pub struct ShardOptions {
    /// Where outputs go.
    pub out_dir: PathBuf,
    /// How many shards.
    pub shards: u32,
    /// What to hash.
    pub key: ShardKeySource,
    /// The hash seed.
    pub seed: u64,
    /// What to do when the key is absent.
    pub missing: MissingPolicy,
    /// The stem of the no-key output.
    pub missing_name: String,
    /// Zero-padding width.
    pub shard_width: Option<usize>,
}

impl ShardOptions {
    /// Options writing `shards` shards into `out_dir`.
    pub fn new(out_dir: impl Into<PathBuf>, shards: u32) -> Self {
        Self {
            out_dir: out_dir.into(),
            shards,
            key: ShardKeySource::QName,
            seed: crate::routing::shard::DEFAULT_SEED,
            missing: MissingPolicy::File,
            missing_name: "no-key".to_string(),
            shard_width: None,
        }
    }
}

/// `bamsplit region` options.
#[cfg(feature = "annotation")]
#[derive(Debug, Clone)]
pub struct RegionOptions {
    /// Where outputs go.
    pub out_dir: PathBuf,
    /// A BED, GTF, or GFF annotation. Mutually exclusive with `window_size`.
    pub regions: Option<PathBuf>,
    /// Generate fixed windows instead of reading an annotation.
    pub window_size: Option<String>,
    /// Override the detected annotation format.
    pub format: Option<crate::annotation::AnnotationFormat>,
    /// Override the detected BED width.
    pub bed_type: Option<crate::annotation::BedType>,
    /// Which segments of each record participate.
    pub feature: crate::annotation::FeatureType,
    /// How an alignment is matched against them.
    pub assignment: crate::routing::region::AssignmentMode,
    /// Compare against aligned blocks or the outer span.
    pub geometry: crate::bam::AlignmentGeometry,
    /// An annotation attribute to prefer as the logical name.
    pub name_field: Option<String>,
    /// The prefix for generated names.
    pub unnamed_prefix: String,
    /// What to do with records lacking the requested feature.
    pub missing_feature: crate::annotation::MissingFeaturePolicy,
    /// Turn the BED span-as-single-exon fallback into an error.
    pub require_blocks: bool,
    /// Write header-only BAMs for regions with no records.
    pub emit_empty: bool,
}

#[cfg(feature = "annotation")]
impl RegionOptions {
    /// Options reading `regions` into `out_dir`.
    pub fn from_annotation(out_dir: impl Into<PathBuf>, regions: impl Into<PathBuf>) -> Self {
        Self {
            out_dir: out_dir.into(),
            regions: Some(regions.into()),
            window_size: None,
            format: None,
            bed_type: None,
            feature: crate::annotation::FeatureType::Span,
            assignment: crate::routing::region::AssignmentMode::Start,
            geometry: crate::bam::AlignmentGeometry::Blocks,
            name_field: None,
            unnamed_prefix: "region".to_string(),
            missing_feature: crate::annotation::MissingFeaturePolicy::Skip,
            require_blocks: false,
            emit_empty: false,
        }
    }

    /// Options generating fixed windows into `out_dir`.
    pub fn from_windows(out_dir: impl Into<PathBuf>, window_size: impl Into<String>) -> Self {
        Self {
            regions: None,
            window_size: Some(window_size.into()),
            ..Self::from_annotation(out_dir, PathBuf::new())
        }
    }
}

/// `bamsplit tag` options.
#[derive(Debug, Clone)]
pub struct TagOptions {
    /// Where outputs go.
    pub out_dir: PathBuf,
    /// The auxiliary tag to route on.
    pub tag: Option<String>,
    /// The `@RG` field to route on.
    pub field: Option<ReadGroupField>,
    /// What to do when the value is absent.
    pub missing: MissingPolicy,
    /// The stem of the no-value output.
    pub missing_name: String,
    /// What to do when a read group is unknown.
    pub unknown_read_group: MissingPolicy,
    /// The cardinality limit.
    pub max_outputs: usize,
    /// Whether to lift the limit.
    pub allow_high_cardinality: bool,
    /// Whether `B` arrays may be keys.
    pub allow_array_tags: bool,
}

impl TagOptions {
    /// Options routing on an auxiliary tag.
    pub fn by_tag(out_dir: impl Into<PathBuf>, tag: impl Into<String>) -> Self {
        Self {
            out_dir: out_dir.into(),
            tag: Some(tag.into()),
            field: None,
            missing: MissingPolicy::File,
            missing_name: "no-tag".to_string(),
            unknown_read_group: MissingPolicy::File,
            max_outputs: 1_000,
            allow_high_cardinality: false,
            allow_array_tags: false,
        }
    }

    /// Options routing on a `@RG` field.
    pub fn by_field(out_dir: impl Into<PathBuf>, field: ReadGroupField) -> Self {
        Self {
            tag: None,
            field: Some(field),
            ..Self::by_tag(out_dir, "RG")
        }
    }
}

/// What a split produced.
#[derive(Debug)]
pub struct SplitReport {
    /// The validated manifest.
    pub manifest: Manifest,
    /// The manifest files written.
    pub manifest_paths: Vec<PathBuf>,
}

impl SplitReport {
    /// How many outputs were created, excluding skipped ones.
    #[must_use]
    pub fn created_outputs(&self) -> usize {
        self.manifest
            .outputs
            .iter()
            .filter(|output| !output.skipped)
            .count()
    }
}

/// Splits a BAM into one output per reference sequence.
///
/// # Errors
///
/// Returns [`Error`] for invalid options, malformed input, an output conflict,
/// an interruption, or a failed conservation check.
pub fn split_by_chromosome(
    input: impl AsRef<Path>,
    options: &ChromOptions,
    run: &RunOptions,
) -> Result<SplitReport> {
    run.validate()?;
    let source = InputSource::from_arg(input.as_ref());
    source.validate()?;
    let header = read_header(&source, run)?;

    let mut include: Vec<Vec<u8>> = options
        .include
        .iter()
        .map(|name| name.as_bytes().to_vec())
        .collect();
    if let Some(path) = &options.reference_list {
        include.extend(crate::routing::chrom::read_reference_list(path)?);
    }

    let router = ChromRouter::new(
        &header,
        ChromRouterOptions {
            placed_unmapped: options.placed_unmapped,
            unplaced: options.unplaced,
            unmapped_name: options.unmapped_name.clone().into_bytes(),
            include,
            exclude: options
                .exclude
                .iter()
                .map(|name| name.as_bytes().to_vec())
                .collect(),
            ignore_missing_references: options.ignore_missing_references,
        },
    )?;

    run_split(
        &source,
        &header,
        &options.out_dir,
        &router,
        RoutingMode::Chrom,
        run,
        options.emit_empty,
        None,
        None,
    )
}

/// Splits a BAM into a fixed number of deterministic shards.
///
/// # Errors
///
/// As [`split_by_chromosome`].
pub fn split_by_shard(
    input: impl AsRef<Path>,
    options: &ShardOptions,
    run: &RunOptions,
) -> Result<SplitReport> {
    run.validate()?;
    let source = InputSource::from_arg(input.as_ref());
    source.validate()?;
    let header = read_header(&source, run)?;

    let router = ShardRouter::new(ShardRouterOptions {
        shards: options.shards,
        key: options.key.clone(),
        seed: options.seed,
        missing: options.missing,
        missing_name: options.missing_name.clone().into_bytes(),
        shard_width: options.shard_width,
    })?;

    run_split(
        &source,
        &header,
        &options.out_dir,
        &router,
        RoutingMode::Shard,
        run,
        false,
        None,
        None,
    )
}

/// Splits a BAM into one output per tag value.
///
/// # Errors
///
/// As [`split_by_chromosome`], plus
/// [`RoutingError::CardinalityExceeded`](crate::error::RoutingError::CardinalityExceeded)
/// when the split would produce more outputs than `max_outputs`.
pub fn split_by_tag(
    input: impl AsRef<Path>,
    options: &TagOptions,
    run: &RunOptions,
) -> Result<SplitReport> {
    run.validate()?;
    crate::routing::tag::require_exactly_one_source(
        options.tag.as_deref(),
        options.field.map(ReadGroupField::as_str),
    )?;
    let source = InputSource::from_arg(input.as_ref());
    source.validate()?;
    let header = read_header(&source, run)?;

    let tag_source = match (&options.tag, options.field) {
        (Some(tag), _) => TagSource::Tag(crate::routing::parse_tag(tag)?),
        (None, Some(field)) => TagSource::Field(field),
        (None, None) => {
            return Err(Error::from(ConfigError::ExactlyOneRequired {
                options: "`--tag`, `--field`",
            }));
        }
    };

    let router = TagRouter::new(TagRouterOptions {
        source: tag_source,
        missing: options.missing,
        missing_name: options.missing_name.clone().into_bytes(),
        unknown_read_group: options.unknown_read_group,
        max_outputs: options.max_outputs,
        allow_high_cardinality: options.allow_high_cardinality,
        allow_array_tags: options.allow_array_tags,
    })?;

    run_split(
        &source,
        &header,
        &options.out_dir,
        &router,
        RoutingMode::Tag,
        run,
        false,
        None,
        None,
    )
}

/// Splits a BAM into one output per annotated region, or per generated window.
///
/// # Errors
///
/// As [`split_by_chromosome`], plus every [`crate::error::AnnotationError`] the
/// annotation can raise.
#[cfg(feature = "annotation")]
pub fn split_by_region(
    input: impl AsRef<Path>,
    options: &RegionOptions,
    run: &RunOptions,
) -> Result<SplitReport> {
    use crate::annotation::{
        FeatureOptions, IntervalIndex, LoadOptions, extract_regions, generate_windows,
        parse_window_size,
    };
    use crate::routing::region::{RegionRouter, RegionRouterOptions};

    run.validate()?;
    if options.regions.is_some() == options.window_size.is_some() {
        return Err(Error::from(ConfigError::ExactlyOneRequired {
            options: "`--regions`, `--window-size`",
        }));
    }

    let source = InputSource::from_arg(input.as_ref());
    source.validate()?;
    let header = read_header(&source, run)?;

    let feature_options = FeatureOptions {
        feature: options.feature,
        missing: options.missing_feature,
        require_blocks: options.require_blocks,
        name_field: options
            .name_field
            .as_ref()
            .map(|field| field.as_bytes().to_vec()),
        unnamed_prefix: options.unnamed_prefix.clone(),
    };

    let (regions, region_manifest) = if let Some(path) = &options.regions {
        let loaded = crate::annotation::load(
            path,
            &LoadOptions {
                format: options.format,
                bed_type: options.bed_type,
                use_mmap: matches!(run.io, IoBackend::Mmap),
            },
        )?;
        let (regions, stats) =
            extract_regions(&loaded.records, &loaded.detection, &feature_options)?;
        let manifest = crate::manifest::RegionManifest {
            annotation_source: path.display().to_string(),
            annotation_format: loaded.detection.format.to_string(),
            bed_type: loaded.detection.bed_type.map(|bed| bed.to_string()),
            feature_type: options.feature.to_string(),
            assignment_mode: options.assignment.to_string(),
            alignment_geometry: options.geometry.to_string(),
            annotation_records: stats.annotation_records,
            derived_segment_count: stats.derived_segments,
            missing_feature_records: stats.missing_feature,
            generated_names: stats.generated_names,
            duplicate_name_resolutions: stats.duplicate_names,
            ambiguous_assignments: 0,
            span_as_exon_fallbacks: stats.span_as_exon_fallbacks,
            empty_regions: stats.empty_regions,
        };
        (regions, manifest)
    } else {
        let text = options.window_size.as_deref().unwrap_or_default();
        let size = parse_window_size(text)?;
        let references: Vec<(Vec<u8>, u32)> = header
            .references()
            .iter()
            .map(|reference| (reference.name.clone(), reference.length))
            .collect();
        let regions = generate_windows(&references, size);
        let manifest = crate::manifest::RegionManifest {
            annotation_source: format!("generated windows of {size} bases"),
            annotation_format: "windows".to_string(),
            bed_type: None,
            feature_type: crate::annotation::FeatureType::Span.to_string(),
            assignment_mode: options.assignment.to_string(),
            alignment_geometry: options.geometry.to_string(),
            annotation_records: regions.len() as u64,
            derived_segment_count: regions.len() as u64,
            ..crate::manifest::RegionManifest::default()
        };
        (regions, manifest)
    };

    if regions.is_empty() {
        return Err(Error::from(crate::error::AnnotationError::Empty {
            path: options.regions.clone().unwrap_or_default(),
        }));
    }

    // An annotation whose reference names match nothing in the BAM is almost
    // always a genome-build mismatch — `1` versus `chr1`. Every record would be
    // unmatched, so it is worth saying so rather than writing an empty split.
    let index = IntervalIndex::build(regions)?;
    if !index
        .reference_names()
        .any(|name| header.reference_id(name).is_some())
    {
        return Err(Error::from(
            crate::error::AnnotationError::NoReferenceOverlap {
                annotation_examples: sample_names(index.reference_names()),
                bam_examples: sample_names(
                    header
                        .references()
                        .iter()
                        .map(|reference| &reference.name[..]),
                ),
            },
        ));
    }

    let dense = crate::planning::annotation_is_dense(index.covered_bases(), &header);
    let region_count = index.region_count();

    // The intervals the indexed engine would query. Built unconditionally
    // because they are cheap and only that engine consults them; the planner
    // decides separately whether it runs.
    let targets = crate::engine::indexed::QueryTargets::from_intervals(
        index.regions().iter().filter_map(|region| {
            header
                .reference_id(&region.chrom)
                .map(|reference_id| (reference_id, region.envelope))
        }),
    );

    let router = RegionRouter::new(
        index,
        RegionRouterOptions {
            assignment: options.assignment,
            geometry: options.geometry,
            deletions: crate::bam::DeletionPolicy::ExcludeFromBlocks,
        },
    );

    let mut report = run_split(
        &source,
        &header,
        &options.out_dir,
        &router,
        RoutingMode::Region {
            dense,
            regions: region_count,
        },
        run,
        options.emit_empty,
        Some(region_manifest),
        Some(targets),
    )?;

    // The ambiguity count is only known once every record has been routed.
    let stats = router.stats();
    if let Some(region) = report.manifest.region.as_mut() {
        region.ambiguous_assignments = stats.ambiguous_assignments;
    }
    report.manifest.max_emissions_for_one_record = stats.max_matches_for_one_record;
    report.manifest.write(&options.out_dir, run.manifest)?;
    Ok(report)
}

/// A few reference names, for a diagnostic.
#[cfg(feature = "annotation")]
fn sample_names<'a>(names: impl Iterator<Item = &'a [u8]>) -> String {
    let mut sample: Vec<String> = names
        .take(4)
        .map(|name| format!("{:?}", String::from_utf8_lossy(name)))
        .collect();
    sample.sort();
    sample.join(", ")
}

/// Reads and validates the input header without splitting anything.
fn read_header(source: &InputSource, run: &RunOptions) -> Result<BamHeader> {
    let (header, _, _) =
        crate::engine::open_bam(source, IoBackend::Buffered, 1, run.max_record_size)
            .map_err(Error::Engine)?;
    Ok(header)
}

/// The shared body of every split.
#[allow(clippy::too_many_arguments)]
fn run_split<R: Router>(
    source: &InputSource,
    header: &BamHeader,
    out_dir: &Path,
    router: &R,
    mode: RoutingMode,
    run: &RunOptions,
    emit_empty: bool,
    region: Option<crate::manifest::RegionManifest>,
    targets: Option<crate::engine::indexed::QueryTargets>,
) -> Result<SplitReport> {
    let started_at = SystemTime::now();
    let clock = Instant::now();

    let facts = InputFacts::gather(source, header);
    let plan = crate::planning::plan(run.engine, run.io, mode, &facts, run.threads)?;

    let longest = header
        .references()
        .iter()
        .max_by_key(|reference| reference.length)
        .map(|reference| reference.display_name().into_owned())
        .unwrap_or_default();
    let index_plan: IndexPlan = crate::index::plan(
        run.index,
        u64::from(header.max_reference_length()),
        &longest,
        out_dir,
    )?;

    // Every output gets the input's complete reference dictionary plus one new
    // `@PG` record; see `bam::header` for why the dictionary is never reduced.
    let mut output_header = header.clone();
    if !run.no_pg {
        output_header.append_program_record(crate::VERSION, &run.command_line)?;
    }

    let template = match &run.filename_template {
        Some(text) => FilenameTemplate::parse(text)?,
        None => FilenameTemplate::default(),
    };

    let budget = match plan.engine {
        EngineKind::Stream => ThreadBudget::for_total(run.threads),
        // The indexed and spool engines both parallelize across *outputs*, so
        // their budget is split by task rather than concentrated on one writer.
        EngineKind::Indexed | EngineKind::Spool | EngineKind::Auto => {
            ThreadBudget::for_tasks(run.threads, run.threads)
        }
    };
    let pool = budget.build_pool().map_err(Error::Engine)?;

    let context = ExecutionContext {
        input: source.clone(),
        header: &output_header,
        output_directory: out_dir.to_path_buf(),
        manager_options: OutputManagerOptions {
            max_open_files: run.max_open_files,
            settings: OutputSettings {
                compression_level: run.compression_level,
                force: run.force,
                index_kind: index_plan.kind,
                compression_workers: budget.compression,
                pool,
                write_buffer_bytes: 1 << 20,
            },
            template,
        },
        threads: budget,
        emit_empty,
        temp_dir: run.temp_dir.clone(),
        keep_temp: run.keep_temp,
        io: plan.io,
        interrupt: crate::output::interrupt_flag(),
        max_record_size: run.max_record_size,
    };

    let report = execute(&context, router, plan.engine, source, header, targets)?;

    let manifest = build_manifest(
        source,
        header,
        &facts,
        &plan,
        &index_plan,
        router,
        run,
        &report,
        region,
        started_at,
        clock.elapsed().as_secs_f64(),
    );
    let manifest_paths = manifest.write(out_dir, run.manifest)?;

    Ok(SplitReport {
        manifest,
        manifest_paths,
    })
}

fn execute<R: Router>(
    context: &ExecutionContext<'_>,
    router: &R,
    engine: EngineKind,
    source: &InputSource,
    header: &BamHeader,
    targets: Option<crate::engine::indexed::QueryTargets>,
) -> Result<ExecutionReport> {
    use crate::engine::indexed::{IndexedEngine, find_index, load_index};
    use crate::engine::spool::SpoolEngine;
    use crate::engine::stream::StreamEngine;

    let outcome = match engine {
        EngineKind::Stream | EngineKind::Auto => {
            let stream = StreamEngine {
                verify_coordinate_order: router.is_grouped_by_coordinate()
                    && header.sort_order() == SortOrder::Coordinate,
                finalize_on_key_change: router.is_grouped_by_coordinate(),
            };
            stream.execute(context, router)
        }
        EngineKind::Spool => SpoolEngine.execute(context, router),
        EngineKind::Indexed => {
            let Some(path) = source.path().and_then(find_index) else {
                return Err(Error::from(ConfigError::IncompatibleEngine {
                    engine: "indexed",
                    reason: "no `.bai` or `.csi` index was found next to the input".to_string(),
                }));
            };
            let index = load_index(&path, header)?;
            let engine = match targets {
                Some(targets) => IndexedEngine::per_region(index, targets),
                None => IndexedEngine::per_reference(index),
            };
            engine.execute(context, router)
        }
    };

    match outcome {
        Ok(report) => Ok(report),
        Err(
            crate::error::EngineError::SortOrderViolation { .. }
            | crate::error::EngineError::UngroupedInput { .. },
        ) if engine == EngineKind::Stream && source.is_seekable() => {
            // The header claimed an ordering the records do not have. Every
            // partial output has already been removed, so a clean retry with the
            // spool engine is safe — and is what `--engine auto` promises.
            tracing::warn!(
                "the input is not grouped by routing key despite its header; retrying with the \
                 spool engine"
            );
            let mut report = SpoolEngine
                .execute(context, router)
                .map_err(Error::Engine)?;
            report.notes.push(
                "the streaming engine detected an ordering violation, so the spool engine was \
                 used instead"
                    .to_string(),
            );
            Ok(report)
        }
        Err(error) => Err(Error::Engine(error)),
    }
}

#[allow(clippy::too_many_arguments)]
fn build_manifest<R: Router>(
    source: &InputSource,
    header: &BamHeader,
    facts: &InputFacts,
    plan: &crate::planning::Plan,
    index_plan: &IndexPlan,
    router: &R,
    run: &RunOptions,
    report: &ExecutionReport,
    region: Option<crate::manifest::RegionManifest>,
    started_at: SystemTime,
    elapsed: f64,
) -> Manifest {
    let mut outputs: Vec<OutputEntry> = report
        .outputs
        .iter()
        .map(|output| OutputEntry {
            logical_key: String::from_utf8_lossy(&output.paths.logical_key).into_owned(),
            logical_key_hex: hex(&output.paths.logical_key),
            resolved_key: output.paths.stem.clone(),
            encoded_key: output.paths.encoded_key.clone(),
            bam_path: file_name(&output.paths.bam),
            index_path: output.index_path.as_deref().map(file_name),
            index_type: output
                .index_path
                .as_ref()
                .and_then(|_| index_plan.kind.map(|kind| kind.extension().to_string())),
            index_note: output.index_note.clone(),
            skipped: false,
            stats: output.stats.clone(),
            region: None,
        })
        .collect();

    outputs.extend(report.skipped_keys.iter().map(|key| OutputEntry {
        logical_key: String::from_utf8_lossy(key).into_owned(),
        logical_key_hex: hex(key),
        resolved_key: crate::output::filename::encode(key),
        encoded_key: crate::output::filename::encode(key),
        bam_path: format!("{}.bam", crate::output::filename::encode(key)),
        index_path: None,
        index_type: None,
        index_note: Some("no record routed here and `--emit-empty` was not set".to_string()),
        skipped: true,
        stats: crate::stats::OutputStats::new(),
        region: None,
    }));

    let mut notes = report.notes.clone();
    notes.push(plan.io_reason.clone());

    let mut manifest = Manifest {
        program: crate::PROGRAM_NAME.to_string(),
        version: crate::VERSION.to_string(),
        command: run.command_line.clone(),
        start_time: format_timestamp(started_at),
        end_time: format_timestamp(SystemTime::now()),
        elapsed_seconds: elapsed,
        input_path: source.display(),
        input_size: source.size(),
        input_header_checksum: format!("{:016x}", header.checksum()),
        input_index_type: facts.indexed.then(|| {
            source
                .path()
                .and_then(crate::engine::indexed::find_index)
                .and_then(|path| {
                    path.extension()
                        .and_then(std::ffi::OsStr::to_str)
                        .map(ToString::to_string)
                })
                .unwrap_or_else(|| "bai".to_string())
        }),
        input_sort_order: header.sort_order().to_string(),
        selected_engine: report.engine.to_string(),
        engine_selection_reason: plan.engine_reason.clone(),
        io_backend: report.io.to_string(),
        threads: run.threads,
        compression_level: run.compression_level,
        routing_mode: router.mode().to_string(),
        index_mode: run.index.to_string(),
        index_selection_reason: index_plan.reason.clone(),
        temporary_bytes: 0,
        input_records: 0,
        unique_emitted_records: 0,
        total_output_emissions: 0,
        dropped_records: 0,
        unmatched_records: 0,
        duplicate_emissions: 0,
        max_emissions_for_one_record: 0,
        may_duplicate: router.may_duplicate(),
        notes,
        region,
        outputs,
    };
    manifest.apply_run_stats(&report.stats);
    manifest
}

fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write as _;
    bytes.iter().fold(String::new(), |mut out, byte| {
        let _ = write!(out, "{byte:02x}");
        out
    })
}

fn file_name(path: &Path) -> String {
    path.file_name().map_or_else(
        || path.display().to_string(),
        |name| name.to_string_lossy().into_owned(),
    )
}

/// `bamsplit inspect` options.
#[derive(Debug, Clone, Default)]
pub struct InspectOptions {
    /// Whether to scan every record, not just the header.
    pub full: bool,
    /// A tag whose cardinality to measure during a full scan.
    pub tag: Option<String>,
    /// Which routing mode to predict output counts for.
    pub by: Option<RoutingMode>,
    /// An annotation to inspect alongside the BAM.
    pub regions: Option<PathBuf>,
    /// The feature type to check the annotation for.
    #[cfg(feature = "annotation")]
    pub feature: crate::annotation::FeatureType,
}

/// What `bamsplit inspect --regions` found.
#[cfg(feature = "annotation")]
#[derive(Debug, Clone, serde::Serialize)]
pub struct AnnotationReport {
    /// The annotation file.
    pub path: String,
    /// The detected format.
    pub format: String,
    /// The detected BED width, for BED input.
    pub bed_type: Option<String>,
    /// How the file is compressed.
    pub compression: String,
    /// How detection reached its conclusion.
    pub detection_reason: String,
    /// Columns beyond the chosen BED width.
    pub trailing_columns: usize,
    /// Records read.
    pub annotation_records: u64,
    /// The feature that was checked.
    pub feature: String,
    /// Regions the feature produced.
    pub usable_regions: usize,
    /// Segments across those regions.
    pub derived_segments: u64,
    /// Records lacking the information the feature needs.
    pub missing_feature_records: u64,
    /// Records with no `thickStart`/`thickEnd`.
    pub missing_coding_bounds: u64,
    /// Records with no usable strand.
    pub missing_strand: u64,
    /// Records with no block structure.
    pub missing_blocks: u64,
    /// Names that had to be invented.
    pub generated_names: u64,
    /// Duplicate names that were disambiguated.
    pub duplicate_names: u64,
    /// BED records whose whole span stood in for missing blocks.
    pub span_as_exon_fallbacks: u64,
    /// Regions that produced no segments.
    pub empty_regions: u64,
    /// Pairs of regions whose envelopes overlap on the same reference.
    pub overlapping_regions: u64,
    /// Total reference bases the regions cover.
    pub covered_bases: u64,
    /// Annotation reference names absent from the BAM header.
    pub references_absent_from_bam: Vec<String>,
    /// BAM reference names the annotation never mentions.
    pub references_absent_from_annotation: Vec<String>,
}

/// What `bamsplit inspect` found.
#[derive(Debug, Clone, serde::Serialize)]
pub struct InspectReport {
    /// The input, as given.
    pub input_path: String,
    /// The input size, when known.
    pub input_size: Option<u64>,
    /// Whether the file is a valid BAM with a parseable header.
    pub valid_bam: bool,
    /// The serialized header size in bytes.
    pub header_bytes: usize,
    /// The declared sort order.
    pub sort_order: String,
    /// How many references the header declares.
    pub reference_count: usize,
    /// The reference names and lengths.
    pub references: Vec<(String, u32)>,
    /// The read-group identifiers.
    pub read_groups: Vec<String>,
    /// Whether the input can be seeked.
    pub seekable: bool,
    /// Whether the BGZF end-of-file marker is present.
    pub bgzf_eof_present: bool,
    /// The index format found next to the input, if any.
    pub index_type: Option<String>,
    /// Whether that index agrees with the header's reference count.
    pub index_consistent: Option<bool>,
    /// References longer than BAI can address.
    pub references_exceeding_bai: Vec<String>,
    /// Reference names whose encoded filename differs from the name itself.
    pub unsafe_output_names: Vec<(String, String)>,
    /// How many outputs a split would produce.
    pub predicted_output_count: usize,
    /// Which engine the planner would choose.
    pub recommended_engine: String,
    /// Why.
    pub engine_reason: String,
    /// Whether memory mapping is possible for this input.
    pub mmap_eligible: bool,
    /// The full-scan results, when `--full` was given.
    pub scan: Option<ScanReport>,
    /// The annotation results, when `--regions` was given.
    #[cfg(feature = "annotation")]
    pub annotation: Option<AnnotationReport>,
}

/// The results of a `--full` record scan.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct ScanReport {
    /// Records read.
    pub total_records: u64,
    /// Neither secondary nor supplementary.
    pub primary_records: u64,
    /// `SECONDARY` set.
    pub secondary_records: u64,
    /// `SUPPLEMENTARY` set.
    pub supplementary_records: u64,
    /// Placed and mapped.
    pub mapped_records: u64,
    /// Placed but unmapped.
    pub placed_unmapped_records: u64,
    /// Unplaced.
    pub unplaced_unmapped_records: u64,
    /// `DUPLICATE` set.
    pub duplicate_records: u64,
    /// `QC_FAIL` set.
    pub qc_fail_records: u64,
    /// Per-reference record counts, by name.
    pub per_reference_counts: Vec<(String, u64)>,
    /// How many records broke coordinate ordering.
    pub coordinate_order_violations: u64,
    /// How many distinct values the requested tag had.
    pub tag_cardinality: Option<usize>,
    /// How many records had a malformed auxiliary section.
    pub malformed_tags: u64,
    /// Total uncompressed record bytes, which bounds spool usage.
    pub estimated_temporary_bytes: u64,
    /// Whether every output could carry an index.
    pub estimated_indexable: bool,
}

/// Inspects a BAM without writing any output.
///
/// # Errors
///
/// Returns [`Error`] if the input cannot be opened or its header is invalid. A
/// malformed *record* during a full scan is reported in the result rather than
/// raised, because the point of `inspect` is to describe a file that may be
/// broken.
pub fn inspect(
    input: impl AsRef<Path>,
    options: &InspectOptions,
    run: &RunOptions,
) -> Result<InspectReport> {
    let source = InputSource::from_arg(input.as_ref());
    source.validate()?;
    let (header, mut reader, _) =
        crate::engine::open_bam(&source, IoBackend::Buffered, 1, run.max_record_size)
            .map_err(Error::Engine)?;

    let facts = InputFacts::gather(&source, &header);
    let mode = options.by.unwrap_or(RoutingMode::Chrom);
    let plan = crate::planning::plan(EngineKind::Auto, IoBackend::Auto, mode, &facts, run.threads)?;

    let index_path = source.path().and_then(crate::engine::indexed::find_index);
    let index_consistent = index_path
        .as_deref()
        .map(|path| crate::engine::indexed::load_index(path, &header).is_ok());

    let unsafe_output_names = header
        .references()
        .iter()
        .filter_map(|reference| {
            let encoded = crate::output::filename::encode(&reference.name);
            let logical = reference.display_name().into_owned();
            (encoded != logical).then_some((logical, encoded))
        })
        .collect();

    let mut report = InspectReport {
        input_path: source.display(),
        input_size: source.size(),
        valid_bam: true,
        header_bytes: header.serialized_len().unwrap_or(0),
        sort_order: header.sort_order().to_string(),
        reference_count: header.reference_count(),
        references: header
            .references()
            .iter()
            .map(|reference| (reference.display_name().into_owned(), reference.length))
            .collect(),
        read_groups: header
            .read_groups()
            .iter()
            .map(|group| String::from_utf8_lossy(&group.id).into_owned())
            .collect(),
        seekable: facts.seekable,
        bgzf_eof_present: source
            .path()
            .is_some_and(|path| has_bgzf_eof(path).unwrap_or(false)),
        index_type: index_path
            .as_deref()
            .and_then(Path::extension)
            .and_then(std::ffi::OsStr::to_str)
            .map(ToString::to_string),
        index_consistent,
        references_exceeding_bai: header
            .references()
            .iter()
            .filter(|reference| !crate::index::bai::can_represent(u64::from(reference.length)))
            .map(|reference| reference.display_name().into_owned())
            .collect(),
        unsafe_output_names,
        predicted_output_count: header.reference_count() + 1,
        recommended_engine: plan.engine.to_string(),
        engine_reason: plan.engine_reason,
        mmap_eligible: facts.local_file,
        scan: None,
        #[cfg(feature = "annotation")]
        annotation: None,
    };

    #[cfg(feature = "annotation")]
    if let Some(annotation) = &options.regions {
        report.annotation = Some(inspect_annotation(annotation, options, &header)?);
    }

    if options.full {
        let tag = options
            .tag
            .as_deref()
            .map(crate::routing::parse_tag)
            .transpose()?;
        let mut scan = ScanReport {
            estimated_indexable: true,
            ..ScanReport::default()
        };
        let mut counters = crate::bam::validation::RecordCounters::default();
        let mut order = crate::bam::validation::CoordinateOrderTracker::new();
        let mut per_reference = vec![0u64; header.reference_count()];
        let mut tag_values: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();

        loop {
            match reader.next_record() {
                Ok(None) => break,
                Ok(Some(record)) => {
                    let ordinal = counters.total + 1;
                    if counters.observe(&record).is_err() {
                        scan.malformed_tags += 1;
                        continue;
                    }
                    let _ = order.observe(&record, ordinal);
                    if let Ok(Some(id)) = record.reference_sequence_id()
                        && let Ok(index) = usize::try_from(id)
                        && let Some(slot) = per_reference.get_mut(index)
                    {
                        *slot += 1;
                    }
                    scan.estimated_temporary_bytes += record.len() as u64;
                    if let Some(tag) = tag {
                        match record.tag(tag) {
                            Ok(Some(value)) => {
                                tag_values.insert(value.to_routing_bytes());
                            }
                            Ok(None) => {}
                            Err(_) => scan.malformed_tags += 1,
                        }
                    }
                }
                Err(_) => {
                    // A malformed record ends the scan: everything after it is
                    // unreliable. The counts up to that point are still useful,
                    // which is why this is reported rather than raised.
                    report.valid_bam = false;
                    break;
                }
            }
        }

        scan.total_records = counters.total;
        scan.primary_records = counters.primary;
        scan.secondary_records = counters.secondary;
        scan.supplementary_records = counters.supplementary;
        scan.mapped_records = counters.mapped;
        scan.placed_unmapped_records = counters.placed_unmapped;
        scan.unplaced_unmapped_records = counters.unplaced_unmapped;
        scan.duplicate_records = counters.duplicate;
        scan.qc_fail_records = counters.qc_fail;
        scan.coordinate_order_violations = order.violations();
        scan.estimated_indexable = order.is_sorted();
        scan.tag_cardinality = tag.map(|_| tag_values.len());
        scan.per_reference_counts = header
            .references()
            .iter()
            .zip(&per_reference)
            .filter(|(_, count)| **count > 0)
            .map(|(reference, count)| (reference.display_name().into_owned(), *count))
            .collect();

        report.predicted_output_count =
            scan.per_reference_counts.len() + usize::from(scan.unplaced_unmapped_records > 0);
        report.scan = Some(scan);
    }

    Ok(report)
}

/// Inspects an annotation alongside the BAM.
///
/// Runs the same detection, loading, and extraction a split would, so what it
/// reports is what a split would actually do — including which records the
/// requested feature cannot be derived from.
#[cfg(feature = "annotation")]
fn inspect_annotation(
    path: &Path,
    options: &InspectOptions,
    header: &BamHeader,
) -> Result<AnnotationReport> {
    use crate::annotation::{FeatureOptions, LoadOptions, MissingFeaturePolicy, extract_regions};
    use std::collections::BTreeSet;

    let loaded = crate::annotation::load(path, &LoadOptions::default())?;

    // `Empty` keeps the records that a split would skip, so the report can
    // count *why* each one is unusable rather than just how many vanished.
    let (regions, stats) = extract_regions(
        &loaded.records,
        &loaded.detection,
        &FeatureOptions {
            feature: options.feature,
            missing: MissingFeaturePolicy::Empty,
            ..FeatureOptions::default()
        },
    )?;

    // Attribute each missing record to a cause, which is the actionable part:
    // "3 000 records have no CDS" tells a user to pick a different feature.
    let mut missing_coding = 0u64;
    let mut missing_strand = 0u64;
    let mut missing_blocks = 0u64;
    for record in &loaded.records {
        let has_blocks = record
            .block_count
            .is_some_and(|count| count > 0 && record.block_starts.is_some());
        let has_coding = matches!(
            (record.thick_start, record.thick_end),
            (Some(start), Some(end)) if start < end
        );
        let has_strand = matches!(
            record.strand,
            Some(genepred::Strand::Forward | genepred::Strand::Reverse)
        );
        if options.feature.needs_blocks() && !has_blocks {
            missing_blocks += 1;
        }
        if options.feature.needs_coding_bounds() && !has_coding {
            missing_coding += 1;
        }
        if options.feature.needs_strand() && !has_strand {
            missing_strand += 1;
        }
    }

    // Overlapping envelopes on one reference, counted pairwise by a sweep.
    let mut by_reference: std::collections::HashMap<&[u8], Vec<crate::Interval>> =
        std::collections::HashMap::new();
    for region in &regions {
        if !region.segments.is_empty() {
            by_reference
                .entry(&region.chrom)
                .or_default()
                .push(region.envelope);
        }
    }
    let mut overlapping = 0u64;
    for envelopes in by_reference.values_mut() {
        envelopes.sort_unstable();
        for (index, envelope) in envelopes.iter().enumerate() {
            overlapping += envelopes[index + 1..]
                .iter()
                .take_while(|later| later.start < envelope.end)
                .count() as u64;
        }
    }

    let annotation_references: BTreeSet<Vec<u8>> =
        regions.iter().map(|region| region.chrom.clone()).collect();
    let bam_references: BTreeSet<Vec<u8>> = header
        .references()
        .iter()
        .map(|reference| reference.name.clone())
        .collect();
    let render = |names: Vec<&Vec<u8>>| -> Vec<String> {
        names
            .into_iter()
            .take(20)
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect()
    };

    let covered_bases: u64 = regions
        .iter()
        .map(|region| region.covered_bases().max(0) as u64)
        .sum();

    Ok(AnnotationReport {
        path: path.display().to_string(),
        format: loaded.detection.format.to_string(),
        bed_type: loaded.detection.bed_type.map(|bed| bed.to_string()),
        compression: loaded.detection.compression.to_string(),
        detection_reason: loaded.detection.reason.clone(),
        trailing_columns: loaded.detection.trailing_columns,
        annotation_records: stats.annotation_records,
        feature: options.feature.to_string(),
        usable_regions: regions.iter().filter(|region| !region.is_empty()).count(),
        derived_segments: stats.derived_segments,
        missing_feature_records: stats.missing_feature,
        missing_coding_bounds: missing_coding,
        missing_strand,
        missing_blocks,
        generated_names: stats.generated_names,
        duplicate_names: stats.duplicate_names,
        span_as_exon_fallbacks: stats.span_as_exon_fallbacks,
        empty_regions: stats.empty_regions,
        overlapping_regions: overlapping,
        covered_bases,
        references_absent_from_bam: render(
            annotation_references.difference(&bam_references).collect(),
        ),
        references_absent_from_annotation: render(
            bam_references.difference(&annotation_references).collect(),
        ),
    })
}

/// Whether a file ends with the 28-byte BGZF end-of-file marker.
///
/// A missing marker means the file was truncated — or written by a tool that
/// omits it — and every downstream reader will warn about it, so `inspect`
/// checks explicitly.
fn has_bgzf_eof(path: &Path) -> std::io::Result<bool> {
    use std::io::{Read as _, Seek as _, SeekFrom};

    let mut file = std::fs::File::open(path)?;
    let length = file.metadata()?.len();
    let marker = crate::output::BGZF_EOF;
    if length < marker.len() as u64 {
        return Ok(false);
    }
    file.seek(SeekFrom::End(-(marker.len() as i64)))?;
    let mut tail = [0u8; 28];
    file.read_exact(&mut tail)?;
    Ok(tail == marker)
}

/// Run-wide statistics, re-exported for callers that drive engines directly.
pub type Stats = RunStats;
