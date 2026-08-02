# bamsplit

Lossless, high-throughput partitioning of BAM files — by chromosome, by
deterministic shard, by auxiliary tag, or by annotated region.

```console
bamsplit chrom sample.bam --out-dir chromosomes --threads 16 --index auto
```

```text
chromosomes/
├── chr1.bam
├── chr1.bam.bai
├── chr2.bam
├── chr2.bam.bai
├── chrX.bam
├── chrX.bam.bai
├── unmapped.bam
├── unmapped.bam.bai
└── bamsplit.manifest.json
```

## Purpose

Scattering a BAM across chromosomes is the first step of almost every parallel
genomics pipeline, and the usual way to do it is a loop:

```bash
for chr in $(samtools idxstats in.bam | cut -f1); do
    samtools view -b in.bam "$chr" > "$chr.bam"
done
```

That works, and it re-reads and re-decompresses the whole input once per
chromosome. On a 3 000-contig assembly it reads the file 3 000 times.

`bamsplit chrom` reads it **once**. A coordinate-sorted BAM already visits its
chromosomes in contiguous runs, so a single sequential pass can open `chr1.bam`,
write every `chr1` record, close it, and move on — with one output open at a
time, constant memory, and the index built as it writes.

What `bamsplit` adds beyond speed:

* **Losslessness you can check.** Record bodies are copied byte for byte, and
  every run writes a manifest carrying per-output digests and the conservation
  equation `input = emitted + dropped + unmatched`. A run whose accounting does
  not add up **fails**, even if every byte was written.
* **Safe filenames.** A reference called `chr1/alternate`, `HLA-A*01:01`, or `..`
  cannot write outside the output directory.
* **Transactional outputs.** An output appears at its final path only once it is
  complete, indexed, and validated. An interrupted run leaves nothing behind.
* **Bounded resources.** Memory and file descriptors stay inside configured
  limits no matter how many outputs a split produces.

## Installation

Requires Rust 1.85 or newer.

```console
git clone https://github.com/alejandrogzi/bamsplit
cd bamsplit
cargo build --release
./target/release/bamsplit --help
```

`samtools` is **not** required to build, test, or run `bamsplit`. It is used only
as an external oracle by the optional differential test-suite.

## Quick start

```console
# One BAM per reference sequence, plus unmapped.bam, plus indexes.
bamsplit chrom sample.bam --out-dir chromosomes --threads 16

# Only the primary assembly.
bamsplit chrom sample.bam --out-dir main --include chr1,chr2,chr3,chrX,chrY

# 32 shards that keep every read pair together.
bamsplit shard sample.bam --shards 32 --key qname --out-dir shards

# One BAM per sample, resolved through the @RG dictionary.
bamsplit tag sample.bam --field sample --out-dir by-sample

# One BAM per transcript, by exonic overlap.
bamsplit region sample.bam --regions genes.gtf.gz --feature exon \
    --assignment overlap --out-dir by-transcript

# What would a split look like, and what is in this file?
bamsplit inspect sample.bam --full
```

## Commands

### Global options

| option | default | meaning |
| --- | --- | --- |
| `--threads <INT>` | `1` | total thread budget, not per-output |
| `--compression-level <INT>` | `6` | DEFLATE level, `0`–`9` |
| `--engine <auto\|stream\|indexed\|spool>` | `auto` | execution engine |
| `--io <auto\|buffered\|mmap>` | `auto` | how input bytes are read |
| `--max-open-files <INT>` | `64` | concurrently open output descriptors |
| `--temp-dir <PATH>` | system temp | where spool files go |
| `--keep-temp` | off | keep spool files, for debugging |
| `--force` | off | replace existing outputs |
| `--log-level <error\|warn\|info\|debug\|trace>` | `info` | verbosity (stderr) |
| `--quiet` | off | suppress all logging |
| `--no-progress` | off | no progress indicator |
| `--no-pg` | off | do not append an `@PG` record |
| `--manifest <json\|tsv\|both\|none>` | `json` | which manifests to write |
| `--max-record-size <BYTES>` | 64 MiB | reject an implausible `block_size` |

`--help`, `--version`, and `<COMMAND> --help` all work.

### Exit codes

| code | meaning |
| --- | --- |
| `0` | success |
| `1` | general runtime failure |
| `2` | invalid command-line arguments |
| `3` | invalid or malformed input |
| `4` | output conflict — an output exists, or two keys collide |
| `5` | validation failure — conservation or digest check failed |
| `6` | interrupted, and rolled back |

### `bamsplit chrom`

```console
bamsplit chrom input.bam --out-dir chromosomes --threads 16 --index auto
```

| option | default | meaning |
| --- | --- | --- |
| `--out-dir <PATH>` | `chromosomes` | where outputs go |
| `--placed-unmapped <by-reference\|unmapped>` | `by-reference` | where an unmapped record *with* a reference goes |
| `--unplaced <keep\|drop\|error>` | `keep` | what to do with `refID == -1` |
| `--unmapped-name <STRING>` | `unmapped` | stem of the unplaced output |
| `--emit-empty` | off | write header-only BAMs for references with no records |
| `--include <NAME,...>` | all | only these references |
| `--exclude <NAME,...>` | none | never these references |
| `--reference-list <FILE>` | none | names from a file, merged into `--include` |
| `--ignore-missing-references` | off | tolerate a requested name the header lacks |
| `--index <auto\|bai\|csi\|none>` | `auto` | which output index to build |
| `--filename-template <TEMPLATE>` | `{key}` | output naming, e.g. `sample1_{key}` |

`--include`, `--exclude`, and `--reference-list` match the **original logical
reference names**, never encoded filenames: write `chr1/alt`, not
`chr1%2Falt`. An unknown name is an error unless
`--ignore-missing-references` is given.

### `bamsplit shard`

```console
bamsplit shard input.bam --shards 32 --key qname --out-dir shards
```

| option | default | meaning |
| --- | --- | --- |
| `--shards <INT>` | `8` | how many shards |
| `--key <qname\|record\|tag:XX\|read-group\|sample\|library\|platform-unit>` | `qname` | what to hash |
| `--seed <U64>` | fixed, documented | hash seed |
| `--missing <file\|drop\|error>` | `file` | when the key is absent |
| `--missing-name <STRING>` | `no-key` | stem of the no-key output |
| `--shard-width <INT>` | from `--shards` | zero-padding width |

The assignment is

```text
shard = xxh3_64_with_seed(key_bytes, seed) % shard_count
```

with a **fixed default seed**, so the same input and shard count always produce
the same layout — across runs, machines, and `bamsplit` versions. Outputs are
`shard-0000.bam`, `shard-0001.bam`, …

`--key qname` keeps every record of a template in one shard: both mates, every
secondary alignment, every supplementary alignment. No QNAME table is kept, so
memory does not grow with the number of templates.

`--key record` deliberately gives that up — it hashes the record's own bytes and
spreads records evenly regardless of template structure. Use it only when
downstream work is strictly per-record.

### `bamsplit tag`

```console
bamsplit tag input.bam --tag RG --out-dir by-read-group
bamsplit tag input.bam --field sample --out-dir by-sample
```

Exactly one of `--tag` and `--field` is required.

| option | default | meaning |
| --- | --- | --- |
| `--tag <XX>` | — | route on this auxiliary tag |
| `--field <read-group\|sample\|library\|platform\|platform-unit\|sequencing-center>` | — | route on a `@RG`-derived field |
| `--missing <file\|drop\|error>` | `file` | when the value is absent |
| `--missing-name <STRING>` | `no-tag` | stem of the no-value output |
| `--unknown-read-group <file\|drop\|error>` | `file` | when `RG` names a group the header lacks |
| `--max-outputs <INT>` | `1000` | refuse to produce more outputs than this |
| `--allow-high-cardinality` | off | lift the limit entirely |
| `--allow-array-tags` | off | permit `B`-array values as keys |

**Cardinality is the hazard.** `--tag RG` gives a handful of outputs; `--tag CB`
on a single-cell BAM gives hundreds of thousands. `bamsplit` counts distinct keys
and stops before exceeding `--max-outputs`, removes every partial file, and points
you at the bounded alternative:

```console
bamsplit shard input.bam --shards 64 --key tag:CB --out-dir by-barcode
```

Array-valued (`B`) tags are rejected by default: `CB:B:i,1,2,3` has no obvious
filename and no obvious equality. `--allow-array-tags` opts into the canonical
`<subtype>,<v0>,<v1>,…` rendering, which is order-sensitive. Floating-point
values are never formatted through anything locale-sensitive.

### `bamsplit region`

```console
bamsplit region input.bam --regions annotation.bed12 --feature span   --assignment start        --out-dir regions
bamsplit region input.bam --regions annotation.bed12 --feature exon   --assignment overlap      --out-dir transcripts
bamsplit region input.bam --regions annotation.gtf.gz --feature cds   --assignment best-overlap --out-dir cds
bamsplit region input.bam --regions annotation.gff3  --feature intron --assignment overlap      --out-dir introns
bamsplit region input.bam --window-size 50M          --assignment start --out-dir windows
```

Exactly one of `--regions` and `--window-size` is required.

| option | default | meaning |
| --- | --- | --- |
| `--regions <PATH>` | — | a BED, GTF, or GFF annotation |
| `--window-size <SIZE>` | — | generate fixed windows instead, e.g. `50M` |
| `--type <auto\|3\|4\|5\|6\|8\|9\|12>` | `auto` | BED width; **BED input only** |
| `--format <auto\|bed\|gtf\|gff>` | `auto` | annotation format |
| `--feature <span\|exon\|intron\|cds\|utr\|five-utr\|three-utr>` | `span` | which segments participate |
| `--assignment <start\|midpoint\|contained\|best-overlap\|overlap>` | `start` | how a read is matched |
| `--name-field <STRING>` | — | annotation attribute to use as the key |
| `--unnamed-prefix <STRING>` | `region` | prefix for generated names |
| `--missing-feature <skip\|empty\|error>` | `skip` | records lacking the feature |
| `--require-blocks` | off | turn the BED span-as-exon fallback into an error |
| `--alignment-geometry <blocks\|span>` | `blocks` | what to compare against |
| `--emit-empty` | off | header-only BAMs for regions with no records |

Annotation parsing goes entirely through
[`genepred`](https://crates.io/crates/genepred) — `bamsplit` has no BED, GTF, or
GFF parser of its own. Format and BED width are detected from **content**, not
just the extension, and `.gz`, `.bgz`, `.zst`, and `.bz2` are recognized from
magic bytes so a gzip stream called `genes.bed` still works.

Each annotation record is **one** output key. `--feature` picks which of its
segments participate; it never writes one BAM per exon.

`--assignment overlap` is the only mode that can send a record to several
outputs, and the manifest's conservation equation switches accordingly.

Full semantics — feature derivation, tie-breaking, the spliced-read containment
rule, and the interval index — are in
[`docs/region-routing.md`](docs/region-routing.md).

### `bamsplit inspect`

```console
bamsplit inspect input.bam
bamsplit inspect input.bam --full
bamsplit inspect input.bam --by chrom
bamsplit inspect input.bam --tag CB --full
bamsplit inspect input.bam --regions genes.gtf --feature cds
bamsplit inspect input.bam --full --json
```

Writes no BAMs. Header-only inspection reports validity, header size, sort
order, references and lengths, read groups, seekability, BGZF end-of-file marker,
index availability and consistency, references beyond BAI's limits, reference
names that need filename encoding, the predicted output count, the engine the
planner would choose and why, and mmap eligibility.

`--full` scans every record and adds total/primary/secondary/supplementary,
mapped/placed-unmapped/unplaced-unmapped, duplicate and QC-fail counts,
per-reference counts, coordinate-order violations, requested tag cardinality,
malformed auxiliary sections, estimated temporary storage, and estimated
indexability.

With `--regions` it additionally runs the same detection, loading, and feature
extraction a split would, and reports the detected format and BED width, how
many records lack the information `--feature` needs and why, duplicate and
generated names, span-as-exon fallbacks, overlapping region pairs, and the
reference names each side has that the other does not — which is how a
`1`-versus-`chr1` build mismatch shows up before you run the split.

The report goes to **stdout**; logs go to stderr, so `--json | jq` works.

## Output layout

For every logical key the split produces:

```text
<out-dir>/<encoded-key>.bam
<out-dir>/<encoded-key>.bam.bai    (or .csi, when an index is built)
<out-dir>/bamsplit.manifest.json   (and/or .tsv)
```

Logical keys are **not** required to be filesystem-safe. The encoder maps them
into a reversible, printable-ASCII form:

| logical key | filename |
| --- | --- |
| `chr1` | `chr1.bam` |
| `chr1/alternate` | `chr1%2Falternate.bam` |
| `HLA-A*01:01` | `HLA-A%2A01%3A01.bam` |
| `..` | `%2E%2E.bam` |
| `CON` (a Windows device name) | `%43ON.bam` |
| a 4 000-byte key | `<prefix>--<16 hex digits>.bam` |

Bytes outside `A-Za-z0-9-_.+~` become `%XX`. `.`, `..`, Windows device names, and
trailing dots get their first (or last) byte escaped as well, so no output can be
`.`, `..`, `CON`, or `chr1.`. Only the length-capped case is lossy, and the
manifest always carries the complete logical value plus its hex encoding.

## Semantics

### Placed versus unplaced unmapped records

BAM has two different kinds of "unmapped":

* **Placed unmapped** — the `UNMAPPED` flag is set *and* `refID`/`pos` are valid.
  A coordinate-sorted BAM parks an unmapped read right next to its mapped mate
  this way.
* **Unplaced** — `refID == -1`. No coordinate at all.

By default placed-unmapped records **stay with their reference**
(`--placed-unmapped by-reference`) and only genuinely unplaced records go to
`unmapped.bam` (`--unplaced keep`). Sending placed-unmapped records to
`unmapped.bam` would split pairs across files for no benefit; `--placed-unmapped
unmapped` is available when you want it.

### Cross-chromosome mates

Mate fields are **never rewritten**. A read on `chr2` whose mate is on `chrX`
ends up in `chr2.bam` with `next_ref_id` still pointing at `chrX` — and that
reference still resolves, because every output keeps the full dictionary.

Two reads of one pair can therefore end up in different files. That is inherent
to splitting by chromosome, not a `bamsplit` choice. If you need pairs together,
use `bamsplit shard --key qname`.

### Why every output keeps the full header

`bamsplit` does not reduce an output header to the one chromosome it holds.
Doing so would force it to renumber every record's `ref_id` *and* `next_ref_id`,
which would:

* break raw-record transfer, the whole point of the fast path;
* silently invalidate mate coordinates for cross-chromosome pairs;
* surprise downstream tools that expect the original dictionary.

The cost is a few kilobytes per output. `samtools idxstats chr1.bam` therefore
lists every reference, most with zero records — that is expected.

### The `@PG` record

Unless `--no-pg` is given, each output gains

```text
@PG	ID:bamsplit	PN:bamsplit	VN:<version>	CL:<command>	PP:<previous>
```

`PP` links to the last leaf of the existing chain, so no existing `@PG` record is
rewritten. If `ID:bamsplit` is already present, `bamsplit.1`, `bamsplit.2`, … are
tried in order — deterministically for a given input header.

## Index behaviour

`--index auto` decides **before writing**, from the reference dictionary:

| longest reference | choice |
| --- | --- |
| ≤ 2<sup>29</sup> − 1 (536 870 911) | BAI |
| larger | CSI, at the smallest depth that covers it |

`--index bai` on an input with a longer reference is an **error**, not a silent
downgrade. `--index csi` always produces CSI. `--index none` produces none.

Indexes are built **while writing**, from the writer's own virtual offsets, so
there is no second pass over the finished file.

An output that is not coordinate-sorted gets **no index**, and the manifest
records why. A binning index over unsorted records would be actively misleading,
so `bamsplit` refuses to write one rather than producing something `samtools`
would silently mis-query.

## Engine selection

Three engines, one interface. The result does not depend on which one runs — the
test-suite asserts byte-identical outputs and equal digests across all three.

| engine | when | passes | memory | descriptors |
| --- | --- | --- | --- | --- |
| `stream` | output keys arrive grouped (coordinate-sorted `chrom`) | 1 | ~constant | 1 |
| `indexed` | seekable, indexed, few large references, threads > 1 | 1 logical | ~constant × workers | workers |
| `spool` | unsorted input, interleaved keys, high cardinality, stdin | 2 | ~constant | bounded |

`--engine auto` picks:

* **`chrom`**, coordinate-sorted → `stream`; additionally `indexed` when the input
  is a local indexed file with at most 64 references, is at least 8 MiB, and more
  than one thread is available.
* **`chrom`**, not coordinate-sorted → `spool`.
* **`shard`**, **`tag`** → `spool`, always: hashed and tag keys are interleaved by
  construction.
* **`region`** → `stream` when the annotation covers most of the genome in at
  most 256 regions and the input is coordinate-sorted; `spool` otherwise.
  Regions overlap and nest, so their keys are never contiguous — the streaming
  engine therefore holds every live output open for them, which only stays
  bounded while there are few.

The chosen engine and the reason are logged at `info` and written into the
manifest. An explicit `--engine` is **validated, never silently overridden**:
asking for something impossible is an error.

The streaming engine does not trust `@HD SO:coordinate`. It verifies monotonic
`(refID, pos)` and tracks the keys it has finalized; a reappearance means the
input is interleaved. It then removes every partial output and — under
`--engine auto`, or when the input is seekable — retries with the spool engine.
So a mislabelled BAM produces correct output and a note in the manifest, not a
directory of half-written files.

## Threads

`--threads` is a **total budget**, not a per-output allowance. A 3 000-contig
split with `--threads 16` uses 16 threads, not 48 000.

* **Streaming engine** — one routing thread, roughly a third of the remainder on
  BGZF inflation, the rest on compressing the single active output. Because only
  one output is open, nearly the whole budget lands on the writer, which is where
  a `chrom` split spends its time.
* **Indexed engine** — bounded per-reference tasks drawn from one shared pool.
* **Spool finalization** — bounded output tasks, lightly threaded compression
  each.

`--threads 1` is correct and deterministic. So is every other value: output is
byte-identical regardless of thread count, because block boundaries and
compression settings do not depend on scheduling. The one exception is the `@PG`
record, which deliberately embeds the command line — `--threads` included — so
byte-comparing two runs at different thread counts needs `--no-pg` on both
sides. Measured in [`docs/performance.md`](docs/performance.md), along with the
budget's poor rounding at `--threads 2`.

## Memory mapping

`--io mmap` is **not** a way to avoid decompression. BAM records inside a mapping
are still BGZF-compressed and still have to be inflated. What a mapping buys is
one fewer copy per block and cheap independent cursors, which is why `--io auto`
selects it only for the indexed engine on a local file.

Never used for stdin, pipes, non-regular files, or remote sources. A mapping
failure falls back to buffered I/O unless `--io mmap` was requested explicitly.
Performance depends on the filesystem and the access pattern: on network storage
a sequential buffered pass is usually faster.

## Temporary storage

The spool engine writes one plain, uncompressed, append-only file per routing key
under `--temp-dir`:

```text
magic "BSPL" | version | (u32 length, record body)* | sentinel | count | checksum
```

No repeated BAM headers, so temporary usage is close to the uncompressed record
volume — roughly what `bamsplit inspect --full` reports as *estimated temporary
storage*. Handles are cached with an LRU bounded by `--max-open-files`, so
descriptor use stays flat regardless of key count.

Spools are removed on success, on error, and on interruption, unless
`--keep-temp` is given. Every replay validates the framing, the record count,
and the checksum, so a spool corrupted between the two phases is caught rather
than silently dropping records.

## Manifest

```json
{
  "program": "bamsplit",
  "version": "0.1.0",
  "command": "bamsplit chrom sample.bam --out-dir chromosomes --threads 16",
  "start_time": "2025-08-01T09:14:02Z",
  "end_time": "2025-08-01T09:14:37Z",
  "elapsed_seconds": 35.41,
  "input_path": "sample.bam",
  "input_size": 4823901184,
  "input_header_checksum": "3f2a91c04d8e1b77",
  "input_index_type": "bai",
  "input_sort_order": "coordinate",
  "selected_engine": "stream",
  "engine_selection_reason": "the input is coordinate-sorted, so one sequential pass suffices (25 reference sequences, index present)",
  "io_backend": "buffered",
  "threads": 16,
  "compression_level": 6,
  "routing_mode": "chrom",
  "index_mode": "auto",
  "index_selection_reason": "the longest reference is 248956422 bases, within BAI's 536870911 limit",
  "temporary_bytes": 0,
  "input_records": 812445990,
  "unique_emitted_records": 812445990,
  "total_output_emissions": 812445990,
  "dropped_records": 0,
  "unmatched_records": 0,
  "duplicate_emissions": 0,
  "max_emissions_for_one_record": 1,
  "may_duplicate": false,
  "outputs": [
    {
      "logical_key": "chr1",
      "logical_key_hex": "63687231",
      "resolved_key": "chr1",
      "encoded_key": "chr1",
      "bam_path": "chr1.bam",
      "index_path": "chr1.bam.bai",
      "index_type": "bai",
      "skipped": false,
      "record_count": 63182044,
      "mapped_count": 62901337,
      "placed_unmapped_count": 280707,
      "unplaced_unmapped_count": 0,
      "primary_count": 62118904,
      "secondary_count": 402118,
      "supplementary_count": 661022,
      "compressed_bytes": 391048221,
      "first_coordinate": { "reference_id": 0, "position": 9999 },
      "last_coordinate": { "reference_id": 0, "position": 248946421 },
      "coordinate_sorted": true,
      "raw_record_digest": "b71e4c0a9f3d2856"
    }
  ]
}
```

Every run checks, before writing anything:

```text
input_records = unique_emitted_records + dropped_records + unmatched_records
total_output_emissions = unique_emitted_records + duplicate_emissions
```

A failure is exit code `5`. The manifest itself is written atomically, and an
invalid one is never written at all.

The TSV manifest is one row per output with the columns a workflow manager
needs, so a `[meta, chromosome, bam, index]` tuple can be built without parsing
JSON.

## Failure and cleanup guarantees

* Every output is written to `<name>.part.<pid>.<nonce>` and renamed only after
  the header, all records, the BGZF end-of-file marker, the index, the statistics,
  and validation are all done.
* An index is renamed **before** its BAM, so a crash between the two renames
  leaves an orphan index — harmless and obvious — rather than a BAM that looks
  finished but has no index.
* Any failure removes every temporary file *and* every output the run had already
  committed. You get a complete output directory or an empty one.
* `SIGINT`/`SIGTERM` set a flag the record loop polls; the run unwinds through
  ordinary cleanup and exits with code `6`. A process-wide sweep catches anything
  whose owner never got to run.
* An existing output is refused with exit code `4`. `--force` replaces it, still
  transactionally: the old file survives until the new one is complete.

## Pipeline examples

### Nextflow

```groovy
process SPLIT_BY_CHROM {
    input:  tuple val(meta), path(bam)
    output: tuple val(meta), path('*.chromosomes/bamsplit.manifest.tsv'), path('*.chromosomes/*')

    script:
    """
    bamsplit chrom ${bam} \\
        --out-dir ${meta.id}.chromosomes \\
        --emit-empty \\
        --index auto \\
        --threads ${task.cpus} \\
        --manifest tsv
    """
}

workflow {
    SPLIT_BY_CHROM(samples)
        // The TSV is directly a [meta, chromosome, bam, index] channel.
        .flatMap { meta, manifest, _files ->
            manifest.splitCsv(sep: '\t', header: true)
                    .findAll { !it.skipped.toBoolean() }
                    .collect { [meta, it.logical_key, file(it.bam), file(it.index)] }
        }
        | CALL_VARIANTS
}
```

`--emit-empty` matters here: it guarantees one output per reference, so the
number of downstream tasks does not depend on which chromosomes happened to have
reads.

### Snakemake

```python
import csv

checkpoint split_by_chrom:
    input:  "aligned/{sample}.bam"
    output: directory("chromosomes/{sample}")
    threads: 16
    shell:
        "bamsplit chrom {input} --out-dir {output} --emit-empty "
        "--threads {threads} --manifest tsv"

def chromosome_bams(wildcards):
    directory = checkpoints.split_by_chrom.get(**wildcards).output[0]
    with open(f"{directory}/bamsplit.manifest.tsv") as handle:
        return [
            f"{directory}/{row['bam']}"
            for row in csv.DictReader(handle, delimiter="\t")
            if row["skipped"] == "false"
        ]

rule call_all:
    input: chromosome_bams
```

### Shell scatter/gather

```bash
set -euo pipefail

bamsplit chrom sample.bam --out-dir chromosomes --threads 16 --manifest tsv

# Scatter: one job per non-empty chromosome.
tail -n +2 chromosomes/bamsplit.manifest.tsv \
  | awk -F'\t' '$17 == "false" { print $1, $4 }' \
  | while read -r chromosome bam; do
        call_variants "chromosomes/$bam" > "calls/$chromosome.vcf" &
    done
wait

# Gather.
bcftools concat calls/*.vcf -Oz -o merged.vcf.gz
```

### Rust library

```rust
use bamsplit_core::{ChromOptions, RunOptions, split_by_chromosome};

let report = split_by_chromosome(
    "sample.bam",
    &ChromOptions {
        emit_empty: true,
        ..ChromOptions::new("chromosomes")
    },
    &RunOptions::default().with_threads(16),
)?;

for output in &report.manifest.outputs {
    println!(
        "{}: {} records, digest {}",
        output.logical_key, output.stats.record_count, output.stats.raw_record_digest
    );
}
# Ok::<_, bamsplit_core::error::Error>(())
```

A runnable version is in [`examples/library_split.rs`](examples/library_split.rs).

## Documentation

| document | contents |
| --- | --- |
| [`docs/architecture.md`](docs/architecture.md) | layers, data flow, and the reasoning behind each boundary |
| [`docs/bam-semantics.md`](docs/bam-semantics.md) | exact routing rules, header handling, mate behaviour |
| [`docs/region-routing.md`](docs/region-routing.md) | region semantics in full: features, assignment, the interval index |
| [`docs/performance.md`](docs/performance.md) | what is fast, what is not, and how to measure it |
| [`docs/pipeline-examples.md`](docs/pipeline-examples.md) | longer workflow-manager recipes |
| [`CHANGELOG.md`](CHANGELOG.md) | release history |
| [`LICENSING.md`](LICENSING.md) | licensing, and a known header discrepancy |

## Development

```console
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
cargo test --doc --workspace --all-features
cargo build --workspace --all-features --release

# Differential tests against samtools; skipped if samtools is absent.
cargo test --workspace --features differential

# Fixtures on disk, for reproducing a bug by hand.
cargo run --bin bamsplit-generate-fixtures -- tests/data/generated

# Benchmarks.
cargo bench
scripts/benchmark.sh sample.bam
scripts/differential-test.sh sample.bam
```

There is one `unsafe` block in the workspace — the `memmap2` call behind
`--io mmap` — isolated in `engine::unsafe_map` with a `SAFETY` comment. The rest
is `deny(unsafe_code)`.

## License

GPL-3.0-only. See [`LICENSE`](LICENSE) and the discrepancy note in
[`LICENSING.md`](LICENSING.md).
