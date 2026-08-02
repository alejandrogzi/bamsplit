# Region routing

One output per annotated region, or per generated window.

All four engines serve region routing. The indexed engine runs a *per-region*
query plan: each region's envelope becomes an index query, the resulting chunks
are merged per reference, and only those chunks are read. A sparse annotation
over a large sorted BAM therefore never touches most of the input. See
[Engine selection](#engine-selection).

## Interface

```console
bamsplit region input.bam --regions annotation.bed12 --feature span     --assignment start        --out-dir regions
bamsplit region input.bam --regions annotation.bed12 --feature exon     --assignment overlap      --out-dir transcripts
bamsplit region input.bam --regions annotation.gtf.gz --feature cds     --assignment best-overlap --out-dir cds
bamsplit region input.bam --regions annotation.gff3  --feature intron   --assignment overlap      --out-dir introns
bamsplit region input.bam --window-size 50M          --assignment start --out-dir windows
```

Exactly one of `--regions` and `--window-size` is required.

| option | default | meaning |
| --- | --- | --- |
| `--regions <PATH>` | — | a BED, GTF, or GFF annotation |
| `--type <auto\|3\|4\|5\|6\|8\|9\|12>` | `auto` | BED width; **BED input only** |
| `--format <auto\|bed\|gtf\|gff>` | `auto` | annotation format |
| `--feature <span\|exon\|intron\|cds\|utr\|five-utr\|three-utr>` | `span` | which segments participate |
| `--assignment <start\|midpoint\|contained\|best-overlap\|overlap>` | `start` | how a read is assigned |
| `--name-field <STRING>` | — | annotation attribute to use as the key |
| `--unnamed-prefix <STRING>` | `region` | prefix for generated names |
| `--missing-feature <skip\|empty\|error>` | `skip` | records lacking the requested feature |
| `--require-blocks` | off | turn the BED span-as-exon fallback into an error |
| `--alignment-geometry <blocks\|span>` | `blocks` | what to compare against |
| `--window-size <SIZE>` | — | generate windows instead, e.g. `50M` |

## Annotation parsing

All parsing goes through the [`genepred`](https://crates.io/crates/genepred)
crate. `bamsplit` does **not** contain a BED, GTF, or GFF parser, and should not
grow one: those formats have enough dialects that a second implementation is a
second set of bugs.

### Format detection

1. an explicit `--format`;
2. the extension, after stripping `.gz`, `.bgz`, `.zst`, `.bz2` — recognizing
   `.bed`, `.gtf`, `.gff`, `.gff3`;
3. content inspection of the first non-empty, non-comment records;
4. otherwise a clear error naming what was seen.

Extension alone is never sufficient. When the extension and the content disagree,
**content wins** and the manifest records that it did — a `.bed` file containing
GTF is more common than anyone would like.

### BED width detection

With `--type auto`: count the columns of the first valid data row, choose the
highest supported width consistent with it (3, 4, 5, 6, 8, 9, 12), and confirm
every later row matches. Extra columns beyond the chosen width become additional
fields. A differing count is an error naming the line.

The chosen width dispatches to a real typed reader — `Reader<Bed3>` through
`Reader<Bed12>` — in a `match`. No type erasure, no `Box<dyn Any>`.

GTF and GFF go through the crate's aggregating path, which produces
transcript-like records with derived exon structure.

## Features

Each annotation record is **one** output key. `--feature` selects which of its
genomic segments participate in overlap testing; it never creates one BAM per
exon.

| feature | segments | requires |
| --- | --- | --- |
| `span` | `start..end` | nothing |
| `exon` | block intervals | blocks, or the BED fallback |
| `intron` | gaps between consecutive exons | ≥ 2 ordered, non-overlapping blocks |
| `cds` | the coding portions of exons | coding bounds **and** blocks |
| `utr` | all untranslated exonic segments | coding bounds and blocks |
| `five-utr` | strand-aware 5′ exonic UTR | strand, coding bounds, blocks |
| `three-utr` | strand-aware 3′ exonic UTR | strand, coding bounds, blocks |

Notes that matter:

* `cds` uses the coding portions of **exons**, not the whole thick span. The
  difference is every intron inside the coding region.
* `utr` is exonic only. Intronic sequence outside the CDS is not UTR.
* For a BED3–BED9 record with no block structure, `exon` falls back to treating
  the whole span as one exon, and the manifest counts the fallback.
  `--require-blocks` makes it an error. For aggregated GTF/GFF, a record with no
  derived exon structure is always an error — the aggregation was supposed to
  supply one.
* A record with no introns produces an empty logical region: no output by
  default, a header-only BAM under `--emit-empty`.
* `--missing-feature` decides what happens to a record lacking the requested
  information: `skip` (default, counted), `empty`, or `error`.

## Naming

In order: `--name-field` when present, then the record's name, then a generated

```text
<unnamed-prefix>-<chrom>-<start>-<end>-<ordinal>
```

Logical names need not be filesystem-safe — the filename encoder handles that.
Duplicates do not overwrite: they resolve to `<name>`, `<name>.2`, `<name>.3`, …
in annotation input order, and the manifest records both the original name and
the resolved key.

## Assignment modes

Feature selection and assignment are **independent**. `--feature exon
--assignment overlap` compares alignment blocks against the union of a
transcript's exons.

| mode | rule | emissions |
| --- | --- | --- |
| `start` | the alignment's leftmost reference position | ≤ 1 |
| `midpoint` | the midpoint of the alignment span | ≤ 1 |
| `contained` | the whole alignment is inside the feature union | ≤ 1 |
| `best-overlap` | most shared reference bases | 1 |
| `overlap` | every feature with a qualifying overlap | ≥ 0, **may duplicate** |

Tie-breaking, applied in order:

* `start` and `midpoint` — the smallest containing feature, then the earliest
  annotation input order.
* `best-overlap` — greatest overlap, then greatest overlap *fraction* relative to
  the aligned reference-consuming length, then the smallest feature span, then
  the earliest input order.

Every tie-break is counted as an ambiguous assignment in the manifest.

`contained`, for a discontinuous feature union such as exons:

* every reference-consuming **aligned block** must be inside the selected
  segments;
* skipped-reference (`N`) operations need **not** be — a spliced read whose exonic
  blocks all fall inside a transcript's exons is contained, even though its
  introns are not.

`overlap` is the only mode that duplicates records, and the manifest's
conservation equation switches accordingly:

```text
total_output_emissions = unique_emitted_records + duplicate_emissions
```

## Geometry

`--alignment-geometry blocks` (the default) compares annotation segments against
the alignment's `M`/`=`/`X` runs. `span` uses the outer footprint instead.

`N` separates blocks and never counts as exonic overlap. A deletion contributes
to the reference span but **not** to aligned-base overlap, and `best-overlap`
uses that consistently. See `docs/bam-semantics.md` for the diagram.

## Engine selection

`--engine auto` picks, in order:

* **`stream`** when the annotation covers at least a quarter of the reference
  dictionary in at most 256 regions *and* the input is coordinate-sorted.
* **`indexed`** when the annotation is *sparse* — fewer than a quarter of the
  references, at most 8192 regions — and the input is an indexed, seekable local
  file. This is the case where skipping matters: the queries reach only the
  chunks the regions fall in.
* **`spool`** otherwise.

### The per-region query plan

The indexed engine's default decomposition is per *reference*, which a region
split cannot use directly: a region is a sub-interval, and several regions share
a reference. So `region` builds a `QueryTargets` map — reference id to the list
of region envelopes on it — and the engine plans from that instead:

1. Each envelope becomes one index query.
2. The resulting chunks are collected per reference, sorted, and merged, so a
   chunk shared by two overlapping regions is read once.
3. One task per reference reads only its merged chunks. Because a region names
   exactly one reference, every output is written by exactly one task, and no
   cross-task coordination is needed.

Querying the *envelope* is sufficient for every assignment mode. A mode can only
assign a record to a region if the record touches one of that region's segments,
and every segment lies inside the envelope; a binning index returns every record
intersecting the queried interval. The router then applies the mode as usual, so
the outputs are byte-identical to the other engines'.

The unplaced tail is not scanned under this plan: an unplaced record has no
position, so no region can claim it.

One consequence shows up in the manifest. `input_records` counts the records
*examined*, not the records in the file, so it is smaller than the other engines
report on the same input. Conservation still holds among what was read, and the
report carries a note saying so.

`shard` and `tag` remain refused under `--engine indexed`: their keys are a
function of the record, not of its position, so they cannot be decomposed into
index queries at all.

Measured on a 2 M-record BAM with 200 regions of 10 kb — 7.9% of the records —
the plan read 157,098 records instead of 2,000,000 and finished in 0.15 s
against the streaming engine's 1.61 s, with byte-identical output. Full numbers
and caveats in [performance.md](performance.md).

### Why `stream` is capped at 256 regions

Regions overlap and nest, so a coordinate-sorted BAM does not visit their keys
in contiguous runs. The streaming engine therefore does **not** finalize a
region output when its key stops appearing — it holds every live output open,
bounded by `--max-open-files`, and finalizes at the end. That stays bounded only
while there are few regions, which is what the 256 limit encodes. A
whole-transcriptome annotation goes to the spool engine, where each output's
index is built and released one at a time.

## Interval index

Each annotation record becomes:

```rust
struct LogicalRegion {
    key: RegionKey,
    chrom: Vec<u8>,
    segments: Vec<Interval>,   // merged where adjacent or overlapping
    envelope: Interval,
    source_ordinal: u64,
    metadata: RegionMetadata,
}
```

Per reference, a sorted array of envelopes with a running prefix-maximum of
envelope ends. A query binary-searches the start bound and walks back only while
`max_end > query.start`. Candidates are then verified against the actual
segments.

Chosen over a centered interval tree for build cost, cache locality, and zero
pointer chasing. It degrades when one region's envelope spans the whole
reference — a full-length transcript makes every query walk past it — which is
the honest trade-off, and it is why the module documents the weakness rather
than hiding it. Correctness is unaffected: a randomized test compares every
query against a full linear scan.

An envelope covers a transcript's introns as well as its exons, so a hit is only
a *candidate*. The router then verifies against the actual segments.

## Manifest fields

Global: detected format, detected or explicit BED type, feature type, assignment
mode, annotation record count, derived segment count, records lacking the
requested feature, generated names, duplicate-name resolutions, ambiguous
assignments, unmatched BAM records, duplicate emissions, span-as-exon fallbacks.

Per output: annotation ordinal, derived segment count, span-as-exon fallback,
original logical name.
