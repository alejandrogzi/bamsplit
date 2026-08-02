# Architecture

## Layers

```text
┌──────────────────────────────────────────────────────────────────────┐
│ bamsplit-cli        argument parsing, logging, exit codes, rendering │
└───────────────────────────────┬──────────────────────────────────────┘
                                │  option structs only
┌───────────────────────────────▼──────────────────────────────────────┐
│ bamsplit-core::api  one function per subcommand                      │
├──────────────────────────────────────────────────────────────────────┤
│ planning            which engine, which I/O backend, and why         │
├──────────────────────────────────────────────────────────────────────┤
│ engine              stream · indexed · spool                         │
│   ├── routing       where a record goes (never touches a file)       │
│   ├── output        naming, BGZF, index, transaction                 │
│   └── manifest      what happened, and proof nothing was lost        │
├──────────────────────────────────────────────────────────────────────┤
│ bam                 header · raw records · CIGAR · tags              │
│ annotation          detect · load · features · interval index        │
└──────────────────────────────────────────────────────────────────────┘
```

Every boundary exists for a reason, and the reasons are worth stating because
they are what keep the fast path fast.

### `bam` — bytes to borrowed views

`RawRecord<'a>` is a `&[u8]` and nothing else. Fields are parsed on demand, so a
router that only needs `refID` costs one bounds-checked 4-byte read per record
rather than a full decode into an owned structure.

The reader owns one growable buffer and reuses it, so a pass over a billion-record
BAM performs a constant number of heap allocations.

No `unsafe`, no unchecked indexing on input, no unaligned reads: multi-byte
scalars are rebuilt with `from_le_bytes`. Every length is validated before use and
every offset computation is checked.

### `routing` — where, and only where

A `Router` borrows a record and returns a `Route`. It cannot open a file because
it is never given one. That is what lets the same router run under three engines
unchanged, and what makes the routing tests pure functions over synthetic
records.

`Route::Many` exists for exactly one documented case — region routing in
`overlap` mode. Everything else returns `One` or `Drop`, and the manifest's
conservation check enforces it.

### `annotation` — files to logical regions

Four stages, none of which touches a BAM: detect the format, load through
`genepred`, derive segments for the requested feature, and index them per
reference. The output is an interval index the region router queries. Keeping it
BAM-free is what lets the feature semantics be tested as pure functions over
`GenePred` values.

### `output` — naming, writing, committing

Three concerns, deliberately separated:

* `filename` turns arbitrary bytes into a reversible, traversal-safe stem.
* `bgzf` compresses blocks, in parallel, while still reporting exact virtual
  offsets.
* `transaction` writes to a temporary path and renames only on success.

`writer` joins them into one output; `manager` bounds how many are live.

### `engine` — how the pass runs

The engines differ only in *scheduling*: which records are read when, and how
many outputs are open at once. They share the router, the output manager, and the
manifest, so a run's result is identical whichever one runs. `tests/integration`
asserts that directly, on the bytes.

## Data flow, `bamsplit chrom` on a coordinate-sorted BAM

```text
   file ──▶ BGZF inflate ──▶ RawRecordReader ──▶ RawRecord (borrowed)
                                                     │
                                          ChromRouter │ reads refID only
                                                     ▼
                                            Route::One("chr1")
                                                     │
                            key changed? ────────────┤
                                 │ yes               │ no
                                 ▼                   ▼
                    finish_one("chr1")      OutputManager::write
                    ├─ flush BGZF                    │
                    ├─ EOF marker           ┌────────▼────────┐
                    ├─ build BAI            │ BgzfBlockWriter │
                    ├─ validate             │  block N ──┐    │
                    └─ rename into place    └────────────┼────┘
                                                         ▼
                                            resolved virtual offsets
                                                         │
                                                         ▼
                                                  IndexBuilder
```

One pass, one active output, constant memory.

## The parallel-BGZF-with-offsets problem

Streaming index construction needs a virtual offset per record. Parallel
compression makes a block's compressed size unknown until it comes back from a
worker. Those two requirements conflict, and `noodles` resolves it by offering
two writers — one that reports offsets and one that is parallel.

`bamsplit` needs both, so `output::bgzf` splits the record's offset in half:

* the **offset within the block** is known the instant the record is staged;
* the **block's compressed base** is known once every preceding block is written.

Records are queued with their `(block index, offset)` slot, and resolved in write
order as batches complete. The queue holds at most one batch, so memory is bounded
at roughly `workers × 64 KiB` regardless of file size. With one worker, resolution
happens immediately and the machinery costs nothing.

A subtlety worth recording: the queue entry is reserved *before* the record's
bytes are appended. Appending can seal a block, and sealing can discard the
block-offset table entry the record's own start needs. Reserving first is what
keeps a record that straddles a block boundary resolvable.

## Parking

When more outputs are live than `--max-open-files` allows, the least recently used
one is *parked*: its current BGZF block is flushed, its handle is closed, and its
compressed offset is remembered. Resuming reopens in append mode with that base
offset, so the index continues without a discontinuity.

This works only because the end-of-file marker is written at finish rather than on
drop. A BGZF stream is a concatenation of independent blocks, so appending to an
unterminated file is well-defined.

## Errors

Leaf error per subsystem, aggregated by `EngineError`, unified by `Error`, each
carrying an `ExitCode`. The CLI never re-classifies: by the time a failure reaches
`main`, it already knows what to exit with.

Malformed record errors carry a `RecordLocation` — the ordinal a user counts to in
`samtools view`, and the BGZF virtual offset a BGZF-aware tool seeks to. Because
`RawRecord` holds only a slice, it cannot know where it came from; the reader and
the engine fill that in with `relocate`.

## What is not here

* **A routing-expression language.** Deliberately out of scope for the first
  release. The four routers cover the cases that matter, and a mini-language
  would need its own parser, its own errors, and its own test-suite.
* **A BED/GTF/GFF parser.** Everything goes through `genepred`. Those formats
  have enough dialects that a second implementation is a second set of bugs.
* **CRAM.** Different container, different reference-resolution problem.
* **Sorting.** `bamsplit` never reorders records. `samtools sort` exists.
* **Index queries for `shard` and `tag`.** Their keys are a function of the
  record, not of its position, so they cannot be decomposed into index queries.
  `--engine indexed` with either is refused rather than approximated.
