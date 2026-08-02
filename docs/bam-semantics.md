# BAM semantics

Exactly what `bamsplit` does with each record, and why.

## Routing rules for `chrom`

| record | destination |
| --- | --- |
| valid `refID` | that reference |
| `refID == -1` | the unplaced output, or dropped, or an error |
| `UNMAPPED` set **and** valid `refID` | that reference (default), or the unplaced output |
| `refID` outside the dictionary | **input error**, exit code 3 |
| `next_ref_id` outside the dictionary | **input error**, exit code 3 |
| secondary (`0x100`) | its own `refID`, never its primary's |
| supplementary (`0x800`) | its own `refID`, never its primary's |
| duplicate (`0x400`), QC-fail (`0x200`) | routed normally; flags are never filters |

`bamsplit` is a splitter, not a filter. A duplicate-flagged record is written
where it belongs, exactly as it arrived.

## Placed unmapped versus unplaced

BAM has two distinct notions of "unmapped", and conflating them loses data:

```text
FLAG 0x4 set, refID = 0, pos = 12345    placed unmapped   → chr1.bam
FLAG 0x4 set, refID = -1, pos = -1      unplaced          → unmapped.bam
```

A coordinate-sorted BAM stores the unmapped mate of a mapped read the first way,
so it sorts next to its mate. Defaults:

```text
--placed-unmapped by-reference
--unplaced        keep
```

so pairs stay together and only records with no coordinate at all go to
`unmapped.bam`.

`--placed-unmapped unmapped` moves them, at the cost of separating those pairs.

## Mate fields are never rewritten

`next_ref_id` and `next_pos` are copied verbatim. A `chr2` read whose mate is on
`chrX` keeps pointing at `chrX`, and that pointer stays valid because the output
retains the whole reference dictionary.

The consequence is that a pair can be split across two files. That is inherent to
splitting by chromosome. `bamsplit shard --key qname` keeps templates together.

## Header handling

### Parsed and preserved

* BAM magic `BAM\1`
* the textual SAM header (`@HD`, `@SQ`, `@RG`, `@PG`, `@CO`)
* the binary reference dictionary

The **binary dictionary is authoritative** for reference identifiers: a record's
`ref_id` indexes it. The text is what humans and most SAM tooling read.

### Validated

* reference names are non-empty and free of tabs, newlines, and NUL
* reference names are unique
* reference lengths are positive
* `@SQ` lines, when present, agree with the binary dictionary in count, name, and
  length
* read-group identifiers are unique
* program identifiers are unique and every `PP` resolves
* `@HD SO` is one of `unknown`, `unsorted`, `queryname`, `coordinate`

A header with **no** `@SQ` lines is legal and common — the binary dictionary
stands alone — so the cross-check is skipped in that case rather than failing.

### Why outputs keep the full dictionary

Reducing an output header to its one chromosome would require renumbering
`ref_id` *and* `next_ref_id` in every record, which would:

* defeat raw-record transfer, since every record body would change;
* invalidate mate coordinates for cross-chromosome pairs;
* break tools that expect the original dictionary.

`samtools idxstats chr1.bam` therefore lists every reference, most with zero
records. That is correct.

### `@PG`

```text
@PG	ID:bamsplit	PN:bamsplit	VN:<version>	CL:<command>	PP:<last leaf>
```

`PP` points at the last `@PG` record that no other record names as its `PP`, so
the chain extends rather than being rewritten. A collision on `ID:bamsplit`
produces `bamsplit.1`, `bamsplit.2`, … deterministically.

`--no-pg` suppresses the record entirely.

## Record layout

```text
block_size  u32 LE               framing, not part of the body
├─ 0..4     ref_id      i32 LE   -1 = unplaced
├─ 4..8     pos         i32 LE   -1 = absent, otherwise 0-based
├─ 8        l_read_name u8       includes the NUL terminator
├─ 9        mapq        u8       255 = unavailable
├─ 10..12   bin         u16 LE   advisory; bamsplit never reads it
├─ 12..14   n_cigar_op  u16 LE
├─ 14..16   flag        u16 LE
├─ 16..20   l_seq       i32 LE
├─ 20..24   next_ref_id i32 LE
├─ 24..28   next_pos    i32 LE
├─ 28..32   tlen        i32 LE
├─ read_name   l_read_name bytes, NUL-terminated
├─ cigar       n_cigar_op * 4 bytes
├─ seq         (l_seq + 1) / 2 bytes, 4-bit packed
├─ qual        l_seq bytes
└─ data        the remainder
```

`bamsplit` preserves the **body** — everything after `block_size` — byte for
byte. `block_size` itself is recomputed from the body length, because a prefix
that disagreed with its body would be malformed no matter what the input said.

`bin` is never read. It is advisory, frequently wrong in the wild, and every
reader that cares recomputes it.

## CIGAR geometry

Two views, both derived in one pass:

* **reference span** — every reference-consuming operation: `M`, `D`, `N`, `=`,
  `X`. The outer footprint.
* **aligned blocks** — maximal runs of `M`, `=`, `X`. Where read bases actually
  sit on the reference.

```text
CIGAR   10S 20M 1000N 30M 4D 10M 5S
        ····▓▓▓▓─────▓▓▓▓▓▓···▓▓▓▓····
            └──┘     └────┘   └──┘      aligned blocks
            └──────────────────────┘    reference span (20+1000+30+4+10)
```

`N` separates blocks: it is the spliced-out intron of an RNA-seq alignment and
must never count as exonic overlap.

**A deletion contributes to the reference span but not to aligned-base overlap.**
So `20M 4D 10M` yields two blocks and a span of 34. This keeps overlap scoring
honest — a read is credited only for bases it actually aligns. The alternative
(htslib-style block merging) exists as `DeletionPolicy::IncludeInBlocks` in the
library but is not selected by the CLI.

`I`, `S`, `H`, and `P` consume no reference bases and appear in neither view.

An alignment whose CIGAR consumes no reference bases still occupies one position,
matching `bam_endpos` in htslib, so `end > start` holds for every placed record.

## Long CIGARs

SAM §4.2.2: a record with more than 65 535 CIGAR operations stores a
two-operation placeholder inline —

```text
<l_seq>S <reference span>N
```

— and the real operations in a `CG:B:I` tag. `bamsplit` resolves this
transparently, and distinguishes it from a genuine two-operation CIGAR by
requiring that the soft clip's length equals `l_seq` *and* that a `CG:B:I` tag is
present. A record that signals the convention without the tag is an error, not a
silent mis-parse.

## Auxiliary data

Every type is decoded: `A c C s S i I f Z H B`.

For `B` arrays the subtype is validated, the element count is checked for
multiplication overflow, and the whole payload is bounds-checked before any
element is read.

Tag values used as routing keys are rendered deterministically:

| type | rendering |
| --- | --- |
| `A` | the character |
| `c C s S i I` | decimal, `-` for negatives |
| `f` | shortest round-tripping decimal; `NaN`, `inf`, `-inf` |
| `Z` `H` | the raw bytes |
| `B` | `<subtype>,<v0>,<v1>,…`, order-sensitive |

Floating-point values never pass through locale-sensitive formatting.

## Malformed input

Rejected with exit code 3, always naming the record ordinal and BGZF virtual
offset:

* negative, undersized, or oversized `block_size`
* a truncated record or a truncated length prefix
* a zero-length or unterminated read name
* a negative `l_seq`
* a CIGAR, sequence, or quality section that overruns the body
* a `refID` or `next_ref_id` outside the dictionary
* a malformed auxiliary section — but only when something actually reads it

That last point is deliberate: `chrom` routing never looks at auxiliary data, so
a record with a broken tag section is transferred untouched. `bamsplit tag --tag
XX` does read it, and then it is an error. Failing on data nobody consulted would
reject files that are perfectly usable for the requested operation.

A **missing BGZF end-of-file marker** is not fatal. Every record is still there,
so `bamsplit` reads them all — and writes outputs that do carry the marker.
`bamsplit inspect` reports the absence.
