# Pipeline examples

Longer recipes than the README's. All use the TSV manifest, because parsing JSON
inside a workflow manager is more trouble than it is worth.

## The TSV manifest

```console
bamsplit chrom sample.bam --out-dir chromosomes --manifest tsv
```

```text
logical_key  resolved_key  encoded_key  bam        index          index_type  records  ...  skipped
chr1         chr1          chr1         chr1.bam   chr1.bam.bai   bai         63182044      false
chr2         chr2          chr2         chr2.bam   chr2.bam.bai   bai         61402118      false
...
unmapped     unmapped      unmapped     unmapped.bam unmapped.bam.bai bai      1200481      false
```

Columns: `logical_key`, `resolved_key`, `encoded_key`, `bam`, `index`,
`index_type`, `records`, `mapped`, `placed_unmapped`, `unplaced_unmapped`,
`primary`, `secondary`, `supplementary`, `compressed_bytes`, `coordinate_sorted`,
`raw_record_digest`, `skipped`.

`bam` and `index` are **relative to the output directory**, so a workflow can
stage the directory and resolve them without absolute paths.

`skipped` marks a key that was listed but not created — a reference with no
records, when `--emit-empty` was not given. Filter on it.

## Nextflow

### Scatter–gather over chromosomes

```groovy
process BAMSPLIT_CHROM {
    tag "${meta.id}"
    cpus 16

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.chromosomes"), emit: directory
    path "${meta.id}.chromosomes/bamsplit.manifest.tsv",  emit: manifest

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

workflow SCATTER_BY_CHROM {
    take: samples
    main:
        BAMSPLIT_CHROM(samples)

        shards = BAMSPLIT_CHROM.out.directory
            .flatMap { meta, directory ->
                file("${directory}/bamsplit.manifest.tsv")
                    .splitCsv(sep: '\t', header: true)
                    .findAll { row -> row.skipped == 'false' }
                    .collect { row ->
                        [ meta + [chromosome: row.logical_key],
                          file("${directory}/${row.bam}"),
                          file("${directory}/${row.index}") ]
                    }
            }

        CALL_VARIANTS(shards)
        MERGE(CALL_VARIANTS.out.groupTuple(by: 0))
    emit:
        MERGE.out
}
```

`--emit-empty` guarantees one output per reference, so the number of downstream
tasks is fixed by the reference dictionary rather than by which chromosomes
happened to have reads. Without it, a sample with no `chrY` coverage silently
produces one fewer task, which is exactly the kind of thing that goes unnoticed
until a merge step has a hole in it.

### Load-balanced shards instead of chromosomes

Chromosomes are wildly unequal — `chr1` is 25× `chr21`. When the work is per-read
rather than per-locus, shards balance better:

```groovy
process BAMSPLIT_SHARD {
    cpus 8
    input:  tuple val(meta), path(bam)
    output: tuple val(meta), path("shards/shard-*.bam")
    script:
    """
    bamsplit shard ${bam} \\
        --shards ${params.shards} \\
        --key qname \\
        --out-dir shards \\
        --threads ${task.cpus} \\
        --index none
    """
}
```

`--key qname` keeps both mates and every secondary and supplementary alignment
of a template in one shard, so per-template work is still correct.
`--index none` because shards are not coordinate-sorted and an index would be
refused anyway.

### Per-sample demultiplexing

```groovy
process BAMSPLIT_BY_SAMPLE {
    input:  path(merged_bam)
    output: path("by-sample/*.bam")
    script:
    """
    bamsplit tag ${merged_bam} \\
        --field sample \\
        --out-dir by-sample \\
        --unknown-read-group error \\
        --threads ${task.cpus}
    """
}
```

`--unknown-read-group error` is the right choice here: a record whose `RG` is not
in the header means the merge was wrong, and silently filing it under `no-tag`
would hide that.

## Snakemake

```python
import csv
from pathlib import Path

checkpoint bamsplit_chrom:
    input:
        bam="aligned/{sample}.bam",
    output:
        directory("chromosomes/{sample}"),
    threads: 16
    shell:
        "bamsplit chrom {input.bam} "
        "--out-dir {output} "
        "--emit-empty --index auto --manifest tsv "
        "--threads {threads}"


def manifest_rows(wildcards):
    directory = Path(checkpoints.bamsplit_chrom.get(**wildcards).output[0])
    with open(directory / "bamsplit.manifest.tsv") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["skipped"] == "false":
                yield directory, row


def chromosome_calls(wildcards):
    return [
        f"calls/{wildcards.sample}/{row['logical_key']}.vcf.gz"
        for _, row in manifest_rows(wildcards)
    ]


rule call_one_chromosome:
    input:
        bam="chromosomes/{sample}/{chromosome}.bam",
        bai="chromosomes/{sample}/{chromosome}.bam.bai",
    output:
        "calls/{sample}/{chromosome}.vcf.gz",
    shell:
        "call_variants {input.bam} | bgzip > {output}"


rule merge_calls:
    input:
        chromosome_calls,
    output:
        "calls/{sample}.vcf.gz",
    shell:
        "bcftools concat {input} -Oz -o {output}"
```

The manifest is what makes the checkpoint work: Snakemake needs to know the
output set after the fact, and reading a TSV is more robust than globbing a
directory whose names have been percent-encoded.

## Shell

### Scatter with a job limit

```bash
#!/usr/bin/env bash
set -euo pipefail

bamsplit chrom sample.bam \
    --out-dir chromosomes \
    --threads 16 \
    --manifest tsv

mkdir -p calls
tail -n +2 chromosomes/bamsplit.manifest.tsv \
  | awk -F'\t' '$17 == "false" { print $1 "\t" $4 }' \
  | xargs -P 8 -n 2 bash -c '
        call_variants "chromosomes/$1" > "calls/$0.vcf"
    '

bcftools concat calls/*.vcf -Oz -o merged.vcf.gz
```

### Verifying a split before trusting it

```bash
bamsplit chrom sample.bam --out-dir chromosomes --manifest json

python3 - <<'PY'
import json
m = json.load(open("chromosomes/bamsplit.manifest.json"))
lhs = m["input_records"]
rhs = m["unique_emitted_records"] + m["dropped_records"] + m["unmatched_records"]
assert lhs == rhs, f"records lost: {lhs} != {rhs}"
print(f"{lhs} records conserved across {len(m['outputs'])} outputs")
for output in m["outputs"]:
    if not output["skipped"]:
        print(f"  {output['logical_key']:<20} {output['record_count']:>12}  {output['raw_record_digest']}")
PY
```

`bamsplit` already performs this check and exits with code 5 on failure, so the
script is belt-and-braces — but the digests are useful to record alongside the
outputs for later comparison.

### Working from a pipe

```bash
samtools view -b -q 30 sample.bam | bamsplit chrom - --out-dir filtered
```

stdin is not seekable, so the planner picks the spool engine and `--io mmap` is
refused. Output is identical either way.

## Region routing

### Per-transcript coverage

```bash
bamsplit region rnaseq.bam \
    --regions gencode.v44.annotation.gtf.gz \
    --feature exon \
    --assignment overlap \
    --out-dir by-transcript \
    --manifest tsv \
    --threads 8
```

`--assignment overlap` is the only mode that duplicates records, which is what
you want for coverage: a read spanning two transcripts is evidence for both. The
manifest's `duplicate_emissions` says how much duplication that produced.

### Per-CDS, without double-counting

```bash
bamsplit region rnaseq.bam \
    --regions gencode.v44.annotation.gtf.gz \
    --feature cds \
    --assignment best-overlap \
    --out-dir by-cds \
    --threads 8
```

`best-overlap` assigns each read exactly once, to the CDS it shares the most
aligned bases with — so summing the outputs gives back the input.

### Fixed windows for a scatter

```bash
bamsplit region sample.bam \
    --window-size 10M \
    --assignment start \
    --emit-empty \
    --out-dir windows \
    --manifest tsv
```

Windows beat chromosomes when the work is per-locus and `chr1` would otherwise
be 25× the size of `chr21`. `--emit-empty` fixes the task count.

### Checking the annotation first

```bash
bamsplit inspect sample.bam --regions gencode.v44.annotation.gtf.gz --feature cds
```

This runs the same detection and extraction a split would, so it tells you
before you commit an hour of compute: how many records have no CDS, whether the
reference names match (`1` versus `chr1` is the classic failure), and how many
regions overlap.

## Rust

See [`examples/library_split.rs`](../examples/library_split.rs) for a runnable
program. The short version:

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

// The manifest is already validated; a failed conservation check would have
// been returned as an error rather than reaching here.
for output in report.manifest.outputs.iter().filter(|o| !o.skipped) {
    println!("{} -> {}", output.logical_key, output.bam_path);
}
# Ok::<_, bamsplit_core::error::Error>(())
```
