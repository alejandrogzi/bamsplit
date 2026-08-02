#!/usr/bin/env bash
# Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
# Distributed under the terms of the Apache License, Version 2.0.
#
# Compares `bamsplit chrom` against `samtools view` on a real BAM.
#
# The Rust differential suite (`cargo test --features differential`) covers the
# synthetic fixtures. This script exists for the other half of the job: pointing
# the same comparison at an actual sequencing BAM, which is the only way to find
# the record shapes nobody thought to synthesize.
#
#   scripts/differential-test.sh sample.bam [threads]
set -euo pipefail

bam=${1:?usage: differential-test.sh <input.bam> [threads]}
threads=${2:-4}

command -v samtools >/dev/null || { echo "samtools is required" >&2; exit 1; }
[[ -f $bam ]] || { echo "no such file: $bam" >&2; exit 1; }

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
bamsplit=${BAMSPLIT:-$root/target/release/bamsplit}
[[ -x $bamsplit ]] || { echo "build first: cargo build --release" >&2; exit 1; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
echo "workspace: $work"

echo "==> bamsplit chrom --threads $threads"
"$bamsplit" chrom "$bam" --out-dir "$work/bamsplit" --threads "$threads" --index auto

[[ -f "$bam.bai" || -f "$bam.csi" || -f "${bam%.bam}.bai" ]] || samtools index -@ "$threads" "$bam"

failures=0
total_split=0
while read -r reference _; do
  [[ $reference == "*" ]] && continue

  expected_count=$(samtools view -c "$bam" "$reference")
  output=$work/bamsplit/$reference.bam
  if [[ ! -f $output ]]; then
    if (( expected_count > 0 )); then
      echo "MISSING  $reference ($expected_count records expected)"
      failures=$((failures + 1))
    fi
    continue
  fi

  actual_count=$(samtools view -c "$output")
  total_split=$((total_split + actual_count))

  # Compare the eleven mandatory fields. Auxiliary data is excluded because
  # `samtools view` renders it in its own order; the Rust suite compares full
  # record bytes.
  expected_digest=$(samtools view "$bam" "$reference" | cut -f1-11 | sha256sum | cut -d' ' -f1)
  actual_digest=$(samtools view "$output" | cut -f1-11 | sha256sum | cut -d' ' -f1)

  if [[ $expected_count == "$actual_count" && $expected_digest == "$actual_digest" ]]; then
    if samtools quickcheck -v "$output" 2>/dev/null; then
      printf 'OK       %-24s %s records\n' "$reference" "$actual_count"
    else
      echo "BADFILE  $reference (quickcheck rejected the output)"
      failures=$((failures + 1))
    fi
  else
    echo "MISMATCH $reference: count $expected_count/$actual_count digest ${expected_digest:0:12}/${actual_digest:0:12}"
    failures=$((failures + 1))
  fi
done < <(samtools idxstats "$bam")

unmapped=$work/bamsplit/unmapped.bam
if [[ -f $unmapped ]]; then
  count=$(samtools view -c "$unmapped")
  total_split=$((total_split + count))
  printf 'OK       %-24s %s records\n' "(unplaced)" "$count"
fi

echo "==> conservation"
input_total=$(samtools view -c "$bam")
if [[ $input_total == "$total_split" ]]; then
  echo "OK       $input_total records in, $total_split out"
else
  echo "MISMATCH $input_total records in, $total_split out"
  failures=$((failures + 1))
fi

echo "==> manifest"
python3 - "$work/bamsplit/bamsplit.manifest.json" <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
lhs = m["input_records"]
rhs = m["unique_emitted_records"] + m["dropped_records"] + m["unmatched_records"]
print(f"OK       manifest conserves {lhs} records" if lhs == rhs
      else f"MISMATCH manifest: {lhs} != {rhs}")
sys.exit(0 if lhs == rhs else 1)
PY

if (( failures > 0 )); then
  echo "FAILED: $failures reference(s) differ"
  exit 1
fi
echo "PASSED"
