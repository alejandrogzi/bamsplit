#!/usr/bin/env bash
# Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
# Distributed under the terms of the Apache License, Version 2.0.
#
# Times `bamsplit chrom` against the equivalent `samtools view` loop.
#
# The samtools loop is the baseline every bioinformatician already uses:
#
#   for chr in $(samtools idxstats in.bam | cut -f1); do
#       samtools view -b in.bam "$chr" > "$chr.bam"
#   done
#
# It opens the input once per chromosome. With an index each open seeks to that
# reference's chunks rather than re-decompressing the whole file, so the honest
# comparison is not "N passes versus one" — it is per-invocation overhead, index
# dependence, and what each tool does with a thread budget. Note also that
# samtools `-@N` means N *additional* threads while `--threads N` is bamsplit's
# total budget, and that these bamsplit runs also build an output index, which
# the samtools loop does not. `--index none` below isolates that cost.
#
#   scripts/benchmark.sh sample.bam [thread-list] [annotation.bed]
#
# With a third argument the script also times `bamsplit region` across the
# engines, which is where the indexed engine's per-region query plan shows up.
set -euo pipefail

bam=${1:?usage: benchmark.sh <input.bam> [thread-list] [annotation.bed]}
threads_list=${2:-"1 2 4 8 16"}
annotation=${3:-}

command -v samtools >/dev/null || { echo "samtools is required" >&2; exit 1; }
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
bamsplit=${BAMSPLIT:-$root/target/release/bamsplit}
[[ -x $bamsplit ]] || { echo "build first: cargo build --release" >&2; exit 1; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

input_bytes=$(stat -c %s "$bam" 2>/dev/null || stat -f %z "$bam")
records=$(samtools view -c "$bam")
references=$(samtools idxstats "$bam" | grep -cv '^\*')

echo "input: $bam"
echo "  $input_bytes bytes, $records records, $references references"
echo

# `%e` wall seconds, `%U`+`%S` CPU seconds, `%M` peak RSS in KiB.
timefmt='%e %U %S %M'
measure() {
  local label=$1 outdir=$2
  shift 2
  rm -rf "$outdir"
  local stats
  # `env time` execs a program, so a shell function has to be reached through
  # an explicit `bash -c`; see `samtools_script` below.
  if ! /usr/bin/env time -f "$timefmt" -o "$work/time" "$@" >/dev/null 2>"$work/err"; then
    echo "FAILED: $label" >&2
    sed -n '1,5p' "$work/err" >&2
    return 1
  fi
  stats=$(cat "$work/time")
  read -r wall user sys rss <<<"$stats"

  local out_bytes opened opens
  out_bytes=$(find "$outdir" -name '*.bam' -printf '%s\n' 2>/dev/null | paste -sd+ | bc || echo 0)
  opened=$(find "$outdir" -name '*.bam' 2>/dev/null | wc -l)
  # The samtools loop opens the input once per reference. With an index each
  # open only *seeks* to that reference's chunks, so the aggregate bytes read
  # are close to one pass; the cost is per-invocation startup and index load,
  # not eight full decompressions.
  if [[ $label == samtools* ]]; then opens=$references; else opens=1; fi

  printf '%-22s wall %7ss  cpu %7ss  rss %8s KiB  out %12s B  files %5s  input-opens %s\n' \
    "$label" "$wall" "$(echo "$user + $sys" | bc)" "$rss" "${out_bytes:-0}" "$opened" "$opens"
  if [[ -n ${records:-} && $wall != 0 ]]; then
    printf '%-22s %s records/s, %s input B/s, compression %.3f\n' "" \
      "$(echo "$records / $wall" | bc)" "$(echo "$input_bytes / $wall" | bc)" \
      "$(echo "scale=3; ${out_bytes:-1} / $input_bytes" | bc)"
  fi
}

echo "==> samtools baseline (one pass per reference)"
samtools_script=$work/samtools-loop.sh
cat >"$samtools_script" <<'LOOP'
#!/usr/bin/env bash
set -euo pipefail
bam=$1 outdir=$2 threads=$3
mkdir -p "$outdir"
while read -r reference _; do
  [[ $reference == "*" ]] && continue
  samtools view -@ "$threads" -b "$bam" "$reference" > "$outdir/$reference.bam"
done < <(samtools idxstats "$bam")
LOOP
chmod +x "$samtools_script"
for threads in $threads_list; do
  measure "samtools -@$threads" "$work/samtools" \
    "$samtools_script" "$bam" "$work/samtools" "$threads"
done

echo
echo "==> bamsplit chrom, stream engine"
for threads in $threads_list; do
  measure "bamsplit stream -t$threads" "$work/stream" \
    "$bamsplit" chrom "$bam" --out-dir "$work/stream" --threads "$threads" \
    --engine stream --index auto --force --quiet
done

echo
echo "==> bamsplit chrom, other configurations (4 threads)"
# Every row below pins one variable against the `stream -t4` row above. Leaving
# `--engine` on `auto` would move two variables at once, since `auto` may pick a
# different engine and swamp whatever the row is meant to isolate.
measure "bamsplit no-index" "$work/noindex" \
  "$bamsplit" chrom "$bam" --out-dir "$work/noindex" --threads 4 \
  --engine stream --index none --force --quiet
measure "bamsplit mmap" "$work/mmap" \
  "$bamsplit" chrom "$bam" --out-dir "$work/mmap" --threads 4 \
  --engine stream --index auto --io mmap --force --quiet
measure "bamsplit level 1" "$work/level1" \
  "$bamsplit" chrom "$bam" --out-dir "$work/level1" --threads 4 \
  --engine stream --index auto --compression-level 1 --force --quiet
measure "bamsplit spool" "$work/spool" \
  "$bamsplit" chrom "$bam" --out-dir "$work/spool" --threads 4 \
  --engine spool --index auto --force --quiet
if [[ -f "$bam.bai" || -f "$bam.csi" ]]; then
  measure "bamsplit indexed" "$work/indexed" \
    "$bamsplit" chrom "$bam" --out-dir "$work/indexed" --threads 4 \
    --engine indexed --index auto --force --quiet
fi
measure "bamsplit auto" "$work/auto" \
  "$bamsplit" chrom "$bam" --out-dir "$work/auto" --threads 4 --index auto --force --quiet

if [[ -n $annotation ]]; then
  echo
  echo "==> bamsplit region on $annotation (4 threads)"
  for engine in stream spool indexed auto; do
    measure "bamsplit region $engine" "$work/region-$engine" \
      "$bamsplit" region "$bam" --regions "$annotation" --out-dir "$work/region-$engine" \
      --threads 4 --engine "$engine" --index auto --assignment overlap --force --quiet
  done
fi

echo
echo "Numbers above are from this machine and this input only. Do not quote them"
echo "as a general claim: engine choice depends on storage, reference count, and"
echo "record size. See docs/performance.md."
