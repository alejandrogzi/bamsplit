# Release checklist

Run in order. Every step is a gate, not a suggestion.

## 1. The build gates

```console
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo build --workspace --all-targets --all-features --locked
cargo build --workspace --all-features --release --locked
cargo test --workspace --all-features --locked
cargo test --doc --workspace --all-features --locked
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps
```

Feature combinations, each on its own:

```console
for features in "--no-default-features" \
                "--no-default-features --features mmap" \
                "--no-default-features --features annotation" \
                "" "--all-features"; do
    cargo test -p bamsplit-core $features --locked
done
```

MSRV, from `rust-version` in the workspace manifest:

```console
rustup toolchain install 1.85
cargo +1.85 build --workspace --all-features --locked
cargo +1.85 test -p bamsplit-core --all-features --locked
```

## 2. Differential correctness

```console
cargo test --workspace --features differential --locked -- --nocapture
```

`--features differential` makes a missing `samtools` a **failure**, so a green
run here means the comparisons actually ran.

Then on real data, not just fixtures — this is the step that finds record shapes
nobody thought to synthesize:

```console
cargo build --release
scripts/differential-test.sh /path/to/real-sample.bam 8
```

Repeat for at least: a short-read WGS BAM, an RNA-seq BAM with spliced CIGARs, a
long-read BAM, and something with many small contigs.

## 3. Benchmarks

```console
cargo bench --bench core
cargo bench --bench end_to_end
scripts/benchmark.sh /path/to/real-sample.bam "1 2 4 8 16" /path/to/annotation.bed
```

Record the numbers in the release notes **with the machine and the input**. Do
not carry a claim forward from a previous release; storage and reference count
change the answer.

## 4. Behaviour under stress

* **Interruption.** Start a split of a large BAM, `Ctrl-C` partway, and confirm
  exit code 6 and an empty output directory:

  ```console
  bamsplit chrom big.bam --out-dir out & sleep 3; kill -INT %1; wait
  echo "exit $?"; ls -la out/
  ```

  Repeat with `--engine spool` and confirm the temporary directory is gone too.

* **Low file-descriptor limit.**

  ```console
  (ulimit -n 64; bamsplit chrom many-contigs.bam --out-dir out --max-open-files 16)
  ```

* **Disk full.** Point `--out-dir` at a small `tmpfs` and confirm the failure is
  a typed I/O error, exit code 1, and that no partial output survives:

  ```console
  sudo mount -t tmpfs -o size=1M tmpfs /mnt/tiny
  bamsplit chrom big.bam --out-dir /mnt/tiny/out; echo "exit $?"; ls -la /mnt/tiny/out
  ```

* **Existing outputs.** Confirm exit code 4 without `--force`, and that `--force`
  leaves the old file intact until the new one is complete.

* **stdin.**

  ```console
  cat sample.bam | bamsplit chrom - --out-dir out
  ```

## 5. Source audit

```console
# Every hand-written Rust module carries the watermark.
git ls-files '*.rs' | while read -r f; do
    head -2 "$f" | grep -q 'Copyright (c) 2025 Alejandro Gonzales-Irribarren' \
        || echo "MISSING: $f"
done

# Every unsafe block has a SAFETY comment.
git ls-files '*.rs' | xargs grep -n '\bunsafe\b *{'

# No stray panics on input paths.
git ls-files 'crates/*/src/*.rs' | xargs grep -n 'unwrap()\|expect(\|panic!\|unreachable!\|todo!' \
    | grep -v '#\[cfg(test)\]' | grep -v '/tests/'
```

Each surviving `unwrap`/`expect` outside tests must be justified by a comment
explaining the invariant that makes it unreachable.

CI runs the first two automatically; the third needs a human.

## 6. Licensing

* `LICENSE` contains the complete GPLv3 text.
* Every manifest declares `GPL-3.0-only`, directly or via
  `license.workspace = true`.
* `LICENSING.md` still describes the Apache-2.0 watermark discrepancy, and the
  project is **not** described anywhere as dual-licensed.
* Optionally: `cargo deny check licenses`.

**The discrepancy in `LICENSING.md` should be resolved before a public release.**
It is a known defect, not a decision.

## 7. Documentation

* `README.md` — every command has an example; the status table matches reality.
* `docs/` — architecture, BAM semantics, region routing, performance, pipeline
  examples all current.
* `CHANGELOG.md` — the release has an entry, and the "not yet implemented"
  section is accurate.
* `--help` for every subcommand mentions its options and their defaults.
* A new user can install and use every implemented command without reading
  source.

## 8. Version and tag

* Bump `version` in the workspace manifest.
* `cargo update --workspace` and commit `Cargo.lock`.
* Re-run section 1 with `--locked`.
* Add the CI badge to `README.md` once the workflow path and repository name are
  final.
* Tag, and attach the benchmark numbers from section 3.

## Completion criteria

Do not release unless every one of these is true.

- [ ] `chrom`, `shard`, `tag`, and `inspect` all work
- [ ] `region` fails loudly and its state is documented in `README.md` and `CHANGELOG.md`
- [ ] chromosome splitting is lossless — digests match `samtools`
- [ ] a coordinate-sorted chromosome split makes exactly one input pass
- [ ] unsorted input is handled through the spool engine
- [ ] every output header is valid and retains the full reference dictionary
- [ ] every output carries a BGZF end-of-file marker
- [ ] BAI and CSI generation both work and are accepted by `samtools`
- [ ] forced BAI beyond its coordinate limit is rejected
- [ ] output creation is transactional and interruption cleans up
- [ ] memory, file descriptors, and threads all stay inside their budgets
- [ ] QNAME sharding is deterministic and template-coherent
- [ ] tag cardinality is guarded and fails atomically
- [ ] manifests satisfy their conservation checks
- [ ] malformed input is rejected safely, with the record ordinal and offset
- [ ] differential `samtools` tests pass
- [ ] benchmarks are reproducible and reported with their context
- [ ] every hand-written Rust module carries the watermark
- [ ] GPL-3.0-only metadata and the GPLv3 text are present and consistent
- [ ] CI checks formatting, linting, debug and release builds, tests, doc tests,
      the feature matrix, MSRV, and the differential suite
