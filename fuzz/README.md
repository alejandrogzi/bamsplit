# Fuzz targets

Every target here feeds arbitrary bytes into a parser and asserts one thing: it
returns, either with a value or with a typed error. A panic, an abort, or a
hang is a bug — parsers in this crate must never do any of those on untrusted
input, because "untrusted input" here means "a BAM someone else produced".

```console
cargo install cargo-fuzz
cargo fuzz list
cargo fuzz run raw_record -- -max_total_time=300
cargo fuzz run tags       -- -max_total_time=300
```

The fuzz crate is deliberately outside the workspace: `cargo-fuzz` needs a
nightly toolchain and `-Z` flags, and requiring that of everyone who runs
`cargo test` would be a poor trade.

| target | what it attacks |
| --- | --- |
| `bam_header` | magic, `l_text`, header text, the reference dictionary |
| `raw_record` | `block_size`, the fixed core, section bounds, framing |
| `tags` | auxiliary type codes, `B` arrays, `Z`/`H` termination |
| `cigar` | operation codes, span overflow, block derivation |
| `filename` | encode/decode round-tripping and containment |
| `spool` | spool framing, length prefixes, checksum, footer |
