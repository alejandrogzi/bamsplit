<p align="center">
  <p align="center">
    <img width=200 align="center" src="./assets/logo.png" >
  </p>

  <span>
    <h1 align="center">
        bamsplit
    </h1>
  </span>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.0.1-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.0.1-green">
    </a>
    <a href="https://crates.io/crates/bamsplit" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/bamsplit">
    </a>
    <a href="https://github.com/alejandrogzi/bamsplit" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/bamsplit?color=blue">
    </a>
    <a href="https://crates.io/crates/bamsplit" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/bamsplit">
    </a>
  </p>

  <p align="center">
    <samp>
        <span>lossless, high-throughput partitioning of BAM files</span>
        <br>
        <br>
        <a href="https://docs.rs/bamsplit/0.0.1/bamsplit/">docs</a> .
        <a href="https://github.com/alejandrogzi/bamsplit/tree/master/docs/architecture">architecture</a> .
        <a href="https://github.com/alejandrogzi/bamsplit/tree/master/docs/usage.md">usage</a> 
    </samp>
  </p>

</p>


## Installation
### Binary
```bash
cargo install --all-features bamsplit
```

### Docker
```bash
docker pull ghcr.io/alejandrogzi/bamsplit:latest
```

### Conda
```bash
conda install -c bioconda bamsplit
```

### Library
Add this to your `Cargo.toml`:

```toml
[dependencies]
bamsplit = { version = "0.0.1" }
```

## Benchmarks

### `chrom`
<div align="center">

| threads | samtools loop | bamsplit `stream` | bamsplit `stream` CPU |
| ---: | ---: | ---: | ---: |
| 1 | 9.44 s | **7.69 s** | 7.65 s |
| 2 | 4.82 s | 7.54 s | 7.50 s |
| 4 | 2.58 s | 5.42 s | 8.61 s |
| 8 | 1.66 s | 2.28 s | 9.73 s |
| 16 | 1.57 s | 1.78 s | 11.96 s |

* **The samtools loop parallelizes across processes.** Eight independent
  `samtools view` invocations each get `-@N`, so the loop's effective width is
  larger than the number the column says.

</div>

---

### `region`

200 regions of 10 kb each — 2 Mb of a 2 Gb declared genome, holding 7.9% of the
records — with `--assignment overlap`, 4 threads:

<div align="center">

| engine | wall | records read | outputs |
| --- | ---: | ---: | ---: |
| `stream` | 1.61 s | 2,000,000 | 200 |
| `spool` | 1.49 s | 2,000,000 | 200 |
| `indexed` | **0.15 s** | **157,098** | 200 |
| `auto` | 0.16 s | 157,098 | 200 |

</div>

*ran using hyperfine 1.18.0 on an AMD Ryzen 7 5700X with 128 GB of RAM and 16 cores with a 262 MB chain.gz file as input
