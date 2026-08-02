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
