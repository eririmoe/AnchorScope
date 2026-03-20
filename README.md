# scfastq-qc

A structure-aware FASTQ QC tool for single-cell long-read libraries.

## Features

- Total reads / total bases / mean / median / N50 read length
- Read length histogram and cumulative curve
- Per-read mean Q distribution and length-vs-Q plot
- Configurable anchor detection using fixed sequences or regex patterns
- Anchor order validation and read structure classification
- HTML report with figures and summary tables

## Quick start

```bash
python -m scfastq_qc.cli run \
  --fastq examples/example.fastq \
  --config examples/config.json \
  --outdir out
```

## Config format

See `examples/config.json`.

## Build Rust accelerator

If `cargo` is available, the Python package now attempts to build the Rust hotspot module automatically on first use, so most users do not need to run a separate setup step.

You can still build it manually if you want:

```bash
cargo build --release --manifest-path rust/anchor_engine/Cargo.toml
```

If the shared library is present under the Cargo `target/release/` directory, the Python code will load it automatically and use Rust for fixed-anchor scans; otherwise it falls back to the pure-Python implementation.
