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

After the command finishes, inspect:

- `out/report.html` for the interactive summary report
- `out/summary.json` for machine-readable metrics
- `out/figures/` for the generated SVG figures

## Typical workflow

1. Prepare a FASTQ or FASTQ.GZ file.
2. Define anchor rules and QC thresholds in a JSON config file.
3. Run the CLI with the FASTQ, config, and an output directory.
4. Review the HTML report for plots and high-level interpretation.
5. Parse `summary.json` in downstream automation if needed.

## CLI usage

```bash
python -m scfastq_qc.cli run --fastq <reads.fastq.gz> --config <config.json> --outdir <outdir>
```

Arguments:

- `--fastq`: input FASTQ/FASTQ.GZ file
- `--config`: JSON config file
- `--outdir`: output directory; created automatically if it does not exist

## Config format

The config file contains four top-level sections:

```json
{
  "sample_name": "example_sample",
  "anchors": [
    {
      "name": "adapter_5p",
      "type": "fixed",
      "sequence": "ACGTACGT",
      "max_mismatches": 1
    },
    {
      "name": "polyT",
      "type": "regex",
      "pattern": "T{6,}"
    }
  ],
  "structure": {
    "expected_order": ["adapter_5p", "polyT"]
  },
  "thresholds": {
    "long_read_min_bp": 20,
    "long_read_min_q": 20,
    "heatmap_max_reads": 100
  }
}
```

Field reference:

- `sample_name`: sample label used in output summaries
- `anchors`: anchor definitions
  - `name`: label used in summaries and structure calls
  - `type`: `fixed` or `regex`
  - `sequence`: required for `fixed` anchors
  - `pattern`: required for `regex` anchors
  - `max_mismatches`: optional for `fixed` anchors; must be a non-negative integer
- `structure.expected_order`: expected left-to-right anchor order used to classify read structure
- `thresholds.long_read_min_bp`: minimum read length used in the long/high-quality ratio
- `thresholds.long_read_min_q`: minimum mean quality used in the long/high-quality ratio
- `thresholds.heatmap_max_reads`: maximum reads rendered in the heatmap

See `examples/config.json` for a working example.

## Output summary

The report currently includes:

- total reads, total bases, mean/median read length, and N50
- read length histogram and cumulative distribution
- mean read quality distribution and length-vs-quality scatter plot
- anchor detection ratios and structure classification counts
- anchor occupancy heatmap across normalized read positions

## Build Rust accelerator

The Rust module is an optional accelerator for fixed-anchor scans. If it is unavailable or fails to build/load, the package automatically falls back to the pure-Python implementation.

If `cargo` is available, the Python package attempts to build the Rust hotspot module automatically on first use, so most users do not need to run a separate setup step.

You can still build it manually if you want:

```bash
cargo build --release --manifest-path rust/anchor_engine/Cargo.toml
```

If the shared library is present under the Cargo `target/release/` directory, the Python code will load it automatically and use Rust for fixed-anchor scans; otherwise it falls back to the pure-Python implementation.

## Troubleshooting

- **Config load fails on another machine**: ensure the JSON file is UTF-8 encoded.
- **Anchor preparation raises a validation error**: verify that `fixed` anchors include `sequence`, `regex` anchors include `pattern`, and `max_mismatches` is non-negative.
- **Rust accelerator is not used**: this is not fatal. Check that `cargo` is installed and a compatible Rust toolchain/linker is available if you want the optional speedup.
