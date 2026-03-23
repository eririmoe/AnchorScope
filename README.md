# scfastq-qc

A structure-aware FASTQ QC tool for single-cell long-read libraries.

## Features

- Total reads / total bases / mean / median / N50 read length
- Read length histogram and cumulative curve
- Per-read mean Q distribution and length-vs-Q plot
- Configurable anchor detection using fixed sequences or regex patterns
- Anchor order validation and read structure classification
- HTML report with figures and summary tables
- Optional Rust accelerator for fixed-anchor scans with automatic Python fallback

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

### Field reference

#### `sample_name`

A free-form label written into the report and summary JSON. Use a value that identifies the library or sequencing run.

#### `anchors`

Each anchor definition describes one motif that should be detected in reads.

- `name`: required label used in summaries, heatmaps, and structure classification.
- `type`: required anchor type. Supported values are `fixed` and `regex`.
- `sequence`: required when `type` is `fixed`. The exact motif to scan for.
- `pattern`: required when `type` is `regex`. Any Python regular expression accepted by `re.compile`.
- `max_mismatches`: optional for `fixed` anchors. Defaults to `0` and must be a non-negative integer.

#### `structure.expected_order`

A left-to-right list of anchor names. Reads are classified against this expected order, so keep the names aligned with the `anchors` section.

#### `thresholds`

Threshold settings control how the QC summary is interpreted.

- `long_read_min_bp`: minimum read length used in the long/high-quality ratio. Default: `1000`.
- `long_read_min_q`: minimum mean quality used in the long/high-quality ratio. Default: `10.0`.
- `heatmap_max_reads`: maximum reads rendered in the heatmap. Default: `200`.

See `examples/config.json` for a minimal working example.

## Config examples

### Example 1: Adapter + polyT

This is the simplest structure model: a fixed adapter followed by a polyT tail.

```json
{
  "sample_name": "adapter_polyT_demo",
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
      "pattern": "T{8,}"
    }
  ],
  "structure": {
    "expected_order": ["adapter_5p", "polyT"]
  },
  "thresholds": {
    "long_read_min_bp": 1000,
    "long_read_min_q": 10.0,
    "heatmap_max_reads": 200
  }
}
```

### Example 2: Adapter + TSO + polyT

For many single-cell protocols, you may want to explicitly track a template-switch oligo (TSO) in addition to the 5' adapter and polyT region.

```json
{
  "sample_name": "adapter_tso_polyT_demo",
  "anchors": [
    {
      "name": "adapter_5p",
      "type": "fixed",
      "sequence": "CTACACGACGCTCTTCCGATCT",
      "max_mismatches": 1
    },
    {
      "name": "tso",
      "type": "fixed",
      "sequence": "AAGCAGTGGTATCAACGCAGAGTACATGGG",
      "max_mismatches": 2
    },
    {
      "name": "polyT",
      "type": "regex",
      "pattern": "T{10,}"
    }
  ],
  "structure": {
    "expected_order": ["adapter_5p", "tso", "polyT"]
  },
  "thresholds": {
    "long_read_min_bp": 1500,
    "long_read_min_q": 12.0,
    "heatmap_max_reads": 300
  }
}
```

### Tips for choosing anchors

- Use `fixed` anchors for known constant sequences such as adapters, barcodes, or TSO motifs.
- Use `regex` anchors for variable-length motifs such as polyT or polyA stretches.
- Start with `max_mismatches: 0` for fixed anchors, then relax only if real data shows expected sequencing noise.
- Keep `expected_order` short and biologically meaningful; the order is used for structure classification, not exhaustive annotation.

## Output summary

The report currently includes:

- total reads, total bases, mean/median read length, and N50
- read length histogram and cumulative distribution
- mean read quality distribution and length-vs-quality scatter plot
- anchor detection ratios and structure classification counts
- anchor occupancy heatmap across normalized read positions

## Build Rust accelerator

The Rust module is an optional accelerator for fixed-anchor scans. If it is unavailable, incompatible, or fails to build/load, the package automatically falls back to the pure-Python implementation.

If `cargo` is available, the Python package attempts to build the Rust hotspot module automatically on first use, so most users do not need to run a separate setup step.

You can still build it manually if you want:

```bash
cargo build --release --manifest-path rust/anchor_engine/Cargo.toml
```

If the shared library is present under the Cargo `target/release/` directory, the Python code will load it automatically and use Rust for fixed-anchor scans; otherwise it falls back to the pure-Python implementation.

## Troubleshooting

- **Config load fails on another machine**: ensure the JSON file is UTF-8 encoded.
- **Anchor preparation raises a validation error**: verify that `fixed` anchors include `sequence`, `regex` anchors include `pattern`, and `max_mismatches` is non-negative.
- **Structure calls look wrong**: make sure every name in `structure.expected_order` exactly matches an anchor `name`.
- **Rust accelerator is not used**: this is not fatal. Check that `cargo` is installed and a compatible Rust toolchain/linker is available if you want the optional speedup.
