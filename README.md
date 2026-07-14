# AnchorScope

`AnchorScope` is a high-performance, anchor-aware Quality Control (QC) and filtering tool designed specifically for single-cell long-read sequencing libraries.

Unlike traditional bulk FASTQ QC tools, `AnchorScope` understands the molecular structure of single-cell libraries. It actively searches for expected sequences (like adapters, polyA tails, cell barcodes) to classify each read, providing deep insights into library construction success, structural integrity, and overall sequencing quality.

## Key Features

*   **Structure-Aware QC**: Identifies expected motifs (Anchors) using exact match with mismatches or regex patterns.
*   **Orientation Agnostic**: Automatically checks both forward and reverse-complement strands to determine the true orientation of each read.
*   **Self-contained HTML Reports**: Generates accessible SVG visualizations without external web dependencies.
*   **Batch Processing**: Natively supports processing multiple samples concurrently.
*   **Parallel Acceleration**: Utilizes multi-processing to significantly speed up single-file processing (`--threads`).
*   **FASTQ Filtering**: Optionally splits reads into `passed.fastq` and `failed.fastq` based on comprehensive QC verdicts.
*   **Indel-aware Alignment**: Semi-global Levenshtein matching reports substitutions, insertions, deletions, and CIGAR operations.
*   **Protocol Structure Grammar**: Anchor count, terminal position, distance, orientation, segment, and concatemer-cycle rules.
*   **Segment-level QC**: Length, Q-score, GC, N-content, observability, and configured-bound conformance between anchors.
*   **Barcode/UMI Recoverability QC**: Diagnostic-only extraction, quality filtering, whitelist distance, and ambiguity metrics without sequence correction.
*   **Workflow Interoperability**: MultiQC custom content, gzip FASTQ output, and BAM tag auditing for common single-cell tags and Dorado poly(A).

---

## Installation

### Prerequisites

*   Python 3.10+
*   Git

### 1. Clone the repository

```bash
git clone https://github.com/eririmoe/AnchorScope.git
cd AnchorScope
```

### 2. Create and activate a virtual environment

Using a dedicated environment avoids conflicts with other Python tools.

Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

Windows PowerShell:

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

### 3. Install AnchorScope

For normal use, install the checked-out source tree into the active environment:

```bash
python -m pip install .
```

For development, use an editable installation so that local code changes take effect immediately:

```bash
python -m pip install -e .
```

To enable optional BAM tag auditing, install the `bam` extra:

```bash
python -m pip install ".[bam]"
```

### 4. Verify the installation

```bash
anchorscope --help
python -c "import anchorscope; print('AnchorScope is installed')"
```

To update a cloned copy later:

```bash
git pull
python -m pip install .
```

---

## Core Concepts

Understanding how `AnchorScope` processes reads is key to configuring it correctly.

### 1. Anchors
An **Anchor** is a known sequence motif expected to be present in your library (e.g., a 5' adapter, a 3' adapter, or a polyA tail).
*   **Fixed Anchors**: Defined by a specific DNA sequence and a maximum number of allowed mismatches.
*   **Regex Anchors**: Defined by a regular expression (e.g., `A{8,}` for a polyA tail).

**Anchor Matching Principles:**
1.  **Fuzzy Searching**: `max_mismatches` selects the fast, fixed-length Hamming matcher and permits substitutions only. Setting `max_edits` selects semi-global Levenshtein alignment and permits substitutions, insertions, and deletions. The latter also reports operation counts and a CIGAR string.
2.  **Best Hit Selection**: If a motif appears multiple times in one read, the tool selects the hit with the lowest edit burden and then the earliest coordinates.
3.  **Strand Agnostic**: Because single-cell long-read libraries often sequence both the forward and reverse-complement strands randomly, the anchor matching runs twice for every read: once on the raw sequence, and once on its reverse-complement. The strand that yields the most complete and correctly ordered set of anchors is determined to be the true biological orientation.

### 2. Structure Classification
`AnchorScope` scans every read (and its reverse complement) for all defined anchors. Based on what it finds, it assigns a **Structure Label**:
*   `full_structure`: All expected anchors were found.
*   `missing_5p_anchor` / `missing_3p_anchor`: A required terminal anchor was not found.
*   `anchor_order_invalid`: Anchors were found in an order inconsistent with the configured grammar.
*   `structure_rule_violation`: Order was present, but a count, terminal-position, or distance rule failed.
*   `concatemer_candidate`: At least the configured number of complete anchor cycles was detected.
*   `no_anchor_detected`: None of the anchors were found.

### 3. QC Buckets
Beyond just finding anchors, the tool evaluates read length, quality scores, truncation, and valid base content to assign each read into a mutually exclusive **QC Bucket**:
*   `pass`: The read meets all quality criteria, has valid anchor order, and is not severely truncated.
*   `too_short` / `low_read_q`: Fails length or mean quality score thresholds.
*   `order_invalid` / `structure_rule_violation`: Anchors conflict with order or protocol rules.
*   `no_anchor`: Fails to detect any anchors.
*   `five_prime_truncated` / `three_prime_truncated`: Terminal anchors are too far from their expected read ends.
*   `high_n` / `invalid_bases`: The read contains too many `N` bases or non-standard characters.

---

## Configuration (`config.json`)

`AnchorScope` is heavily driven by a JSON configuration file. Here is a detailed breakdown of all parameters:

```json
{
  "sample_name": "My_Experiment",
  "qscore_method": "conservative",
  "export_non_full_structure_fastq": false,

  "anchors": [
    {
      "name": "adapter_5p",
      "type": "fixed",
      "sequence": "CGACATGGCTACGATCCGACTT",
      "max_edits": 2,
      "search_region": "5p",
      "search_window_bp": 150
    },
    {
      "name": "polyA",
      "type": "regex",
      "pattern": "A{8,}"
    }
  ],

  "structure": {
    "expected_order": ["adapter_5p", "polyA"],
    "expected_orientation": "either",
    "concatemer_min_cycles": 2,
    "split_concatemers": true,
    "anchor_rules": {
      "adapter_5p": {"terminal": "5p", "max_terminal_offset": 0.15},
      "polyA": {"terminal": "3p", "max_terminal_offset": 0.15}
    },
    "segments": [
      {"name": "transcript", "start_anchor": "adapter_5p", "end_anchor": "polyA", "role": "insert"}
    ]
  },

  "barcode": {
    "enabled": false,
    "left_anchor": "adapter_5p",
    "right_anchor": "polyA",
    "barcode_length": 16,
    "umi_length": 12,
    "min_base_quality": 15,
    "whitelist_path": "3M-february-2018.txt",
    "max_edit_distance": 2,
    "min_distance_margin": 1
  },

  "thresholds": {
    "long_read_min_bp": 500,
    "long_read_min_q": 10.0,
    "terminal_anchor_max_offset": 0.15,
    "high_n_fraction": 0.1,

    "warn_long_high_quality_ratio": 0.7,
    "fail_long_high_quality_ratio": 0.5,
    "warn_correct_anchor_order_ratio": 0.7,
    "fail_correct_anchor_order_ratio": 0.5,
    "warn_no_anchor_ratio": 0.2,
    "fail_no_anchor_ratio": 0.35,
    "warn_reversed_read_ratio": 0.3,
    "fail_reversed_read_ratio": 0.5,
    "warn_high_n_ratio": 0.1,
    "fail_high_n_ratio": 0.2
  },

  "samples": [
    { "sample_name": "Sample_A", "path": "data/sampleA.fastq.gz" },
    { "sample_name": "Sample_B", "path": "data/sampleB.fastq.gz" }
  ]
}
```

### Parameter Details

#### Global Settings
*   `sample_name` (string): Default name for the sample (used in reports).
*   `qscore_method` (string): Method for calculating mean read Q-score. `"conservative"` (default, converts to probabilities first) or `"arithmetic_mean"`.
*   `export_non_full_structure_fastq` (bool): If true, creates a separate FASTQ file for *each* structure class (e.g., `missing_polyA.fastq`) containing the reads that fell into that class.

#### `anchors` (List of Objects)
Defines the motifs to search for.
*   `name` (string): Unique identifier for the anchor.
*   `type` (string): `"fixed"` or `"regex"`.
*   `sequence` (string): Required if type is `"fixed"`.
*   `max_mismatches` (int): Substitution-only threshold for the Hamming matcher; omit or leave at zero when `max_edits` is used.
*   `max_edits` (int): Optional indel-aware alternative to `max_mismatches`; uses Levenshtein distance.
*   `search_region` / `search_window_bp`: Optionally restrict matching to a 5-prime or 3-prime window.
*   `pattern` (string): Required if type is `"regex"`.

#### `structure`
*   `expected_order` (List of strings): The order in which the defined anchors should appear from 5' to 3' on the read.
*   `expected_orientation`: `"either"`, `"forward"`, or `"reverse"`. Bidirectional protocols should use `"either"`.
*   `anchor_rules`: Per-anchor count, terminal-position, and inter-anchor distance constraints.
*   `segments`: Named molecular regions bounded by two anchors for segment-level QC.
*   `split_concatemers`: Export each complete repeated structure cycle as an oriented subread.

#### `thresholds`
Controls the logic for QC Buckets and the final pass/warn/fail status of the entire run.
*   `long_read_min_bp`: Minimum length for a read to be considered "long/valid".
*   `long_read_min_q`: Minimum mean Q-score for a read to be considered "high quality".
*   `terminal_anchor_max_offset`: Used to detect truncation. If the first/last expected anchors are found, but their distance from the end of the read exceeds this fraction (e.g., `0.15` = 15% of read length), the read is flagged as truncated.
*   `high_n_fraction`: Maximum allowed fraction of 'N' bases in a read.
*   `warn_*` / `fail_*`: Thresholds for the overall sample HTML report verdicts. For example, if the ratio of reads with no anchors exceeds `fail_no_anchor_ratio`, the report will flag the library construction as "Failed".

#### `samples` (Optional)
A list of sample dictionaries (`sample_name`, `path`). If provided, you can run batch jobs entirely driven by the config file without specifying paths on the command line.

---

## Usage & CLI Reference

The CLI has two main subcommands: `run` (for single files or config-driven batches) and `batch` (for directory/file-list inputs).

### Global Options
*   `--log-file <path>`: Write logs to a file.
*   `--log-level <DEBUG|INFO|WARNING|ERROR|CRITICAL>`: Default is `INFO`.

### Subcommand: `run`
Used for processing a single FASTQ file, or triggering a batch if the `--config` file contains a `"samples"` list.

```bash
anchorscope run \
  --fastq input.fastq.gz \
  --config config.json \
  --outdir ./results \
  --threads 8 \
  --output-passed-fastq \
  --output-failed-fastq \
  --gzip-output \
  --export-csv
```

**Parameters for `run`:**
*   `--fastq`: Path to input FASTQ/FASTQ.GZ. *Required unless your config file has a `samples` list.*
*   `--config`: Path to `config.json`.
*   `--outdir`: Where to save the output files.
*   `--threads N`: Number of CPU processes to use for parallelizing the anchor search within the single FASTQ file. Dramatically speeds up processing (Default: 1).
*   `--output-passed-fastq`: Export reads that landed in the `pass` QC Bucket to `passed.fastq`.
*   `--output-failed-fastq`: Export reads that landed in any non-pass QC Bucket to `failed.fastq`.
*   `--export-csv`: Export detailed per-read and per-anchor statistics to CSV files.
*   `--gzip-output`: Compress all emitted FASTQ files as `.fastq.gz`.

### Subcommand: `audit-bam`

Audit tags generated by an external barcode/UMI workflow without rerunning correction:

```bash
anchorscope audit-bam --bam tagged.bam --outdir tag_audit --sample-name sample_A
```

The command reads `CR`, `CB`, `CY`, `UR`, `UB`, `UY`, and Dorado `pt` tags. It writes `tag_audit_summary.json` and `anchorscope_tag_audit_mqc.json`. Text SAM works without optional dependencies; binary BAM support is installed with `pip install -e ".[bam]"`.

### Subcommand: `batch`
Used for processing multiple FASTQ files by scanning a directory or reading a list of files.

```bash
# Scan a directory for FASTQ files
anchorscope batch \
  --input /path/to/data_dir \
  --pattern "*.fastq.gz" \
  --config config.json \
  --outdir ./batch_results \
  --parallel 4 \
  --output-passed-fastq
```

**Parameters for `batch`:**
*   `--input`: Path to a directory containing FASTQ files, or a text file containing one FASTQ path per line. *Required unless config has a `samples` list.*
*   `--pattern`: If `--input` is a directory, glob pattern to match files (Default: `*.fastq*`).
*   `--config`: Path to `config.json`.
*   `--outdir`: Base output directory. A subdirectory will be created for each sample.
*   `--parallel N`: Number of *samples* to process concurrently (Sample-level parallelism).
*   `--continue-on-error`: If one sample fails, log the error and continue with the rest.
*   `--output-passed-fastq` / `--output-failed-fastq`: Export filtered FASTQs inside each sample's subdirectory.
*   `--export-csv`: Export CSVs inside each sample's subdirectory.

> **Parallelism Note:**
> *   Use `batch --parallel N` to process *N different samples* at the same time.
> *   Use `run --threads N` to process *1 sample* using N CPU cores.

---

## Output Files

Depending on the mode and flags used, `AnchorScope` generates the following inside the `--outdir`:

### Standard Outputs (Always Generated)
*   **`report.html`**: A highly visual, interactive HTML report containing length distributions, Q-score plots, anchor detection heatmaps, and overall sample verdicts. It is completely self-contained (no internet required to view).
*   **`summary.json`**: A machine-readable JSON file containing all the raw metrics and calculated ratios presented in the HTML report.
*   **`anchorscope_mqc.json`**: MultiQC custom-content general statistics, emitted automatically.

### Filtered FASTQs (When requested)
*   **`passed.fastq`**: Generated if `--output-passed-fastq` is provided. Contains reads deemed "good" by the QC bucket logic.
*   **`failed.fastq`**: Generated if `--output-failed-fastq` is provided. Contains the reads that failed QC.
*   **`<structure_label>.fastq`**: Generated if `export_non_full_structure_fastq` is `true` in the config. Groups reads strictly by missing anchors or unexpected order.

### CSV Exports (When `--export-csv` is requested)
*   **`summary.csv`**: Flattened version of the high-level metrics.
*   **`structure_classification.csv`**: Counts and ratios of reads falling into each structure category.
*   **`read_details.csv`**: A large file containing 1 row per read, detailing its length, q-score, offsets, structure, and QC bucket.
*   **`anchor_hits.csv`**: One row per anchor hit with coordinates, total edits, substitutions, insertions, deletions, CIGAR, and matched sequence.

### Batch Outputs (Only in `batch` mode)
In addition to a subdirectory for each sample containing the files above, batch mode generates at the root of `--outdir`:
*   **`batch_report.html`**: A high-level summary table comparing key metrics (Total reads, High-quality %, Valid order %, Overall status) across all processed samples.
*   **`batch_summary.json`**: JSON representation of the batch results.

---

## Reproducibility and publication validation

AnchorScope uses deterministic algorithms throughout the QC path. The synthetic truth benchmark fixes its random seed and records the Git commit, Python version, platform, matcher backend, confusion matrix, accuracy, and throughput.

```bash
python -B -m unittest discover -s tests -v
python -B benchmarks/benchmark_synthetic.py \
  --reads-per-class 100 \
  --seed 20260710 \
  --output benchmarks/results/synthetic_benchmark.json
python -B benchmarks/validate_comprehensive_synthetic.py \
  --output-dir validation/synthetic_comprehensive \
  --replicates 5 \
  --seed 20260713 \
  --threads 2
```

Formal metric definitions, algorithm details, limitations, benchmarking guidance, and validation design are documented in `docs/METHODS.md`, `docs/BENCHMARKING.md`, and `docs/SYNTHETIC_VALIDATION.md`. Configuration keys are machine-described by `src/anchorscope/config.schema.json`.
