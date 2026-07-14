# Comprehensive synthetic truth validation

`benchmarks/validate_comprehensive_synthetic.py` generates deterministic FASTQ
inputs with known truth and runs the complete AnchorScope pipeline. It covers:

- exact, reverse-complement, lowercase, substituted, inserted, and deleted anchors;
- the edit-distance boundary and missing, reordered, duplicated, and concatemer structures;
- terminal offsets, segment-length failures, low read quality, high `N` content,
  invalid bases, short reads, and empty reads;
- exact, uniquely recoverable, ambiguous, unrecoverable, low-quality, incomplete,
  missing-anchor, and reverse-oriented barcode/UMI cases;
- malformed FASTQ rejection;
- serial/plain-FASTQ versus multiprocessing/gzip-FASTQ equivalence;
- gzip passed, failed, and concatemer-split exports.

Run the full reproducible validation from the repository root:

```bash
python -B benchmarks/validate_comprehensive_synthetic.py \
  --output-dir validation/synthetic_comprehensive \
  --replicates 5 \
  --seed 20260713 \
  --threads 2
```

The output directory contains the simulated FASTQ files, `truth.tsv`, configs,
per-run AnchorScope reports, `validation_results.json`, and a concise
`validation_report.md`. A passing result supports correctness for the encoded
synthetic cases; it does not replace validation on independent experimental
libraries or platform-specific error models.
