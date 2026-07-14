# Changelog

## 0.2.0 — 2026-07-10

- Added semi-global Levenshtein anchor alignment with substitution, insertion, deletion, and CIGAR reporting.
- Added strict configuration validation, terminal/count/distance rules, and protocol-aware orientation handling.
- Added anchor-relative segment length, quality, GC, N-content, and bound-conformance statistics.
- Added complete-cycle concatemer detection and optional oriented FASTQ splitting.
- Added diagnostic-only barcode/UMI extraction, whitelist distance, ambiguity, and recoverability metrics.
- Added external BAM tag auditing for `CR`, `CB`, `CY`, `UR`, `UB`, `UY`, and Dorado `pt` tags.
- Added automatic MultiQC custom-content and optional gzip FASTQ output.
- Added deterministic synthetic benchmarking and serial/parallel consistency tests.

## 0.1.0

- Initial anchor-aware FASTQ QC, filtering, batch processing, HTML reporting, and optional Rust Hamming matcher.
