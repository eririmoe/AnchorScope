# Publication and release checklist

## Completed in 0.2.0

- [x] Versioned package metadata, license, and citation file.
- [x] Machine-readable configuration schema and semantic validation.
- [x] Explicit mathematical and algorithmic method definitions.
- [x] Benchmark design and interpretation document.
- [x] Deterministic synthetic truth generator with fixed seed.
- [x] Confusion matrix, accuracy, throughput, environment, backend, and Git commit capture.
- [x] Tests for substitutions, insertions, deletions, structure rules, concatemer cycles, segments, barcode ambiguity, external tags, gzip, and MultiQC.
- [x] Serial versus multiprocessing end-to-end equivalence test.
- [x] End-to-end SAM tag audit plus optional Linux/Python 3.12 BAM CI job.
- [x] Self-contained HTML, JSON, CSV, FASTQ, and MultiQC outputs.
- [x] Limitations stated without claiming barcode correction or biological poly(A) measurement.

## Required before manuscript submission

- [ ] Replace the collective author in `CITATION.cff` with final authors and ORCIDs.
- [ ] Add public repository, issue tracker, archived release DOI, and accessioned datasets.
- [ ] Run benchmarks on the publication compute environment and retain unedited JSON.
- [ ] Benchmark at least three real datasets spanning protocol, chemistry, instrument, and read quality.
- [ ] Compare anchor detection with an independent aligner and manually curated truth subset.
- [ ] Compare structure calls with at least one protocol-specific pipeline on identical reads.
- [ ] Report wall time, CPU time, peak RSS, disk, throughput, and input size.
- [ ] Add confidence intervals where scientifically appropriate.
- [ ] Freeze dependencies and containers; archive configs, commands, logs, and checksums.
- [x] Add continuous integration on Linux, macOS, Windows, and supported Python versions.
- [ ] Obtain independent review of terminology such as truncation and full-length.
- [ ] Perform manual visual regression checks in current Chrome, Firefox, and Safari; automated local-browser inspection was unavailable in the present sandbox.

Passing automated tests demonstrates software consistency, not biological validity. Manuscript claims must remain limited to the datasets and truth definitions actually benchmarked.
