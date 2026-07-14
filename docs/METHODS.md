# AnchorScope computational methods

## Scope

AnchorScope 0.2.0 performs pre-alignment quality control of long-read libraries whose molecular architecture can be described by constant sequence anchors. It does not infer biological cell identities, correct barcodes or UMIs, call transcripts, or replace a single-cell quantification workflow.

## FASTQ parsing and read quality

FASTQ records are processed as a stream and sequence/quality lengths must agree. By default, per-read quality is calculated from mean error probability:

`Q_read = -10 log10(mean(10^(-Q_i/10)))`.

The optional arithmetic method reports the mean of base-level Phred scores. Length, base-Q, read-Q, length/Q scatter, N content, invalid bases, and N50 are accumulated without retaining complete reads in memory.

## Anchor alignment

Fixed anchors have two explicit modes.

- `max_mismatches`: substitution-only Hamming distance, preserving the accelerated legacy implementation.
- `max_edits`: semi-global Levenshtein alignment of a motif against any substring of the configured search region. Dynamic-programming row zero is initialized to zero to permit a free read prefix. Substitution, insertion, and deletion have unit cost.

Qualifying endpoints are traced back to obtain coordinates, operation counts, and compact CIGAR. Overlapping endpoints representing one occurrence are reduced to the lowest-edit, closest-length alignment; distinct non-overlapping occurrences are retained. `search_region` and `search_window_bp` can constrain an anchor to a 5-prime or 3-prime window.

Regex anchors remain available for motifs such as sequence-level poly(A) proxies. A regex homopolymer length must not be interpreted as a signal-derived biological poly(A) tail length. Dorado `pt` tags can instead be audited from BAM.

## Orientation and structure classification

Each read is evaluated in its input orientation and as its reverse complement. Orientations are ranked by valid order, structure-class priority, detected-anchor count, hit count, and total edit burden. The selected orientation becomes a failure criterion only when `expected_orientation` is `forward` or `reverse`; `either` treats both as valid.

Rules can specify anchor requirement, count, terminal role, normalized terminal offset, and distance from the previous anchor. Segments specify two bounding anchors and optional length limits.

A concatemer cycle is a contiguous occurrence of the complete expected anchor order. Reads with at least `concatemer_min_cycles` cycles are labeled `concatemer_candidate`. If splitting is enabled, each complete interval from first to last anchor is emitted in biological orientation with `|cycle=N` appended. Partial cycles are not emitted as complete molecules.

## Anchor-relative segment QC

The segment starts at the end of its left anchor and ends at the start of its right anchor. Reversed or absent boundaries produce `missing`. Observed segments report length, conservative mean Q-score, GC fraction, N fraction, and configured-bound status. Aggregates include observed/missing ratios, median length and Q-score, mean GC/N fractions, and bound-conformance ratio.

## Barcode and UMI recoverability

Barcode diagnostics are deliberately non-corrective. A fixed-length barcode and UMI are extracted relative to a configured left anchor, optionally constrained by a right anchor. Reads below the minimum barcode base quality are `low_quality`.

With a whitelist, exact members are `exact`. Other sequences are queried in a BK-tree by Levenshtein distance. A sequence is `uniquely_recoverable` only when its nearest candidate is within `max_edit_distance` and its distance advantage over the runner-up reaches `min_distance_margin`; otherwise it is `ambiguous` or `unrecoverable`. The candidate is reported for audit but is never written back.

## External tag audit

`audit-bam` reads raw/corrected cell barcode (`CR`, `CB`), barcode quality (`CY`), raw/corrected UMI (`UR`, `UB`), UMI quality (`UY`), and Dorado poly(A) (`pt`) tags from BAM/uBAM or text SAM. It reports tag coverage, raw-to-corrected edit-distance distributions, and non-negative poly(A) length statistics. Numeric BAM barcode suffixes such as `-1` are excluded before comparison. Binary BAM uses optional `pysam`; SAM parsing uses the standard library and shares the same accumulator and outputs.

## QC verdicts

Higher-is-better metrics fail below the fail threshold and warn below the warn threshold. Lower-is-better metrics fail above the fail threshold and warn above the warn threshold. Validation rejects inverted thresholds, out-of-range fractions, unknown anchors, duplicate names, and invalid structural bounds.

The mutually exclusive primary bucket is selected from an ordered read-level flag list. Detailed CSV output retains the complete flag list.

## Determinism and parallelism

The QC path contains no random sampling. Multiprocessing preserves input order. Automated tests compare serial and two-process summaries for total reads, structures, QC buckets, segments, barcode diagnostics, and order-valid ratio.

The synthetic benchmark uses seeded random payloads but rejects any payload containing a sequence within one Levenshtein edit of either anchor or either reverse-complement anchor. This prevents unlabelled approximate anchors from contaminating the structure truth set.

## Limitations

- Sequence-only anchor and poly(A) calls depend on basecalling accuracy.
- The pure-Python edit aligner is slower than the Rust Hamming path; benchmarks must report mode and backend.
- Barcode recoverability is not validated barcode correction or cell calling.
- Fused-read splitting emits only complete configured cycles and does not generate consensus.
- BAM auditing requires optional `pysam` and evaluates existing tags rather than reproducing the upstream algorithm.
