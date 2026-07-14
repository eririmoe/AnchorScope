"""Parallel (multiprocessing) single-sample QC pipeline.

Splits the per-read anchor search / classification work across multiple
worker processes while keeping I/O (FASTQ reading, CSV writing, FASTQ
filtering, report generation) in the main process.

Design notes:
- Uses ``multiprocessing.Pool.imap`` to preserve original FASTQ read order.
- ``PreparedAnchor`` objects contain ``re.Pattern`` and simple dataclasses,
  both of which are pickle-safe.
"""

from __future__ import annotations

import csv
import json
import logging
from collections import Counter, defaultdict
from itertools import islice
from multiprocessing import Pool
from pathlib import Path
from typing import Any

from .anchors import AnchorHit, find_anchor_hits, get_matcher_status, prepare_anchors, reverse_complement
from .barcode import BarcodeAccumulator, evaluate_barcode_recoverability, prepare_barcode_matcher
from .classify import classify_best_orientation
from .config import AppConfig
from .concatemers import split_complete_cycles
from .export import export_structure_classification_to_csv, export_summary_to_csv
from .fastq import FastqRead, fastq_output_path, open_output_text, read_fastq
from .segments import SegmentAccumulator, observe_segments
from .multiqc import write_multiqc_summary
from .report import (
    HEATMAP_BINS,
    READ_QSCORE_SCALE,
    SCATTER_LENGTH_BIN_BP,
    SCATTER_QSCORE_BIN,
    _best_hits_by_anchor,
    _build_qc_verdicts,
    _compute_terminal_offsets,
    _derive_qc_flags,
    _format_num,
    _primary_qc_bucket,
    _weighted_median,
    _write_html_report,
    _write_text,
    compute_n50_counts,
)

logger = logging.getLogger(__name__)

# Number of reads per chunk sent to a worker process.
_CHUNK_SIZE = 512


class _SerialPool:
    """Pool-compatible executor used for deterministic single-process runs."""

    @staticmethod
    def imap(function: Any, iterable: Any) -> Any:
        return map(function, iterable)

    @staticmethod
    def close() -> None:
        return None

    @staticmethod
    def terminate() -> None:
        return None

    @staticmethod
    def join() -> None:
        return None


def _process_chunk(
    args: tuple[list[tuple[str, str, str, str]], list[Any], AppConfig],
) -> list[dict[str, Any]]:
    """Worker function: process a chunk of reads.

    Receives a list of (name, sequence, quality, qscore_method) tuples,
    prepared anchors, and config.  Returns per-read result dicts.
    """
    raw_reads, prepared_anchors, config = args
    results: list[dict[str, Any]] = []
    for name, sequence, quality, qscore_method in raw_reads:
        read = FastqRead(name=name, sequence=sequence, quality=quality, qscore_method=qscore_method)
        forward_hits = find_anchor_hits(read.sequence, prepared_anchors)
        reverse_hits = find_anchor_hits(reverse_complement(read.sequence), prepared_anchors)
        structure, hits = classify_best_orientation(
            forward_hits, reverse_hits, config.structure, read_length=read.length
        )
        best_hits = _best_hits_by_anchor(hits)
        five_prime_offset, three_prime_offset = _compute_terminal_offsets(
            best_hits, config.structure.expected_order, read.length
        )
        oriented_sequence = reverse_complement(read.sequence) if structure.is_reversed else read.sequence
        oriented_quality = read.quality[::-1] if structure.is_reversed else read.quality
        segment_observations = observe_segments(
            oriented_sequence, oriented_quality, best_hits, config.structure.segments
        )
        qc_flags = _derive_qc_flags(
            read,
            structure.label,
            five_prime_offset,
            three_prime_offset,
            config,
            hits,
            segment_observations,
        )
        qc_bucket = _primary_qc_bucket(qc_flags)
        orientation = "reversed" if structure.is_reversed else "forward"
        results.append({
            "name": read.name,
            "sequence": read.sequence,
            "quality": read.quality,
            "length": read.length,
            "phred_scores": read.phred_scores,
            "read_qscore": read.read_qscore,
            "n_fraction": read.n_fraction,
            "invalid_base_fraction": read.invalid_base_fraction,
            "structure_label": structure.label,
            "structure_order_valid": structure.order_valid,
            "structure_detected_anchors": structure.detected_anchors,
            "structure_is_reversed": structure.is_reversed,
            "hits": hits,
            "best_hits": best_hits,
            "five_prime_offset": five_prime_offset,
            "three_prime_offset": three_prime_offset,
            "qc_flags": qc_flags,
            "qc_bucket": qc_bucket,
            "orientation": orientation,
            "segment_observations": segment_observations,
            "oriented_sequence": oriented_sequence,
            "oriented_quality": oriented_quality,
        })
    return results


def _chunked_reads(
    fastq_path: str, qscore_method: str, chunk_size: int,
) -> Any:
    """Yield chunks of (name, sequence, quality, qscore_method) tuples from a FASTQ file."""
    reader = read_fastq(fastq_path, qscore_method=qscore_method)
    while True:
        chunk = list(islice(reader, chunk_size))
        if not chunk:
            break
        yield [(r.name, r.sequence, r.quality, r.qscore_method) for r in chunk]


def run_qc_parallel(
    fastq_path: str,
    config: AppConfig,
    outdir: str,
    export_csv: bool = False,
    output_passed_fastq: bool = False,
    output_failed_fastq: bool = False,
    threads: int = 2,
    gzip_output: bool = False,
) -> dict[str, Any]:
    """Parallel version of run_qc for a single sample.

    Distributes per-read anchor search / classification across *threads*
    worker processes while the main process handles aggregation and I/O.
    """
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)
    prepared_anchors = prepare_anchors(config.anchors)
    barcode_matcher = prepare_barcode_matcher(config.barcode)
    barcode_accumulator = BarcodeAccumulator() if config.barcode.enabled else None

    # ---- accumulators (same as run_qc) ----
    length_counts: Counter[int] = Counter()
    base_qscore_counts: Counter[int] = Counter()
    read_qscore_counts: Counter[int] = Counter()
    scatter_counts: Counter[tuple[int, int]] = Counter()
    anchor_counts: Counter[str] = Counter()
    anchor_best_mismatches: dict[str, list[int]] = defaultdict(list)
    anchor_best_substitutions: dict[str, list[int]] = defaultdict(list)
    anchor_best_insertions: dict[str, list[int]] = defaultdict(list)
    anchor_best_deletions: dict[str, list[int]] = defaultdict(list)
    anchor_hit_counts: dict[str, list[int]] = defaultdict(list)
    anchor_positions: dict[str, dict[str, float]] = defaultdict(lambda: {"sum": 0.0, "count": 0})
    heatmap_groups: dict[tuple[str, str], dict[str, Any]] = {}
    segment_accumulators = {
        segment.name: SegmentAccumulator() for segment in config.structure.segments
    }
    structure_counts: Counter[str] = Counter()
    structure_orientation_counts: Counter[tuple[str, str]] = Counter()
    qc_bucket_counts: Counter[str] = Counter()
    total_reads = total_bases = long_high_quality = order_valid_count = reversed_read_count = 0
    five_prime_truncation_count = three_prime_truncation_count = high_n_read_count = invalid_base_read_count = 0
    passed_reads_exported = failed_reads_exported = 0
    concatemer_parent_reads = split_molecules_exported = 0

    # ---- optional CSV / FASTQ handles ----
    read_details_handle = anchor_hits_handle = None
    read_details_writer = anchor_hits_writer = None
    passed_fastq_handle = failed_fastq_handle = None
    split_fastq_handle = None

    if export_csv:
        read_details_handle = open(out_path / "read_details.csv", "w", newline="", encoding="utf-8")
        read_details_writer = csv.writer(read_details_handle)
        read_details_writer.writerow([
            "Read ID", "Read Length", "Read Qscore", "Structure Label", "QC Bucket",
            "QC Flags", "Is Reversed", "Order Valid", "5p Offset", "3p Offset",
            "N Fraction", "Invalid Base Fraction", "Anchor Count", "Segment Statuses",
            "Barcode Status", "Raw Barcode", "Raw UMI", "Barcode Nearest Distance",
            "Barcode Second Distance",
        ])
        anchor_hits_handle = open(out_path / "anchor_hits.csv", "w", newline="", encoding="utf-8")
        anchor_hits_writer = csv.writer(anchor_hits_handle)
        anchor_hits_writer.writerow([
            "Read ID", "Read Length", "Read Qscore", "Structure Label", "QC Bucket",
            "Is Reversed", "Anchor Name", "Start", "End", "Mismatches",
            "Substitutions", "Insertions", "Deletions", "CIGAR", "Matched Sequence",
        ])

    fastq_export_handles: dict[str, Any] = {}

    try:
        if output_passed_fastq:
            passed_fastq_handle = open_output_text(fastq_output_path(out_path, "passed", gzip_output))
        if output_failed_fastq:
            failed_fastq_handle = open_output_text(fastq_output_path(out_path, "failed", gzip_output))
        if config.structure.split_concatemers:
            split_fastq_handle = open_output_text(
                fastq_output_path(out_path, "split_concatemers", gzip_output)
            )

        logger.info("Starting parallel QC with %d worker processes", threads)

        pool: Any = Pool(processes=threads) if threads > 1 else _SerialPool()
        try:
            chunk_iter = _chunked_reads(fastq_path, config.qscore_method, _CHUNK_SIZE)
            task_iter = ((chunk, prepared_anchors, config) for chunk in chunk_iter)
            for chunk_results in pool.imap(_process_chunk, task_iter):
                for r in chunk_results:
                    total_reads += 1
                    total_bases += r["length"]
                    length_counts[r["length"]] += 1
                    base_qscore_counts.update(r["phred_scores"])
                    read_qscore_key = int(round(r["read_qscore"] * READ_QSCORE_SCALE))
                    read_qscore_counts[read_qscore_key] += 1
                    scatter_length = int(round(r["length"] / SCATTER_LENGTH_BIN_BP) * SCATTER_LENGTH_BIN_BP)
                    scatter_qscore = round(round(r["read_qscore"] / SCATTER_QSCORE_BIN) * SCATTER_QSCORE_BIN, 2)
                    scatter_counts[(scatter_length, scatter_qscore)] += 1

                    structure_label = r["structure_label"]
                    orientation = r["orientation"]
                    qc_bucket = r["qc_bucket"]
                    qc_flags = r["qc_flags"]
                    segment_observations = r["segment_observations"]
                    hits: list[AnchorHit] = r["hits"]
                    best_hits: dict[str, AnchorHit] = r["best_hits"]
                    barcode_observation = None
                    if barcode_accumulator is not None:
                        barcode_observation = evaluate_barcode_recoverability(
                            r["oriented_sequence"],
                            r["oriented_quality"],
                            best_hits,
                            config.barcode,
                            barcode_matcher,
                        )
                        barcode_accumulator.update(barcode_observation)
                    five_prime_offset = r["five_prime_offset"]
                    three_prime_offset = r["three_prime_offset"]

                    structure_counts[structure_label] += 1
                    structure_orientation_counts[(structure_label, orientation)] += 1
                    qc_bucket_counts[qc_bucket] += 1
                    for name, observation in segment_observations.items():
                        segment_accumulators[name].update(observation)
                    if r["structure_order_valid"]:
                        order_valid_count += 1
                    if r["structure_is_reversed"]:
                        reversed_read_count += 1
                    if r["length"] >= config.thresholds.long_read_min_bp and r["read_qscore"] >= config.thresholds.long_read_min_q:
                        long_high_quality += 1
                    if five_prime_offset is not None and five_prime_offset > config.thresholds.terminal_anchor_max_offset:
                        five_prime_truncation_count += 1
                    if three_prime_offset is not None and three_prime_offset > config.thresholds.terminal_anchor_max_offset:
                        three_prime_truncation_count += 1
                    if r["n_fraction"] >= config.thresholds.high_n_fraction:
                        high_n_read_count += 1
                    if r["invalid_base_fraction"] > 0:
                        invalid_base_read_count += 1

                    # heatmap
                    heatmap_key = (qc_bucket, structure_label)
                    heatmap_group = heatmap_groups.setdefault(
                        heatmap_key,
                        {"count": 0, "anchor_cover": [0] * HEATMAP_BINS},
                    )
                    heatmap_group["count"] += 1
                    if r["length"] > 0 and hits:
                        covered = [False] * HEATMAP_BINS
                        for hit in hits:
                            start = min(HEATMAP_BINS - 1, max(0, int((hit.start / r["length"]) * HEATMAP_BINS)))
                            end = min(HEATMAP_BINS, max(start + 1, int((hit.end / r["length"]) * HEATMAP_BINS)))
                            for idx in range(start, end):
                                covered[idx] = True
                        for idx, has_cover in enumerate(covered):
                            if has_cover:
                                heatmap_group["anchor_cover"][idx] += 1

                    # anchor stats
                    for name, hit in best_hits.items():
                        anchor_counts[name] += 1
                        anchor_best_mismatches[name].append(hit.mismatches)
                        anchor_best_substitutions[name].append(hit.substitutions)
                        anchor_best_insertions[name].append(hit.insertions)
                        anchor_best_deletions[name].append(hit.deletions)
                        anchor_positions[name]["sum"] += hit.start / r["length"] if r["length"] > 0 else 0.0
                        anchor_positions[name]["count"] += 1
                        anchor_hit_counts[name].append(sum(1 for item in hits if item.anchor_name == name))

                    # CSV export
                    if read_details_writer is not None:
                        read_details_writer.writerow([
                            r["name"], r["length"], f"{r['read_qscore']:.2f}",
                            structure_label, qc_bucket, ";".join(qc_flags),
                            r["structure_is_reversed"], r["structure_order_valid"],
                            f"{five_prime_offset:.4f}" if five_prime_offset is not None else "",
                            f"{three_prime_offset:.4f}" if three_prime_offset is not None else "",
                            f"{r['n_fraction']:.4f}", f"{r['invalid_base_fraction']:.4f}",
                            len(hits),
                            ";".join(
                                f"{name}:{observation.status}"
                                for name, observation in segment_observations.items()
                            ),
                            barcode_observation.status if barcode_observation else "",
                            barcode_observation.raw_barcode if barcode_observation and barcode_observation.raw_barcode else "",
                            barcode_observation.raw_umi if barcode_observation and barcode_observation.raw_umi else "",
                            barcode_observation.nearest_distance if barcode_observation and barcode_observation.nearest_distance is not None else "",
                            barcode_observation.second_distance if barcode_observation and barcode_observation.second_distance is not None else "",
                        ])
                    if anchor_hits_writer is not None:
                        for hit in hits:
                            anchor_hits_writer.writerow([
                                r["name"], r["length"], f"{r['read_qscore']:.2f}",
                                structure_label, qc_bucket, r["structure_is_reversed"],
                                hit.anchor_name, hit.start, hit.end, hit.mismatches,
                                hit.substitutions, hit.insertions, hit.deletions, hit.cigar,
                                hit.matched_sequence,
                            ])

                    # structure-based FASTQ export
                    if config.export_non_full_structure_fastq and structure_label != "full_structure":
                        handle = fastq_export_handles.get(structure_label)
                        if handle is None:
                            handle = open_output_text(
                                fastq_output_path(out_path, structure_label, gzip_output)
                            )
                            fastq_export_handles[structure_label] = handle
                        handle.write(f"@{r['name']}\n{r['sequence']}\n+\n{r['quality']}\n")

                    # passed / failed FASTQ export
                    if qc_bucket == "pass" and passed_fastq_handle is not None:
                        passed_fastq_handle.write(f"@{r['name']}\n{r['sequence']}\n+\n{r['quality']}\n")
                        passed_reads_exported += 1
                    if qc_bucket != "pass" and failed_fastq_handle is not None:
                        failed_fastq_handle.write(f"@{r['name']}\n{r['sequence']}\n+\n{r['quality']}\n")
                        failed_reads_exported += 1
                    if structure_label == "concatemer_candidate" and split_fastq_handle is not None:
                        molecules = split_complete_cycles(
                            r["name"],
                            r["oriented_sequence"],
                            r["oriented_quality"],
                            hits,
                            config.structure.expected_order,
                        )
                        if molecules:
                            concatemer_parent_reads += 1
                        for molecule in molecules:
                            split_fastq_handle.write(
                                f"@{molecule.read_id}\n{molecule.sequence}\n+\n{molecule.quality}\n"
                            )
                            split_molecules_exported += 1

            pool.close()
        except Exception:
            pool.terminate()
            raise
        finally:
            pool.join()

    finally:
        if read_details_handle is not None:
            read_details_handle.close()
        if anchor_hits_handle is not None:
            anchor_hits_handle.close()
        for handle in fastq_export_handles.values():
            handle.close()
        if passed_fastq_handle is not None:
            passed_fastq_handle.close()
        if failed_fastq_handle is not None:
            failed_fastq_handle.close()
        if split_fastq_handle is not None:
            split_fastq_handle.close()

    # ---- summary (same as run_qc) ----
    anchor_quality_stats: dict[str, dict[str, float]] = {}
    for anchor in config.anchors:
        mismatch_values = anchor_best_mismatches.get(anchor.name, [])
        substitution_values = anchor_best_substitutions.get(anchor.name, [])
        insertion_values = anchor_best_insertions.get(anchor.name, [])
        deletion_values = anchor_best_deletions.get(anchor.name, [])
        count_values = anchor_hit_counts.get(anchor.name, [])
        anchor_quality_stats[anchor.name] = {
            "detection_ratio": (anchor_counts.get(anchor.name, 0) / total_reads) if total_reads else 0.0,
            "mean_best_mismatches": (sum(mismatch_values) / len(mismatch_values)) if mismatch_values else 0.0,
            "mean_best_substitutions": (sum(substitution_values) / len(substitution_values)) if substitution_values else 0.0,
            "mean_best_insertions": (sum(insertion_values) / len(insertion_values)) if insertion_values else 0.0,
            "mean_best_deletions": (sum(deletion_values) / len(deletion_values)) if deletion_values else 0.0,
            "multi_hit_ratio": (sum(1 for value in count_values if value > 1) / len(count_values)) if count_values else 0.0,
            "mean_hits_per_positive_read": (sum(count_values) / len(count_values)) if count_values else 0.0,
        }

    qc_bucket_counts_dict = dict(qc_bucket_counts.most_common())
    summary = {
        "sample_name": config.sample_name,
        "fastq_path": str(fastq_path),
        "config_summary": {
            "anchors": [anchor.name for anchor in config.anchors],
            "expected_order": list(config.structure.expected_order),
            "expected_orientation": config.structure.expected_orientation,
            "segments": [segment.name for segment in config.structure.segments],
            "barcode_qc_enabled": config.barcode.enabled,
            "export_non_full_structure_fastq": config.export_non_full_structure_fastq,
            "long_read_min_bp": config.thresholds.long_read_min_bp,
            "long_read_min_q": config.thresholds.long_read_min_q,
            "terminal_anchor_max_offset": config.thresholds.terminal_anchor_max_offset,
            "high_n_fraction": config.thresholds.high_n_fraction,
        },
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mean_read_length": (total_bases / total_reads) if total_reads else 0.0,
        "median_read_length": _weighted_median(length_counts),
        "n50": compute_n50_counts(length_counts),
        "min_read_length": min(length_counts) if length_counts else 0,
        "max_read_length": max(length_counts) if length_counts else 0,
        "median_read_qscore": _weighted_median(read_qscore_counts, READ_QSCORE_SCALE),
        "long_high_quality_ratio": (long_high_quality / total_reads) if total_reads else 0.0,
        "reversed_read_ratio": (reversed_read_count / total_reads) if total_reads else 0.0,
        "five_prime_truncation_ratio": (five_prime_truncation_count / total_reads) if total_reads else 0.0,
        "three_prime_truncation_ratio": (three_prime_truncation_count / total_reads) if total_reads else 0.0,
        "high_n_read_ratio": (high_n_read_count / total_reads) if total_reads else 0.0,
        "invalid_base_read_ratio": (invalid_base_read_count / total_reads) if total_reads else 0.0,
        "anchor_detection_ratio": {anchor.name: anchor_quality_stats[anchor.name]["detection_ratio"] for anchor in config.anchors},
        "anchor_quality_stats": anchor_quality_stats,
        "segment_qc": {
            name: accumulator.summary() for name, accumulator in segment_accumulators.items()
        },
        "barcode_qc": barcode_accumulator.summary() if barcode_accumulator is not None else None,
        "structure_counts": dict(structure_counts),
        "structure_orientation_counts": {label: {orientation: structure_orientation_counts.get((label, orientation), 0) for orientation in ("forward", "reversed") if structure_orientation_counts.get((label, orientation), 0) > 0} for label in structure_counts},
        "correct_anchor_order_ratio": (order_valid_count / total_reads) if total_reads else 0.0,
        "qc_bucket_counts": qc_bucket_counts_dict,
        "qc_bucket_ratios": {bucket: (count / total_reads) if total_reads else 0.0 for bucket, count in qc_bucket_counts_dict.items()},
        "matcher_backend": get_matcher_status(config.anchors),
        "passed_reads_exported": passed_reads_exported if output_passed_fastq else None,
        "failed_reads_exported": failed_reads_exported if output_failed_fastq else None,
        "concatemer_parent_reads": concatemer_parent_reads,
        "split_molecules_exported": (
            split_molecules_exported if config.structure.split_concatemers else None
        ),
    }
    summary["qc_verdicts"] = _build_qc_verdicts(
        summary, config.thresholds, config.structure.expected_orientation
    )
    _write_text(out_path / "summary.json", json.dumps(summary, indent=2))
    write_multiqc_summary(summary, out_path / "anchorscope_mqc.json")
    if export_csv:
        export_summary_to_csv(summary, out_path / "summary.csv")
        export_structure_classification_to_csv(summary, out_path / "structure_classification.csv")

    from .report import _value_counts_to_float_map
    _write_html_report(
        out_path / "report.html",
        summary,
        length_counts,
        base_qscore_counts,
        read_qscore_counts,
        scatter_counts,
        anchor_positions,
        heatmap_groups,
    )
    logger.info("Parallel QC completed: %d reads processed with %d workers", total_reads, threads)
    return summary
