from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def export_summary_to_csv(summary: dict[str, Any], output_path: Path) -> None:
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Sample Name", summary.get("sample_name", "")])
        writer.writerow(["Total Reads", summary.get("total_reads", 0)])
        writer.writerow(["Total Bases", summary.get("total_bases", 0)])
        writer.writerow(["Mean Read Length", f"{summary.get('mean_read_length', 0):.2f}"])
        writer.writerow(["Median Read Length", summary.get("median_read_length", 0)])
        writer.writerow(["N50", summary.get("n50", 0)])
        writer.writerow(["Min Read Length", summary.get("min_read_length", 0)])
        writer.writerow(["Max Read Length", summary.get("max_read_length", 0)])
        writer.writerow(["Median Read Qscore", f"{summary.get('median_read_qscore', 0):.2f}"])
        writer.writerow(["Long High-Quality Ratio", f"{summary.get('long_high_quality_ratio', 0):.4f}"])
        writer.writerow(["Reversed Read Ratio", f"{summary.get('reversed_read_ratio', 0):.4f}"])
        writer.writerow(["Correct Anchor Order Ratio", f"{summary.get('correct_anchor_order_ratio', 0):.4f}"])
        writer.writerow(["5' Truncation Ratio", f"{summary.get('five_prime_truncation_ratio', 0):.4f}"])
        writer.writerow(["3' Truncation Ratio", f"{summary.get('three_prime_truncation_ratio', 0):.4f}"])
        writer.writerow(["High-N Read Ratio", f"{summary.get('high_n_read_ratio', 0):.4f}"])
        writer.writerow(["Matcher Backend", summary.get("matcher_backend", summary.get("rust_accelerator", {})).get("mode", "unknown")])

        writer.writerow([])
        writer.writerow(["Anchor Detection Ratios"])
        for anchor_name, ratio in summary.get("anchor_detection_ratio", {}).items():
            writer.writerow([anchor_name, f"{ratio:.4f}"])

        writer.writerow([])
        writer.writerow(["Anchor Edit Statistics"])
        writer.writerow(["Anchor", "Mean Edits", "Mean Substitutions", "Mean Insertions", "Mean Deletions"])
        for anchor_name, values in summary.get("anchor_quality_stats", {}).items():
            writer.writerow([
                anchor_name,
                f"{values.get('mean_best_mismatches', 0):.4f}",
                f"{values.get('mean_best_substitutions', 0):.4f}",
                f"{values.get('mean_best_insertions', 0):.4f}",
                f"{values.get('mean_best_deletions', 0):.4f}",
            ])

        if summary.get("segment_qc"):
            writer.writerow([])
            writer.writerow(["Segment QC"])
            writer.writerow(["Segment", "Observed Ratio", "Median Length", "Median Qscore", "Mean GC", "Within Bounds"])
            for name, values in summary["segment_qc"].items():
                writer.writerow([
                    name,
                    f"{values['observed_ratio']:.4f}",
                    f"{values['median_length']:.2f}",
                    f"{values['median_qscore']:.2f}",
                    f"{values['mean_gc_fraction']:.4f}",
                    f"{values['within_length_ratio']:.4f}",
                ])

        if summary.get("barcode_qc"):
            writer.writerow([])
            writer.writerow(["Barcode/UMI Recoverability"])
            for key in (
                "observed_ratio",
                "recoverable_ratio",
                "ambiguous_ratio",
                "low_quality_ratio",
                "unique_raw_barcodes",
                "unique_raw_umis",
                "umi_singleton_ratio",
            ):
                writer.writerow([key, summary["barcode_qc"].get(key, "")])

        writer.writerow([])
        writer.writerow(["QC Bucket Counts"])
        for reason, count in summary.get("qc_bucket_counts", {}).items():
            writer.writerow([reason, count])

        writer.writerow([])
        writer.writerow(["Structure Classification Counts"])
        for structure, count in summary.get("structure_counts", {}).items():
            writer.writerow([structure, count])

        writer.writerow([])
        writer.writerow(["Structure Orientation Counts"])
        for structure, orientations in summary.get("structure_orientation_counts", {}).items():
            for orientation, count in orientations.items():
                writer.writerow([f"{structure} ({orientation})", count])


def export_structure_classification_to_csv(summary: dict[str, Any], output_path: Path) -> None:
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Structure Class", "Orientation", "Count", "Percentage"])

        total_reads = summary.get("total_reads", 0)
        for structure, orientations in summary.get("structure_orientation_counts", {}).items():
            for orientation, count in orientations.items():
                percentage = (count / total_reads * 100) if total_reads > 0 else 0
                writer.writerow([structure, orientation, count, f"{percentage:.2f}%"])
