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
        writer.writerow(["Rust Accelerator", summary.get("rust_accelerator", {}).get("mode", "python-fallback")])

        writer.writerow([])
        writer.writerow(["Anchor Detection Ratios"])
        for anchor_name, ratio in summary.get("anchor_detection_ratio", {}).items():
            writer.writerow([anchor_name, f"{ratio:.4f}"])

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


def export_anchor_hits_to_csv(read_results: list[dict[str, Any]], output_path: Path) -> None:
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Read ID",
                "Read Length",
                "Read Qscore",
                "Structure Label",
                "QC Bucket",
                "Is Reversed",
                "Anchor Name",
                "Start",
                "End",
                "Mismatches",
                "Matched Sequence",
            ]
        )

        for result in read_results:
            for hit in result.get("hits", []):
                writer.writerow(
                    [
                        result.get("read_id", ""),
                        result.get("length", 0),
                        f"{result.get('read_qscore', 0):.2f}",
                        result.get("structure_label", ""),
                        result.get("qc_bucket", ""),
                        result.get("is_reversed", False),
                        hit.get("anchor_name", ""),
                        hit.get("start", 0),
                        hit.get("end", 0),
                        hit.get("mismatches", 0),
                        hit.get("matched_sequence", ""),
                    ]
                )


def export_read_results_to_csv(read_results: list[dict[str, Any]], output_path: Path) -> None:
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Read ID",
                "Read Length",
                "Read Qscore",
                "Structure Label",
                "QC Bucket",
                "QC Flags",
                "Is Reversed",
                "Order Valid",
                "5p Offset",
                "3p Offset",
                "N Fraction",
                "Invalid Base Fraction",
                "Anchor Count",
            ]
        )
        for result in read_results:
            writer.writerow(
                [
                    result.get("read_id", ""),
                    result.get("length", 0),
                    f"{result.get('read_qscore', 0):.2f}",
                    result.get("structure_label", ""),
                    result.get("qc_bucket", ""),
                    ";".join(result.get("qc_flags", [])),
                    result.get("is_reversed", False),
                    result.get("order_valid", False),
                    f"{result.get('five_prime_offset', 0):.4f}" if result.get("five_prime_offset") is not None else "",
                    f"{result.get('three_prime_offset', 0):.4f}" if result.get("three_prime_offset") is not None else "",
                    f"{result.get('n_fraction', 0):.4f}",
                    f"{result.get('invalid_base_fraction', 0):.4f}",
                    result.get("anchor_count", 0),
                ]
            )


def export_structure_classification_to_csv(summary: dict[str, Any], output_path: Path) -> None:
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Structure Class", "Orientation", "Count", "Percentage"])

        total_reads = summary.get("total_reads", 0)
        for structure, orientations in summary.get("structure_orientation_counts", {}).items():
            for orientation, count in orientations.items():
                percentage = (count / total_reads * 100) if total_reads > 0 else 0
                writer.writerow([structure, orientation, count, f"{percentage:.2f}%"])
