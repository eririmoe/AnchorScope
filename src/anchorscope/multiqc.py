from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _metric(value: Any, suffix: str = "", fmt: str = "{:.3g}") -> dict[str, Any]:
    return {"value": value, "suffix": suffix, "format": fmt}


def write_multiqc_summary(summary: dict[str, Any], output_path: Path) -> None:
    sample = summary["sample_name"]
    metrics: dict[str, Any] = {
        "total_reads": summary["total_reads"],
        "median_length": summary["median_read_length"],
        "median_qscore": summary["median_read_qscore"],
        "full_structure": summary["correct_anchor_order_ratio"],
        "pass_reads": summary["qc_bucket_ratios"].get("pass", 0.0),
        "no_anchor": summary["qc_bucket_ratios"].get("no_anchor", 0.0),
        "concatemer": summary["qc_bucket_ratios"].get("concatemer_candidate", 0.0),
    }
    barcode = summary.get("barcode_qc")
    if barcode:
        metrics["barcode_recoverable"] = barcode["recoverable_ratio"]
        metrics["barcode_ambiguous"] = barcode["ambiguous_ratio"]
    payload = {
        "id": "anchorscope_general_stats",
        "section_name": "AnchorScope",
        "description": "Anchor-aware long-read library structure quality control.",
        "plot_type": "generalstats",
        "pconfig": {"id": "anchorscope_general_stats", "title": "AnchorScope"},
        "headers": {
            "total_reads": _metric(None, fmt="{:.0f}"),
            "median_length": _metric(None, suffix=" bp", fmt="{:.0f}"),
            "median_qscore": _metric(None, fmt="{:.2f}"),
            "full_structure": _metric(None, suffix="%", fmt="{:.1%}"),
            "pass_reads": _metric(None, suffix="%", fmt="{:.1%}"),
            "no_anchor": _metric(None, suffix="%", fmt="{:.1%}"),
            "concatemer": _metric(None, suffix="%", fmt="{:.1%}"),
            "barcode_recoverable": _metric(None, suffix="%", fmt="{:.1%}"),
            "barcode_ambiguous": _metric(None, suffix="%", fmt="{:.1%}"),
        },
        "data": {sample: metrics},
    }
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_tag_audit_multiqc(summary: dict[str, Any], output_path: Path) -> None:
    sample = summary["sample_name"]
    payload = {
        "id": "anchorscope_tag_audit",
        "section_name": "AnchorScope barcode-tag audit",
        "description": "Audit of raw/corrected barcode and UMI tags in an aligned or unaligned BAM.",
        "plot_type": "generalstats",
        "pconfig": {"id": "anchorscope_tag_audit", "title": "AnchorScope tag audit"},
        "data": {
            sample: {
                "reads": summary["total_reads"],
                "corrected_barcode_rate": summary["corrected_barcode_ratio"],
                "corrected_umi_rate": summary["corrected_umi_ratio"],
                "ambiguous_barcode_rate": summary["ambiguous_barcode_ratio"],
                "poly_a_tag_rate": summary["poly_a_tag_ratio"],
            }
        },
    }
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
