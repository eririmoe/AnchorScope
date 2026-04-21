from __future__ import annotations

import json
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from html import escape
from math import ceil, floor
from pathlib import Path
from statistics import median
from typing import Any

from .anchors import AnchorHit, find_anchor_hits, get_rust_status, prepare_anchors, reverse_complement
from .classify import classify_best_orientation
from .config import AppConfig, ThresholdConfig
from .export import export_anchor_hits_to_csv, export_read_results_to_csv, export_structure_classification_to_csv, export_summary_to_csv
from .fastq import FastqRead, read_fastq

logger = logging.getLogger(__name__)


@dataclass
class ReadResult:
    read_id: str
    length: int
    read_qscore: float
    structure_label: str
    is_reversed: bool
    order_valid: bool
    hits: list[AnchorHit]
    qc_bucket: str
    qc_flags: list[str]
    five_prime_offset: float | None
    three_prime_offset: float | None
    n_fraction: float
    invalid_base_fraction: float


def compute_n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= total / 2:
            return length
    return 0


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _format_num(value: float, decimals: int = 1) -> str:
    if float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:,.{decimals}f}"


def _best_hits_by_anchor(hits: list[AnchorHit]) -> dict[str, AnchorHit]:
    best_hits: dict[str, AnchorHit] = {}
    for hit in hits:
        current = best_hits.get(hit.anchor_name)
        if current is None or (hit.mismatches, hit.start, hit.end) < (current.mismatches, current.start, current.end):
            best_hits[hit.anchor_name] = hit
    return best_hits


def _compute_terminal_offsets(best_hits: dict[str, AnchorHit], expected_order: list[str], read_length: int) -> tuple[float | None, float | None]:
    if read_length <= 0 or not expected_order:
        return None, None
    first = best_hits.get(expected_order[0])
    last = best_hits.get(expected_order[-1])
    five_prime = first.start / read_length if first else None
    three_prime = (read_length - last.end) / read_length if last else None
    return five_prime, three_prime


def _derive_qc_flags(
    read: FastqRead,
    structure_label: str,
    five_prime_offset: float | None,
    three_prime_offset: float | None,
    config: AppConfig,
    hits: list[AnchorHit],
) -> list[str]:
    flags: list[str] = []
    t = config.thresholds
    if read.length == 0:
        flags.append("empty_read")
    if read.length < t.long_read_min_bp:
        flags.append("too_short")
    if read.read_qscore < t.long_read_min_q:
        flags.append("low_read_q")
    if read.n_fraction >= t.high_n_fraction:
        flags.append("high_n")
    if read.invalid_base_fraction > 0:
        flags.append("invalid_bases")

    structure_map = {
        "no_anchor_detected": "no_anchor",
        "missing_5p_anchor": "missing_5p",
        "missing_3p_anchor": "missing_3p",
        "anchor_order_invalid": "order_invalid",
        "duplicated_anchor": "multi_anchor",
        "internal_5p_anchor": "internal_adapter",
        "internal_3p_anchor": "internal_adapter",
        "concatemer_candidate": "concatemer_candidate",
    }
    mapped = structure_map.get(structure_label)
    if mapped:
        flags.append(mapped)
    if five_prime_offset is not None and five_prime_offset > t.terminal_anchor_max_offset:
        flags.append("five_prime_truncated")
    if three_prime_offset is not None and three_prime_offset > t.terminal_anchor_max_offset:
        flags.append("three_prime_truncated")
    if len(hits) > max(len(config.structure.expected_order), 1) * 2:
        flags.append("dense_anchor_hits")
    return list(dict.fromkeys(flags))


def _primary_qc_bucket(flags: list[str]) -> str:
    priority = [
        "empty_read", "too_short", "low_read_q", "invalid_bases", "high_n", "no_anchor",
        "internal_adapter", "concatemer_candidate", "multi_anchor", "order_invalid",
        "missing_5p", "missing_3p", "five_prime_truncated", "three_prime_truncated", "dense_anchor_hits",
    ]
    for item in priority:
        if item in flags:
            return item
    return "pass"


def _serialize_read_result(result: ReadResult) -> dict[str, Any]:
    return {
        "read_id": result.read_id,
        "length": result.length,
        "read_qscore": result.read_qscore,
        "structure_label": result.structure_label,
        "qc_bucket": result.qc_bucket,
        "qc_flags": result.qc_flags,
        "is_reversed": result.is_reversed,
        "order_valid": result.order_valid,
        "anchor_count": len(result.hits),
        "five_prime_offset": result.five_prime_offset,
        "three_prime_offset": result.three_prime_offset,
        "n_fraction": result.n_fraction,
        "invalid_base_fraction": result.invalid_base_fraction,
        "hits": [{"anchor_name": hit.anchor_name, "start": hit.start, "end": hit.end, "mismatches": hit.mismatches, "matched_sequence": hit.matched_sequence} for hit in result.hits],
    }


def _verdict(value: float, warn: float, fail: float, higher_is_better: bool) -> str:
    if higher_is_better:
        if value < fail:
            return "fail"
        if value < warn:
            return "warn"
        return "pass"
    if value > fail:
        return "fail"
    if value > warn:
        return "warn"
    return "pass"


def _build_qc_verdicts(summary: dict[str, Any], t: ThresholdConfig) -> dict[str, dict[str, Any]]:
    verdicts = {
        "long_high_quality_ratio": {"label": "Long high-quality ratio", "value": summary["long_high_quality_ratio"], "status": _verdict(summary["long_high_quality_ratio"], t.warn_long_high_quality_ratio, t.fail_long_high_quality_ratio, True)},
        "correct_anchor_order_ratio": {"label": "Correct anchor order ratio", "value": summary["correct_anchor_order_ratio"], "status": _verdict(summary["correct_anchor_order_ratio"], t.warn_correct_anchor_order_ratio, t.fail_correct_anchor_order_ratio, True)},
        "no_anchor_ratio": {"label": "No-anchor ratio", "value": summary["qc_bucket_ratios"].get("no_anchor", 0.0), "status": _verdict(summary["qc_bucket_ratios"].get("no_anchor", 0.0), t.warn_no_anchor_ratio, t.fail_no_anchor_ratio, False)},
        "reversed_read_ratio": {"label": "Reversed read ratio", "value": summary["reversed_read_ratio"], "status": _verdict(summary["reversed_read_ratio"], t.warn_reversed_read_ratio, t.fail_reversed_read_ratio, False)},
        "high_n_read_ratio": {"label": "High-N read ratio", "value": summary["high_n_read_ratio"], "status": _verdict(summary["high_n_read_ratio"], t.warn_high_n_ratio, t.fail_high_n_ratio, False)},
    }
    overall = "pass"
    if any(item["status"] == "fail" for item in verdicts.values()):
        overall = "fail"
    elif any(item["status"] == "warn" for item in verdicts.values()):
        overall = "warn"
    verdicts["overall"] = {"label": "Overall QC", "value": None, "status": overall}
    return verdicts


def _scale(value: float, lower: float, upper: float, span: float) -> float:
    if upper <= lower:
        return 0.0
    return ((value - lower) / (upper - lower)) * span


def _svg_wrapper(title: str, body: str, width: int = 780, height: int = 320) -> str:
    return f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'><style>text{{font-family:IBM Plex Sans,Segoe UI,sans-serif;fill:#12232f}}.grid{{stroke:#d7e0e6;stroke-width:1}}.axis{{stroke:#304252;stroke-width:1}}.bar{{rx:5;ry:5}}</style><text x='18' y='28' font-size='20' font-weight='700'>{escape(title)}</text>{body}</svg>"


def _histogram(values: list[float], title: str, xlabel: str, color: str = "#2364aa", bins: int = 20) -> str:
    if not values:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>")
    width, height = 780, 320
    left, top, right, bottom = 60, 48, 20, 45
    plot_w, plot_h = width - left - right, height - top - bottom
    lower, upper = min(values), max(values)
    if lower == upper:
        lower -= 0.5
        upper += 0.5
    step = (upper - lower) / bins
    counts = [0] * bins
    for value in values:
        counts[min(bins - 1, int((value - lower) / step))] += 1
    max_count = max(counts) or 1
    parts: list[str] = []
    for idx in range(6):
        y = top + plot_h - (idx / 5) * plot_h
        parts.append(f"<line class='grid' x1='{left}' y1='{y:.2f}' x2='{left+plot_w}' y2='{y:.2f}'/>")
        parts.append(f"<text x='{left-8}' y='{y+4:.2f}' text-anchor='end' font-size='11'>{_format_num((idx/5)*max_count, 0)}</text>")
    bar_w = plot_w / bins
    for idx, count in enumerate(counts):
        bar_h = (count / max_count) * plot_h
        x = left + idx * bar_w
        y = top + plot_h - bar_h
        parts.append(f"<rect class='bar' x='{x:.2f}' y='{y:.2f}' width='{max(bar_w-1,1):.2f}' height='{bar_h:.2f}' fill='{color}'/>")
    parts.append(f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>")
    parts.append(f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-12}' text-anchor='middle' font-size='12'>{escape(xlabel)}</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _line_density(values: list[float], title: str, xlabel: str, color: str = "#2d8f85", bin_width: float = 0.2) -> str:
    if not values:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>")
    width, height = 780, 320
    left, top, right, bottom = 60, 48, 20, 45
    plot_w, plot_h = width - left - right, height - top - bottom
    lower, upper = min(values), max(values)
    if lower == upper:
        lower -= 0.5
        upper += 0.5
    lower = floor(lower) - 1
    upper = ceil(upper) + 1
    bin_count = max(2, int(ceil((upper - lower) / bin_width)))
    counts = [0.0] * bin_count
    for value in values:
        counts[min(bin_count - 1, max(0, int((value - lower) / bin_width)))] += 1
    densities = [count / (len(values) * bin_width) for count in counts]
    y_max = max(densities) * 1.1 if densities else 1.0
    points = []
    for idx, density in enumerate(densities):
        x_value = lower + (idx + 0.5) * bin_width
        x = left + _scale(x_value, lower, upper, plot_w)
        y = top + plot_h - _scale(density, 0, y_max, plot_h)
        points.append(f"{x:.2f},{y:.2f}")
    body = f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/><line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/><polygon points='{left},{top+plot_h} {' '.join(points)} {left+plot_w},{top+plot_h}' fill='{color}' fill-opacity='0.16'/><polyline fill='none' stroke='{color}' stroke-width='2.5' points='{' '.join(points)}'/><text x='{width/2:.0f}' y='{height-12}' text-anchor='middle' font-size='12'>{escape(xlabel)}</text>"
    return _svg_wrapper(title, body, width, height)


def _bar_chart(items: list[tuple[str, float]], title: str, ylabel: str, color: str = "#d1495b", percent: bool = False) -> str:
    if not items:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>", 820, 360)
    width, height = 820, 360
    left, top, right, bottom = 70, 52, 24, 96
    plot_w, plot_h = width - left - right, height - top - bottom
    y_max = 1.0 if percent else max(value for _, value in items) or 1.0
    parts = [f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>", f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>"]
    slot_w = plot_w / len(items)
    bar_w = min(slot_w * 0.58, 120)
    for idx, (label, value) in enumerate(items):
        bar_h = (value / y_max) * plot_h if y_max else 0
        x = left + idx * slot_w + (slot_w - bar_w) / 2
        y = top + plot_h - bar_h
        disp = f"{value:.1%}" if percent else _format_num(value, 0)
        tx = left + idx * slot_w + slot_w / 2
        parts.append(f"<rect class='bar' x='{x:.2f}' y='{y:.2f}' width='{bar_w:.2f}' height='{bar_h:.2f}' fill='{color}'/>")
        parts.append(f"<text x='{tx:.2f}' y='{y-8:.2f}' text-anchor='middle' font-size='11'>{escape(disp)}</text>")
        parts.append(f"<text x='{tx:.2f}' y='{top+plot_h+22:.2f}' text-anchor='end' transform='rotate(-24 {tx:.2f},{top+plot_h+22:.2f})' font-size='11'>{escape(label)}</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle' font-size='12'>{escape(ylabel)}</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _scatter(lengths: list[int], read_qscores: list[float], title: str) -> str:
    width, height = 780, 320
    left, top, right, bottom = 60, 48, 20, 45
    plot_w, plot_h = width - left - right, height - top - bottom
    if not lengths or not read_qscores:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>")
    min_x, max_x = min(lengths), max(lengths)
    min_y, max_y = min(read_qscores), max(read_qscores)
    if min_x == max_x:
        max_x += 1
    if min_y == max_y:
        max_y += 1
    parts = [f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>", f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>"]
    for x_value, y_value in zip(lengths, read_qscores):
        x = left + _scale(x_value, min_x, max_x, plot_w)
        y = top + plot_h - _scale(y_value, min_y, max_y, plot_h)
        parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='3' fill='#2d8f85' fill-opacity='0.58'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-12}' text-anchor='middle' font-size='12'>Read length (bp)</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle' font-size='12'>Read Qscore</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _heatmap(read_results: list[ReadResult], title: str, max_reads: int) -> str:
    width, height = 900, 360
    left, top, right, bottom = 60, 48, 20, 36
    plot_w, plot_h = width - left - right, height - top - bottom
    sampled = read_results[:max_reads]
    if not sampled:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>", width, height)
    row_h = max(1, plot_h / len(sampled))
    col_w = plot_w / 100
    palette = {"pass": "#2d8f85", "too_short": "#d1495b", "low_read_q": "#d97706", "no_anchor": "#8d99ae"}
    parts = [f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>"]
    for row_idx, result in enumerate(sampled):
        y = top + row_idx * row_h
        base = palette.get(result.qc_bucket, "#e6eef3")
        parts.append(f"<rect x='{left}' y='{y:.2f}' width='{plot_w}' height='{row_h:.2f}' fill='{base}' fill-opacity='0.12'/>")
        for hit in result.hits:
            if result.length <= 0:
                continue
            start = min(99, int((hit.start / result.length) * 100))
            end = min(100, max(start + 1, int((hit.end / result.length) * 100)))
            parts.append(f"<rect x='{left + start*col_w:.2f}' y='{y:.2f}' width='{(end-start)*col_w:.2f}' height='{row_h:.2f}' fill='#2364aa'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-10}' text-anchor='middle' font-size='12'>Normalized read position</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _status_pill(status: str) -> str:
    return f"<span class='status-pill status-{escape(status)}'>{escape(status.upper())}</span>"


def run_qc(fastq_path: str, config: AppConfig, outdir: str, export_csv: bool = False) -> dict[str, Any]:
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)
    prepared_anchors = prepare_anchors(config.anchors)

    lengths: list[int] = []
    base_qscores: list[float] = []
    read_qscores: list[float] = []
    read_results: list[ReadResult] = []
    anchor_counts: Counter[str] = Counter()
    anchor_best_mismatches: dict[str, list[int]] = defaultdict(list)
    anchor_hit_counts: dict[str, list[int]] = defaultdict(list)
    anchor_positions: dict[str, list[float]] = defaultdict(list)
    structure_counts: Counter[str] = Counter()
    structure_orientation_counts: Counter[tuple[str, str]] = Counter()
    qc_bucket_counts: Counter[str] = Counter()
    total_reads = total_bases = long_high_quality = order_valid_count = reversed_read_count = 0
    five_prime_truncation_count = three_prime_truncation_count = high_n_read_count = invalid_base_read_count = 0

    for read in read_fastq(fastq_path, qscore_method=config.qscore_method):
        total_reads += 1
        total_bases += read.length
        lengths.append(read.length)
        base_qscores.extend(float(q) for q in read.phred_scores)
        read_qscores.append(read.read_qscore)
        forward_hits = find_anchor_hits(read.sequence, prepared_anchors)
        reverse_hits = find_anchor_hits(reverse_complement(read.sequence), prepared_anchors)
        structure, hits = classify_best_orientation(forward_hits, reverse_hits, config.structure.expected_order)
        best_hits = _best_hits_by_anchor(hits)
        five_prime_offset, three_prime_offset = _compute_terminal_offsets(best_hits, config.structure.expected_order, read.length)
        qc_flags = _derive_qc_flags(read, structure.label, five_prime_offset, three_prime_offset, config, hits)
        qc_bucket = _primary_qc_bucket(qc_flags)
        result = ReadResult(read.name, read.length, read.read_qscore, structure.label, structure.is_reversed, structure.order_valid, hits, qc_bucket, qc_flags, five_prime_offset, three_prime_offset, read.n_fraction, read.invalid_base_fraction)
        read_results.append(result)
        structure_counts[structure.label] += 1
        structure_orientation_counts[(structure.label, "reversed" if structure.is_reversed else "forward")] += 1
        qc_bucket_counts[qc_bucket] += 1
        if structure.order_valid:
            order_valid_count += 1
        if structure.is_reversed:
            reversed_read_count += 1
        if read.length >= config.thresholds.long_read_min_bp and read.read_qscore >= config.thresholds.long_read_min_q:
            long_high_quality += 1
        if five_prime_offset is not None and five_prime_offset > config.thresholds.terminal_anchor_max_offset:
            five_prime_truncation_count += 1
        if three_prime_offset is not None and three_prime_offset > config.thresholds.terminal_anchor_max_offset:
            three_prime_truncation_count += 1
        if read.n_fraction >= config.thresholds.high_n_fraction:
            high_n_read_count += 1
        if read.invalid_base_fraction > 0:
            invalid_base_read_count += 1
        for name, hit in best_hits.items():
            anchor_counts[name] += 1
            anchor_best_mismatches[name].append(hit.mismatches)
            anchor_positions[name].append(hit.start / read.length if read.length > 0 else 0.0)
            anchor_hit_counts[name].append(sum(1 for item in hits if item.anchor_name == name))

    anchor_quality_stats: dict[str, dict[str, float]] = {}
    for anchor in config.anchors:
        mismatch_values = anchor_best_mismatches.get(anchor.name, [])
        count_values = anchor_hit_counts.get(anchor.name, [])
        anchor_quality_stats[anchor.name] = {
            "detection_ratio": (anchor_counts.get(anchor.name, 0) / total_reads) if total_reads else 0.0,
            "mean_best_mismatches": (sum(mismatch_values) / len(mismatch_values)) if mismatch_values else 0.0,
            "multi_hit_ratio": (sum(1 for value in count_values if value > 1) / len(count_values)) if count_values else 0.0,
            "mean_hits_per_positive_read": (sum(count_values) / len(count_values)) if count_values else 0.0,
        }

    qc_bucket_counts_dict = dict(qc_bucket_counts.most_common())
    summary = {
        "sample_name": config.sample_name,
        "fastq_path": str(fastq_path),
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mean_read_length": (total_bases / total_reads) if total_reads else 0.0,
        "median_read_length": median(lengths) if lengths else 0,
        "n50": compute_n50(lengths),
        "min_read_length": min(lengths) if lengths else 0,
        "max_read_length": max(lengths) if lengths else 0,
        "median_read_qscore": median(read_qscores) if read_qscores else 0.0,
        "long_high_quality_ratio": (long_high_quality / total_reads) if total_reads else 0.0,
        "reversed_read_ratio": (reversed_read_count / total_reads) if total_reads else 0.0,
        "five_prime_truncation_ratio": (five_prime_truncation_count / total_reads) if total_reads else 0.0,
        "three_prime_truncation_ratio": (three_prime_truncation_count / total_reads) if total_reads else 0.0,
        "high_n_read_ratio": (high_n_read_count / total_reads) if total_reads else 0.0,
        "invalid_base_read_ratio": (invalid_base_read_count / total_reads) if total_reads else 0.0,
        "anchor_detection_ratio": {anchor.name: anchor_quality_stats[anchor.name]["detection_ratio"] for anchor in config.anchors},
        "anchor_quality_stats": anchor_quality_stats,
        "structure_counts": dict(structure_counts),
        "structure_orientation_counts": {label: {orientation: structure_orientation_counts.get((label, orientation), 0) for orientation in ("forward", "reversed") if structure_orientation_counts.get((label, orientation), 0) > 0} for label in structure_counts},
        "correct_anchor_order_ratio": (order_valid_count / total_reads) if total_reads else 0.0,
        "qc_bucket_counts": qc_bucket_counts_dict,
        "qc_bucket_ratios": {bucket: (count / total_reads) if total_reads else 0.0 for bucket, count in qc_bucket_counts_dict.items()},
        "rust_accelerator": get_rust_status(),
    }
    summary["qc_verdicts"] = _build_qc_verdicts(summary, config.thresholds)
    _write_text(out_path / "summary.json", json.dumps(summary, indent=2))
    if export_csv:
        serialized = [_serialize_read_result(result) for result in read_results]
        export_summary_to_csv(summary, out_path / "summary.csv")
        export_structure_classification_to_csv(summary, out_path / "structure_classification.csv")
        export_read_results_to_csv(serialized, out_path / "read_details.csv")
        export_anchor_hits_to_csv(serialized, out_path / "anchor_hits.csv")
    _write_html_report(out_path / "report.html", summary, read_results, lengths, base_qscores, read_qscores, anchor_positions, config.thresholds.heatmap_max_reads)
    return summary


def _write_html_report(path: Path, summary: dict[str, Any], read_results: list[ReadResult], lengths: list[int], base_qscores: list[float], read_qscores: list[float], anchor_positions: dict[str, list[float]], heatmap_max_reads: int) -> None:
    overall = summary["qc_verdicts"]["overall"]["status"]
    overview_cards = "".join(f"<div class='metric-card'><div class='metric-label'>{escape(label)}</div><div class='metric-value'>{escape(value)}</div></div>" for label, value in [("Total reads", f"{summary['total_reads']:,}"), ("Total bases", f"{summary['total_bases']:,}"), ("Median read length", f"{summary['median_read_length']:,} bp"), ("N50", f"{summary['n50']:,} bp"), ("Median read Qscore", f"{summary['median_read_qscore']:.2f}"), ("5' truncation", f"{summary['five_prime_truncation_ratio']:.1%}"), ("3' truncation", f"{summary['three_prime_truncation_ratio']:.1%}"), ("High-N reads", f"{summary['high_n_read_ratio']:.1%}")])
    verdict_cards = "".join(f"<div class='metric-card verdict-{escape(item['status'])}'><div class='metric-top'><div class='metric-label'>{escape(item['label'])}</div>{_status_pill(item['status'])}</div><div class='metric-value'>{item['value']:.1%}</div></div>" for key, item in summary["qc_verdicts"].items() if key != "overall")
    anchor_rows = "".join(f"<tr><td>{escape(name)}</td><td>{values['detection_ratio']:.1%}</td><td>{values['mean_best_mismatches']:.2f}</td><td>{values['multi_hit_ratio']:.1%}</td><td>{values['mean_hits_per_positive_read']:.2f}</td></tr>" for name, values in summary["anchor_quality_stats"].items())
    bucket_rows = "".join(f"<tr><td>{escape(bucket)}</td><td>{count:,}</td><td>{summary['qc_bucket_ratios'][bucket]:.1%}</td></tr>" for bucket, count in summary["qc_bucket_counts"].items())
    structure_rows = "".join(f"<tr><td>{escape(label)}</td><td>{escape(orientation)}</td><td>{count:,}</td></tr>" for label, counts in summary["structure_orientation_counts"].items() for orientation, count in counts.items())
    anchor_position_svg = _bar_chart([(name, sum(values) / len(values)) for name, values in anchor_positions.items() if values], "Anchor mean relative position", "Relative position", color="#2364aa", percent=True)
    html = f"""<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>scfastq-qc report - {escape(summary['sample_name'])}</title><style>:root{{--bg:#f7f5ef;--paper:rgba(255,255,255,.84);--ink:#12232f;--muted:#5b6b78;--line:rgba(18,35,47,.1);--teal:#2d8f85;--amber:#d97706;--red:#d1495b;--shadow:0 18px 48px rgba(21,40,54,.1)}}*{{box-sizing:border-box}}body{{margin:0;font-family:'IBM Plex Sans','Segoe UI',sans-serif;color:var(--ink);background:radial-gradient(circle at top left, rgba(35,100,170,.12), transparent 26%),radial-gradient(circle at top right, rgba(45,143,133,.14), transparent 22%),linear-gradient(180deg,#f2efe5 0%,#fbfaf7 100%)}}.page{{max-width:1280px;margin:0 auto;padding:30px 18px 56px}}.hero,section{{background:var(--paper);border:1px solid var(--line);border-radius:24px;box-shadow:var(--shadow)}}.hero{{padding:26px}}section{{margin-top:22px;padding:22px}}.eyebrow{{letter-spacing:.18em;text-transform:uppercase;color:var(--muted);font-size:.74rem;margin-bottom:10px}}h1{{font-family:Georgia,'Times New Roman',serif;font-size:clamp(2.1rem,4vw,3.4rem);margin:0 0 10px}}h2,h3{{margin:0}}p{{margin:0;color:var(--muted)}}.hero-grid{{display:grid;grid-template-columns:2.2fr 1fr;gap:18px;align-items:end}}.hero-chip{{border:1px solid var(--line);border-radius:18px;background:rgba(255,255,255,.72);padding:14px 16px;margin-top:12px}}.cards{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px}}.metric-card{{background:rgba(255,255,255,.74);border:1px solid rgba(18,35,47,.08);border-radius:18px;padding:16px;display:grid;gap:10px}}.metric-top{{display:flex;justify-content:space-between;gap:10px;align-items:center}}.metric-label{{color:var(--muted);font-size:.82rem;letter-spacing:.06em;text-transform:uppercase}}.metric-value{{font-size:clamp(1.35rem,2.8vw,1.95rem);font-weight:700}}.status-pill{{display:inline-flex;border-radius:999px;padding:7px 12px;font-size:.76rem;font-weight:700;letter-spacing:.08em}}.status-pass{{background:rgba(45,143,133,.12);color:var(--teal)}}.status-warn{{background:rgba(217,119,6,.12);color:var(--amber)}}.status-fail{{background:rgba(209,73,91,.12);color:var(--red)}}.verdict-pass{{border-color:rgba(45,143,133,.28)}}.verdict-warn{{border-color:rgba(217,119,6,.28)}}.verdict-fail{{border-color:rgba(209,73,91,.28)}}.section-head{{display:flex;justify-content:space-between;gap:16px;align-items:end;margin-bottom:16px}}.media-grid,.two-col{{display:grid;gap:16px}}.media-grid{{grid-template-columns:repeat(2,minmax(0,1fr))}}.two-col{{grid-template-columns:1.2fr 1fr}}.panel{{background:rgba(255,255,255,.74);border:1px solid rgba(18,35,47,.08);border-radius:18px;padding:14px;overflow:auto}}table{{width:100%;border-collapse:collapse;font-size:.95rem}}th,td{{padding:12px 14px;border-bottom:1px solid rgba(18,35,47,.08);text-align:left}}th{{background:rgba(18,35,47,.04);color:var(--muted);font-size:.78rem;text-transform:uppercase;letter-spacing:.06em}}tr:last-child td{{border-bottom:none}}.svg-wrap svg{{width:100%;height:auto;display:block}}code{{word-break:break-all}}@media (max-width:980px){{.hero-grid,.cards,.media-grid,.two-col{{grid-template-columns:1fr}}}}</style></head><body><div class='page'><header class='hero'><div class='eyebrow'>Structure-aware fastq QC</div><div class='hero-grid'><div><h1>{escape(summary['sample_name'])}</h1><p>Single-file HTML report for long-read single-cell basic QC: read yield, quality, anchor integrity, structure compliance, truncation, and orientation.</p><p style='margin-top:10px;'>Input FASTQ: <code>{escape(summary['fastq_path'])}</code></p><p style='margin-top:10px;'>Overall assessment: {_status_pill(overall)}</p></div><div><div class='hero-chip'><strong>Rust accelerator</strong><div>{escape(summary['rust_accelerator'].get('mode','unknown'))}</div><p style='margin-top:6px;'>{escape(summary['rust_accelerator'].get('detail',''))}</p></div><div class='hero-chip'><strong>Anchor order</strong><div>{summary['correct_anchor_order_ratio']:.1%}</div><p style='margin-top:6px;'>Reads whose detected anchors match the configured order.</p></div></div></div></header><section><div class='section-head'><div><h2>QC Verdicts</h2><p>Threshold-based pass, warn, fail calls for the core basic-QC indicators.</p></div></div><div class='cards'>{verdict_cards}</div></section><section><div class='section-head'><div><h2>Overview</h2><p>Yield, read quality, truncation, and contamination proxies at a glance.</p></div></div><div class='cards'>{overview_cards}</div></section><section><div class='section-head'><div><h2>Yield And Quality</h2><p>Baseline library quality by read length and Q-score distributions.</p></div></div><div class='media-grid'><div class='panel svg-wrap'>{_histogram([float(v) for v in lengths], "Read length distribution", "Read length (bp)")}</div><div class='panel svg-wrap'>{_histogram(base_qscores, "Overall base quality distribution", "Phred quality score", color="#2d8f85")}</div><div class='panel svg-wrap'>{_line_density(read_qscores, "Per-read Qscore density", "Read Qscore")}</div><div class='panel svg-wrap'>{_scatter(lengths, read_qscores, "Read length vs read Qscore")}</div></div></section><section><div class='section-head'><div><h2>Anchors</h2><p>Detection rate alone is not enough; mismatch burden and hit multiplicity expose weak motifs and internal artifacts.</p></div></div><div class='media-grid'><div class='panel svg-wrap'>{_bar_chart(list(summary['anchor_detection_ratio'].items()), "Anchor detection ratio", "Fraction of reads", color="#2364aa", percent=True)}</div><div class='panel svg-wrap'>{_bar_chart([(name, values['mean_best_mismatches']) for name, values in summary['anchor_quality_stats'].items()], "Anchor mismatch burden", "Mean best-hit mismatches", color="#d97706")}</div></div><div class='panel svg-wrap' style='margin-top:16px;'>{anchor_position_svg}</div><div class='panel' style='margin-top:16px;'><table><thead><tr><th>Anchor</th><th>Detection ratio</th><th>Mean mismatches</th><th>Multi-hit ratio</th><th>Mean hits / positive read</th></tr></thead><tbody>{anchor_rows}</tbody></table></div></section><section><div class='section-head'><div><h2>Structure And Failure Modes</h2><p>Basic-QC only: structure integrity, truncation, orientation, and read-level hygiene.</p></div></div><div class='media-grid'><div class='panel svg-wrap'>{_bar_chart([(key, float(value)) for key, value in summary['structure_counts'].items()], "Read structure classification", "Read count", color="#2d8f85")}</div><div class='panel svg-wrap'>{_bar_chart([(key, float(value)) for key, value in summary['qc_bucket_counts'].items()], "QC failure buckets", "Read count", color="#d1495b")}</div></div><div class='panel svg-wrap' style='margin-top:16px;'>{_heatmap(read_results, "Read structure heatmap", heatmap_max_reads)}</div><div class='two-col' style='margin-top:16px;'><div class='panel'><h3>QC Buckets</h3><table><thead><tr><th>Bucket</th><th>Reads</th><th>Ratio</th></tr></thead><tbody>{bucket_rows}</tbody></table></div><div class='panel'><h3>Structure Orientation</h3><table><thead><tr><th>Class</th><th>Orientation</th><th>Reads</th></tr></thead><tbody>{structure_rows}</tbody></table></div></div></section></div></body></html>"""
    _write_text(path, html)
