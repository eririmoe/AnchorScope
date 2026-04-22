from __future__ import annotations

import json
import logging
import csv
from collections import Counter, defaultdict
from html import escape
from math import ceil, floor
from pathlib import Path
from typing import Any

from .anchors import AnchorHit, find_anchor_hits, get_rust_status, prepare_anchors, reverse_complement
from .classify import classify_best_orientation
from .config import AppConfig, ThresholdConfig
from .export import export_structure_classification_to_csv, export_summary_to_csv
from .fastq import FastqRead, read_fastq

logger = logging.getLogger(__name__)

READ_QSCORE_SCALE = 100
SCATTER_LENGTH_BIN_BP = 50
SCATTER_QSCORE_BIN = 0.5
HEATMAP_BINS = 100
COLOR_PRIMARY = "#2364aa"
COLOR_QUALITY = "#2d8f85"
COLOR_WARNING = "#d97706"
COLOR_FAILURE = "#d1495b"
COLOR_MUTED = "#8d99ae"

def compute_n50_counts(length_counts: Counter[int]) -> int:
    if not length_counts:
        return 0
    total = sum(length * count for length, count in length_counts.items())
    running = 0
    for length in sorted(length_counts, reverse=True):
        running += length * length_counts[length]
        if running >= total / 2:
            return length
    return 0


def _weighted_median(value_counts: Counter[int], scale: float = 1.0) -> float:
    if not value_counts:
        return 0.0
    total = sum(value_counts.values())
    midpoint_low = (total - 1) // 2
    midpoint_high = total // 2
    cumulative = 0
    low_value: int | None = None
    high_value: int | None = None
    for value in sorted(value_counts):
        cumulative += value_counts[value]
        if low_value is None and cumulative > midpoint_low:
            low_value = value
        if cumulative > midpoint_high:
            high_value = value
            break
    assert low_value is not None and high_value is not None
    return ((low_value + high_value) / 2) / scale


def _mean_from_counts(value_counts: Counter[int], scale: float = 1.0) -> float:
    total = sum(value_counts.values())
    if total == 0:
        return 0.0
    return sum(value * count for value, count in value_counts.items()) / (total * scale)


def _value_counts_to_float_map(value_counts: Counter[int], scale: float = 1.0) -> dict[float, int]:
    return {value / scale: count for value, count in value_counts.items()}


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


def _format_axis_value(value: float, percent: bool = False) -> str:
    if percent:
        return f"{value:.0%}"
    if abs(value) >= 100 or float(value).is_integer():
        return _format_num(value, 0)
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _linear_ticks(lower: float, upper: float, steps: int = 5) -> list[float]:
    if upper <= lower:
        return [lower]
    return [lower + ((upper - lower) * idx / steps) for idx in range(steps + 1)]


def _histogram(values: list[float] | dict[float, int], title: str, xlabel: str, color: str = COLOR_QUALITY, bins: int = 20) -> str:
    if isinstance(values, dict):
        value_counts = {float(key): count for key, count in values.items() if count > 0}
    else:
        value_counts = dict(Counter(float(value) for value in values))
    if not value_counts:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>")
    width, height = 780, 320
    left, top, right, bottom = 60, 48, 20, 45
    plot_w, plot_h = width - left - right, height - top - bottom
    lower, upper = min(value_counts), max(value_counts)
    if lower == upper:
        lower -= 0.5
        upper += 0.5
    step = (upper - lower) / bins
    counts = [0] * bins
    for value, count in value_counts.items():
        counts[min(bins - 1, int((value - lower) / step))] += count
    max_count = max(counts) or 1
    parts: list[str] = []
    hover_targets: list[str] = []
    for idx in range(6):
        y = top + plot_h - (idx / 5) * plot_h
        parts.append(f"<line class='grid' x1='{left}' y1='{y:.2f}' x2='{left+plot_w}' y2='{y:.2f}'/>")
        parts.append(f"<text x='{left-8}' y='{y+4:.2f}' text-anchor='end' font-size='11'>{_format_num((idx/5)*max_count, 0)}</text>")
    x_ticks = _linear_ticks(lower, upper, 5)
    for tick in x_ticks:
        x = left + _scale(tick, lower, upper, plot_w)
        parts.append(f"<line class='grid' x1='{x:.2f}' y1='{top}' x2='{x:.2f}' y2='{top+plot_h}' stroke-opacity='0.45'/>")
        parts.append(f"<text x='{x:.2f}' y='{top+plot_h+20:.2f}' text-anchor='middle' font-size='11'>{escape(_format_axis_value(tick))}</text>")
    bar_w = plot_w / bins
    for idx, count in enumerate(counts):
        bar_h = (count / max_count) * plot_h
        x = left + idx * bar_w
        y = top + plot_h - bar_h
        width_value = max(bar_w - 1, 1)
        if bar_h > 0:
            parts.append(f"<rect class='bar' x='{x:.2f}' y='{y:.2f}' width='{width_value:.2f}' height='{bar_h:.2f}' fill='{color}'/>")
            bin_start = lower + idx * step
            bin_end = lower + (idx + 1) * step
            hover_targets.append(
                f"<rect class='hist-bar-target' data-range='{bin_start:.2f} to {bin_end:.2f}' data-count='{count}' "
                f"x='{x:.2f}' y='{top:.2f}' width='{width_value:.2f}' height='{plot_h:.2f}' fill='transparent'/>"
            )
    parts.append(f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>")
    parts.append(f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-12}' text-anchor='middle' font-size='12'>{escape(xlabel)}</text>")
    parts.append(f"<text x='18' y='{height/2:.0f}' transform='rotate(-90 18,{height/2:.0f})' text-anchor='middle' font-size='12'>Count</text>")
    chart_id = "hist-" + "".join(ch if ch.isalnum() else "-" for ch in title.lower())
    svg = _svg_wrapper(title, ''.join(parts) + ''.join(hover_targets), width, height).replace(
        "<svg ",
        f"<svg id='{escape(chart_id)}' ",
        1,
    )
    return (
        f"<div class='interactive-plot'><div class='plot-note'>Hover over a bar to inspect its bin range and read count.</div>{svg}"
        f"<div class='plot-tip' id='{escape(chart_id)}-tip' hidden></div></div>"
        f"<script>(function(){{"
        f"const root=document.getElementById('{escape(chart_id)}');"
        f"if(!root) return;"
        f"const tip=document.getElementById('{escape(chart_id)}-tip');"
        f"for(const node of root.querySelectorAll('.hist-bar-target')){{"
        f"node.addEventListener('mouseenter', function(){{"
        f"tip.hidden=false;"
        f"tip.textContent='Range '+this.dataset.range+' | count '+this.dataset.count;"
        f"}});"
        f"node.addEventListener('mousemove', function(event){{"
        f"const box=root.getBoundingClientRect();"
        f"tip.style.left=(event.clientX-box.left+14)+'px';"
        f"tip.style.top=(event.clientY-box.top-10)+'px';"
        f"}});"
        f"node.addEventListener('mouseleave', function(){{ tip.hidden=true; }});"
        f"}}"
        f"}})();</script>"
    )


def _line_density(values: list[float] | dict[float, int], title: str, xlabel: str, color: str = COLOR_QUALITY, bin_width: float = 0.2, chart_id: str = "density") -> str:
    if isinstance(values, dict):
        value_counts = {float(key): count for key, count in values.items() if count > 0}
    else:
        value_counts = dict(Counter(float(value) for value in values))
    if not value_counts:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>")
    width, height = 780, 320
    left, top, right, bottom = 60, 48, 20, 45
    plot_w, plot_h = width - left - right, height - top - bottom
    lower, upper = min(value_counts), max(value_counts)
    if lower == upper:
        lower -= 0.5
        upper += 0.5
    lower = floor(lower) - 1
    upper = ceil(upper) + 1
    bin_count = max(2, int(ceil((upper - lower) / bin_width)))
    counts = [0.0] * bin_count
    total_points = sum(value_counts.values())
    for value, count in value_counts.items():
        counts[min(bin_count - 1, max(0, int((value - lower) / bin_width)))] += count
    densities = [count / (total_points * bin_width) for count in counts]
    y_max = max(densities) * 1.1 if densities else 1.0
    points = []
    hover_targets: list[str] = []
    x_ticks = _linear_ticks(lower, upper, 5)
    y_ticks = _linear_ticks(0, y_max, 5)
    parts: list[str] = []
    for tick in x_ticks:
        x = left + _scale(tick, lower, upper, plot_w)
        parts.append(f"<line class='grid' x1='{x:.2f}' y1='{top}' x2='{x:.2f}' y2='{top+plot_h}' stroke-opacity='0.45'/>")
        parts.append(f"<text x='{x:.2f}' y='{top+plot_h+20:.2f}' text-anchor='middle' font-size='11'>{escape(_format_axis_value(tick))}</text>")
    for tick in y_ticks:
        y = top + plot_h - _scale(tick, 0, y_max, plot_h)
        parts.append(f"<line class='grid' x1='{left}' y1='{y:.2f}' x2='{left+plot_w}' y2='{y:.2f}'/>")
        parts.append(f"<text x='{left-8}' y='{y+4:.2f}' text-anchor='end' font-size='11'>{escape(_format_axis_value(tick))}</text>")
    for idx, density in enumerate(densities):
        x_value = lower + (idx + 0.5) * bin_width
        x = left + _scale(x_value, lower, upper, plot_w)
        y = top + plot_h - _scale(density, 0, y_max, plot_h)
        points.append(f"{x:.2f},{y:.2f}")
        hover_targets.append(
            f"<circle class='density-point' data-q='{x_value:.2f}' data-density='{density:.4f}' cx='{x:.2f}' cy='{y:.2f}' r='5' fill='{color}' fill-opacity='0.01' stroke='transparent'/>"
        )
        parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='2.3' fill='{color}'/>")
    body = (
        f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>"
        f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>"
        f"<polygon points='{left},{top+plot_h} {' '.join(points)} {left+plot_w},{top+plot_h}' fill='{color}' fill-opacity='0.16'/>"
        f"<polyline fill='none' stroke='{color}' stroke-width='2.5' points='{' '.join(points)}'/>"
        f"{''.join(parts)}"
        f"{''.join(hover_targets)}"
        f"<text x='{width/2:.0f}' y='{height-12}' text-anchor='middle' font-size='12'>{escape(xlabel)}</text>"
        f"<text x='18' y='{height/2:.0f}' transform='rotate(-90 18,{height/2:.0f})' text-anchor='middle' font-size='12'>Probability density</text>"
    )
    svg = _svg_wrapper(title, body, width, height).replace(
        "<svg ",
        f"<svg id='{escape(chart_id)}' ",
        1,
    )
    return (
        f"<div class='interactive-plot'><div class='plot-note'>Hover over the curve to inspect local density.</div>{svg}"
        f"<div class='plot-tip' id='{escape(chart_id)}-tip' hidden></div></div>"
        f"<script>(function(){{"
        f"const root=document.getElementById('{escape(chart_id)}');"
        f"if(!root) return;"
        f"const tip=document.getElementById('{escape(chart_id)}-tip');"
        f"for(const node of root.querySelectorAll('.density-point')){{"
        f"node.addEventListener('mouseenter', function(){{"
        f"tip.hidden=false;"
        f"tip.textContent='Qscore '+this.dataset.q+' | density '+this.dataset.density;"
        f"}});"
        f"node.addEventListener('mousemove', function(event){{"
        f"const box=root.getBoundingClientRect();"
        f"tip.style.left=(event.clientX-box.left+14)+'px';"
        f"tip.style.top=(event.clientY-box.top-10)+'px';"
        f"}});"
        f"node.addEventListener('mouseleave', function(){{ tip.hidden=true; }});"
        f"}}"
        f"}})();</script>"
    )


def _bar_chart(items: list[tuple[str, float]], title: str, ylabel: str, color: str = "#d1495b", percent: bool = False) -> str:
    if not items:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>", 820, 360)
    width, height = 820, 360
    left, top, right, bottom = 70, 52, 24, 96
    plot_w, plot_h = width - left - right, height - top - bottom
    y_max = 1.0 if percent else max(value for _, value in items) or 1.0
    parts = [f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>", f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>"]
    for tick in _linear_ticks(0, y_max, 5):
        y = top + plot_h - _scale(tick, 0, y_max, plot_h)
        parts.append(f"<line class='grid' x1='{left}' y1='{y:.2f}' x2='{left+plot_w}' y2='{y:.2f}'/>")
        parts.append(f"<text x='{left-8}' y='{y+4:.2f}' text-anchor='end' font-size='11'>{escape(_format_axis_value(tick, percent=percent))}</text>")
    slot_w = plot_w / len(items)
    bar_w = min(slot_w * 0.58, 120)
    for idx, (label, value) in enumerate(items):
        bar_h = (value / y_max) * plot_h if y_max else 0
        x = left + idx * slot_w + (slot_w - bar_w) / 2
        y = top + plot_h - bar_h
        disp = _format_axis_value(value, percent=percent)
        tx = left + idx * slot_w + slot_w / 2
        if bar_h > 0:
            parts.append(f"<rect class='bar' x='{x:.2f}' y='{y:.2f}' width='{bar_w:.2f}' height='{bar_h:.2f}' fill='{color}'/>")
            parts.append(f"<text x='{tx:.2f}' y='{max(y-8, top+10):.2f}' text-anchor='middle' font-size='11'>{escape(disp)}</text>")
        else:
            parts.append(f"<text x='{tx:.2f}' y='{top+plot_h-6:.2f}' text-anchor='middle' font-size='11'>{escape(disp)}</text>")
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
    for tick in _linear_ticks(min_x, max_x, 5):
        x = left + _scale(tick, min_x, max_x, plot_w)
        parts.append(f"<line class='grid' x1='{x:.2f}' y1='{top}' x2='{x:.2f}' y2='{top+plot_h}' stroke-opacity='0.45'/>")
        parts.append(f"<text x='{x:.2f}' y='{top+plot_h+20:.2f}' text-anchor='middle' font-size='11'>{escape(_format_axis_value(tick))}</text>")
    for tick in _linear_ticks(min_y, max_y, 5):
        y = top + plot_h - _scale(tick, min_y, max_y, plot_h)
        parts.append(f"<line class='grid' x1='{left}' y1='{y:.2f}' x2='{left+plot_w}' y2='{y:.2f}'/>")
        parts.append(f"<text x='{left-8}' y='{y+4:.2f}' text-anchor='end' font-size='11'>{escape(_format_axis_value(tick))}</text>")
    for x_value, y_value in zip(lengths, read_qscores):
        x = left + _scale(x_value, min_x, max_x, plot_w)
        y = top + plot_h - _scale(y_value, min_y, max_y, plot_h)
        parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='3' fill='#2d8f85' fill-opacity='0.58'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-12}' text-anchor='middle' font-size='12'>Read length (bp)</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle' font-size='12'>Read Qscore</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _scatter_binned(scatter_counts: Counter[tuple[int, int]], title: str) -> str:
    width, height = 780, 320
    left, top, right, bottom = 60, 48, 20, 45
    plot_w, plot_h = width - left - right, height - top - bottom
    if not scatter_counts:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>")
    x_values = [point[0] for point in scatter_counts]
    y_values = [point[1] for point in scatter_counts]
    min_x, max_x = min(x_values), max(x_values)
    min_y, max_y = min(y_values), max(y_values)
    if min_x == max_x:
        max_x += SCATTER_LENGTH_BIN_BP
    if min_y == max_y:
        max_y += SCATTER_QSCORE_BIN
    max_count = max(scatter_counts.values()) or 1
    hover_targets: list[str] = []
    parts = [f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>", f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>"]
    for tick in _linear_ticks(min_x, max_x, 5):
        x = left + _scale(tick, min_x, max_x, plot_w)
        parts.append(f"<line class='grid' x1='{x:.2f}' y1='{top}' x2='{x:.2f}' y2='{top+plot_h}' stroke-opacity='0.45'/>")
        parts.append(f"<text x='{x:.2f}' y='{top+plot_h+20:.2f}' text-anchor='middle' font-size='11'>{escape(_format_axis_value(tick))}</text>")
    for tick in _linear_ticks(min_y, max_y, 5):
        y = top + plot_h - _scale(tick, min_y, max_y, plot_h)
        parts.append(f"<line class='grid' x1='{left}' y1='{y:.2f}' x2='{left+plot_w}' y2='{y:.2f}'/>")
        parts.append(f"<text x='{left-8}' y='{y+4:.2f}' text-anchor='end' font-size='11'>{escape(_format_axis_value(tick))}</text>")
    for (x_value, y_value), count in scatter_counts.items():
        x = left + _scale(x_value, min_x, max_x, plot_w)
        y = top + plot_h - _scale(y_value, min_y, max_y, plot_h)
        radius = 2.2 + 1.8 * ((count / max_count) ** 0.35)
        opacity = 0.14 + 0.64 * ((count / max_count) ** 0.6)
        parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='{radius:.2f}' fill='{COLOR_QUALITY}' fill-opacity='{opacity:.2f}'/>")
        length_start = x_value - (SCATTER_LENGTH_BIN_BP / 2)
        length_end = x_value + (SCATTER_LENGTH_BIN_BP / 2)
        qscore_start = y_value - (SCATTER_QSCORE_BIN / 2)
        qscore_end = y_value + (SCATTER_QSCORE_BIN / 2)
        hover_targets.append(
            f"<circle class='scatter-point-target' "
            f"data-length='{length_start:.0f}-{length_end:.0f} bp' "
            f"data-qscore='{qscore_start:.2f}-{qscore_end:.2f}' "
            f"data-count='{count}' "
            f"cx='{x:.2f}' cy='{y:.2f}' r='{max(radius + 4, 7):.2f}' fill='transparent'/>"
        )
    parts.append(''.join(hover_targets))
    parts.append(f"<text x='{width/2:.0f}' y='{height-12}' text-anchor='middle' font-size='12'>Read length (bp)</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle' font-size='12'>Read Qscore</text>")
    chart_id = "scatter-" + "".join(ch if ch.isalnum() else "-" for ch in title.lower())
    svg = _svg_wrapper(title, ''.join(parts), width, height).replace(
        "<svg ",
        f"<svg id='{escape(chart_id)}' ",
        1,
    )
    return (
        f"<div class='interactive-plot'><div class='plot-note'>Hover over a density point to inspect its length bin, Qscore bin, and read count.</div>{svg}"
        f"<div class='plot-tip' id='{escape(chart_id)}-tip' hidden></div></div>"
        f"<script>(function(){{"
        f"const root=document.getElementById('{escape(chart_id)}');"
        f"if(!root) return;"
        f"const tip=document.getElementById('{escape(chart_id)}-tip');"
        f"for(const node of root.querySelectorAll('.scatter-point-target')){{"
        f"node.addEventListener('mouseenter', function(){{"
        f"tip.hidden=false;"
        f"tip.textContent='Length '+this.dataset.length+' | Qscore '+this.dataset.qscore+' | reads '+this.dataset.count;"
        f"}});"
        f"node.addEventListener('mousemove', function(event){{"
        f"const box=root.getBoundingClientRect();"
        f"tip.style.left=(event.clientX-box.left+14)+'px';"
        f"tip.style.top=(event.clientY-box.top-10)+'px';"
        f"}});"
        f"node.addEventListener('mouseleave', function(){{ tip.hidden=true; }});"
        f"}}"
        f"}})();</script>"
    )


def _heatmap(heatmap_groups: dict[tuple[str, str], dict[str, Any]], title: str) -> str:
    bins = HEATMAP_BINS
    palette = {"pass": COLOR_QUALITY, "too_short": COLOR_FAILURE, "low_read_q": COLOR_WARNING, "no_anchor": COLOR_MUTED}
    if not heatmap_groups:
        return _svg_wrapper(title, "<text x='40' y='90'>No data</text>", 980, 360)
    rows = sorted(
        heatmap_groups.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            -item[1]["count"],
        ),
    )
    max_label_len = max(len(f"{qc_bucket} | {structure_label} (n={group['count']:,})") for (qc_bucket, structure_label), group in rows)
    left = min(460, max(280, 36 + max_label_len * 6))
    width = max(980, left + 760)
    top, right, bottom = 52, 26, 64
    row_h = 28
    plot_h = max(180, len(rows) * row_h)
    height = top + plot_h + bottom
    plot_w = width - left - right
    col_w = plot_w / bins

    def _blend_with_white(hex_color: str, ratio: float) -> str:
        ratio = max(0.0, min(1.0, ratio))
        base = hex_color.lstrip("#")
        r = int(base[0:2], 16)
        g = int(base[2:4], 16)
        b = int(base[4:6], 16)
        out_r = round(255 - (255 - r) * ratio)
        out_g = round(255 - (255 - g) * ratio)
        out_b = round(255 - (255 - b) * ratio)
        return f"#{out_r:02x}{out_g:02x}{out_b:02x}"

    parts = [
        f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>",
        f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}'/>",
    ]
    for tick in range(6):
        x = left + (tick / 5) * plot_w
        label = f"{int((tick / 5) * 100)}%"
        parts.append(f"<line class='grid' x1='{x:.2f}' y1='{top}' x2='{x:.2f}' y2='{top+plot_h}' stroke-opacity='0.45'/>")
        parts.append(f"<text x='{x:.2f}' y='{top+plot_h+22:.2f}' text-anchor='middle' font-size='11'>{label}</text>")

    for row_idx, ((qc_bucket, structure_label), group) in enumerate(rows):
        y = top + row_idx * row_h
        count = group["count"]
        label = f"{qc_bucket} | {structure_label} (n={count:,})"
        bucket_color = palette.get(qc_bucket, "#e6eef3")
        parts.append(f"<rect x='{left-18:.2f}' y='{y+6:.2f}' width='10' height='{max(row_h-12, 8):.2f}' rx='3' ry='3' fill='{bucket_color}'/>")
        parts.append(f"<text x='{left-24:.2f}' y='{y + row_h/2 + 4:.2f}' text-anchor='end' font-size='11'>{escape(label)}</text>")
        parts.append(f"<rect x='{left}' y='{y:.2f}' width='{plot_w}' height='{row_h:.2f}' fill='#f6f8fb'/>")
        for idx, covered_count in enumerate(group["anchor_cover"]):
            ratio = covered_count / count if count else 0.0
            intensity = ratio ** 0.65
            fill = _blend_with_white(COLOR_PRIMARY, intensity)
            x = left + idx * col_w
            parts.append(f"<rect x='{x:.2f}' y='{y:.2f}' width='{max(col_w, 1):.2f}' height='{row_h:.2f}' fill='{fill}'/>")
        parts.append(f"<line x1='{left}' y1='{y+row_h:.2f}' x2='{left+plot_w}' y2='{y+row_h:.2f}' stroke='rgba(18,35,47,.08)' stroke-width='1'/>")

    parts.append(f"<text x='{width/2:.0f}' y='{height-16:.2f}' text-anchor='middle' font-size='12'>Normalized read position</text>")
    svg = _svg_wrapper(title, "".join(parts), width, height)
    total_reads = sum(group["count"] for group in heatmap_groups.values())
    note = (
        f"Aggregated across all {total_reads:,} reads into {len(rows):,} QC/structure groups. "
        f"Darker blue means a larger fraction of reads in that group contain an anchor hit at that relative position."
    )
    legend_items = [
        ("QC bucket: pass", "#2d8f85"),
        ("QC bucket: too_short", "#d1495b"),
        ("QC bucket: low_read_q", "#d97706"),
        ("QC bucket: no_anchor", "#8d99ae"),
        ("QC bucket: other", "#e6eef3"),
        ("Anchor occupancy: low", _blend_with_white(COLOR_PRIMARY, 0.2)),
        ("Anchor occupancy: high", _blend_with_white(COLOR_PRIMARY, 1.0)),
    ]
    legend_html = "".join(
        f"<span class='heatmap-legend-item'>"
        f"<span class='heatmap-swatch' style='background:{color}'></span>"
        f"<span>{escape(label)}</span>"
        f"</span>"
        for label, color in legend_items
    )
    return (
        f"<div class='heatmap-frame'>"
        f"<div class='plot-note'>{escape(note)}</div>"
        f"<div class='heatmap-legend'>{legend_html}</div>"
        f"<div class='chart-note'>Each row aggregates reads sharing the same `qc_bucket` and `structure_label`. Forward and reverse-complement reads are merged here to keep the map readable. Cell color encodes the fraction of reads in that group whose anchor hits cover that normalized position bin.</div>"
        f"{svg}"
        f"</div>"
    )


def _anchor_position_tracks(anchor_positions: dict[str, dict[str, float]]) -> str:
    items = [
        (name, stats["sum"] / stats["count"])
        for name, stats in anchor_positions.items()
        if stats.get("count", 0) > 0
    ]
    if not items:
        return _svg_wrapper("Anchor mean relative position", "<text x='40' y='90'>No data</text>", 900, 240)
    width = 900
    height = max(220, 130 + len(items) * 46)
    left, top, right, bottom = 160, 54, 36, 48
    plot_w = width - left - right
    row_gap = 40
    parts = []
    for tick in range(6):
        x = left + (tick / 5) * plot_w
        parts.append(f"<line class='grid' x1='{x:.2f}' y1='{top-6}' x2='{x:.2f}' y2='{height-bottom+6}' stroke-opacity='0.4'/>")
        parts.append(f"<text x='{x:.2f}' y='{height-18:.2f}' text-anchor='middle' font-size='11'>{int((tick / 5) * 100)}%</text>")
    for idx, (name, value) in enumerate(items):
        y = top + idx * row_gap
        x = left + value * plot_w
        parts.append(f"<text x='{left-12:.2f}' y='{y+4:.2f}' text-anchor='end' font-size='12'>{escape(name)}</text>")
        parts.append(f"<line class='axis' x1='{left}' y1='{y:.2f}' x2='{left+plot_w}' y2='{y:.2f}'/>")
        parts.append(f"<polygon points='{x:.2f},{y-11:.2f} {x-7:.2f},{y+1:.2f} {x+7:.2f},{y+1:.2f}' fill='{COLOR_PRIMARY}'/>")
        parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='4.5' fill='{COLOR_PRIMARY}'/>")
        parts.append(f"<text x='{min(x + 10, left + plot_w - 6):.2f}' y='{y-10:.2f}' font-size='11'>{value:.1%}</text>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-4:.2f}' text-anchor='middle' font-size='12'>Relative position from 5' to 3' end of the read</text>")
    return _svg_wrapper("Anchor mean relative position", "".join(parts), width, height)


def _status_pill(status: str) -> str:
    return f"<span class='status-pill status-{escape(status)}'>{escape(status.upper())}</span>"


def _tooltip_label(label: str, body: str) -> str:
    return (
        f"<span class='tooltip metric-tip metric-tip-inline' tabindex='0' role='note' aria-label='{escape(label)} help'>"
        f"<span class='tooltip-trigger'>{escape(label)}</span>"
        f"<span class='tooltip-content'><strong>{escape(label)}</strong><div>{escape(body)}</div></span>"
        f"</span>"
    )


def run_qc(fastq_path: str, config: AppConfig, outdir: str, export_csv: bool = False) -> dict[str, Any]:
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)
    prepared_anchors = prepare_anchors(config.anchors)

    length_counts: Counter[int] = Counter()
    base_qscore_counts: Counter[int] = Counter()
    read_qscore_counts: Counter[int] = Counter()
    scatter_counts: Counter[tuple[int, int]] = Counter()
    anchor_counts: Counter[str] = Counter()
    anchor_best_mismatches: dict[str, list[int]] = defaultdict(list)
    anchor_hit_counts: dict[str, list[int]] = defaultdict(list)
    anchor_positions: dict[str, dict[str, float]] = defaultdict(lambda: {"sum": 0.0, "count": 0})
    heatmap_groups: dict[tuple[str, str], dict[str, Any]] = {}
    structure_counts: Counter[str] = Counter()
    structure_orientation_counts: Counter[tuple[str, str]] = Counter()
    qc_bucket_counts: Counter[str] = Counter()
    total_reads = total_bases = long_high_quality = order_valid_count = reversed_read_count = 0
    five_prime_truncation_count = three_prime_truncation_count = high_n_read_count = invalid_base_read_count = 0

    read_details_handle = anchor_hits_handle = None
    read_details_writer = anchor_hits_writer = None
    if export_csv:
        read_details_handle = open(out_path / "read_details.csv", "w", newline="", encoding="utf-8")
        read_details_writer = csv.writer(read_details_handle)
        read_details_writer.writerow(
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
        anchor_hits_handle = open(out_path / "anchor_hits.csv", "w", newline="", encoding="utf-8")
        anchor_hits_writer = csv.writer(anchor_hits_handle)
        anchor_hits_writer.writerow(
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

    fastq_export_handles: dict[str, Any] = {}

    try:
        for read in read_fastq(fastq_path, qscore_method=config.qscore_method):
            total_reads += 1
            total_bases += read.length
            length_counts[read.length] += 1
            base_qscore_counts.update(read.phred_scores)
            read_qscore_key = int(round(read.read_qscore * READ_QSCORE_SCALE))
            read_qscore_counts[read_qscore_key] += 1
            scatter_length = int(round(read.length / SCATTER_LENGTH_BIN_BP) * SCATTER_LENGTH_BIN_BP)
            scatter_qscore = round(round(read.read_qscore / SCATTER_QSCORE_BIN) * SCATTER_QSCORE_BIN, 2)
            scatter_counts[(scatter_length, scatter_qscore)] += 1

            forward_hits = find_anchor_hits(read.sequence, prepared_anchors)
            reverse_hits = find_anchor_hits(reverse_complement(read.sequence), prepared_anchors)
            structure, hits = classify_best_orientation(forward_hits, reverse_hits, config.structure.expected_order)
            best_hits = _best_hits_by_anchor(hits)
            five_prime_offset, three_prime_offset = _compute_terminal_offsets(best_hits, config.structure.expected_order, read.length)
            qc_flags = _derive_qc_flags(read, structure.label, five_prime_offset, three_prime_offset, config, hits)
            qc_bucket = _primary_qc_bucket(qc_flags)
            orientation = "reversed" if structure.is_reversed else "forward"

            structure_counts[structure.label] += 1
            structure_orientation_counts[(structure.label, orientation)] += 1
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

            heatmap_key = (qc_bucket, structure.label)
            heatmap_group = heatmap_groups.setdefault(
                heatmap_key,
                {"count": 0, "anchor_cover": [0] * HEATMAP_BINS},
            )
            heatmap_group["count"] += 1
            if read.length > 0 and hits:
                covered = [False] * HEATMAP_BINS
                for hit in hits:
                    start = min(HEATMAP_BINS - 1, max(0, int((hit.start / read.length) * HEATMAP_BINS)))
                    end = min(HEATMAP_BINS, max(start + 1, int((hit.end / read.length) * HEATMAP_BINS)))
                    for idx in range(start, end):
                        covered[idx] = True
                for idx, has_cover in enumerate(covered):
                    if has_cover:
                        heatmap_group["anchor_cover"][idx] += 1

            for name, hit in best_hits.items():
                anchor_counts[name] += 1
                anchor_best_mismatches[name].append(hit.mismatches)
                anchor_positions[name]["sum"] += hit.start / read.length if read.length > 0 else 0.0
                anchor_positions[name]["count"] += 1
                anchor_hit_counts[name].append(sum(1 for item in hits if item.anchor_name == name))

            if read_details_writer is not None:
                read_details_writer.writerow(
                    [
                        read.name,
                        read.length,
                        f"{read.read_qscore:.2f}",
                        structure.label,
                        qc_bucket,
                        ";".join(qc_flags),
                        structure.is_reversed,
                        structure.order_valid,
                        f"{five_prime_offset:.4f}" if five_prime_offset is not None else "",
                        f"{three_prime_offset:.4f}" if three_prime_offset is not None else "",
                        f"{read.n_fraction:.4f}",
                        f"{read.invalid_base_fraction:.4f}",
                        len(hits),
                    ]
                )
            if anchor_hits_writer is not None:
                for hit in hits:
                    anchor_hits_writer.writerow(
                        [
                            read.name,
                            read.length,
                            f"{read.read_qscore:.2f}",
                            structure.label,
                            qc_bucket,
                            structure.is_reversed,
                            hit.anchor_name,
                            hit.start,
                            hit.end,
                            hit.mismatches,
                            hit.matched_sequence,
                        ]
                    )

            if config.export_non_full_structure_fastq and structure.label != "full_structure":
                handle = fastq_export_handles.get(structure.label)
                if handle is None:
                    handle = open(out_path / f"{structure.label}.fastq", "w", encoding="utf-8")
                    fastq_export_handles[structure.label] = handle
                handle.write(f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n")
    finally:
        if read_details_handle is not None:
            read_details_handle.close()
        if anchor_hits_handle is not None:
            anchor_hits_handle.close()
        for handle in fastq_export_handles.values():
            handle.close()

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
        "config_summary": {
            "anchors": [anchor.name for anchor in config.anchors],
            "expected_order": list(config.structure.expected_order),
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
        export_summary_to_csv(summary, out_path / "summary.csv")
        export_structure_classification_to_csv(summary, out_path / "structure_classification.csv")
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
    return summary


def _write_html_report(
    path: Path,
    summary: dict[str, Any],
    length_counts: Counter[int],
    base_qscore_counts: Counter[int],
    read_qscore_counts: Counter[int],
    scatter_counts: Counter[tuple[int, int]],
    anchor_positions: dict[str, dict[str, float]],
    heatmap_groups: dict[tuple[str, str], dict[str, Any]],
) -> None:
    overall = summary["qc_verdicts"]["overall"]["status"]
    config_summary = summary["config_summary"]
    metric_help = {
        "Long high-quality ratio": "Fraction of reads that meet both the configured minimum read length and minimum read Q-score thresholds.",
        "Correct anchor order ratio": "Fraction of reads whose detected anchors follow the configured expected order.",
        "No-anchor ratio": "Fraction of reads where none of the configured anchors were detected.",
        "Reversed read ratio": "Fraction of reads whose best structure classification was obtained after reverse-complementing the read.",
        "High-N read ratio": "Fraction of reads whose N content reaches or exceeds the configured high-N cutoff.",
        "Total reads": "Total number of FASTQ records processed in this run.",
        "Total bases": "Sum of base counts across all processed reads.",
        "Median read length": "Median read length across all processed reads.",
        "N50": "Read length such that reads of this length or longer contribute at least half of all sequenced bases.",
        "Median read Qscore": "Median per-read Q-score across all reads, using the configured read-level Q-score method.",
        "5' truncation": "Fraction of reads whose first expected anchor starts farther from the 5' end than the configured terminal offset threshold.",
        "3' truncation": "Fraction of reads whose last expected anchor ends farther from the 3' end than the configured terminal offset threshold.",
        "High-N reads": "Fraction of reads with N content at or above the configured cutoff.",
        "Anchor detection ratio": "Fraction of reads in which each anchor is detected at least once. High detection with low mismatch burden usually indicates a strong and well-positioned motif.",
        "Anchor mismatch burden": "For each anchor, this is the average mismatch count of the best-scoring hit among reads where that anchor was actually detected. A value near 0 means the anchor is typically recovered with an exact or near-exact sequence match; larger values suggest sequence drift, weak specificity, or false-positive pressure.",
        "Anchor mean relative position": "Average normalized start position of the best anchor hit along the read length, where 0% is the 5' end and 100% is the 3' end. This helps confirm whether anchors appear where the library design expects them to appear.",
        "Detection ratio": "Fraction of reads in which the anchor was detected at least once.",
        "Mean mismatches": "Average mismatch count of the best-scoring hit for this anchor among anchor-positive reads only. Reads where the anchor is absent are not included in this average.",
        "Multi-hit ratio": "Fraction of anchor-positive reads that contain more than one hit for the same anchor. Elevated values can indicate repeated motifs, internal adapters, or overly permissive matching.",
        "Mean hits / positive read": "Average number of hits for this anchor among reads where the anchor was detected. Values materially above 1 suggest duplicated or internal occurrences.",
        "Read structure classification": "Counts of reads grouped by the structure class inferred from anchor presence and order.",
        "QC failure buckets": "Primary read-level QC failure reason assigned to each read after applying priority ordering.",
        "Read structure heatmap": "Fixed-size aggregated heatmap of anchor occupancy across normalized read position. Each row pools reads sharing the same QC bucket and structure class, with forward and reverse-complement reads merged so the plot stays readable even for very large libraries.",
        "QC Buckets": "Counts and fractions for the primary QC failure bucket assigned to each read.",
        "Structure Orientation": "Counts of each structure class stratified by forward versus reverse-complement orientation. This is useful for spotting strand inversions or reverse-complement enrichment outside the expected full-structure population.",
    }
    config_items = "".join(
        f"<div class='config-row'><span>{escape(label)}</span><strong>{escape(value)}</strong></div>"
        for label, value in [
            ("Rust", f"{summary['rust_accelerator'].get('mode', 'unknown')}"),
            ("Anchors", ", ".join(config_summary["anchors"]) or "none"),
            ("Expected order", " -> ".join(config_summary["expected_order"]) or "none"),
            ("Min long read", f"{config_summary['long_read_min_bp']:,} bp"),
            ("Min read Q", f"{config_summary['long_read_min_q']:.2f}"),
            ("Terminal offset", f"{config_summary['terminal_anchor_max_offset']:.0%}"),
            ("High-N cutoff", f"{config_summary['high_n_fraction']:.0%}"),
            ("Export non-full FASTQ", "enabled" if config_summary["export_non_full_structure_fastq"] else "disabled"),
        ]
    )
    overview_cards = "".join(f"<div class='metric-card'><div class='metric-label'>{_tooltip_label(label, metric_help[label])}</div><div class='metric-value'>{escape(value)}</div></div>" for label, value in [("Total reads", f"{summary['total_reads']:,}"), ("Total bases", f"{summary['total_bases']:,}"), ("Median read length", f"{summary['median_read_length']:,} bp"), ("N50", f"{summary['n50']:,} bp"), ("Median read Qscore", f"{summary['median_read_qscore']:.2f}"), ("5' truncation", f"{summary['five_prime_truncation_ratio']:.1%}"), ("3' truncation", f"{summary['three_prime_truncation_ratio']:.1%}"), ("High-N reads", f"{summary['high_n_read_ratio']:.1%}")])
    verdict_cards = "".join(f"<div class='metric-card verdict-{escape(item['status'])}'><div class='metric-top'><div class='metric-label'>{_tooltip_label(item['label'], metric_help[item['label']])}</div>{_status_pill(item['status'])}</div><div class='metric-value'>{item['value']:.1%}</div></div>" for key, item in summary["qc_verdicts"].items() if key != "overall")
    anchor_rows = "".join(f"<tr><td>{escape(name)}</td><td>{values['detection_ratio']:.1%}</td><td>{values['mean_best_mismatches']:.2f}</td><td>{values['multi_hit_ratio']:.1%}</td><td>{values['mean_hits_per_positive_read']:.2f}</td></tr>" for name, values in summary["anchor_quality_stats"].items())
    bucket_rows = "".join(f"<tr><td>{escape(bucket)}</td><td>{count:,}</td><td>{summary['qc_bucket_ratios'][bucket]:.1%}</td></tr>" for bucket, count in summary["qc_bucket_counts"].items())
    structure_rows = "".join(f"<tr><td>{escape(label)}</td><td>{escape(orientation)}</td><td>{count:,}</td></tr>" for label, counts in summary["structure_orientation_counts"].items() for orientation, count in counts.items())
    anchor_position_svg = _anchor_position_tracks(anchor_positions)
    density_chart_id = "density-" + "".join(
        ch if ch.isalnum() else "-" for ch in summary["sample_name"]
    )
    length_histogram = _histogram(_value_counts_to_float_map(length_counts), "Read length distribution", "Read length (bp)")
    base_qscore_histogram = _histogram(
        _value_counts_to_float_map(base_qscore_counts),
        "Overall base quality distribution",
        "Phred quality score",
        color=COLOR_QUALITY,
    )
    read_qscore_density = _line_density(
        _value_counts_to_float_map(read_qscore_counts, READ_QSCORE_SCALE),
        "Per-read Qscore density",
        "Read Qscore",
        chart_id=density_chart_id,
    )
    length_vs_qscore = _scatter_binned(scatter_counts, "Read length vs read Qscore")
    structure_heatmap = _heatmap(heatmap_groups, "Read structure heatmap")
    html = f"""<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>scfastq-qc report - {escape(summary['sample_name'])}</title><style>:root{{--bg:#f7f5ef;--paper:rgba(255,255,255,.84);--ink:#12232f;--muted:#5b6b78;--line:rgba(18,35,47,.1);--teal:{COLOR_QUALITY};--blue:{COLOR_PRIMARY};--amber:{COLOR_WARNING};--red:{COLOR_FAILURE};--shadow:0 18px 48px rgba(21,40,54,.1)}}*{{box-sizing:border-box}}body{{margin:0;font-family:'IBM Plex Sans','Segoe UI',sans-serif;color:var(--ink);background:radial-gradient(circle at top left, rgba(35,100,170,.12), transparent 26%),radial-gradient(circle at top right, rgba(45,143,133,.14), transparent 22%),linear-gradient(180deg,#f2efe5 0%,#fbfaf7 100%)}}.page{{max-width:1280px;margin:0 auto;padding:30px 18px 56px}}.hero,section{{background:var(--paper);border:1px solid var(--line);border-radius:24px;box-shadow:var(--shadow)}}.hero{{padding:26px}}section{{margin-top:22px;padding:22px;overflow:visible}}.eyebrow{{letter-spacing:.18em;text-transform:uppercase;color:var(--muted);font-size:.74rem;margin-bottom:10px}}h1{{font-family:Georgia,'Times New Roman',serif;font-size:clamp(2.1rem,4vw,3.4rem);margin:0 0 10px}}h2,h3{{margin:0}}p{{margin:0;color:var(--muted)}}.hero-meta{{display:grid;gap:10px}}.assessment-line{{display:flex;align-items:center;gap:10px;flex-wrap:wrap}}.tooltip{{position:relative;display:inline-flex;cursor:help;vertical-align:middle;z-index:2}}.metric-tip{{margin:0}}.metric-tip-inline{{display:inline-flex;align-items:center;justify-content:flex-start}}.tooltip-trigger{{display:inline;border-bottom:1px dashed rgba(35,100,170,.45);color:inherit;transition:border-color .16s ease,color .16s ease}}.tooltip:hover .tooltip-trigger,.tooltip:focus-within .tooltip-trigger{{border-bottom-color:var(--blue);color:var(--blue)}}.tooltip-content{{position:absolute;left:50%;bottom:calc(100% + 12px);transform:translateX(-50%) translateY(4px);width:min(360px,80vw);padding:12px 14px;border-radius:16px;background:rgba(255,255,255,.98);border:1px solid rgba(18,35,47,.12);box-shadow:0 18px 48px rgba(21,40,54,.16);font-size:.86rem;line-height:1.5;color:var(--ink);opacity:0;pointer-events:none;transition:opacity .16s ease,transform .16s ease;z-index:50;text-align:left;text-transform:none;letter-spacing:normal;font-weight:400}}.tooltip-content::after{{content:'';position:absolute;left:50%;top:100%;width:12px;height:12px;background:rgba(255,255,255,.98);border-right:1px solid rgba(18,35,47,.12);border-bottom:1px solid rgba(18,35,47,.12);transform:translateX(-50%) rotate(45deg)}}.tooltip-below .tooltip-content{{top:calc(100% + 12px);bottom:auto;transform:translateX(-50%) translateY(-4px)}}.tooltip-below .tooltip-content::after{{top:auto;bottom:100%;border-right:none;border-bottom:none;border-left:1px solid rgba(18,35,47,.12);border-top:1px solid rgba(18,35,47,.12)}}.tooltip:hover .tooltip-content,.tooltip:focus-within .tooltip-content{{opacity:1;transform:translateX(-50%) translateY(0)}}.tooltip-below:hover .tooltip-content,.tooltip-below:focus-within .tooltip-content{{transform:translateX(-50%) translateY(0)}}.tooltip-content strong{{display:block;margin-bottom:6px}}.hero-chip{{border:1px solid var(--line);border-radius:18px;background:rgba(255,255,255,.72);padding:14px 16px}}.hero-chip strong{{display:block;margin-bottom:10px}}.config-row{{display:flex;justify-content:space-between;gap:12px;padding:8px 0;border-bottom:1px solid rgba(18,35,47,.08);font-size:.93rem}}.config-row:last-of-type{{border-bottom:none}}.config-row span{{color:var(--muted)}}.config-note{{margin-top:10px;font-size:.88rem;line-height:1.45}}.cards{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px}}.metric-card{{background:rgba(255,255,255,.74);border:1px solid rgba(18,35,47,.08);border-radius:18px;padding:16px;display:grid;gap:10px}}.metric-top{{display:flex;justify-content:space-between;gap:10px;align-items:center}}.metric-label{{color:var(--muted);font-size:.82rem;letter-spacing:.06em;text-transform:uppercase}}.metric-value{{font-size:clamp(1.35rem,2.8vw,1.95rem);font-weight:700}}.status-pill{{display:inline-flex;border-radius:999px;padding:7px 12px;font-size:.76rem;font-weight:700;letter-spacing:.08em}}.status-pass{{background:rgba(45,143,133,.12);color:var(--teal)}}.status-warn{{background:rgba(217,119,6,.12);color:var(--amber)}}.status-fail{{background:rgba(209,73,91,.12);color:var(--red)}}.verdict-pass{{border-color:rgba(45,143,133,.28)}}.verdict-warn{{border-color:rgba(217,119,6,.28)}}.verdict-fail{{border-color:rgba(209,73,91,.28)}}.section-head{{display:flex;justify-content:space-between;gap:16px;align-items:end;margin-bottom:16px}}.media-grid,.two-col{{display:grid;gap:16px}}.media-grid{{grid-template-columns:repeat(2,minmax(0,1fr))}}.two-col{{grid-template-columns:1.2fr 1fr}}.panel{{background:rgba(255,255,255,.74);border:1px solid rgba(18,35,47,.08);border-radius:18px;padding:14px;overflow:visible;position:relative}}.chart-note{{margin-top:10px;font-size:.88rem;line-height:1.5;color:var(--muted)}}table{{width:100%;border-collapse:collapse;font-size:.95rem}}th,td{{padding:12px 14px;border-bottom:1px solid rgba(18,35,47,.08);text-align:left}}th{{background:rgba(18,35,47,.04);color:var(--muted);font-size:.78rem;text-transform:uppercase;letter-spacing:.06em}}tr:last-child td{{border-bottom:none}}.svg-wrap svg{{width:100%;height:auto;display:block}}.interactive-plot{{position:relative}}.plot-note{{margin:0 0 8px;color:var(--muted);font-size:.84rem}}.plot-tip{{position:absolute;z-index:5;padding:6px 8px;border-radius:10px;background:rgba(18,35,47,.92);color:#fff;font-size:.78rem;pointer-events:none;white-space:nowrap}}.heatmap-legend{{display:flex;flex-wrap:wrap;gap:10px 14px;margin:10px 0 4px}}.heatmap-legend-item{{display:inline-flex;align-items:center;gap:8px;font-size:.84rem;color:var(--muted)}}.heatmap-swatch{{width:14px;height:14px;border-radius:4px;border:1px solid rgba(18,35,47,.12);flex:0 0 auto}}code{{word-break:break-all}}@media (max-width:980px){{.cards,.media-grid,.two-col{{grid-template-columns:1fr}}.tooltip-content{{left:auto;right:0;bottom:calc(100% + 10px);transform:translateY(4px)}}.tooltip-content::after{{left:auto;right:10px;transform:rotate(45deg)}}.tooltip-below .tooltip-content{{top:calc(100% + 10px);bottom:auto;transform:translateY(-4px)}}.tooltip-below .tooltip-content::after{{left:auto;right:10px}}.tooltip:hover .tooltip-content,.tooltip:focus-within .tooltip-content,.tooltip-below:hover .tooltip-content,.tooltip-below:focus-within .tooltip-content{{transform:translateY(0)}}}}</style></head><body><div class='page'><header class='hero'><div class='eyebrow'>Structure-aware fastq QC</div><div class='hero-meta'><h1>{escape(summary['sample_name'])}</h1><p>Input FASTQ: <code>{escape(summary['fastq_path'])}</code></p><div class='assessment-line'><span>Overall assessment: {_status_pill(overall)}</span><span class='tooltip metric-tip-inline tooltip-below' tabindex='0' role='note' aria-label='Assessment guide help'><span class='tooltip-trigger'>Assessment guide</span><span class='tooltip-content'><strong>Assessment guide</strong><div>`PASS`: all tracked verdicts are within configured acceptable bounds.</div><div>`WARN`: at least one tracked verdict crossed a warning threshold, but none crossed a fail threshold.</div><div>`FAIL`: at least one tracked verdict crossed a fail threshold.</div><div style='margin-top:8px;'>Current verdicts: Long high-quality {summary['qc_verdicts']['long_high_quality_ratio']['status'].upper()}, Anchor order {summary['qc_verdicts']['correct_anchor_order_ratio']['status'].upper()}, No-anchor {summary['qc_verdicts']['no_anchor_ratio']['status'].upper()}, Reversed reads {summary['qc_verdicts']['reversed_read_ratio']['status'].upper()}, High-N {summary['qc_verdicts']['high_n_read_ratio']['status'].upper()}.</div></span></span></div></div></header><section><div class='section-head'><div><h2>Run Configuration</h2><p>Key user-provided parameters and runtime backend details.</p></div></div><div class='hero-chip'><strong>Configuration</strong>{config_items}<p class='config-note'>{escape(summary['rust_accelerator'].get('detail',''))}</p></div></section><section><div class='section-head'><div><h2>QC Verdicts</h2><p>Threshold-based pass, warn, fail calls for the core basic-QC indicators.</p></div></div><div class='cards'>{verdict_cards}</div></section><section><div class='section-head'><div><h2>Overview</h2><p>Yield, read quality, truncation, and contamination proxies at a glance.</p></div></div><div class='cards'>{overview_cards}</div></section><section><div class='section-head'><div><h2>Yield And Quality</h2><p>Baseline library quality by read length and Q-score distributions.</p></div></div><div class='media-grid'><div class='panel svg-wrap'>{length_histogram}<p class='chart-note'>Histogram of per-read lengths summarized from streaming length counts. Hover over each bar to inspect the bin span and how many reads fall into that interval.</p></div><div class='panel svg-wrap'>{base_qscore_histogram}<p class='chart-note'>Histogram of base-level Phred scores pooled across all bases from all reads. Hover over each bar to inspect the score interval and base count.</p></div><div class='panel svg-wrap'>{read_qscore_density}</div><div class='panel svg-wrap'>{length_vs_qscore}</div></div></section><section><div class='section-head'><div><h2>Anchors</h2><p>Detection rate alone is not enough; mismatch burden and hit multiplicity expose weak motifs and internal artifacts.</p></div></div><div class='media-grid'><div class='panel svg-wrap'>{_bar_chart(list(summary['anchor_detection_ratio'].items()), "Anchor detection ratio", "Fraction of reads", color=COLOR_PRIMARY, percent=True)}<p class='chart-note'>{escape(metric_help['Anchor detection ratio'])}</p></div><div class='panel svg-wrap'>{_bar_chart([(name, values['mean_best_mismatches']) for name, values in summary['anchor_quality_stats'].items()], "Anchor mismatch burden", "Mean best-hit mismatches", color=COLOR_WARNING)}<p class='chart-note'>{escape(metric_help['Anchor mismatch burden'])}</p></div></div><div class='panel svg-wrap' style='margin-top:16px;'>{anchor_position_svg}<p class='chart-note'>{escape(metric_help['Anchor mean relative position'])}</p></div><div class='panel' style='margin-top:16px;'><table><thead><tr><th>Anchor</th><th>{_tooltip_label('Detection ratio', metric_help['Detection ratio'])}</th><th>{_tooltip_label('Mean mismatches', metric_help['Mean mismatches'])}</th><th>{_tooltip_label('Multi-hit ratio', metric_help['Multi-hit ratio'])}</th><th>{_tooltip_label('Mean hits / positive read', metric_help['Mean hits / positive read'])}</th></tr></thead><tbody>{anchor_rows}</tbody></table></div></section><section><div class='section-head'><div><h2>Structure And Failure Modes</h2><p>Basic-QC only: structure integrity, truncation, orientation, and read-level hygiene.</p></div></div><div class='media-grid'><div class='panel svg-wrap'>{_bar_chart([(key, float(value)) for key, value in summary['structure_counts'].items()], "Read structure classification", "Read count", color=COLOR_PRIMARY)}</div><div class='panel svg-wrap'>{_bar_chart([(key, float(value)) for key, value in summary['qc_bucket_counts'].items()], "QC failure buckets", "Read count", color=COLOR_FAILURE)}</div></div><div class='panel svg-wrap' style='margin-top:16px;'>{structure_heatmap}</div><div class='two-col' style='margin-top:16px;'><div class='panel'><h3>{_tooltip_label('QC Buckets', metric_help['QC Buckets'])}</h3><table><thead><tr><th>Bucket</th><th>Reads</th><th>Ratio</th></tr></thead><tbody>{bucket_rows}</tbody></table></div><div class='panel'><h3>{_tooltip_label('Structure Orientation', metric_help['Structure Orientation'])}</h3><table><thead><tr><th>Class</th><th>Orientation</th><th>Reads</th></tr></thead><tbody>{structure_rows}</tbody></table></div></div></section></div></body></html>"""
    _write_text(path, html)
