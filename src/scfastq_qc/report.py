from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from html import escape
from math import ceil, exp, floor
from pathlib import Path
from statistics import median
from typing import Any, Callable

from .anchors import AnchorHit, find_anchor_hits, prepare_anchors, reverse_complement
from .classify import classify_best_orientation
from .config import AppConfig
from .fastq import read_fastq


@dataclass
class ReadResult:
    read_id: str
    length: int
    read_qscore: float
    structure_label: str
    is_reversed: bool
    order_valid: bool
    hits: list[AnchorHit]


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


def _scale(value: float, lower: float, upper: float, span: float) -> float:
    if upper <= lower:
        return 0.0
    return ((value - lower) / (upper - lower)) * span


def _format_tick(value: float) -> str:
    if abs(value) >= 1000:
        return f"{value:,.0f}"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _format_int_tick(value: float) -> str:
    return f"{int(round(value)):,}"


def _format_percent_tick(value: float) -> str:
    return f"{value * 100:.0f}%"


def _integer_tick_count(max_value: float, cap: int = 5) -> int:
    rounded = int(round(max_value))
    return max(1, min(rounded, cap))


def _axis_ticks(
    left: int,
    top: int,
    width: int,
    height: int,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    x_ticks: int = 5,
    y_ticks: int = 5,
    x_formatter: Callable[[float], str] | None = None,
    y_formatter: Callable[[float], str] | None = None,
    show_x_grid: bool = True,
    show_x_ticks: bool = True,
    show_x_labels: bool = True,
    show_y_grid: bool = True,
    show_y_ticks: bool = True,
    show_y_labels: bool = True,
) -> list[str]:
    parts: list[str] = []
    bottom = top + height
    right = left + width
    x_span = x_max - x_min if x_max != x_min else 1
    y_span = y_max - y_min if y_max != y_min else 1
    x_formatter = x_formatter or _format_tick
    y_formatter = y_formatter or _format_tick

    for idx in range(x_ticks + 1):
        fraction = idx / x_ticks if x_ticks else 0
        x = left + fraction * width
        value = x_min + fraction * x_span
        if show_x_grid:
            parts.append(f"<line class='grid' x1='{x:.2f}' y1='{top}' x2='{x:.2f}' y2='{bottom}'/>")
        if show_x_ticks:
            parts.append(f"<line class='axis-tick' x1='{x:.2f}' y1='{bottom}' x2='{x:.2f}' y2='{bottom + 6}'/>")
        if show_x_labels:
            parts.append(
                f"<text x='{x:.2f}' y='{bottom + 24}' text-anchor='middle' font-size='12'>{escape(x_formatter(value))}</text>"
            )

    for idx in range(y_ticks + 1):
        fraction = idx / y_ticks if y_ticks else 0
        y = bottom - fraction * height
        value = y_min + fraction * y_span
        if show_y_grid:
            parts.append(f"<line class='grid' x1='{left}' y1='{y:.2f}' x2='{right}' y2='{y:.2f}'/>")
        if show_y_ticks:
            parts.append(f"<line class='axis-tick' x1='{left - 6}' y1='{y:.2f}' x2='{left}' y2='{y:.2f}'/>")
        if show_y_labels:
            parts.append(
                f"<text x='{left - 10}' y='{y + 4:.2f}' text-anchor='end' font-size='12'>{escape(y_formatter(value))}</text>"
            )
    return parts


def _centered_histogram_x_ticks(
    left: int,
    top: int,
    width: int,
    height: int,
    lower: float,
    upper: float,
    bins: int,
    formatter: Callable[[float], str],
    label_count: int = 6,
) -> list[str]:
    parts: list[str] = []
    bottom = top + height
    if bins <= 0:
        return parts
    bar_w = width / bins
    step = (upper - lower) / bins
    tick_bins: list[int] = []
    target_count = min(label_count, bins)
    for idx in range(target_count):
        if target_count == 1:
            bin_idx = bins // 2
        else:
            bin_idx = round(idx * (bins - 1) / (target_count - 1))
        if bin_idx not in tick_bins:
            tick_bins.append(bin_idx)
    for idx in tick_bins:
        x = left + (idx + 0.5) * bar_w
        value = lower + (idx + 0.5) * step
        parts.append(f"<line class='axis-tick' x1='{x:.2f}' y1='{bottom}' x2='{x:.2f}' y2='{bottom + 6}'/>")
        parts.append(
            f"<text x='{x:.2f}' y='{bottom + 24}' text-anchor='middle' font-size='12'>{escape(formatter(value))}</text>"
        )
    return parts


def _svg_wrapper(title: str, body: str, width: int = 800, height: int = 400) -> str:
    return f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
<style>
text {{ font-family: Arial, sans-serif; fill: #222; }}
.axis {{ stroke: #333; stroke-width: 1; }}
.axis-tick {{ stroke: #333; stroke-width: 1; }}
.grid {{ stroke: #ddd; stroke-width: 1; }}
</style>
<text x='20' y='30' font-size='22' font-weight='bold'>{escape(title)}</text>
{body}
</svg>"""


def _histogram(values: list[float], title: str, xlabel: str, bins: int = 20, color: str = "#4C78A8") -> str:
    width, height = 800, 400
    left, right, top, bottom = 70, 30, 60, 60
    plot_w, plot_h = width - left - right, height - top - bottom
    if not values:
        return _svg_wrapper(title, "<text x='70' y='120'>No data</text>")
    lower, upper = min(values), max(values)
    if lower == upper:
        lower -= 0.5
        upper += 0.5
    step = (upper - lower) / bins
    counts = [0] * bins
    x_formatter = _format_int_tick if all(float(value).is_integer() for value in values) else _format_tick
    for value in values:
        idx = min(bins - 1, int((value - lower) / step))
        counts[idx] += 1
    max_count = max(counts) or 1
    parts = _axis_ticks(
        left,
        top,
        plot_w,
        plot_h,
        lower,
        upper,
        0,
        max_count,
        y_ticks=_integer_tick_count(max_count),
        y_formatter=_format_int_tick,
        show_x_grid=False,
        show_x_ticks=False,
        show_x_labels=False,
    )
    parts.extend(_centered_histogram_x_ticks(left, top, plot_w, plot_h, lower, upper, bins, x_formatter))
    parts.extend(
        [
            f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>",
            f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{height-bottom}'/>",
        ]
    )
    bar_w = plot_w / bins
    for idx, count in enumerate(counts):
        bar_h = (count / max_count) * plot_h
        x = left + idx * bar_w
        y = height - bottom - bar_h
        parts.append(f"<rect x='{x:.2f}' y='{y:.2f}' width='{max(bar_w-1,1):.2f}' height='{bar_h:.2f}' fill='{color}'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-15}' text-anchor='middle'>{escape(xlabel)}</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle'>Count</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _density_plot(
    values: list[float],
    title: str,
    xlabel: str,
    color: str = "#54A24B",
    bin_width: float = 0.2,
    smoothing_sigma: float = 0.8,
) -> str:
    width, height = 800, 400
    left, right, top, bottom = 70, 30, 60, 60
    plot_w, plot_h = width - left - right, height - top - bottom
    if not values:
        return _svg_wrapper(title, "<text x='70' y='120'>No data</text>")

    lower, upper = min(values), max(values)
    if lower == upper:
        lower -= 0.5
        upper += 0.5
    lower = floor(lower) - 1
    upper = ceil(upper) + 1

    x_formatter = _format_int_tick if all(float(value).is_integer() for value in values) else _format_tick
    bin_count = max(2, int(ceil((upper - lower) / bin_width)))
    counts = [0.0] * bin_count
    for value in values:
        idx = min(bin_count - 1, max(0, int((value - lower) / bin_width)))
        counts[idx] += 1.0
    densities = [count / (len(values) * bin_width) for count in counts]

    sigma_bins = max(smoothing_sigma / bin_width, 1.0)
    kernel_radius = max(1, int(round(3 * sigma_bins)))
    kernel = [exp(-0.5 * (offset / sigma_bins) ** 2) for offset in range(-kernel_radius, kernel_radius + 1)]
    kernel_sum = sum(kernel) or 1.0
    kernel = [weight / kernel_sum for weight in kernel]

    points: list[tuple[float, float]] = []
    peak_density = 0.0
    for idx in range(bin_count):
        density = 0.0
        for kernel_offset, weight in enumerate(kernel, start=-kernel_radius):
            source_idx = idx + kernel_offset
            if 0 <= source_idx < bin_count:
                density += densities[source_idx] * weight
        x_value = lower + (idx + 0.5) * bin_width
        points.append((x_value, density))
        peak_density = max(peak_density, density)

    y_max = peak_density * 1.05 if peak_density > 0 else 1.0
    parts = _axis_ticks(
        left,
        top,
        plot_w,
        plot_h,
        lower,
        upper,
        0,
        y_max,
        x_formatter=x_formatter,
    )
    parts.extend(
        [
            f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>",
            f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{height-bottom}'/>",
        ]
    )
    polyline_points = []
    for x_value, density in points:
        x = left + _scale(x_value, lower, upper, plot_w)
        y = height - bottom - _scale(density, 0, y_max, plot_h)
        polyline_points.append(f"{x:.2f},{y:.2f}")
    baseline_y = height - bottom
    area_points = [f"{left:.2f},{baseline_y:.2f}", *polyline_points, f"{width-right:.2f},{baseline_y:.2f}"]
    parts.append(f"<polygon points='{' '.join(area_points)}' fill='{color}' fill-opacity='0.18'/>")
    parts.append(f"<polyline fill='none' stroke='{color}' stroke-width='2.5' points='{' '.join(polyline_points)}'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-15}' text-anchor='middle'>{escape(xlabel)}</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle'>Probability density</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _cdf(values: list[int], title: str, xlabel: str, color: str = "#F58518") -> str:
    width, height = 800, 400
    left, right, top, bottom = 70, 30, 60, 60
    plot_w, plot_h = width - left - right, height - top - bottom
    if not values:
        return _svg_wrapper(title, "<text x='70' y='120'>No data</text>")
    sorted_values = sorted(values)
    min_v, max_v = sorted_values[0], sorted_values[-1]
    parts = _axis_ticks(
        left,
        top,
        plot_w,
        plot_h,
        min_v,
        max_v if max_v != min_v else min_v + 1,
        0,
        1,
        x_formatter=_format_int_tick,
    )
    parts.extend(
        [
            f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>",
            f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{height-bottom}'/>",
        ]
    )
    pts = []
    for idx, value in enumerate(sorted_values):
        x = left + _scale(value, min_v, max_v if max_v != min_v else min_v + 1, plot_w)
        y = height - bottom - ((idx + 1) / len(sorted_values)) * plot_h
        pts.append(f"{x:.2f},{y:.2f}")
    parts.append(f"<polyline fill='none' stroke='{color}' stroke-width='2' points='{' '.join(pts)}'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-15}' text-anchor='middle'>{escape(xlabel)}</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle'>Cumulative fraction</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _bar_chart(
    items: list[tuple[str, float]],
    title: str,
    ylabel: str,
    color: str = "#E45756",
    y_max: float | None = None,
    y_formatter: Callable[[float], str] | None = None,
    value_formatter: Callable[[float], str] | None = None,
) -> str:
    width, height = 860, 500
    left, right, top, bottom = 80, 40, 70, 130
    plot_w, plot_h = width - left - right, height - top - bottom
    if not items:
        return _svg_wrapper(title, "<text x='70' y='120'>No data</text>", width, height)
    if y_max is None:
        y_max = max(value for _, value in items) or 1
    y_formatter = y_formatter or _format_int_tick
    value_formatter = value_formatter or y_formatter
    parts = _axis_ticks(
        left,
        top,
        plot_w,
        plot_h,
        0,
        max(len(items), 1),
        0,
        y_max,
        y_ticks=_integer_tick_count(y_max) if y_formatter is _format_int_tick else 5,
        y_formatter=y_formatter,
        show_x_grid=False,
        show_x_ticks=False,
        show_x_labels=False,
    )
    parts.extend(
        [
            f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>",
            f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{height-bottom}'/>",
        ]
    )
    slot_w = plot_w / len(items)
    bar_w = min(slot_w * 0.55, 140)
    for idx, (label, value) in enumerate(items):
        height_px = (value / y_max) * plot_h if y_max else 0
        x = left + idx * slot_w + (slot_w - bar_w) / 2
        y = height - bottom - height_px
        parts.append(f"<rect x='{x:.2f}' y='{y:.2f}' width='{bar_w:.2f}' height='{height_px:.2f}' fill='{color}' rx='4' ry='4'/>")
        tx = left + idx * slot_w + slot_w / 2
        parts.append(f"<text x='{tx:.2f}' y='{height-bottom+28}' text-anchor='end' transform='rotate(-28 {tx:.2f},{height-bottom+28})'>{escape(label)}</text>")
        parts.append(f"<text x='{tx:.2f}' y='{y-8:.2f}' text-anchor='middle' font-size='12'>{escape(value_formatter(value))}</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle'>{escape(ylabel)}</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _scatter(lengths: list[int], read_qscores: list[float], title: str) -> str:
    width, height = 800, 420
    left, right, top, bottom = 70, 30, 60, 60
    plot_w, plot_h = width - left - right, height - top - bottom
    if lengths and read_qscores:
        min_x, max_x = min(lengths), max(lengths)
        min_y, max_y = min(read_qscores), max(read_qscores)
        parts = _axis_ticks(
            left,
            top,
            plot_w,
            plot_h,
            min_x,
            max_x if max_x != min_x else min_x + 1,
            min_y,
            max_y if max_y != min_y else min_y + 1,
            x_formatter=_format_int_tick,
        )
        parts.extend(
            [
                f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>",
                f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{height-bottom}'/>",
            ]
        )
        for xval, yval in zip(lengths, read_qscores):
            x = left + _scale(xval, min_x, max_x if max_x != min_x else min_x + 1, plot_w)
            y = height - bottom - _scale(yval, min_y, max_y if max_y != min_y else min_y + 1, plot_h)
            parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='3' fill='#54A24B' fill-opacity='0.55'/>")
    else:
        parts = [
            f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>",
            f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{height-bottom}'/>",
        ]
    parts.append(f"<text x='{width/2:.0f}' y='{height-15}' text-anchor='middle'>Read length (bp)</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle'>Read Qscore</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _multi_hist(distributions: dict[str, list[float]], title: str) -> str:
    width, height = 800, 420
    left, right, top, bottom = 70, 30, 60, 60
    plot_w, plot_h = width - left - right, height - top - bottom
    colors = ["#E45756", "#4C78A8", "#54A24B", "#F58518", "#B279A2"]
    all_values = [value for values in distributions.values() for value in values]
    if not all_values:
        return _svg_wrapper(title, "<text x='70' y='120'>No data</text>", width, height)
    lower, upper = min(all_values), max(all_values)
    bins = 25
    step = (upper - lower) / bins if upper != lower else 1
    max_count = 1
    counts_by_series: list[tuple[str, list[int]]] = []
    for name, values in distributions.items():
        counts = [0] * bins
        for value in values:
            bin_idx = min(bins - 1, int((value - lower) / step))
            counts[bin_idx] += 1
        max_count = max(max_count, max(counts) if counts else 0)
        counts_by_series.append((name, counts))
    parts = _axis_ticks(
        left,
        top,
        plot_w,
        plot_h,
        lower,
        upper if upper != lower else lower + 1,
        0,
        max_count,
        y_ticks=_integer_tick_count(max_count),
        y_formatter=_format_int_tick,
    )
    parts.extend(
        [
            f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>",
            f"<line class='axis' x1='{left}' y1='{top}' x2='{left}' y2='{height-bottom}'/>",
        ]
    )
    series_points: list[tuple[str, str]] = []
    for idx, (name, counts) in enumerate(counts_by_series):
        points = []
        for b_idx, count in enumerate(counts):
            x = left + (b_idx / (bins - 1 if bins > 1 else 1)) * plot_w
            y = height - bottom - (count / max_count) * plot_h
            points.append(f"{x:.2f},{y:.2f}")
        series_points.append((name, f"<polyline fill='none' stroke='{colors[idx % len(colors)]}' stroke-width='2' points='{' '.join(points)}'/>") )
    legend_y = top
    for idx, (name, polyline) in enumerate(series_points):
        parts.append(polyline)
        parts.append(f"<rect x='{width-right-180}' y='{legend_y + idx*22 - 10}' width='14' height='14' fill='{colors[idx % len(colors)]}'/>")
        parts.append(f"<text x='{width-right-160}' y='{legend_y + idx*22 + 2}'>{escape(name)}</text>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-15}' text-anchor='middle'>Relative read position</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle'>Hit count</text>")
    return _svg_wrapper(title, ''.join(parts), width, height)


def _heatmap(read_results: list[ReadResult], max_reads: int) -> str:
    width, height = 950, 500
    left, right, top, bottom = 70, 30, 60, 50
    plot_w, plot_h = width - left - right, height - top - bottom
    sampled = read_results[:max_reads]
    row_h = max(1, plot_h / max(len(sampled), 1))
    col_w = plot_w / 100
    parts = _axis_ticks(
        left,
        top,
        plot_w,
        plot_h,
        0,
        100,
        0,
        len(sampled) if sampled else 1,
        y_ticks=_integer_tick_count(len(sampled) if sampled else 1),
        y_formatter=_format_int_tick,
    )
    for row_idx, result in enumerate(sampled):
        y = top + row_idx * row_h
        for col_idx in range(100):
            parts.append(f"<rect x='{left + col_idx * col_w:.2f}' y='{y:.2f}' width='{col_w:.2f}' height='{row_h:.2f}' fill='#f5f5f5'/>")
        if result.length > 0:
            for hit in result.hits:
                start = min(99, int((hit.start / result.length) * 100))
                end = min(100, max(start + 1, int((hit.end / result.length) * 100)))
                for col_idx in range(start, end):
                    parts.append(f"<rect x='{left + col_idx * col_w:.2f}' y='{y:.2f}' width='{col_w:.2f}' height='{row_h:.2f}' fill='#4C78A8'/>")
    parts.append(f"<line class='axis' x1='{left}' y1='{height-bottom}' x2='{width-right}' y2='{height-bottom}'/>")
    parts.append(f"<text x='{width/2:.0f}' y='{height-15}' text-anchor='middle'>Normalized read position</text>")
    parts.append(f"<text x='20' y='{height/2:.0f}' transform='rotate(-90 20,{height/2:.0f})' text-anchor='middle'>Reads</text>")
    return _svg_wrapper("Read structure heatmap (anchor occupancy)", ''.join(parts), width, height)


def run_qc(fastq_path: str, config: AppConfig, outdir: str) -> dict[str, Any]:
    out_path = Path(outdir)
    fig_dir = out_path / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    prepared_anchors = prepare_anchors(config.anchors)
    heatmap_reads: list[ReadResult] = []
    base_qscores: list[float] = []
    lengths: list[int] = []
    read_qscores: list[float] = []
    anchor_counts: Counter[str] = Counter()
    anchor_positions: dict[str, list[float]] = defaultdict(list)
    structure_counts: Counter[str] = Counter()
    structure_orientation_counts: Counter[tuple[str, str]] = Counter()
    long_high_quality = 0
    order_valid_count = 0
    reversed_read_count = 0
    total_reads = 0
    total_bases = 0

    for read in read_fastq(fastq_path):
        total_reads += 1
        total_bases += read.length
        forward_hits = find_anchor_hits(read.sequence, prepared_anchors)
        reverse_hits = find_anchor_hits(reverse_complement(read.sequence), prepared_anchors)
        structure, hits = classify_best_orientation(forward_hits, reverse_hits, config.structure.expected_order)
        if len(heatmap_reads) < config.thresholds.heatmap_max_reads:
            heatmap_reads.append(
                ReadResult(read.name, read.length, read.read_qscore, structure.label, structure.is_reversed, structure.order_valid, hits)
            )
        lengths.append(read.length)
        base_qscores.extend(float(q) for q in read.phred_scores)
        read_qscores.append(read.read_qscore)
        structure_counts[structure.label] += 1
        structure_orientation_counts[(structure.label, "reversed" if structure.is_reversed else "forward")] += 1
        if structure.order_valid:
            order_valid_count += 1
        if structure.is_reversed:
            reversed_read_count += 1
        if read.length >= config.thresholds.long_read_min_bp and read.read_qscore >= config.thresholds.long_read_min_q:
            long_high_quality += 1
        seen = set()
        for hit in hits:
            if hit.anchor_name not in seen:
                anchor_counts[hit.anchor_name] += 1
                seen.add(hit.anchor_name)
            if read.length > 0:
                anchor_positions[hit.anchor_name].append(hit.start / read.length)

    summary = {
        "sample_name": config.sample_name,
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mean_read_length": (total_bases / total_reads) if total_reads else 0,
        "median_read_length": median(lengths) if lengths else 0,
        "n50": compute_n50(lengths),
        "min_read_length": min(lengths) if lengths else 0,
        "max_read_length": max(lengths) if lengths else 0,
        "median_read_qscore": median(read_qscores) if read_qscores else 0,
        "long_high_quality_ratio": (long_high_quality / total_reads) if total_reads else 0.0,
        "reversed_read_ratio": (reversed_read_count / total_reads) if total_reads else 0.0,
        "anchor_detection_ratio": {anchor.name: (anchor_counts.get(anchor.name, 0) / total_reads) if total_reads else 0.0 for anchor in config.anchors},
        "structure_counts": dict(structure_counts),
        "structure_orientation_counts": {
            label: {
                orientation: structure_orientation_counts.get((label, orientation), 0)
                for orientation in ("forward", "reversed")
                if structure_orientation_counts.get((label, orientation), 0) > 0
            }
            for label in structure_counts
        },
        "correct_anchor_order_ratio": (order_valid_count / total_reads) if total_reads else 0.0,
    }
    _write_text(out_path / "summary.json", json.dumps(summary, indent=2))

    _write_text(fig_dir / "read_length_hist.svg", _histogram([float(v) for v in lengths], "Read length distribution", "Read length (bp)"))
    _write_text(fig_dir / "read_length_cdf.svg", _cdf(lengths, "Read length cumulative distribution", "Read length (bp)"))
    _write_text(
        fig_dir / "mean_q_hist.svg",
        _density_plot(base_qscores, "Per-base quality score density", "Quality score", color="#54A24B"),
    )
    _write_text(fig_dir / "length_q_scatter.svg", _scatter(lengths, read_qscores, "Read length vs read Qscore"))
    _write_text(
        fig_dir / "anchor_detect_bar.svg",
        _bar_chart(
            list(summary["anchor_detection_ratio"].items()),
            "Anchor detection ratio",
            "Fraction of reads",
            y_max=1.0,
            y_formatter=_format_percent_tick,
            value_formatter=_format_percent_tick,
        ),
    )
    _write_text(fig_dir / "anchor_position_density.svg", _multi_hist(anchor_positions, "Anchor relative position distribution"))
    _write_text(
        fig_dir / "structure_class_bar.svg",
        _bar_chart([(k, float(v)) for k, v in structure_counts.items()], "Read structure classification", "Read count", color="#72B7B2"),
    )
    _write_text(fig_dir / "structure_heatmap.svg", _heatmap(heatmap_reads, config.thresholds.heatmap_max_reads))
    _write_html_report(out_path / "report.html", summary)
    return summary


def _write_html_report(path: Path, summary: dict[str, Any]) -> None:
    rows = "".join(
        f"<tr><td>{escape(label)}</td><td>{escape(orientation)}</td><td>{count}</td></tr>"
        for label, counts in sorted(
            summary["structure_orientation_counts"].items(),
            key=lambda item: sum(item[1].values()),
            reverse=True,
        )
        for orientation, count in sorted(counts.items(), key=lambda item: item[1], reverse=True)
    )
    anchor_rows = "".join(
        f"<tr><td>{escape(name)}</td><td>{ratio:.2%}</td></tr>"
        for name, ratio in summary["anchor_detection_ratio"].items()
    )
    read_length_stats = [
        ("N50", f"{summary['n50']:,} bp"),
        ("Median length", f"{summary['median_read_length']:,} bp"),
        ("Mean length", f"{summary['mean_read_length']:.1f} bp"),
        ("Longest read", f"{summary['max_read_length']:,} bp"),
        ("Shortest read", f"{summary['min_read_length']:,} bp"),
    ]
    read_length_stat_html = "".join(
        f"<div class='stat-item'><div class='stat-label'>{escape(label)}</div><div class='stat-value'>{escape(value)}</div></div>"
        for label, value in read_length_stats
    )
    html = f"""<!DOCTYPE html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <title>scfastq-qc report - {escape(summary['sample_name'])}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem auto; max-width: 1150px; color: #222; }}
    h1, h2 {{ color: #1f4e79; }}
    .cards {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 1rem; }}
    .card {{ border: 1px solid #ddd; border-radius: 8px; padding: 1rem; background: #fafafa; }}
    .label {{ font-size: 0.85rem; color: #666; }}
    .value {{ font-size: 1.5rem; font-weight: bold; }}
    object {{ width: 100%; min-height: 420px; border: 1px solid #ddd; margin-bottom: 1rem; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 1rem; margin-bottom: 1rem; }}
    th, td {{ border: 1px solid #ddd; padding: 0.5rem; text-align: left; }}
    .split {{ display: grid; grid-template-columns: 2fr 1fr; gap: 1rem; }}
    .feature-panel {{ display: grid; grid-template-columns: 2fr 1fr; gap: 1rem; align-items: stretch; margin-bottom: 1rem; }}
    .feature-panel object {{ margin-bottom: 0; min-height: 460px; }}
    .stat-panel {{ border: 1px solid #ddd; border-radius: 8px; padding: 1.25rem; background: #fafafa; display: flex; flex-direction: column; gap: 0.9rem; }}
    .stat-panel h3 {{ margin: 0; font-size: 1.05rem; color: #1f4e79; }}
    .stat-item {{ padding-bottom: 0.75rem; border-bottom: 1px solid #e5e5e5; }}
    .stat-item:last-child {{ padding-bottom: 0; border-bottom: none; }}
    .stat-label {{ font-size: 0.85rem; color: #666; }}
    .stat-value {{ font-size: 1.35rem; font-weight: bold; color: #222; }}
    @media (max-width: 900px) {{
      .cards, .split, .feature-panel {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <h1>scfastq-qc report</h1>
  <p><strong>Sample:</strong> {escape(summary['sample_name'])}</p>
  <div class='cards'>
    <div class='card'><div class='label'>Total reads</div><div class='value'>{summary['total_reads']}</div></div>
    <div class='card'><div class='label'>Total bases</div><div class='value'>{summary['total_bases']}</div></div>
    <div class='card'><div class='label'>Median read length</div><div class='value'>{summary['median_read_length']}</div></div>
    <div class='card'><div class='label'>N50</div><div class='value'>{summary['n50']}</div></div>
    <div class='card'><div class='label'>Median read Qscore</div><div class='value'>{summary['median_read_qscore']:.2f}</div></div>
    <div class='card'><div class='label'>Long high-quality ratio</div><div class='value'>{summary['long_high_quality_ratio']:.2%}</div></div>
    <div class='card'><div class='label'>Reversed read ratio</div><div class='value'>{summary['reversed_read_ratio']:.2%}</div></div>
    <div class='card'><div class='label'>Correct anchor order ratio</div><div class='value'>{summary['correct_anchor_order_ratio']:.2%}</div></div>
    <div class='card'><div class='label'>Anchor types</div><div class='value'>{len(summary['anchor_detection_ratio'])}</div></div>
  </div>

  <h2>Yield and quality</h2>
  <div class='feature-panel'>
    <object type='image/svg+xml' data='figures/read_length_hist.svg'></object>
    <div class='stat-panel'>
      <h3>Read length summary</h3>
      {read_length_stat_html}
    </div>
  </div>
  <object type='image/svg+xml' data='figures/read_length_cdf.svg'></object>
  <object type='image/svg+xml' data='figures/mean_q_hist.svg'></object>
  <object type='image/svg+xml' data='figures/length_q_scatter.svg'></object>

  <h2>Anchor detection</h2>
  <div class='split'>
    <div>
      <object type='image/svg+xml' data='figures/anchor_detect_bar.svg'></object>
      <object type='image/svg+xml' data='figures/anchor_position_density.svg'></object>
    </div>
    <div>
      <table>
        <thead><tr><th>Anchor</th><th>Detection ratio</th></tr></thead>
        <tbody>{anchor_rows}</tbody>
      </table>
    </div>
  </div>

  <h2>Structure classification</h2>
  <object type='image/svg+xml' data='figures/structure_class_bar.svg'></object>
  <object type='image/svg+xml' data='figures/structure_heatmap.svg'></object>
  <table>
    <thead><tr><th>Structure class</th><th>Orientation</th><th>Read count</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</body>
</html>
"""
    _write_text(path, html)
