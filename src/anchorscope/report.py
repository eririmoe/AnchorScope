from __future__ import annotations

import json
import logging
from collections import Counter
from html import escape
from pathlib import Path
from typing import Any

from .anchors import AnchorHit
from .config import AppConfig, ThresholdConfig
from .fastq import FastqRead
from .segments import SegmentObservation

logger = logging.getLogger(__name__)

READ_QSCORE_SCALE = 100
SCATTER_LENGTH_BIN_BP = 50
SCATTER_QSCORE_BIN = 0.5
HEATMAP_BINS = 100

COLOR_PRIMARY = "#2364aa"
COLOR_QUALITY = "#2d8f85"
COLOR_WARNING = "#d97706"
COLOR_FAILURE = "#d1495b"


def compute_n50_counts(length_counts: Counter[int]) -> int:
    if not length_counts:
        return 0
    target = sum(length * count for length, count in length_counts.items()) / 2
    running = 0
    for length in sorted(length_counts, reverse=True):
        running += length * length_counts[length]
        if running >= target:
            return length
    return 0


def _weighted_median(value_counts: Counter[int], scale: float = 1.0) -> float:
    if not value_counts:
        return 0.0
    total = sum(value_counts.values())
    targets = ((total - 1) // 2, total // 2)
    found: list[int] = []
    cumulative = 0
    for value in sorted(value_counts):
        next_cumulative = cumulative + value_counts[value]
        for target in targets:
            if cumulative <= target < next_cumulative:
                found.append(value)
        cumulative = next_cumulative
    return (sum(found) / len(found) / scale) if found else 0.0


def _value_counts_to_float_map(
    value_counts: Counter[int], scale: float = 1.0
) -> dict[float, int]:
    return {value / scale: count for value, count in value_counts.items()}


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _format_num(value: float, decimals: int = 1) -> str:
    if float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:,.{decimals}f}"


def _best_hits_by_anchor(hits: list[AnchorHit]) -> dict[str, AnchorHit]:
    best: dict[str, AnchorHit] = {}
    for hit in hits:
        current = best.get(hit.anchor_name)
        score = (hit.mismatches, hit.start, hit.end)
        if current is None or score < (current.mismatches, current.start, current.end):
            best[hit.anchor_name] = hit
    return best


def _compute_terminal_offsets(
    best_hits: dict[str, AnchorHit],
    expected_order: list[str],
    read_length: int,
) -> tuple[float | None, float | None]:
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
    segment_observations: dict[str, SegmentObservation] | None = None,
) -> list[str]:
    flags: list[str] = []
    thresholds = config.thresholds
    if read.length == 0:
        flags.append("empty_read")
    if read.length < thresholds.long_read_min_bp:
        flags.append("too_short")
    if read.read_qscore < thresholds.long_read_min_q:
        flags.append("low_read_q")
    if read.n_fraction >= thresholds.high_n_fraction:
        flags.append("high_n")
    if read.invalid_base_fraction > 0:
        flags.append("invalid_bases")

    structure_flags = {
        "no_anchor_detected": "no_anchor",
        "missing_5p_anchor": "missing_5p",
        "missing_3p_anchor": "missing_3p",
        "anchor_order_invalid": "order_invalid",
        "duplicated_anchor": "multi_anchor",
        "internal_5p_anchor": "internal_adapter",
        "internal_3p_anchor": "internal_adapter",
        "concatemer_candidate": "concatemer_candidate",
        "structure_rule_violation": "structure_rule_violation",
        "orientation_unexpected": "orientation_unexpected",
        "partial_structure": "partial_structure",
    }
    mapped = structure_flags.get(structure_label)
    if mapped:
        flags.append(mapped)
    if five_prime_offset is not None and five_prime_offset > thresholds.terminal_anchor_max_offset:
        flags.append("five_prime_truncated")
    if three_prime_offset is not None and three_prime_offset > thresholds.terminal_anchor_max_offset:
        flags.append("three_prime_truncated")
    if len(hits) > max(len(config.structure.expected_order), 1) * 2:
        flags.append("dense_anchor_hits")
    for observation in (segment_observations or {}).values():
        if observation.status not in {"pass", "missing"}:
            flags.append(f"segment_{observation.status}")
    return list(dict.fromkeys(flags))


def _primary_qc_bucket(flags: list[str]) -> str:
    priority = [
        "empty_read",
        "too_short",
        "low_read_q",
        "invalid_bases",
        "high_n",
        "no_anchor",
        "internal_adapter",
        "concatemer_candidate",
        "structure_rule_violation",
        "orientation_unexpected",
        "multi_anchor",
        "order_invalid",
        "partial_structure",
        "missing_5p",
        "missing_3p",
        "five_prime_truncated",
        "three_prime_truncated",
        "segment_too_short",
        "segment_too_long",
        "dense_anchor_hits",
    ]
    return next((name for name in priority if name in flags), "pass")


def _verdict(value: float, warn: float, fail: float, higher_is_better: bool) -> str:
    if higher_is_better:
        return "fail" if value < fail else "warn" if value < warn else "pass"
    return "fail" if value > fail else "warn" if value > warn else "pass"


def _build_qc_verdicts(
    summary: dict[str, Any],
    thresholds: ThresholdConfig,
    expected_orientation: str = "either",
) -> dict[str, dict[str, Any]]:
    no_anchor = summary["qc_bucket_ratios"].get("no_anchor", 0.0)
    reversed_ratio = summary["reversed_read_ratio"]
    reversed_status = (
        "pass"
        if expected_orientation == "either"
        else _verdict(
            reversed_ratio if expected_orientation == "forward" else 1.0 - reversed_ratio,
            thresholds.warn_reversed_read_ratio,
            thresholds.fail_reversed_read_ratio,
            False,
        )
    )
    reversed_label = (
        "Reversed read ratio (informational)"
        if expected_orientation == "either"
        else "Unexpected orientation ratio"
    )
    reversed_value = (
        reversed_ratio
        if expected_orientation in {"either", "forward"}
        else 1.0 - reversed_ratio
    )
    verdicts: dict[str, dict[str, Any]] = {
        "long_high_quality_ratio": {
            "label": "Long high-quality ratio",
            "value": summary["long_high_quality_ratio"],
            "status": _verdict(
                summary["long_high_quality_ratio"],
                thresholds.warn_long_high_quality_ratio,
                thresholds.fail_long_high_quality_ratio,
                True,
            ),
        },
        "correct_anchor_order_ratio": {
            "label": "Correct anchor order ratio",
            "value": summary["correct_anchor_order_ratio"],
            "status": _verdict(
                summary["correct_anchor_order_ratio"],
                thresholds.warn_correct_anchor_order_ratio,
                thresholds.fail_correct_anchor_order_ratio,
                True,
            ),
        },
        "no_anchor_ratio": {
            "label": "No-anchor ratio",
            "value": no_anchor,
            "status": _verdict(
                no_anchor,
                thresholds.warn_no_anchor_ratio,
                thresholds.fail_no_anchor_ratio,
                False,
            ),
        },
        "reversed_read_ratio": {
            "label": reversed_label,
            "value": reversed_value,
            "status": reversed_status,
            "informational": expected_orientation == "either",
        },
        "high_n_read_ratio": {
            "label": "High-N read ratio",
            "value": summary["high_n_read_ratio"],
            "status": _verdict(
                summary["high_n_read_ratio"],
                thresholds.warn_high_n_ratio,
                thresholds.fail_high_n_ratio,
                False,
            ),
        },
    }
    scored = [
        item
        for item in verdicts.values()
        if not item.get("informational", False)
    ]
    overall = (
        "fail"
        if any(item["status"] == "fail" for item in scored)
        else "warn"
        if any(item["status"] == "warn" for item in scored)
        else "pass"
    )
    verdicts["overall"] = {"label": "Overall QC", "value": None, "status": overall}
    return verdicts


def _status_pill(status: str) -> str:
    return f"<span class='pill {escape(status)}'>{escape(status.upper())}</span>"


def _svg(title: str, body: str, width: int = 760, height: int = 300) -> str:
    return (
        f"<svg xmlns='http://www.w3.org/2000/svg' role='img' "
        f"aria-label='{escape(title)}' viewBox='0 0 {width} {height}'>"
        "<style>text{font-family:Segoe UI,sans-serif;fill:#18313d}"
        ".grid{stroke:#dce5e8;stroke-width:1}.axis{stroke:#60727d;stroke-width:1}</style>"
        f"<text x='18' y='27' font-size='18' font-weight='700'>{escape(title)}</text>{body}</svg>"
    )


def _histogram(
    values: dict[float, int],
    title: str,
    xlabel: str,
    color: str = COLOR_QUALITY,
    bins: int = 24,
) -> str:
    clean = {float(key): count for key, count in values.items() if count > 0}
    if not clean:
        return _svg(title, "<text x='30' y='80'>No data</text>")
    low, high = min(clean), max(clean)
    if low == high:
        low -= 0.5
        high += 0.5
    step = (high - low) / bins
    counts = [0] * bins
    for value, count in clean.items():
        index = min(bins - 1, max(0, int((value - low) / step)))
        counts[index] += count
    maximum = max(counts) or 1
    left, top, plot_w, plot_h = 55, 45, 680, 205
    bar_w = plot_w / bins
    parts = [
        f"<line class='axis' x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}'/>"
    ]
    for index, count in enumerate(counts):
        height = plot_h * count / maximum
        x = left + index * bar_w + 1
        y = top + plot_h - height
        parts.append(
            f"<rect x='{x:.1f}' y='{y:.1f}' width='{max(bar_w-2, 1):.1f}' "
            f"height='{height:.1f}' fill='{color}'><title>{count:,} reads</title></rect>"
        )
    parts.extend(
        [
            f"<text x='{left}' y='274' font-size='11'>{low:.2g}</text>",
            f"<text x='{left+plot_w}' y='274' text-anchor='end' font-size='11'>{high:.2g}</text>",
            f"<text x='{left+plot_w/2}' y='292' text-anchor='middle' font-size='12'>{escape(xlabel)}</text>",
        ]
    )
    return _svg(title, "".join(parts))


def _scatter_binned(
    counts: Counter[tuple[int, int]], title: str
) -> str:
    if not counts:
        return _svg(title, "<text x='30' y='80'>No data</text>")
    points = list(counts.items())
    lengths = [point[0][0] for point in points]
    qscores = [point[0][1] for point in points]
    max_count = max(counts.values())
    low_x, high_x = min(lengths), max(lengths)
    low_y, high_y = min(qscores), max(qscores)
    x_span = max(high_x - low_x, 1)
    y_span = max(high_y - low_y, 1)
    circles = []
    for (length, qscore), count in points:
        x = 55 + ((length - low_x) / x_span) * 680
        y = 250 - ((qscore - low_y) / y_span) * 200
        radius = 2.5 + 6 * (count / max_count) ** 0.5
        circles.append(
            f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{radius:.1f}' fill='{COLOR_PRIMARY}' "
            f"fill-opacity='.5'><title>{length} bp, Q{qscore:.1f}: {count} reads</title></circle>"
        )
    axes = (
        "<line class='axis' x1='55' y1='250' x2='735' y2='250'/>"
        "<line class='axis' x1='55' y1='45' x2='55' y2='250'/>"
        "<text x='395' y='292' text-anchor='middle' font-size='12'>Read length (bp)</text>"
        "<text x='15' y='150' transform='rotate(-90 15 150)' text-anchor='middle' font-size='12'>Read Qscore</text>"
    )
    return _svg(title, axes + "".join(circles))


def _anchor_position_tracks(anchor_positions: dict[str, dict[str, float]]) -> str:
    if not anchor_positions:
        return _svg("Mean anchor position", "<text x='30' y='80'>No anchor observations</text>")
    rows = []
    for index, (name, values) in enumerate(sorted(anchor_positions.items())):
        ratio = values["sum"] / values["count"] if values["count"] else 0.0
        y = 62 + index * 30
        x = 145 + ratio * 570
        rows.append(
            f"<text x='135' y='{y+4}' text-anchor='end' font-size='12'>{escape(name)}</text>"
            f"<line x1='145' y1='{y}' x2='715' y2='{y}' stroke='#dce5e8'/>"
            f"<circle cx='{x:.1f}' cy='{y}' r='6' fill='{COLOR_WARNING}'>"
            f"<title>{ratio:.1%} of read length</title></circle>"
        )
    height = max(150, 95 + len(rows) * 30)
    return _svg("Mean anchor position", "".join(rows), height=height)


def _heatmap(
    groups: dict[tuple[str, str], dict[str, Any]], title: str
) -> str:
    if not groups:
        return _svg(title, "<text x='30' y='80'>No data</text>")
    rows = sorted(groups.items(), key=lambda item: (-item[1]["count"], item[0]))
    rows = rows[:18]
    parts: list[str] = []
    for row_index, ((bucket, structure), values) in enumerate(rows):
        y = 50 + row_index * 22
        count = max(values["count"], 1)
        parts.append(
            f"<text x='205' y='{y+12}' text-anchor='end' font-size='10'>"
            f"{escape(bucket)} / {escape(structure)} ({count:,})</text>"
        )
        covers = values["anchor_cover"]
        for index, cover in enumerate(covers):
            opacity = min(1.0, cover / count)
            parts.append(
                f"<rect x='{215 + index*5.15:.1f}' y='{y}' width='5.2' height='16' "
                f"fill='{COLOR_PRIMARY}' fill-opacity='{opacity:.3f}'/>"
            )
    height = max(150, 75 + len(rows) * 22)
    return _svg(title, "".join(parts), height=height)


def _table_rows(rows: list[list[str]]) -> str:
    return "".join(
        "<tr>" + "".join(f"<td>{value}</td>" for value in row) + "</tr>"
        for row in rows
    )


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
    verdict_cards = "".join(
        "<article class='card'>"
        f"<span class='label'>{escape(item['label'])}</span>"
        f"<strong>{item['value']:.1%}</strong>{_status_pill(item['status'])}</article>"
        for key, item in summary["qc_verdicts"].items()
        if key != "overall"
    )
    overview = [
        ("Total reads", f"{summary['total_reads']:,}"),
        ("Total bases", f"{summary['total_bases']:,}"),
        ("Median length", f"{summary['median_read_length']:,.0f} bp"),
        ("N50", f"{summary['n50']:,} bp"),
        ("Median read Q", f"{summary['median_read_qscore']:.2f}"),
        ("Pass reads", f"{summary['qc_bucket_ratios'].get('pass', 0.0):.1%}"),
        ("Concatemer parents", f"{summary.get('concatemer_parent_reads', 0):,}"),
        ("Matcher", escape(summary["matcher_backend"].get("mode", "unknown"))),
    ]
    overview_cards = "".join(
        f"<article class='card'><span class='label'>{escape(label)}</span><strong>{value}</strong></article>"
        for label, value in overview
    )
    config = summary["config_summary"]
    config_rows = _table_rows(
        [
            [escape(label), escape(str(value))]
            for label, value in [
                ("Anchors", ", ".join(config["anchors"])),
                ("Expected order", " -> ".join(config["expected_order"])),
                ("Expected orientation", config["expected_orientation"]),
                ("Segments", ", ".join(config["segments"]) or "none"),
                ("Barcode recoverability QC", "enabled" if config["barcode_qc_enabled"] else "disabled"),
                ("Minimum long-read length", f"{config['long_read_min_bp']} bp"),
                ("Minimum read Q", f"{config['long_read_min_q']:.2f}"),
            ]
        ]
    )
    anchor_rows = _table_rows(
        [
            [
                escape(name),
                f"{values['detection_ratio']:.1%}",
                f"{values['mean_best_mismatches']:.2f}",
                f"{values['mean_best_insertions']:.2f}",
                f"{values['mean_best_deletions']:.2f}",
                f"{values['multi_hit_ratio']:.1%}",
            ]
            for name, values in summary["anchor_quality_stats"].items()
        ]
    )
    bucket_rows = _table_rows(
        [
            [escape(name), f"{count:,}", f"{summary['qc_bucket_ratios'][name]:.1%}"]
            for name, count in summary["qc_bucket_counts"].items()
        ]
    )
    structure_rows = _table_rows(
        [
            [escape(label), escape(orientation), f"{count:,}"]
            for label, orientations in summary["structure_orientation_counts"].items()
            for orientation, count in orientations.items()
        ]
    )
    segment_rows = _table_rows(
        [
            [
                escape(name),
                f"{values['observed_ratio']:.1%}",
                f"{values['median_length']:.1f}",
                f"{values['median_qscore']:.2f}",
                f"{values['mean_gc_fraction']:.1%}",
                f"{values['within_length_ratio']:.1%}",
            ]
            for name, values in summary.get("segment_qc", {}).items()
        ]
    )
    segment_section = ""
    if segment_rows:
        segment_section = (
            "<section><h2>Anchor-relative segment QC</h2>"
            "<p class='note'>Quality, composition and length statistics for molecular regions bounded by anchors.</p>"
            "<div class='scroll'><table><thead><tr><th>Segment</th><th>Observed</th>"
            "<th>Median length</th><th>Median Q</th><th>Mean GC</th><th>Within bounds</th>"
            f"</tr></thead><tbody>{segment_rows}</tbody></table></div></section>"
        )
    barcode_section = ""
    barcode = summary.get("barcode_qc")
    if barcode:
        barcode_rows = _table_rows(
            [
                [
                    escape(status),
                    f"{count:,}",
                    f"{count / barcode['total_reads']:.1%}" if barcode["total_reads"] else "0.0%",
                ]
                for status, count in barcode["status_counts"].items()
            ]
        )
        barcode_section = (
            "<section><h2>Barcode and UMI recoverability</h2>"
            "<p class='note'>Diagnostic only: raw sequences are assessed but never corrected.</p>"
            "<div class='cards'>"
            f"<article class='card'><span class='label'>Observed</span><strong>{barcode['observed_ratio']:.1%}</strong></article>"
            f"<article class='card'><span class='label'>Recoverable</span><strong>{barcode['recoverable_ratio']:.1%}</strong></article>"
            f"<article class='card'><span class='label'>Ambiguous</span><strong>{barcode['ambiguous_ratio']:.1%}</strong></article>"
            f"<article class='card'><span class='label'>Unique raw barcodes</span><strong>{barcode['unique_raw_barcodes']:,}</strong></article>"
            "</div><div class='scroll'><table><thead><tr><th>Status</th><th>Reads</th><th>Ratio</th>"
            f"</tr></thead><tbody>{barcode_rows}</tbody></table></div></section>"
        )

    plots = [
        _histogram(_value_counts_to_float_map(length_counts), "Read length distribution", "Read length (bp)", COLOR_PRIMARY),
        _histogram(_value_counts_to_float_map(base_qscore_counts), "Base quality distribution", "Phred score"),
        _histogram(_value_counts_to_float_map(read_qscore_counts, READ_QSCORE_SCALE), "Per-read Qscore distribution", "Read Qscore", COLOR_WARNING),
        _scatter_binned(scatter_counts, "Read length vs read Qscore"),
        _anchor_position_tracks(anchor_positions),
        _heatmap(heatmap_groups, "Read structure heatmap"),
    ]
    plot_html = "".join(f"<div class='panel'>{plot}</div>" for plot in plots)
    html = f"""<!doctype html>
<html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>AnchorScope report - {escape(summary['sample_name'])}</title>
<style>
:root{{--ink:#18313d;--muted:#60727d;--line:#dce5e8;--paper:#fff;--bg:#f4f7f6}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,'Segoe UI',sans-serif}}
main{{max-width:1240px;margin:auto;padding:28px 18px 60px}}header,section{{background:var(--paper);border:1px solid var(--line);border-radius:18px;padding:22px;margin-bottom:18px}}
h1{{margin:0 0 8px;font-family:Georgia,serif}}h2{{margin:0 0 8px}}.note,.sub{{color:var(--muted)}}
.cards{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin-top:14px}}.card,.panel{{border:1px solid var(--line);border-radius:14px;padding:14px;background:#fcfdfd}}
.card{{display:grid;gap:8px}}.label{{color:var(--muted);font-size:.78rem;text-transform:uppercase;letter-spacing:.05em}}.card strong{{font-size:1.4rem}}
.plots{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}}.panel svg{{width:100%;height:auto;display:block}}
table{{width:100%;border-collapse:collapse}}th,td{{padding:10px;border-bottom:1px solid var(--line);text-align:left}}th{{color:var(--muted)}}.scroll{{overflow:auto;margin-top:14px}}
.pill{{display:inline-flex;width:max-content;border-radius:999px;padding:5px 9px;font-size:.72rem;font-weight:700}}.pill.pass{{background:#e4f4ed;color:#176b50}}.pill.warn{{background:#fff1d6;color:#9a5a00}}.pill.fail{{background:#fde4e7;color:#a62f42}}
code{{overflow-wrap:anywhere}}@media(max-width:850px){{.cards,.plots{{grid-template-columns:1fr}}}}
</style></head><body><main>
<header><div class='sub'>Anchor-aware long-read library QC</div><h1>{escape(summary['sample_name'])}</h1>
<p class='sub'><code>{escape(summary['fastq_path'])}</code></p><p>Overall assessment: {_status_pill(overall)}</p></header>
<section><h2>QC verdicts</h2><div class='cards'>{verdict_cards}</div></section>
<section><h2>Overview</h2><div class='cards'>{overview_cards}</div></section>
<section><h2>Run configuration</h2><div class='scroll'><table><tbody>{config_rows}</tbody></table></div></section>
<section><h2>Yield, quality and position</h2><div class='plots'>{plot_html}</div></section>
<section><h2>Anchor performance</h2><div class='scroll'><table><thead><tr><th>Anchor</th><th>Detected</th><th>Mean edits</th><th>Insertions</th><th>Deletions</th><th>Multi-hit</th></tr></thead><tbody>{anchor_rows}</tbody></table></div></section>
{segment_section}{barcode_section}
<section><h2>Structure and failure modes</h2><div class='plots'><div class='panel scroll'><table><thead><tr><th>QC bucket</th><th>Reads</th><th>Ratio</th></tr></thead><tbody>{bucket_rows}</tbody></table></div><div class='panel scroll'><table><thead><tr><th>Structure</th><th>Orientation</th><th>Reads</th></tr></thead><tbody>{structure_rows}</tbody></table></div></div></section>
</main></body></html>"""
    _write_text(path, html)


def run_qc(
    fastq_path: str,
    config: AppConfig,
    outdir: str,
    export_csv: bool = False,
    output_passed_fastq: bool = False,
    output_failed_fastq: bool = False,
    threads: int = 1,
    gzip_output: bool = False,
) -> dict[str, Any]:
    """Run the shared deterministic pipeline with serial or multiprocessing execution."""
    if threads < 1:
        raise ValueError("threads must be at least 1")
    from .parallel import run_qc_parallel

    return run_qc_parallel(
        fastq_path,
        config,
        outdir,
        export_csv=export_csv,
        output_passed_fastq=output_passed_fastq,
        output_failed_fastq=output_failed_fastq,
        threads=threads,
        gzip_output=gzip_output,
    )
