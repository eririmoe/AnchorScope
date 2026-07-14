from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from math import log10

from .anchors import AnchorHit
from .config import SegmentConfig


@dataclass(frozen=True)
class SegmentObservation:
    name: str
    role: str
    start: int | None
    end: int | None
    length: int | None
    mean_qscore: float | None
    gc_fraction: float | None
    n_fraction: float | None
    status: str


def _mean_qscore(quality: str) -> float:
    if not quality:
        return 0.0
    error_rate = sum(10 ** (-(ord(ch) - 33) / 10) for ch in quality) / len(quality)
    return -10 * log10(error_rate) if error_rate > 0 else 0.0


def observe_segments(
    sequence: str,
    quality: str,
    best_hits: dict[str, AnchorHit],
    segments: list[SegmentConfig],
) -> dict[str, SegmentObservation]:
    observations: dict[str, SegmentObservation] = {}
    for segment in segments:
        left = best_hits.get(segment.start_anchor)
        right = best_hits.get(segment.end_anchor)
        if left is None or right is None or right.start < left.end:
            observations[segment.name] = SegmentObservation(
                name=segment.name,
                role=segment.role,
                start=None,
                end=None,
                length=None,
                mean_qscore=None,
                gc_fraction=None,
                n_fraction=None,
                status="missing",
            )
            continue
        start, end = left.end, right.start
        segment_sequence = sequence[start:end].upper()
        segment_quality = quality[start:end]
        length = len(segment_sequence)
        status = "pass"
        if segment.min_length is not None and length < segment.min_length:
            status = "too_short"
        elif segment.max_length is not None and length > segment.max_length:
            status = "too_long"
        gc_count = segment_sequence.count("G") + segment_sequence.count("C")
        n_count = segment_sequence.count("N")
        observations[segment.name] = SegmentObservation(
            name=segment.name,
            role=segment.role,
            start=start,
            end=end,
            length=length,
            mean_qscore=_mean_qscore(segment_quality),
            gc_fraction=(gc_count / length) if length else 0.0,
            n_fraction=(n_count / length) if length else 0.0,
            status=status,
        )
    return observations


def _counter_median(counter: Counter[int], scale: int = 1) -> float:
    total = sum(counter.values())
    if not total:
        return 0.0
    targets = {(total - 1) // 2, total // 2}
    cumulative = 0
    values: list[int] = []
    for value in sorted(counter):
        next_cumulative = cumulative + counter[value]
        for target in sorted(targets):
            if cumulative <= target < next_cumulative:
                values.append(value)
        cumulative = next_cumulative
    return (sum(values) / len(values) / scale) if values else 0.0


@dataclass
class SegmentAccumulator:
    total_reads: int = 0
    observed_reads: int = 0
    status_counts: Counter[str] = field(default_factory=Counter)
    length_counts: Counter[int] = field(default_factory=Counter)
    qscore_counts: Counter[int] = field(default_factory=Counter)
    gc_sum: float = 0.0
    n_sum: float = 0.0

    def update(self, observation: SegmentObservation) -> None:
        self.total_reads += 1
        self.status_counts[observation.status] += 1
        if observation.length is None:
            return
        self.observed_reads += 1
        self.length_counts[observation.length] += 1
        self.qscore_counts[int(round((observation.mean_qscore or 0.0) * 100))] += 1
        self.gc_sum += observation.gc_fraction or 0.0
        self.n_sum += observation.n_fraction or 0.0

    def summary(self) -> dict[str, object]:
        observed = self.observed_reads
        total = self.total_reads
        return {
            "total_reads": total,
            "observed_reads": observed,
            "observed_ratio": observed / total if total else 0.0,
            "missing_ratio": self.status_counts.get("missing", 0) / total if total else 0.0,
            "within_length_ratio": self.status_counts.get("pass", 0) / observed if observed else 0.0,
            "mean_length": (
                sum(length * count for length, count in self.length_counts.items()) / observed
                if observed
                else 0.0
            ),
            "median_length": _counter_median(self.length_counts),
            "median_qscore": _counter_median(self.qscore_counts, scale=100),
            "mean_gc_fraction": self.gc_sum / observed if observed else 0.0,
            "mean_n_fraction": self.n_sum / observed if observed else 0.0,
            "status_counts": dict(self.status_counts),
        }
