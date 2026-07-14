from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .anchors import AnchorHit, levenshtein_distance
from .config import BarcodeConfig


class BKTree:
    """Compact metric tree for bounded edit-distance whitelist queries."""

    def __init__(self, values: list[str]):
        self.root: tuple[str, dict[int, object]] | None = None
        for value in values:
            self.add(value)

    def add(self, value: str) -> None:
        if self.root is None:
            self.root = (value, {})
            return
        node = self.root
        while True:
            node_value, children = node
            distance = levenshtein_distance(value, node_value)
            child = children.get(distance)
            if child is None:
                children[distance] = (value, {})
                return
            node = child  # type: ignore[assignment]

    def query(self, value: str, max_distance: int) -> list[tuple[int, str]]:
        if self.root is None:
            return []
        results: list[tuple[int, str]] = []
        stack = [self.root]
        while stack:
            node_value, children = stack.pop()
            distance = levenshtein_distance(value, node_value)
            if distance <= max_distance:
                results.append((distance, node_value))
            lower = distance - max_distance
            upper = distance + max_distance
            stack.extend(
                child
                for edge, child in children.items()
                if lower <= edge <= upper
            )
        return sorted(results, key=lambda item: (item[0], item[1]))


@dataclass(frozen=True)
class BarcodeObservation:
    status: str
    raw_barcode: str | None = None
    raw_umi: str | None = None
    min_barcode_q: int | None = None
    nearest_distance: int | None = None
    second_distance: int | None = None
    nearest_barcode: str | None = None


@dataclass
class BarcodeMatcher:
    values: set[str]
    tree: BKTree

    @classmethod
    def from_path(cls, path: str) -> "BarcodeMatcher":
        values = {
            line.strip().upper()
            for line in Path(path).read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        }
        if not values:
            raise ValueError(f"Barcode whitelist is empty: {path}")
        return cls(values=values, tree=BKTree(sorted(values)))


def prepare_barcode_matcher(config: BarcodeConfig) -> BarcodeMatcher | None:
    if not config.enabled or not config.whitelist_path:
        return None
    return BarcodeMatcher.from_path(config.whitelist_path)


def evaluate_barcode_recoverability(
    sequence: str,
    quality: str,
    best_hits: dict[str, AnchorHit],
    config: BarcodeConfig,
    matcher: BarcodeMatcher | None,
) -> BarcodeObservation:
    if not config.enabled or not config.left_anchor:
        return BarcodeObservation(status="disabled")
    left = best_hits.get(config.left_anchor)
    right = best_hits.get(config.right_anchor) if config.right_anchor else None
    if left is None or (config.right_anchor and right is None):
        return BarcodeObservation(status="anchor_missing")
    start = left.end + config.offset_from_left
    total_length = config.barcode_length + config.umi_length
    end = start + total_length
    if start < 0 or end > len(sequence) or (right is not None and end > right.start):
        return BarcodeObservation(status="region_incomplete")
    raw_barcode = sequence[start : start + config.barcode_length].upper()
    raw_umi = sequence[start + config.barcode_length : end].upper()
    barcode_quality = quality[start : start + config.barcode_length]
    if len(barcode_quality) != config.barcode_length:
        return BarcodeObservation(status="region_incomplete")
    min_q = min((ord(ch) - 33 for ch in barcode_quality), default=0)
    if min_q < config.min_base_quality:
        return BarcodeObservation(
            status="low_quality",
            raw_barcode=raw_barcode,
            raw_umi=raw_umi,
            min_barcode_q=min_q,
        )
    if matcher is None:
        return BarcodeObservation(
            status="extracted",
            raw_barcode=raw_barcode,
            raw_umi=raw_umi,
            min_barcode_q=min_q,
        )
    if raw_barcode in matcher.values:
        return BarcodeObservation(
            status="exact",
            raw_barcode=raw_barcode,
            raw_umi=raw_umi,
            min_barcode_q=min_q,
            nearest_distance=0,
            nearest_barcode=raw_barcode,
        )
    candidates = matcher.tree.query(
        raw_barcode, config.max_edit_distance + config.min_distance_margin
    )
    if not candidates or candidates[0][0] > config.max_edit_distance:
        return BarcodeObservation(
            status="unrecoverable",
            raw_barcode=raw_barcode,
            raw_umi=raw_umi,
            min_barcode_q=min_q,
            nearest_distance=candidates[0][0] if candidates else None,
        )
    nearest_distance, nearest = candidates[0]
    second_distance = candidates[1][0] if len(candidates) > 1 else None
    if second_distance is not None and second_distance - nearest_distance < config.min_distance_margin:
        status = "ambiguous"
    else:
        status = "uniquely_recoverable"
    return BarcodeObservation(
        status=status,
        raw_barcode=raw_barcode,
        raw_umi=raw_umi,
        min_barcode_q=min_q,
        nearest_distance=nearest_distance,
        second_distance=second_distance,
        nearest_barcode=nearest,
    )


@dataclass
class BarcodeAccumulator:
    total_reads: int = 0
    status_counts: Counter[str] = field(default_factory=Counter)
    barcode_counts: Counter[str] = field(default_factory=Counter)
    umi_counts: Counter[str] = field(default_factory=Counter)
    edit_distance_counts: Counter[int] = field(default_factory=Counter)
    min_q_counts: Counter[int] = field(default_factory=Counter)

    def update(self, observation: BarcodeObservation) -> None:
        self.total_reads += 1
        self.status_counts[observation.status] += 1
        if observation.raw_barcode:
            self.barcode_counts[observation.raw_barcode] += 1
        if observation.raw_umi:
            self.umi_counts[observation.raw_umi] += 1
        if observation.nearest_distance is not None:
            self.edit_distance_counts[observation.nearest_distance] += 1
        if observation.min_barcode_q is not None:
            self.min_q_counts[observation.min_barcode_q] += 1

    def summary(self) -> dict[str, object]:
        total = self.total_reads
        observed = sum(
            count
            for status, count in self.status_counts.items()
            if status not in {"disabled", "anchor_missing", "region_incomplete"}
        )
        recoverable = sum(
            self.status_counts.get(status, 0)
            for status in ("exact", "uniquely_recoverable", "extracted")
        )
        return {
            "total_reads": total,
            "observed_reads": observed,
            "observed_ratio": observed / total if total else 0.0,
            "recoverable_ratio": recoverable / total if total else 0.0,
            "ambiguous_ratio": self.status_counts.get("ambiguous", 0) / total if total else 0.0,
            "low_quality_ratio": self.status_counts.get("low_quality", 0) / total if total else 0.0,
            "status_counts": dict(self.status_counts),
            "edit_distance_counts": dict(sorted(self.edit_distance_counts.items())),
            "unique_raw_barcodes": len(self.barcode_counts),
            "unique_raw_umis": len(self.umi_counts),
            "top_raw_barcodes": self.barcode_counts.most_common(50),
            "umi_singleton_ratio": (
                sum(1 for count in self.umi_counts.values() if count == 1) / len(self.umi_counts)
                if self.umi_counts
                else 0.0
            ),
        }
