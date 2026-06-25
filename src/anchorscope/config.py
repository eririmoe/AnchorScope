from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


@dataclass
class AnchorConfig:
    name: str
    type: str
    sequence: str | None = None
    pattern: str | None = None
    max_mismatches: int = 0


@dataclass
class SampleEntry:
    sample_name: str
    path: str


@dataclass
class StructureConfig:
    expected_order: list[str] = field(default_factory=list)


@dataclass
class ThresholdConfig:
    long_read_min_bp: int = 1000
    long_read_min_q: float = 10.0
    terminal_anchor_max_offset: float = 0.15
    high_n_fraction: float = 0.1
    warn_long_high_quality_ratio: float = 0.7
    fail_long_high_quality_ratio: float = 0.5
    warn_correct_anchor_order_ratio: float = 0.7
    fail_correct_anchor_order_ratio: float = 0.5
    warn_no_anchor_ratio: float = 0.2
    fail_no_anchor_ratio: float = 0.35
    warn_reversed_read_ratio: float = 0.3
    fail_reversed_read_ratio: float = 0.5
    warn_high_n_ratio: float = 0.1
    fail_high_n_ratio: float = 0.2


@dataclass
class AppConfig:
    sample_name: str = "sample"
    samples: list[SampleEntry] | None = None
    export_non_full_structure_fastq: bool = False
    anchors: list[AnchorConfig] = field(default_factory=list)
    structure: StructureConfig = field(default_factory=StructureConfig)
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)
    qscore_method: str = "conservative"  # "conservative" or "arithmetic_mean"


def load_config(path: str | Path) -> AppConfig:
    raw: dict[str, Any] = json.loads(Path(path).read_text(encoding="utf-8"))
    anchors = [AnchorConfig(**anchor) for anchor in raw.get("anchors", [])]
    structure = StructureConfig(**raw.get("structure", {}))
    threshold_values = dict(raw.get("thresholds", {}))
    legacy_threshold_keys = {
        "heatmap_max_reads",
        "end_proximity_bp",
        "end_proximity_fraction",
        "max_n_fraction",
    }
    ignored_legacy_keys = sorted(key for key in legacy_threshold_keys if key in threshold_values)
    for key in ignored_legacy_keys:
        threshold_values.pop(key, None)
    if ignored_legacy_keys:
        logger.warning(
            "Ignoring legacy threshold setting(s) that are no longer used by the report pipeline: %s",
            ", ".join(ignored_legacy_keys),
        )
    thresholds = ThresholdConfig(**threshold_values)
    qscore_method = raw.get("qscore_method", "conservative")
    if qscore_method not in ("conservative", "arithmetic_mean"):
        raise ValueError(f"qscore_method must be 'conservative' or 'arithmetic_mean', got '{qscore_method}'")
    samples: list[SampleEntry] | None = None
    if "samples" in raw:
        samples = [SampleEntry(**entry) for entry in raw["samples"]]
    return AppConfig(
        sample_name=raw.get("sample_name", "sample"),
        samples=samples,
        export_non_full_structure_fastq=raw.get("export_non_full_structure_fastq", False),
        anchors=anchors,
        structure=structure,
        thresholds=thresholds,
        qscore_method=qscore_method,
    )
