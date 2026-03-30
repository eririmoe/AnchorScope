from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class AnchorConfig:
    name: str
    type: str
    sequence: str | None = None
    pattern: str | None = None
    max_mismatches: int = 0


@dataclass
class StructureConfig:
    expected_order: list[str] = field(default_factory=list)


@dataclass
class ThresholdConfig:
    long_read_min_bp: int = 1000
    long_read_min_q: float = 10.0
    heatmap_max_reads: int = 200


@dataclass
class AppConfig:
    sample_name: str = "sample"
    anchors: list[AnchorConfig] = field(default_factory=list)
    structure: StructureConfig = field(default_factory=StructureConfig)
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)
    qscore_method: str = "conservative"  # "conservative" or "arithmetic_mean"


def load_config(path: str | Path) -> AppConfig:
    raw: dict[str, Any] = json.loads(Path(path).read_text(encoding="utf-8"))
    anchors = [AnchorConfig(**anchor) for anchor in raw.get("anchors", [])]
    structure = StructureConfig(**raw.get("structure", {}))
    thresholds = ThresholdConfig(**raw.get("thresholds", {}))
    qscore_method = raw.get("qscore_method", "conservative")
    if qscore_method not in ("conservative", "arithmetic_mean"):
        raise ValueError(f"qscore_method must be 'conservative' or 'arithmetic_mean', got '{qscore_method}'")
    return AppConfig(
        sample_name=raw.get("sample_name", "sample"),
        anchors=anchors,
        structure=structure,
        thresholds=thresholds,
        qscore_method=qscore_method,
    )
