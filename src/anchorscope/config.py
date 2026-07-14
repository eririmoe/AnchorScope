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
    # When set, fixed anchors are aligned with Levenshtein edit distance and
    # therefore tolerate substitutions, insertions, and deletions.  The legacy
    # max_mismatches option remains a Hamming-distance fast path.
    max_edits: int | None = None
    search_region: str = "any"  # any, 5p, or 3p
    search_window_bp: int | None = None


@dataclass
class SampleEntry:
    sample_name: str
    path: str


@dataclass
class AnchorRule:
    required: bool = True
    min_count: int = 1
    max_count: int | None = 1
    terminal: str = "any"  # any, 5p, 3p, or internal
    max_terminal_offset: float | None = None
    min_distance_from_previous: int | None = None
    max_distance_from_previous: int | None = None


@dataclass
class SegmentConfig:
    name: str
    start_anchor: str
    end_anchor: str
    role: str = "insert"
    min_length: int | None = None
    max_length: int | None = None


@dataclass
class StructureConfig:
    expected_order: list[str] = field(default_factory=list)
    anchor_rules: dict[str, AnchorRule] = field(default_factory=dict)
    segments: list[SegmentConfig] = field(default_factory=list)
    expected_orientation: str = "either"  # either, forward, or reverse
    concatemer_min_cycles: int = 2
    split_concatemers: bool = False


@dataclass
class BarcodeConfig:
    enabled: bool = False
    left_anchor: str | None = None
    right_anchor: str | None = None
    barcode_length: int = 16
    umi_length: int = 12
    offset_from_left: int = 0
    min_base_quality: int = 15
    whitelist_path: str | None = None
    max_edit_distance: int = 2
    min_distance_margin: int = 1


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
    barcode: BarcodeConfig = field(default_factory=BarcodeConfig)
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)
    qscore_method: str = "conservative"  # "conservative" or "arithmetic_mean"


def _validate_config(config: AppConfig) -> None:
    names = [anchor.name for anchor in config.anchors]
    if len(names) != len(set(names)):
        raise ValueError("Anchor names must be unique")
    known = set(names)
    for anchor in config.anchors:
        if anchor.type not in {"fixed", "regex"}:
            raise ValueError(f"Unsupported anchor type: {anchor.type}")
        if anchor.search_region not in {"any", "5p", "3p"}:
            raise ValueError(f"Anchor '{anchor.name}' search_region must be any, 5p, or 3p")
        if anchor.search_window_bp is not None and anchor.search_window_bp <= 0:
            raise ValueError(f"Anchor '{anchor.name}' search_window_bp must be positive")
        if anchor.max_edits is not None and (
            isinstance(anchor.max_edits, bool)
            or not isinstance(anchor.max_edits, int)
            or anchor.max_edits < 0
        ):
            raise ValueError(f"Anchor '{anchor.name}' max_edits must be a non-negative integer")

    order = config.structure.expected_order
    unknown_order = [name for name in order if name not in known]
    if unknown_order:
        raise ValueError(f"Unknown anchor(s) in expected_order: {', '.join(unknown_order)}")
    if len(order) != len(set(order)):
        raise ValueError("expected_order must not contain duplicate anchor names")
    if config.structure.expected_orientation not in {"either", "forward", "reverse"}:
        raise ValueError("expected_orientation must be either, forward, or reverse")
    if config.structure.concatemer_min_cycles < 2:
        raise ValueError("concatemer_min_cycles must be at least 2")

    for name, rule in config.structure.anchor_rules.items():
        if name not in known:
            raise ValueError(f"Unknown anchor in anchor_rules: {name}")
        if rule.terminal not in {"any", "5p", "3p", "internal"}:
            raise ValueError(f"Anchor rule '{name}' has invalid terminal value")
        if rule.min_count < 0 or (rule.max_count is not None and rule.max_count < rule.min_count):
            raise ValueError(f"Anchor rule '{name}' has invalid count bounds")
        if rule.max_terminal_offset is not None and not 0 <= rule.max_terminal_offset <= 1:
            raise ValueError(f"Anchor rule '{name}' max_terminal_offset must be between 0 and 1")
        if (
            rule.min_distance_from_previous is not None
            and rule.max_distance_from_previous is not None
            and rule.min_distance_from_previous > rule.max_distance_from_previous
        ):
            raise ValueError(f"Anchor rule '{name}' has invalid distance bounds")

    segment_names: set[str] = set()
    for segment in config.structure.segments:
        if segment.name in segment_names:
            raise ValueError(f"Duplicate segment name: {segment.name}")
        segment_names.add(segment.name)
        if segment.start_anchor not in known or segment.end_anchor not in known:
            raise ValueError(f"Segment '{segment.name}' references an unknown anchor")
        if segment.min_length is not None and segment.min_length < 0:
            raise ValueError(f"Segment '{segment.name}' min_length must be non-negative")
        if segment.max_length is not None and segment.max_length < 0:
            raise ValueError(f"Segment '{segment.name}' max_length must be non-negative")
        if (
            segment.min_length is not None
            and segment.max_length is not None
            and segment.min_length > segment.max_length
        ):
            raise ValueError(f"Segment '{segment.name}' has invalid length bounds")

    barcode = config.barcode
    if barcode.enabled:
        if not barcode.left_anchor or barcode.left_anchor not in known:
            raise ValueError("barcode.left_anchor must reference a configured anchor")
        if barcode.right_anchor is not None and barcode.right_anchor not in known:
            raise ValueError("barcode.right_anchor must reference a configured anchor")
        if barcode.barcode_length <= 0 or barcode.umi_length < 0:
            raise ValueError("barcode and UMI lengths must be positive/non-negative")
        if barcode.max_edit_distance < 0 or barcode.min_distance_margin < 0:
            raise ValueError("barcode edit-distance settings must be non-negative")
        if barcode.whitelist_path and not Path(barcode.whitelist_path).is_file():
            raise ValueError(f"Barcode whitelist not found: {barcode.whitelist_path}")

    t = config.thresholds
    fractional = {
        key: value
        for key, value in vars(t).items()
        if key.endswith("_ratio") or key in {"terminal_anchor_max_offset", "high_n_fraction"}
    }
    for key, value in fractional.items():
        if not 0 <= value <= 1:
            raise ValueError(f"thresholds.{key} must be between 0 and 1")
    if t.fail_long_high_quality_ratio > t.warn_long_high_quality_ratio:
        raise ValueError("Long high-quality fail threshold must not exceed warn threshold")
    if t.fail_correct_anchor_order_ratio > t.warn_correct_anchor_order_ratio:
        raise ValueError("Anchor-order fail threshold must not exceed warn threshold")
    for metric in ("no_anchor", "reversed_read", "high_n"):
        if getattr(t, f"warn_{metric}_ratio") > getattr(t, f"fail_{metric}_ratio"):
            raise ValueError(f"{metric} warn threshold must not exceed fail threshold")


def load_config(path: str | Path) -> AppConfig:
    config_path = Path(path)
    raw: dict[str, Any] = json.loads(config_path.read_text(encoding="utf-8"))
    anchors = [AnchorConfig(**anchor) for anchor in raw.get("anchors", [])]
    structure_raw = dict(raw.get("structure", {}))
    rules = {
        name: AnchorRule(**values)
        for name, values in structure_raw.pop("anchor_rules", {}).items()
    }
    segments = [SegmentConfig(**values) for values in structure_raw.pop("segments", [])]
    structure = StructureConfig(anchor_rules=rules, segments=segments, **structure_raw)
    barcode_raw = dict(raw.get("barcode", {}))
    whitelist_path = barcode_raw.get("whitelist_path")
    if whitelist_path:
        candidate = Path(whitelist_path)
        if not candidate.is_absolute():
            barcode_raw["whitelist_path"] = str((config_path.parent / candidate).resolve())
    barcode = BarcodeConfig(**barcode_raw)
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
    config = AppConfig(
        sample_name=raw.get("sample_name", "sample"),
        samples=samples,
        export_non_full_structure_fastq=raw.get("export_non_full_structure_fastq", False),
        anchors=anchors,
        structure=structure,
        barcode=barcode,
        thresholds=thresholds,
        qscore_method=qscore_method,
    )
    _validate_config(config)
    return config
