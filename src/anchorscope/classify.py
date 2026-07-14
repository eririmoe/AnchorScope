from __future__ import annotations

from dataclasses import dataclass

from .anchors import AnchorHit
from .config import StructureConfig


@dataclass
class ReadStructure:
    label: str
    order_valid: bool
    detected_anchors: list[str]
    is_reversed: bool = False
    violations: list[str] | None = None
    cycle_count: int = 0

    def __post_init__(self) -> None:
        if self.violations is None:
            self.violations = []


def _as_structure_config(value: list[str] | StructureConfig) -> StructureConfig:
    if isinstance(value, StructureConfig):
        return value
    return StructureConfig(expected_order=list(value))


def find_structure_cycles(
    hits: list[AnchorHit], expected_order: list[str]
) -> list[list[AnchorHit]]:
    """Return non-overlapping, contiguous complete anchor cycles."""
    if not expected_order:
        return []
    ordered = sorted(hits, key=lambda hit: (hit.start, hit.end, hit.anchor_name))
    cycles: list[list[AnchorHit]] = []
    current: list[AnchorHit] = []
    expected_index = 0
    for hit in ordered:
        expected_name = expected_order[expected_index]
        if hit.anchor_name == expected_name:
            current.append(hit)
            expected_index += 1
            if expected_index == len(expected_order):
                cycles.append(current)
                current = []
                expected_index = 0
            continue
        if hit.anchor_name == expected_order[0]:
            current = [hit]
            expected_index = 1
        elif hit.anchor_name in expected_order:
            current = []
            expected_index = 0
    return cycles


def _rule_violations(
    hits: list[AnchorHit], config: StructureConfig, read_length: int | None
) -> list[str]:
    violations: list[str] = []
    by_name: dict[str, list[AnchorHit]] = {}
    for hit in hits:
        by_name.setdefault(hit.anchor_name, []).append(hit)
    for name, rule in config.anchor_rules.items():
        anchor_hits = by_name.get(name, [])
        count = len(anchor_hits)
        if rule.required and count < rule.min_count:
            violations.append(f"{name}:below_min_count")
        if rule.max_count is not None and count > rule.max_count:
            violations.append(f"{name}:above_max_count")
        if not anchor_hits or not read_length or rule.terminal == "any":
            continue
        best = min(anchor_hits, key=lambda hit: (hit.mismatches, hit.start, hit.end))
        limit = rule.max_terminal_offset if rule.max_terminal_offset is not None else 0.15
        five_offset = best.start / read_length
        three_offset = (read_length - best.end) / read_length
        if rule.terminal == "5p" and five_offset > limit:
            violations.append(f"{name}:not_5p_terminal")
        elif rule.terminal == "3p" and three_offset > limit:
            violations.append(f"{name}:not_3p_terminal")
        elif rule.terminal == "internal" and (five_offset <= limit or three_offset <= limit):
            violations.append(f"{name}:not_internal")

    best_by_name = {
        name: min(values, key=lambda hit: (hit.mismatches, hit.start, hit.end))
        for name, values in by_name.items()
    }
    for index, name in enumerate(config.expected_order[1:], start=1):
        rule = config.anchor_rules.get(name)
        if rule is None:
            continue
        previous_name = config.expected_order[index - 1]
        current = best_by_name.get(name)
        previous = best_by_name.get(previous_name)
        if current is None or previous is None:
            continue
        distance = current.start - previous.end
        if rule.min_distance_from_previous is not None and distance < rule.min_distance_from_previous:
            violations.append(f"{name}:distance_below_min")
        if rule.max_distance_from_previous is not None and distance > rule.max_distance_from_previous:
            violations.append(f"{name}:distance_above_max")
    return violations


def classify_structure(
    hits: list[AnchorHit],
    expected_order: list[str] | StructureConfig,
    read_length: int | None = None,
) -> ReadStructure:
    structure_config = _as_structure_config(expected_order)
    expected_order = structure_config.expected_order
    detected = [hit.anchor_name for hit in hits]
    unique_in_order: list[str] = []
    for name in detected:
        if not unique_in_order or unique_in_order[-1] != name:
            unique_in_order.append(name)

    if not hits:
        return ReadStructure(label="no_anchor_detected", order_valid=False, detected_anchors=[])

    cycles = find_structure_cycles(hits, expected_order)
    if len(cycles) >= structure_config.concatemer_min_cycles:
        return ReadStructure(
            label="concatemer_candidate",
            order_valid=False,
            detected_anchors=unique_in_order,
            cycle_count=len(cycles),
        )

    repeated = len(set(detected)) != len(detected)
    if repeated:
        repeated_excess = any(
            structure_config.anchor_rules.get(name) is None
            or structure_config.anchor_rules[name].max_count is not None
            and detected.count(name) > structure_config.anchor_rules[name].max_count
            for name in set(detected)
            if detected.count(name) > 1
        )
        if not repeated_excess:
            repeated = False
    if repeated:
        if expected_order:
            expected_count = len(expected_order)
            for start in range(0, max(0, len(unique_in_order) - expected_count + 1)):
                if unique_in_order[start : start + expected_count] == expected_order:
                    trailing = unique_in_order[start + expected_count :]
                    if trailing and trailing[0] == expected_order[0]:
                        return ReadStructure(
                            label="concatemer_candidate",
                            order_valid=False,
                            detected_anchors=unique_in_order,
                        )
            if detected.count(expected_order[0]) > 1:
                return ReadStructure(label="internal_5p_anchor", order_valid=False, detected_anchors=unique_in_order)
            if detected.count(expected_order[-1]) > 1:
                return ReadStructure(label="internal_3p_anchor", order_valid=False, detected_anchors=unique_in_order)
        return ReadStructure(label="duplicated_anchor", order_valid=False, detected_anchors=unique_in_order)

    if not expected_order:
        return ReadStructure(label="anchors_detected", order_valid=True, detected_anchors=unique_in_order)

    missing = [
        name
        for name in expected_order
        if name not in unique_in_order
        and (name not in structure_config.anchor_rules or structure_config.anchor_rules[name].required)
    ]
    if not missing:
        positions = [expected_order.index(name) for name in unique_in_order if name in expected_order]
        if positions == sorted(positions):
            violations = _rule_violations(hits, structure_config, read_length)
            if violations:
                return ReadStructure(
                    label="structure_rule_violation",
                    order_valid=False,
                    detected_anchors=unique_in_order,
                    violations=violations,
                    cycle_count=len(cycles),
                )
            return ReadStructure(label="full_structure", order_valid=True, detected_anchors=unique_in_order)
        return ReadStructure(label="anchor_order_invalid", order_valid=False, detected_anchors=unique_in_order)

    if unique_in_order and expected_order[0] not in unique_in_order:
        return ReadStructure(label="missing_5p_anchor", order_valid=False, detected_anchors=unique_in_order)

    if unique_in_order and expected_order[-1] not in unique_in_order:
        return ReadStructure(label="missing_3p_anchor", order_valid=False, detected_anchors=unique_in_order)

    return ReadStructure(label="partial_structure", order_valid=False, detected_anchors=unique_in_order)


def _structure_score(structure: ReadStructure, hits: list[AnchorHit]) -> tuple[int, int, int, int, int]:
    label_priority = {
        "full_structure": 6,
        "anchors_detected": 5,
        "anchor_order_invalid": 4,
        "partial_structure": 3,
        "missing_5p_anchor": 2,
        "missing_3p_anchor": 2,
        "concatemer_candidate": 1,
        "internal_5p_anchor": 1,
        "internal_3p_anchor": 1,
        "duplicated_anchor": 1,
        "structure_rule_violation": 1,
        "no_anchor_detected": 0,
    }
    return (
        1 if structure.order_valid else 0,
        label_priority.get(structure.label, 0),
        len(structure.detected_anchors),
        len(hits),
        -sum(hit.mismatches for hit in hits),
    )


def classify_best_orientation(
    forward_hits: list[AnchorHit],
    reverse_hits: list[AnchorHit],
    expected_order: list[str] | StructureConfig,
    read_length: int | None = None,
) -> tuple[ReadStructure, list[AnchorHit]]:
    structure_config = _as_structure_config(expected_order)
    forward = classify_structure(forward_hits, structure_config, read_length=read_length)
    reverse = classify_structure(reverse_hits, structure_config, read_length=read_length)
    reverse.is_reversed = True
    if _structure_score(reverse, reverse_hits) > _structure_score(forward, forward_hits):
        selected, selected_hits = reverse, reverse_hits
    else:
        selected, selected_hits = forward, forward_hits
    expected_orientation = structure_config.expected_orientation
    orientation = "reverse" if selected.is_reversed else "forward"
    if expected_orientation != "either" and orientation != expected_orientation:
        selected.violations = [*selected.violations, "orientation_unexpected"]
        if selected.label == "full_structure":
            selected.label = "orientation_unexpected"
            selected.order_valid = False
    return selected, selected_hits
