from __future__ import annotations

from dataclasses import dataclass

from .anchors import AnchorHit


@dataclass
class ReadStructure:
    label: str
    order_valid: bool
    detected_anchors: list[str]
    is_reversed: bool = False


def classify_structure(hits: list[AnchorHit], expected_order: list[str]) -> ReadStructure:
    detected = [hit.anchor_name for hit in hits]
    unique_in_order: list[str] = []
    for name in detected:
        if not unique_in_order or unique_in_order[-1] != name:
            unique_in_order.append(name)

    if not hits:
        return ReadStructure(label="no_anchor_detected", order_valid=False, detected_anchors=[])

    repeated = len(set(detected)) != len(detected)
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

    missing = [name for name in expected_order if name not in unique_in_order]
    if not missing:
        positions = [expected_order.index(name) for name in unique_in_order if name in expected_order]
        if positions == sorted(positions):
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
    forward_hits: list[AnchorHit], reverse_hits: list[AnchorHit], expected_order: list[str]
) -> tuple[ReadStructure, list[AnchorHit]]:
    forward = classify_structure(forward_hits, expected_order)
    reverse = classify_structure(reverse_hits, expected_order)
    reverse.is_reversed = True
    if _structure_score(reverse, reverse_hits) > _structure_score(forward, forward_hits):
        return reverse, reverse_hits
    return forward, forward_hits
