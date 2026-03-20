from __future__ import annotations

from dataclasses import dataclass

from .anchors import AnchorHit


@dataclass
class ReadStructure:
    label: str
    order_valid: bool
    detected_anchors: list[str]


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
