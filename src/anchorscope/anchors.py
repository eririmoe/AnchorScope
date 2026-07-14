from __future__ import annotations

import re
from dataclasses import dataclass

from .config import AnchorConfig


@dataclass(frozen=True)
class AnchorHit:
    anchor_name: str
    start: int
    end: int
    mismatches: int
    matched_sequence: str
    substitutions: int = 0
    insertions: int = 0
    deletions: int = 0
    cigar: str = ""


@dataclass(frozen=True)
class PreparedAnchor:
    config: AnchorConfig
    sequence: str | None = None
    pattern: re.Pattern[str] | None = None
    segments: tuple[tuple[int, str], ...] = ()


_RC_TRANSLATION = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def get_matcher_status(anchors: list[AnchorConfig]) -> dict[str, str]:
    edit_names = [
        anchor.name
        for anchor in anchors
        if anchor.type == "fixed" and anchor.max_edits is not None
    ]
    hamming_names = [
        anchor.name
        for anchor in anchors
        if anchor.type == "fixed" and anchor.max_edits is None
    ]
    if edit_names and hamming_names:
        return {
            "mode": "python-mixed",
            "detail": (
                f"Python Hamming matching: {', '.join(hamming_names)}; "
                f"Python Levenshtein matching: {', '.join(edit_names)}"
            ),
        }
    if edit_names:
        return {
            "mode": "python-edit",
            "detail": "Python Levenshtein matching: " + ", ".join(edit_names),
        }
    if hamming_names:
        return {
            "mode": "python-hamming",
            "detail": "Python Hamming matching: " + ", ".join(hamming_names),
        }
    return {"mode": "regex-only", "detail": "No fixed-sequence anchors configured"}


def _validated_max_mismatches(anchor: AnchorConfig) -> int:
    value = anchor.max_mismatches
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Anchor '{anchor.name}' max_mismatches must be an integer")
    if value < 0:
        raise ValueError(f"Anchor '{anchor.name}' max_mismatches must be non-negative")
    return value


def _split_segments(motif: str, max_mismatches: int) -> tuple[tuple[int, str], ...]:
    segment_count = max_mismatches + 1
    base_size, remainder = divmod(len(motif), segment_count)
    segments: list[tuple[int, str]] = []
    start = 0
    for idx in range(segment_count):
        seg_len = base_size + (1 if idx < remainder else 0)
        if seg_len <= 0:
            continue
        end = start + seg_len
        segments.append((start, motif[start:end]))
        start = end
    return tuple(segments)


def prepare_anchors(anchors: list[AnchorConfig]) -> list[PreparedAnchor]:
    prepared: list[PreparedAnchor] = []
    for anchor in anchors:
        if anchor.type == "fixed":
            if not anchor.sequence:
                raise ValueError(f"Fixed anchor '{anchor.name}' requires a non-empty sequence")
            max_mismatches = _validated_max_mismatches(anchor)
            motif = anchor.sequence.upper()
            prepared.append(
                PreparedAnchor(
                    config=anchor,
                    sequence=motif,
                    segments=_split_segments(motif, max_mismatches),
                )
            )
        elif anchor.type == "regex":
            if not anchor.pattern:
                raise ValueError(f"Regex anchor '{anchor.name}' requires a non-empty pattern")
            prepared.append(
                PreparedAnchor(
                    config=anchor,
                    pattern=re.compile(anchor.pattern),
                )
            )
        else:
            raise ValueError(f"Unsupported anchor type: {anchor.type}")
    return prepared


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_RC_TRANSLATION)[::-1]


def _hamming_distance_bounded(left: str, right: str, max_mismatches: int) -> int | None:
    mismatches = 0
    for lch, rch in zip(left, right):
        if lch != rch:
            mismatches += 1
            if mismatches > max_mismatches:
                return None
    return mismatches


def levenshtein_distance(left: str, right: str, max_distance: int | None = None) -> int:
    """Return Levenshtein distance, optionally stopping once a row exceeds a bound."""
    if left == right:
        return 0
    if not left:
        return len(right)
    if not right:
        return len(left)
    if max_distance is not None and abs(len(left) - len(right)) > max_distance:
        return max_distance + 1
    if len(left) > len(right):
        left, right = right, left
    previous = list(range(len(left) + 1))
    for row_index, rch in enumerate(right, start=1):
        current = [row_index]
        row_min = row_index
        for column_index, lch in enumerate(left, start=1):
            value = min(
                current[-1] + 1,
                previous[column_index] + 1,
                previous[column_index - 1] + (lch != rch),
            )
            current.append(value)
            row_min = min(row_min, value)
        if max_distance is not None and row_min > max_distance:
            return max_distance + 1
        previous = current
    return previous[-1]


def _compress_cigar(operations: list[str]) -> str:
    if not operations:
        return ""
    parts: list[str] = []
    current = operations[0]
    count = 1
    for operation in operations[1:]:
        if operation == current:
            count += 1
        else:
            parts.append(f"{count}{current}")
            current = operation
            count = 1
    parts.append(f"{count}{current}")
    return "".join(parts)


def _traceback_substring(
    matrix: list[list[int]], motif: str, sequence: str, end: int
) -> tuple[int, int, int, int, str]:
    i = len(motif)
    j = end
    substitutions = insertions = deletions = 0
    operations: list[str] = []
    while i > 0:
        if j > 0:
            cost = 0 if motif[i - 1] == sequence[j - 1] else 1
            if matrix[i][j] == matrix[i - 1][j - 1] + cost:
                operations.append("M" if cost == 0 else "X")
                substitutions += cost
                i -= 1
                j -= 1
                continue
        if j > 0 and matrix[i][j] == matrix[i][j - 1] + 1:
            operations.append("I")
            insertions += 1
            j -= 1
            continue
        operations.append("D")
        deletions += 1
        i -= 1
    operations.reverse()
    return j, substitutions, insertions, deletions, _compress_cigar(operations)


def _python_find_edit_hits(
    sequence: str, anchor_name: str, motif: str, max_edits: int
) -> list[AnchorHit]:
    """Find local motif alignments using semi-global Levenshtein distance."""
    motif_len = len(motif)
    seq_len = len(sequence)
    matrix = [[0] * (seq_len + 1) for _ in range(motif_len + 1)]
    for i in range(1, motif_len + 1):
        matrix[i][0] = i
    for i in range(1, motif_len + 1):
        mch = motif[i - 1]
        for j in range(1, seq_len + 1):
            matrix[i][j] = min(
                matrix[i - 1][j] + 1,
                matrix[i][j - 1] + 1,
                matrix[i - 1][j - 1] + (mch != sequence[j - 1]),
            )

    candidates: list[AnchorHit] = []
    for end in range(1, seq_len + 1):
        distance = matrix[motif_len][end]
        if distance > max_edits:
            continue
        start, substitutions, insertions, deletions, cigar = _traceback_substring(
            matrix, motif, sequence, end
        )
        if end <= start:
            continue
        candidates.append(
            AnchorHit(
                anchor_name=anchor_name,
                start=start,
                end=end,
                mismatches=distance,
                matched_sequence=sequence[start:end],
                substitutions=substitutions,
                insertions=insertions,
                deletions=deletions,
                cigar=cigar,
            )
        )

    # A single biological occurrence can generate several overlapping dynamic-
    # programming endpoints. Keep the best representative while preserving
    # distinct non-overlapping motif occurrences.
    selected: list[AnchorHit] = []
    for candidate in sorted(
        candidates,
        key=lambda hit: (hit.mismatches, abs((hit.end - hit.start) - motif_len), hit.start, hit.end),
    ):
        if any(candidate.start < hit.end and hit.start < candidate.end for hit in selected):
            continue
        selected.append(candidate)
    return sorted(selected, key=lambda hit: (hit.start, hit.end, hit.mismatches))


def _search_bounds(sequence_length: int, anchor: PreparedAnchor) -> tuple[int, int]:
    window = anchor.config.search_window_bp
    if anchor.config.search_region == "5p":
        return 0, min(sequence_length, window or sequence_length)
    if anchor.config.search_region == "3p":
        size = min(sequence_length, window or sequence_length)
        return sequence_length - size, sequence_length
    return 0, sequence_length


def _offset_hits(hits: list[AnchorHit], offset: int) -> list[AnchorHit]:
    if offset == 0:
        return hits
    return [
        AnchorHit(
            anchor_name=hit.anchor_name,
            start=hit.start + offset,
            end=hit.end + offset,
            mismatches=hit.mismatches,
            matched_sequence=hit.matched_sequence,
            substitutions=hit.substitutions,
            insertions=hit.insertions,
            deletions=hit.deletions,
            cigar=hit.cigar,
        )
        for hit in hits
    ]


def _python_find_exact_hits(sequence: str, anchor_name: str, motif: str) -> list[AnchorHit]:
    hits: list[AnchorHit] = []
    start = sequence.find(motif)
    while start != -1:
        hits.append(
            AnchorHit(
                anchor_name=anchor_name,
                start=start,
                end=start + len(motif),
                mismatches=0,
                matched_sequence=motif,
            )
        )
        start = sequence.find(motif, start + 1)
    return hits


def _python_find_approximate_hits(sequence: str, anchor_name: str, motif: str, max_mismatches: int, segments: tuple[tuple[int, str], ...]) -> list[AnchorHit]:
    if not segments:
        return []
    size = len(motif)
    seq_len = len(sequence)
    candidate_starts: set[int] = set()
    for segment_offset, segment in segments:
        start = sequence.find(segment)
        while start != -1:
            candidate_start = start - segment_offset
            if 0 <= candidate_start <= seq_len - size:
                candidate_starts.add(candidate_start)
            start = sequence.find(segment, start + 1)
    hits: list[AnchorHit] = []
    for candidate_start in sorted(candidate_starts):
        window = sequence[candidate_start : candidate_start + size]
        mismatches = _hamming_distance_bounded(window, motif, max_mismatches)
        if mismatches is not None:
            hits.append(
                AnchorHit(
                    anchor_name=anchor_name,
                    start=candidate_start,
                    end=candidate_start + size,
                    mismatches=mismatches,
                    matched_sequence=window,
                )
            )
    return hits


def find_fixed_hits(sequence: str, anchor: PreparedAnchor) -> list[AnchorHit]:
    if anchor.sequence is None:
        raise ValueError(f"Fixed anchor '{anchor.config.name}' is missing prepared sequence")
    motif = anchor.sequence
    region_start, region_end = _search_bounds(len(sequence), anchor)
    search_sequence = sequence[region_start:region_end]
    if not motif or not search_sequence:
        return []
    if anchor.config.max_edits is not None:
        hits = _python_find_edit_hits(
            search_sequence, anchor.config.name, motif, anchor.config.max_edits
        )
        return _offset_hits(hits, region_start)
    if len(search_sequence) < len(motif):
        return []
    if anchor.config.max_mismatches <= 0:
        hits = _python_find_exact_hits(search_sequence, anchor.config.name, motif)
    else:
        hits = _python_find_approximate_hits(search_sequence, anchor.config.name, motif, anchor.config.max_mismatches, anchor.segments)
    return _offset_hits(hits, region_start)


def find_regex_hits(sequence: str, anchor: PreparedAnchor) -> list[AnchorHit]:
    if anchor.pattern is None:
        raise ValueError(f"Regex anchor '{anchor.config.name}' is missing compiled pattern")
    region_start, region_end = _search_bounds(len(sequence), anchor)
    search_sequence = sequence[region_start:region_end]
    hits: list[AnchorHit] = []
    for match in anchor.pattern.finditer(search_sequence):
        hits.append(
            AnchorHit(
                anchor_name=anchor.config.name,
                start=match.start() + region_start,
                end=match.end() + region_start,
                mismatches=0,
                matched_sequence=match.group(0),
            )
        )
    return hits


def find_anchor_hits(sequence: str, anchors: list[PreparedAnchor]) -> list[AnchorHit]:
    seq = sequence.upper()
    all_hits: list[AnchorHit] = []
    for anchor in anchors:
        if anchor.config.type == "fixed":
            all_hits.extend(find_fixed_hits(seq, anchor))
        elif anchor.config.type == "regex":
            all_hits.extend(find_regex_hits(seq, anchor))
        else:
            raise ValueError(f"Unsupported anchor type: {anchor.config.type}")
    return sorted(all_hits, key=lambda hit: (hit.start, hit.end, hit.anchor_name))
