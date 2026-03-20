from __future__ import annotations

import ctypes
import re
import shutil
import subprocess
from threading import Lock
from dataclasses import dataclass
from pathlib import Path

from .config import AnchorConfig


@dataclass(frozen=True)
class AnchorHit:
    anchor_name: str
    start: int
    end: int
    mismatches: int
    matched_sequence: str


@dataclass(frozen=True)
class PreparedAnchor:
    config: AnchorConfig
    sequence: str | None = None
    pattern: re.Pattern[str] | None = None
    segments: tuple[tuple[int, str], ...] = ()


_RUST_ENGINE_LOCK = Lock()
_RUST_BUILD_ATTEMPTED = False
_RUST_ENGINE: RustAnchorEngine | None = None


class RustAnchorEngine:
    def __init__(self, library_path: Path):
        self._lib = ctypes.CDLL(str(library_path))
        self._lib.scfastq_find_fixed_hits.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_size_t]
        self._lib.scfastq_find_fixed_hits.restype = ctypes.c_void_p
        self._lib.scfastq_free_string.argtypes = [ctypes.c_void_p]
        self._lib.scfastq_free_string.restype = None

    def find_fixed_hits(self, sequence: str, motif: str, max_mismatches: int) -> list[tuple[int, int]]:
        raw_ptr = self._lib.scfastq_find_fixed_hits(sequence.encode("utf-8"), motif.encode("utf-8"), max_mismatches)
        if not raw_ptr:
            return []
        try:
            payload = ctypes.string_at(raw_ptr).decode("utf-8")
        finally:
            self._lib.scfastq_free_string(raw_ptr)
        if not payload:
            return []
        hits: list[tuple[int, int]] = []
        for item in payload.split(","):
            start, mismatches = item.split(":", 1)
            hits.append((int(start), int(mismatches)))
        return hits


def _rust_library_candidates(repo_root: Path) -> list[Path]:
    return [
        repo_root / "target" / "release" / "libscfastq_qc_rs.so",
        repo_root / "target" / "release" / "scfastq_qc_rs.dll",
        repo_root / "target" / "release" / "libscfastq_qc_rs.dylib",
        repo_root / "rust" / "anchor_engine" / "target" / "release" / "libscfastq_qc_rs.so",
        repo_root / "rust" / "anchor_engine" / "target" / "release" / "scfastq_qc_rs.dll",
        repo_root / "rust" / "anchor_engine" / "target" / "release" / "libscfastq_qc_rs.dylib",
    ]


def _find_rust_library(repo_root: Path) -> Path | None:
    for candidate in _rust_library_candidates(repo_root):
        if candidate.exists():
            return candidate
    return None


def _try_build_rust_library(repo_root: Path) -> None:
    manifest_path = repo_root / "rust" / "anchor_engine" / "Cargo.toml"
    if not manifest_path.exists() or shutil.which("cargo") is None:
        return
    subprocess.run(
        ["cargo", "build", "--release", "--manifest-path", str(manifest_path)],
        check=False,
        cwd=repo_root,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def get_rust_anchor_engine() -> RustAnchorEngine | None:
    global _RUST_BUILD_ATTEMPTED, _RUST_ENGINE
    if _RUST_ENGINE is not None:
        return _RUST_ENGINE
    repo_root = Path(__file__).resolve().parents[2]
    library_path = _find_rust_library(repo_root)
    if library_path is not None:
        _RUST_ENGINE = RustAnchorEngine(library_path)
        return _RUST_ENGINE
    if _RUST_BUILD_ATTEMPTED:
        return None
    with _RUST_ENGINE_LOCK:
        if _RUST_ENGINE is not None:
            return _RUST_ENGINE
        if _RUST_BUILD_ATTEMPTED:
            return None
        _RUST_BUILD_ATTEMPTED = True
        _try_build_rust_library(repo_root)
    library_path = _find_rust_library(repo_root)
    if library_path is not None:
        _RUST_ENGINE = RustAnchorEngine(library_path)
        return _RUST_ENGINE
    return None


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


def _hamming_distance_bounded(left: str, right: str, max_mismatches: int) -> int | None:
    mismatches = 0
    for lch, rch in zip(left, right):
        if lch != rch:
            mismatches += 1
            if mismatches > max_mismatches:
                return None
    return mismatches


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


def _rust_find_fixed_hits(sequence: str, anchor_name: str, motif: str, max_mismatches: int) -> list[AnchorHit] | None:
    if max_mismatches < 0:
        raise ValueError(f"Anchor '{anchor_name}' max_mismatches must be non-negative")
    engine = get_rust_anchor_engine()
    if engine is None:
        return None
    hits: list[AnchorHit] = []
    for start, mismatches in engine.find_fixed_hits(sequence, motif, max_mismatches):
        end = start + len(motif)
        hits.append(
            AnchorHit(
                anchor_name=anchor_name,
                start=start,
                end=end,
                mismatches=mismatches,
                matched_sequence=sequence[start:end],
            )
        )
    return hits


def find_fixed_hits(sequence: str, anchor: PreparedAnchor) -> list[AnchorHit]:
    if anchor.sequence is None:
        raise ValueError(f"Fixed anchor '{anchor.config.name}' is missing prepared sequence")
    motif = anchor.sequence
    if not motif or len(sequence) < len(motif):
        return []
    rust_hits = _rust_find_fixed_hits(sequence, anchor.config.name, motif, anchor.config.max_mismatches)
    if rust_hits is not None:
        return rust_hits
    if anchor.config.max_mismatches <= 0:
        return _python_find_exact_hits(sequence, anchor.config.name, motif)
    return _python_find_approximate_hits(sequence, anchor.config.name, motif, anchor.config.max_mismatches, anchor.segments)


def find_regex_hits(sequence: str, anchor: PreparedAnchor) -> list[AnchorHit]:
    if anchor.pattern is None:
        raise ValueError(f"Regex anchor '{anchor.config.name}' is missing compiled pattern")
    hits: list[AnchorHit] = []
    for match in anchor.pattern.finditer(sequence):
        hits.append(
            AnchorHit(
                anchor_name=anchor.config.name,
                start=match.start(),
                end=match.end(),
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
