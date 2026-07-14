from __future__ import annotations

import argparse
import json
import platform
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from anchorscope.anchors import (
    find_anchor_hits,
    get_rust_status,
    levenshtein_distance,
    prepare_anchors,
    reverse_complement,
)
from anchorscope.classify import classify_best_orientation
from anchorscope.config import AnchorConfig, StructureConfig


# Adapter-scale anchors avoid the poor specificity of unrealistically short
# motifs while retaining deterministic single-edit truth cases.
LEFT = "CGACATGGCTACGATCCGACTT"
RIGHT = "ATGTACTCTGCGTTGATACCACTGCTT"
CLASSES = (
    "full_structure",
    "substitution",
    "insertion",
    "deletion",
    "missing_5p_anchor",
    "anchor_order_invalid",
    "concatemer_candidate",
    "reverse_full_structure",
)


def _random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _contains_near_anchor(sequence: str) -> bool:
    motifs = (LEFT, RIGHT, reverse_complement(LEFT), reverse_complement(RIGHT))
    for motif in motifs:
        for size in (len(motif) - 1, len(motif), len(motif) + 1):
            for start in range(0, len(sequence) - size + 1):
                if levenshtein_distance(sequence[start : start + size], motif, 1) <= 1:
                    return True
    return False


def _clean_random_sequence(rng: random.Random, length: int) -> str:
    while True:
        candidate = _random_sequence(rng, length)
        if not _contains_near_anchor(candidate):
            return candidate


def _mutate_substitution(sequence: str) -> str:
    replacement = "A" if sequence[3] != "A" else "C"
    return sequence[:3] + replacement + sequence[4:]


def generate_truth(reads_per_class: int, seed: int) -> list[tuple[str, str]]:
    rng = random.Random(seed)
    reads: list[tuple[str, str]] = []
    expected_labels = {
        "substitution": "full_structure",
        "insertion": "full_structure",
        "deletion": "full_structure",
        "reverse_full_structure": "full_structure",
    }
    for class_name in CLASSES:
        for _ in range(reads_per_class):
            payload = _clean_random_sequence(rng, 80)
            unit = LEFT + payload + RIGHT
            if class_name == "full_structure":
                sequence = unit
            elif class_name == "substitution":
                sequence = _mutate_substitution(LEFT) + payload + RIGHT
            elif class_name == "insertion":
                sequence = LEFT[:3] + "A" + LEFT[3:] + payload + RIGHT
            elif class_name == "deletion":
                sequence = LEFT[:3] + LEFT[4:] + payload + RIGHT
            elif class_name == "missing_5p_anchor":
                # A deterministic low-complexity negative region avoids an
                # accidental one-edit LEFT motif in the generated truth set.
                sequence = "G" * len(payload) + RIGHT
            elif class_name == "anchor_order_invalid":
                sequence = RIGHT + payload + LEFT
            elif class_name == "concatemer_candidate":
                sequence = unit + unit
            else:
                sequence = reverse_complement(unit)
            reads.append((sequence, expected_labels.get(class_name, class_name)))
    rng.shuffle(reads)
    return reads


def _evaluate(
    reads: list[tuple[str, str]], indel_aware: bool
) -> tuple[dict[str, Any], float]:
    if indel_aware:
        anchors = [
            AnchorConfig("left", "fixed", sequence=LEFT, max_edits=1),
            AnchorConfig("right", "fixed", sequence=RIGHT, max_edits=1),
        ]
    else:
        anchors = [
            AnchorConfig("left", "fixed", sequence=LEFT, max_mismatches=1),
            AnchorConfig("right", "fixed", sequence=RIGHT, max_mismatches=1),
        ]
    prepared = prepare_anchors(anchors)
    structure = StructureConfig(expected_order=["left", "right"])
    confusion: dict[str, Counter[str]] = defaultdict(Counter)
    correct = 0
    start_time = time.perf_counter()
    for sequence, truth in reads:
        forward = find_anchor_hits(sequence, prepared)
        reverse = find_anchor_hits(reverse_complement(sequence), prepared)
        predicted, _ = classify_best_orientation(
            forward, reverse, structure, read_length=len(sequence)
        )
        confusion[truth][predicted.label] += 1
        correct += predicted.label == truth
    elapsed = time.perf_counter() - start_time
    total = len(reads)
    metrics = {
        "total_reads": total,
        "correct_reads": correct,
        "accuracy": correct / total if total else 0.0,
        "reads_per_second": total / elapsed if elapsed else 0.0,
        "confusion_matrix": {
            truth: dict(predictions) for truth, predictions in sorted(confusion.items())
        },
    }
    return metrics, elapsed


def run_benchmark(reads_per_class: int = 100, seed: int = 20260710) -> dict[str, Any]:
    reads = generate_truth(reads_per_class, seed)
    edit_metrics, _ = _evaluate(reads, indel_aware=True)
    hamming_metrics, _ = _evaluate(reads, indel_aware=False)
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        git_dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
        git_dirty = None
    return {
        "benchmark": "AnchorScope synthetic structure truth set",
        "seed": seed,
        "reads_per_class": reads_per_class,
        "classes": list(CLASSES),
        "git_commit": commit,
        "git_dirty": git_dirty,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "rust_backend": get_rust_status(),
        "indel_aware": edit_metrics,
        "legacy_hamming": hamming_metrics,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads-per-class", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260710)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    result = run_benchmark(args.reads_per_class, args.seed)
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
