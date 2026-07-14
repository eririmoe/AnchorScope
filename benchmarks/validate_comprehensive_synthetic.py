#!/usr/bin/env python3
"""End-to-end synthetic truth validation for AnchorScope.

Writes deterministic FASTQ files, runs the complete reporting pipeline, and
compares every read against truth for structure, orientation, QC bucket, edit
type, segment status, and barcode/UMI recoverability.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from anchorscope.anchors import find_anchor_hits, prepare_anchors, reverse_complement
from anchorscope.config import AnchorConfig, load_config
from anchorscope.fastq import read_fastq
from anchorscope.report import run_qc

LEFT = "AACCGGTTACGTCAGT"
RIGHT = "TGCATGCAAGCTTCGA"
BARCODE_WHITELIST = ("AACCGG", "AACCGA", "TTGGCC")
HIGH_Q = "I"
LOW_Q = "&"  # Phred 5


@dataclass(frozen=True)
class TruthRead:
    read_id: str
    suite: str
    scenario: str
    sequence: str
    quality: str
    expected_structure: str
    expected_bucket: str
    expected_orientation: str = "forward"
    required_flags: tuple[str, ...] = ()
    expected_segment_status: str | None = None
    expected_barcode_status: str | None = None
    expected_raw_barcode: str | None = None
    expected_raw_umi: str | None = None
    edit_anchor: str | None = None
    substitutions: int | None = None
    insertions: int | None = None
    deletions: int | None = None


@dataclass
class CheckBook:
    checks: list[dict[str, Any]] = field(default_factory=list)

    def add(self, name: str, passed: bool, detail: str = "") -> None:
        self.checks.append({"name": name, "passed": passed, "detail": detail})

    @property
    def passed(self) -> bool:
        return all(item["passed"] for item in self.checks)

    @property
    def passed_count(self) -> int:
        return sum(1 for item in self.checks if item["passed"])


def _substitute(sequence: str, positions: tuple[int, ...]) -> str:
    values = list(sequence)
    alternatives = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in positions:
        values[position] = alternatives[values[position]]
    return "".join(values)


def _safe_payload(rng: random.Random, length: int) -> str:
    """Generate DNA without an accidental anchor hit in either orientation."""
    prepared = prepare_anchors(
        [
            AnchorConfig("left", "fixed", sequence=LEFT, max_edits=1),
            AnchorConfig("right", "fixed", sequence=RIGHT, max_edits=1),
        ]
    )
    for _ in range(10_000):
        sequence = "".join(rng.choice("ACGT") for _ in range(length))
        if not find_anchor_hits(sequence, prepared) and not find_anchor_hits(
            reverse_complement(sequence), prepared
        ):
            return sequence
    raise RuntimeError("Unable to generate an anchor-free payload")


def _truth(
    suite: str,
    scenario: str,
    replicate: int,
    sequence: str,
    expected_structure: str,
    expected_bucket: str,
    *,
    quality: str | None = None,
    orientation: str = "forward",
    flags: tuple[str, ...] = (),
    segment: str | None = None,
    barcode_status: str | None = None,
    raw_barcode: str | None = None,
    raw_umi: str | None = None,
    edit_anchor: str | None = None,
    substitutions: int | None = None,
    insertions: int | None = None,
    deletions: int | None = None,
) -> TruthRead:
    return TruthRead(
        read_id=f"{suite}__{scenario}__{replicate:03d}",
        suite=suite,
        scenario=scenario,
        sequence=sequence,
        quality=quality if quality is not None else HIGH_Q * len(sequence),
        expected_structure=expected_structure,
        expected_bucket=expected_bucket,
        expected_orientation=orientation,
        required_flags=flags,
        expected_segment_status=segment,
        expected_barcode_status=barcode_status,
        expected_raw_barcode=raw_barcode,
        expected_raw_umi=raw_umi,
        edit_anchor=edit_anchor,
        substitutions=substitutions,
        insertions=insertions,
        deletions=deletions,
    )


def build_structure_truth(replicates: int, seed: int) -> list[TruthRead]:
    rng = random.Random(seed)
    reads: list[TruthRead] = []
    for replicate in range(1, replicates + 1):
        payload = _safe_payload(rng, 60)
        short_payload = _safe_payload(rng, 38)
        very_short_payload = _safe_payload(rng, 20)
        long_payload = _safe_payload(rng, 90)
        prefix = _safe_payload(rng, 10)
        suffix = _safe_payload(rng, 10)
        cycle = LEFT + payload + RIGHT
        reads.extend(
            [
                _truth("structure", "exact_forward", replicate, cycle, "full_structure", "pass", segment="pass"),
                _truth("structure", "exact_reverse", replicate, reverse_complement(cycle), "full_structure", "pass", orientation="reversed", segment="pass"),
                _truth("structure", "left_substitution", replicate, _substitute(LEFT, (5,)) + payload + RIGHT, "full_structure", "pass", segment="pass", edit_anchor="left", substitutions=1, insertions=0, deletions=0),
                _truth("structure", "right_substitution", replicate, LEFT + payload + _substitute(RIGHT, (7,)), "full_structure", "pass", segment="pass", edit_anchor="right", substitutions=1, insertions=0, deletions=0),
                _truth("structure", "left_insertion", replicate, LEFT[:6] + "A" + LEFT[6:] + payload + RIGHT, "full_structure", "pass", segment="pass", edit_anchor="left", substitutions=0, insertions=1, deletions=0),
                _truth("structure", "right_insertion", replicate, LEFT + payload + RIGHT[:8] + "C" + RIGHT[8:], "full_structure", "pass", segment="pass", edit_anchor="right", substitutions=0, insertions=1, deletions=0),
                _truth("structure", "left_deletion", replicate, LEFT[:6] + LEFT[7:] + payload + RIGHT, "full_structure", "pass", segment="pass", edit_anchor="left", substitutions=0, insertions=0, deletions=1),
                _truth("structure", "right_deletion", replicate, LEFT + payload + RIGHT[:8] + RIGHT[9:], "full_structure", "pass", segment="pass", edit_anchor="right", substitutions=0, insertions=0, deletions=1),
                _truth("structure", "over_edit_limit", replicate, _substitute(LEFT, (3, 10)) + payload + RIGHT, "missing_5p_anchor", "missing_5p", flags=("missing_5p",), segment="missing"),
                _truth("structure", "missing_5p", replicate, payload + RIGHT, "missing_5p_anchor", "missing_5p", flags=("missing_5p",), segment="missing"),
                _truth("structure", "missing_3p", replicate, LEFT + payload, "missing_3p_anchor", "missing_3p", flags=("missing_3p",), segment="missing"),
                _truth("structure", "no_anchor", replicate, _safe_payload(rng, 92), "no_anchor_detected", "no_anchor", flags=("no_anchor",), segment="missing"),
                _truth("structure", "order_invalid", replicate, RIGHT + payload + LEFT, "anchor_order_invalid", "order_invalid", flags=("order_invalid",), segment="missing"),
                _truth("structure", "duplicated_5p", replicate, LEFT + payload + LEFT + payload + RIGHT, "internal_5p_anchor", "internal_adapter", flags=("internal_adapter",)),
                _truth("structure", "duplicated_3p", replicate, LEFT + payload + RIGHT + payload + RIGHT, "internal_3p_anchor", "internal_adapter", flags=("internal_adapter",), segment="pass"),
                _truth("structure", "concatemer_two_cycles", replicate, cycle * 2, "concatemer_candidate", "concatemer_candidate", flags=("concatemer_candidate",), segment="pass"),
                _truth("structure", "concatemer_three_cycles", replicate, cycle * 3, "concatemer_candidate", "concatemer_candidate", flags=("concatemer_candidate", "dense_anchor_hits"), segment="pass"),
                _truth("structure", "five_prime_offset", replicate, prefix + cycle, "full_structure", "five_prime_truncated", flags=("five_prime_truncated",), segment="pass"),
                _truth("structure", "three_prime_offset", replicate, cycle + suffix, "full_structure", "three_prime_truncated", flags=("three_prime_truncated",), segment="pass"),
                _truth("structure", "segment_too_short", replicate, LEFT + short_payload + RIGHT, "full_structure", "segment_too_short", flags=("segment_too_short",), segment="too_short"),
                _truth("structure", "segment_too_long", replicate, LEFT + long_payload + RIGHT, "full_structure", "segment_too_long", flags=("segment_too_long",), segment="too_long"),
                _truth("structure", "low_read_quality", replicate, cycle, "full_structure", "low_read_q", quality=LOW_Q * len(cycle), flags=("low_read_q",), segment="pass"),
                _truth("structure", "high_n_fraction", replicate, LEFT + ("N" * 20) + payload[20:] + RIGHT, "full_structure", "high_n", flags=("high_n",), segment="pass"),
                _truth("structure", "invalid_base", replicate, LEFT + "X" + payload[1:] + RIGHT, "full_structure", "invalid_bases", flags=("invalid_bases",), segment="pass"),
                _truth("structure", "read_too_short", replicate, LEFT + very_short_payload + RIGHT, "full_structure", "too_short", flags=("too_short", "segment_too_short"), segment="too_short"),
                _truth("structure", "lowercase_input", replicate, cycle.lower(), "full_structure", "pass", segment="pass"),
                _truth("structure", "empty_read", replicate, "", "no_anchor_detected", "empty_read", quality="", flags=("empty_read", "too_short", "low_read_q", "no_anchor"), segment="missing"),
            ]
        )
    return reads


def build_barcode_truth(replicates: int) -> list[TruthRead]:
    reads: list[TruthRead] = []
    umi = "TGCA"
    cases = (
        ("exact", "AACCGG", "exact"),
        ("uniquely_recoverable", "AACCTG", "uniquely_recoverable"),
        ("ambiguous", "AACCGT", "ambiguous"),
        ("unrecoverable", "GGGGGG", "unrecoverable"),
    )
    for replicate in range(1, replicates + 1):
        for scenario, barcode, status in cases:
            sequence = LEFT + barcode + umi + RIGHT
            reads.append(_truth("barcode", scenario, replicate, sequence, "full_structure", "pass", segment="pass", barcode_status=status, raw_barcode=barcode, raw_umi=umi))
        exact_sequence = LEFT + BARCODE_WHITELIST[0] + umi + RIGHT
        low_quality = list(HIGH_Q * len(exact_sequence))
        for index in range(len(LEFT), len(LEFT) + 6):
            low_quality[index] = "!"
        reads.append(_truth("barcode", "low_barcode_quality", replicate, exact_sequence, "full_structure", "pass", quality="".join(low_quality), segment="pass", barcode_status="low_quality", raw_barcode=BARCODE_WHITELIST[0], raw_umi=umi))
        incomplete = LEFT + "AACCG" + RIGHT
        reads.append(_truth("barcode", "region_incomplete", replicate, incomplete, "full_structure", "segment_too_short", flags=("segment_too_short",), segment="too_short", barcode_status="region_incomplete"))
        anchor_missing = LEFT + BARCODE_WHITELIST[0] + umi
        reads.append(_truth("barcode", "anchor_missing", replicate, anchor_missing, "missing_3p_anchor", "missing_3p", flags=("missing_3p",), segment="missing", barcode_status="anchor_missing"))
        reads.append(_truth("barcode", "reverse_exact", replicate, reverse_complement(exact_sequence), "full_structure", "pass", orientation="reversed", segment="pass", barcode_status="exact", raw_barcode=BARCODE_WHITELIST[0], raw_umi=umi))
    return reads


def _write_fastq(path: Path, reads: list[TruthRead]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="\n") as handle:
        for read in reads:
            handle.write(f"@{read.read_id}\n{read.sequence}\n+\n{read.quality}\n")


def _write_truth(path: Path, reads: list[TruthRead]) -> None:
    fields = ["read_id", "suite", "scenario", "sequence", "quality", "expected_structure", "expected_bucket", "expected_orientation", "required_flags", "expected_segment_status", "expected_barcode_status", "expected_raw_barcode", "expected_raw_umi", "edit_anchor", "substitutions", "insertions", "deletions"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for read in reads:
            row = {name: getattr(read, name) for name in fields}
            row["required_flags"] = ";".join(read.required_flags)
            writer.writerow(row)


def _write_configs(root: Path) -> tuple[Path, Path]:
    structure = {
        "sample_name": "synthetic_structure_truth",
        "export_non_full_structure_fastq": True,
        "anchors": [
            {"name": "left", "type": "fixed", "sequence": LEFT, "max_edits": 1},
            {"name": "right", "type": "fixed", "sequence": RIGHT, "max_edits": 1},
        ],
        "structure": {
            "expected_order": ["left", "right"], "expected_orientation": "either",
            "anchor_rules": {"left": {"required": True, "min_count": 1, "max_count": 1}, "right": {"required": True, "min_count": 1, "max_count": 1}},
            "segments": [{"name": "payload", "start_anchor": "left", "end_anchor": "right", "min_length": 40, "max_length": 80}],
            "concatemer_min_cycles": 2, "split_concatemers": True,
        },
        "barcode": {"enabled": False},
        "thresholds": {"long_read_min_bp": 70, "long_read_min_q": 12.0, "terminal_anchor_max_offset": 0.05, "high_n_fraction": 0.10},
    }
    barcode = {
        "sample_name": "synthetic_barcode_truth",
        "anchors": [
            {"name": "left", "type": "fixed", "sequence": LEFT, "max_edits": 0},
            {"name": "right", "type": "fixed", "sequence": RIGHT, "max_edits": 0},
        ],
        "structure": {
            "expected_order": ["left", "right"], "expected_orientation": "either",
            "anchor_rules": {"left": {"required": True, "min_count": 1, "max_count": 1}, "right": {"required": True, "min_count": 1, "max_count": 1}},
            "segments": [{"name": "barcode_umi", "start_anchor": "left", "end_anchor": "right", "min_length": 10, "max_length": 10}],
        },
        "barcode": {"enabled": True, "left_anchor": "left", "right_anchor": "right", "barcode_length": 6, "umi_length": 4, "min_base_quality": 15, "whitelist_path": "barcode_whitelist.txt", "max_edit_distance": 1, "min_distance_margin": 1},
        "thresholds": {"long_read_min_bp": 1, "long_read_min_q": 0, "terminal_anchor_max_offset": 0.05, "high_n_fraction": 0.10},
    }
    structure_path = root / "structure_config.json"
    barcode_path = root / "barcode_config.json"
    structure_path.write_text(json.dumps(structure, indent=2), encoding="utf-8")
    barcode_path.write_text(json.dumps(barcode, indent=2), encoding="utf-8")
    (root / "barcode_whitelist.txt").write_text("\n".join(BARCODE_WHITELIST) + "\n", encoding="utf-8")
    return structure_path, barcode_path


def _read_csv_by_id(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["Read ID"]: row for row in csv.DictReader(handle)}


def _read_anchor_hits(path: Path) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            grouped.setdefault(row["Read ID"], []).append(row)
    return grouped


def _segment_status(row: dict[str, str], name: str) -> str | None:
    for item in row["Segment Statuses"].split(";"):
        if item.startswith(name + ":"):
            return item.split(":", 1)[1]
    return None


def _validate_reads(reads: list[TruthRead], details_path: Path, hits_path: Path, segment_name: str, checks: CheckBook) -> int:
    details = _read_csv_by_id(details_path)
    hits = _read_anchor_hits(hits_path)
    matched = 0
    for truth in reads:
        row = details.get(truth.read_id)
        errors: list[str] = []
        if row is None:
            errors.append("missing read_details row")
        else:
            observed_flags = {value for value in row["QC Flags"].split(";") if value}
            expected = {"Structure Label": truth.expected_structure, "QC Bucket": truth.expected_bucket, "Is Reversed": str(truth.expected_orientation == "reversed")}
            for field_name, expected_value in expected.items():
                if row[field_name] != expected_value:
                    errors.append(f"{field_name}={row[field_name]!r}, expected {expected_value!r}")
            missing_flags = set(truth.required_flags) - observed_flags
            if missing_flags:
                errors.append(f"missing flags {sorted(missing_flags)}")
            if truth.expected_segment_status is not None:
                observed_segment = _segment_status(row, segment_name)
                if observed_segment != truth.expected_segment_status:
                    errors.append(f"segment={observed_segment!r}, expected {truth.expected_segment_status!r}")
            for field_name, expected_value in (("Barcode Status", truth.expected_barcode_status), ("Raw Barcode", truth.expected_raw_barcode), ("Raw UMI", truth.expected_raw_umi)):
                if expected_value is not None and row[field_name] != expected_value:
                    errors.append(f"{field_name}={row[field_name]!r}, expected {expected_value!r}")
        if truth.edit_anchor:
            candidates = [item for item in hits.get(truth.read_id, []) if item["Anchor Name"] == truth.edit_anchor]
            if not candidates:
                errors.append(f"missing {truth.edit_anchor} edit hit")
            else:
                best = min(candidates, key=lambda item: (int(item["Mismatches"]), int(item["Start"])))
                for field_name, expected_value in (("Substitutions", truth.substitutions), ("Insertions", truth.insertions), ("Deletions", truth.deletions)):
                    if expected_value is not None and int(best[field_name]) != expected_value:
                        errors.append(f"{field_name}={best[field_name]}, expected {expected_value}")
        passed = not errors
        matched += int(passed)
        checks.add(f"read:{truth.read_id}", passed, "; ".join(errors))
    return matched


def _validate_malformed_fastq(root: Path, checks: CheckBook) -> None:
    malformed = {
        "invalid_header": ("r1\nACGT\n+\nIIII\n", "header"),
        "invalid_plus": ("@r1\nACGT\n-\nIIII\n", "Malformed FASTQ record"),
        "length_mismatch": ("@r1\nACGT\n+\nIII\n", "lengths differ"),
        "incomplete_record": ("@r1\nACGT\n+\n", "incomplete 4-line record"),
    }
    malformed_dir = root / "malformed_fastq"
    malformed_dir.mkdir(parents=True, exist_ok=True)
    for name, (content, expected_message) in malformed.items():
        path = malformed_dir / f"{name}.fastq"
        path.write_text(content, encoding="utf-8")
        try:
            list(read_fastq(path))
        except ValueError as error:
            checks.add(f"malformed:{name}", expected_message in str(error), str(error))
        else:
            checks.add(f"malformed:{name}", False, "input was accepted")


def _write_report(root: Path, result: dict[str, Any], reads: list[TruthRead]) -> None:
    scenario_counts = Counter((read.suite, read.scenario) for read in reads)
    failed = [item for item in result["checks"] if not item["passed"]]
    lines = [
        "# AnchorScope comprehensive synthetic validation", "",
        f"- Overall: **{'PASS' if result['overall_passed'] else 'FAIL'}**",
        f"- Seed: `{result['seed']}`", f"- Replicates per scenario: `{result['replicates']}`",
        f"- Simulated FASTQ reads: `{result['simulated_reads']}`", f"- Scenario types: `{result['scenario_types']}`",
        f"- Assertions passed: `{result['checks_passed']}/{result['checks_total']}`",
        f"- Per-read truth accuracy: `{result['per_read_accuracy']:.3%}`", "",
        "## Coverage", "", "| Suite | Scenario | Reads |", "|---|---|---:|",
    ]
    for (suite, scenario), count in sorted(scenario_counts.items()):
        lines.append(f"| {suite} | {scenario} | {count} |")
    lines.extend([
        "", "## Cross-path checks", "",
        "The same structure truth set was processed as plain FASTQ in serial mode and as gzip FASTQ in multiprocessing mode. Read-level CSV outputs and selected aggregate summaries were required to match exactly. Passed/failed and concatemer-split gzip outputs were also checked.",
        "", "## Interpretation boundary", "",
        "A pass demonstrates correctness for the explicitly encoded synthetic truth cases and deterministic execution paths. It does not replace validation on independent experimental libraries, competing tools, or platform-specific error distributions.",
    ])
    if failed:
        lines.extend(["", "## Failures", ""])
        for item in failed:
            lines.append(f"- `{item['name']}`: {item['detail']}")
    (root / "validation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_validation(output_dir: str | Path, *, replicates: int = 5, seed: int = 20260713, threads: int = 2) -> dict[str, Any]:
    if replicates < 1 or threads < 1:
        raise ValueError("replicates and threads must be at least 1")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    checks = CheckBook()
    structure_reads = build_structure_truth(replicates, seed)
    barcode_reads = build_barcode_truth(replicates)
    all_reads = structure_reads + barcode_reads
    structure_fastq = root / "structure_cases.fastq"
    structure_fastq_gz = root / "structure_cases.fastq.gz"
    barcode_fastq = root / "barcode_cases.fastq"
    _write_fastq(structure_fastq, structure_reads)
    _write_fastq(structure_fastq_gz, structure_reads)
    _write_fastq(barcode_fastq, barcode_reads)
    _write_truth(root / "truth.tsv", all_reads)
    structure_config_path, barcode_config_path = _write_configs(root)
    structure_config = load_config(structure_config_path)
    barcode_config = load_config(barcode_config_path)
    serial_dir = root / "structure_serial"
    parallel_dir = root / "structure_parallel_gzip"
    barcode_dir = root / "barcode_serial"
    serial_summary = run_qc(str(structure_fastq), structure_config, str(serial_dir), export_csv=True, output_passed_fastq=True, output_failed_fastq=True, threads=1, gzip_output=True)
    parallel_summary = run_qc(str(structure_fastq_gz), structure_config, str(parallel_dir), export_csv=True, threads=threads)
    barcode_summary = run_qc(str(barcode_fastq), barcode_config, str(barcode_dir), export_csv=True, threads=1)
    structure_matched = _validate_reads(structure_reads, serial_dir / "read_details.csv", serial_dir / "anchor_hits.csv", "payload", checks)
    barcode_matched = _validate_reads(barcode_reads, barcode_dir / "read_details.csv", barcode_dir / "anchor_hits.csv", "barcode_umi", checks)
    for csv_name in ("read_details.csv", "anchor_hits.csv"):
        serial_text = (serial_dir / csv_name).read_text(encoding="utf-8")
        parallel_text = (parallel_dir / csv_name).read_text(encoding="utf-8")
        checks.add(f"serial_parallel:{csv_name}", serial_text == parallel_text, "CSV outputs differ")
    for key in ("total_reads", "structure_counts", "structure_orientation_counts", "qc_bucket_counts", "segment_qc", "correct_anchor_order_ratio", "split_molecules_exported"):
        checks.add(f"serial_parallel:summary:{key}", serial_summary[key] == parallel_summary[key], f"serial={serial_summary[key]!r}; parallel={parallel_summary[key]!r}")
    expected_pass = sum(read.expected_bucket == "pass" for read in structure_reads)
    expected_fail = len(structure_reads) - expected_pass
    checks.add("export:passed_count", serial_summary["passed_reads_exported"] == expected_pass)
    checks.add("export:failed_count", serial_summary["failed_reads_exported"] == expected_fail)
    checks.add("export:concatemer_split_count", serial_summary["split_molecules_exported"] == replicates * 5)
    for name in ("passed.fastq.gz", "failed.fastq.gz", "split_concatemers.fastq.gz"):
        path = serial_dir / name
        try:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                content = handle.read()
            checks.add(f"export:gzip:{name}", path.is_file() and bool(content))
        except (OSError, EOFError) as error:
            checks.add(f"export:gzip:{name}", False, str(error))
    expected_barcode_counts = Counter(read.expected_barcode_status for read in barcode_reads if read.expected_barcode_status)
    observed_barcode_counts = Counter(barcode_summary["barcode_qc"]["status_counts"])
    checks.add("barcode:aggregate_status_counts", expected_barcode_counts == observed_barcode_counts, f"expected={dict(expected_barcode_counts)}; observed={dict(observed_barcode_counts)}")
    with (barcode_dir / "read_details.csv").open(encoding="utf-8") as handle:
        detail_header = next(csv.reader(handle))
    checks.add("barcode:diagnostic_only_schema", "Corrected Barcode" not in detail_header and "Raw Barcode" in detail_header, "barcode validation must remain diagnostic-only")
    _validate_malformed_fastq(root, checks)
    matched = structure_matched + barcode_matched
    result: dict[str, Any] = {
        "overall_passed": checks.passed, "seed": seed, "replicates": replicates,
        "simulated_reads": len(all_reads), "scenario_types": len({(read.suite, read.scenario) for read in all_reads}),
        "per_read_matched": matched, "per_read_accuracy": matched / len(all_reads) if all_reads else 0.0,
        "checks_total": len(checks.checks), "checks_passed": checks.passed_count,
        "structure_summary": {"total_reads": serial_summary["total_reads"], "structure_counts": serial_summary["structure_counts"], "qc_bucket_counts": serial_summary["qc_bucket_counts"], "split_molecules_exported": serial_summary["split_molecules_exported"]},
        "barcode_summary": barcode_summary["barcode_qc"], "checks": checks.checks,
    }
    (root / "validation_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    _write_report(root, result, all_reads)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "validation" / "synthetic_comprehensive"))
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260713)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    result = run_validation(args.output_dir, replicates=args.replicates, seed=args.seed, threads=args.threads)
    print(json.dumps({"overall_passed": result["overall_passed"], "simulated_reads": result["simulated_reads"], "scenario_types": result["scenario_types"], "checks": f"{result['checks_passed']}/{result['checks_total']}", "per_read_accuracy": result["per_read_accuracy"], "report": str(Path(args.output_dir) / "validation_report.md")}, indent=2))
    return 0 if result["overall_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
