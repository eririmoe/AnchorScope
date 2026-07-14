from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def nested_value(document: dict[str, Any], dotted_path: str) -> Any:
    value: Any = document
    for key in dotted_path.split("."):
        if not isinstance(value, dict) or key not in value:
            raise KeyError(dotted_path)
        value = value[key]
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def values_match(observed: Any, expected: Any, tolerance: float) -> bool:
    if isinstance(expected, bool):
        return observed is expected
    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        return math.isclose(float(observed), float(expected), rel_tol=0.0, abs_tol=tolerance)
    return observed == expected


def validate(
    run_dir: Path,
    fastq: Path,
    provenance_path: Path,
    expected_path: Path,
    tolerance: float = 1e-9,
) -> list[str]:
    expected = load_json(expected_path)
    summary = load_json(run_dir / "summary.json")
    provenance = load_json(provenance_path)
    errors: list[str] = []

    expected_digest = str(expected["subset_sha256"])
    observed_digest = sha256(fastq)
    if observed_digest != expected_digest:
        errors.append(f"FASTQ SHA-256: expected {expected_digest}, observed {observed_digest}")
    if provenance.get("subset_sha256") != expected_digest:
        errors.append("provenance.json does not contain the expected subset SHA-256")
    if provenance.get("written_reads") != summary.get("total_reads"):
        errors.append("provenance read count and summary total_reads differ")

    for dotted_path, expected_value in expected["summary_fields"].items():
        try:
            observed_value = nested_value(summary, dotted_path)
        except KeyError:
            errors.append(f"missing summary field: {dotted_path}")
            continue
        if not values_match(observed_value, expected_value, tolerance):
            errors.append(
                f"{dotted_path}: expected {expected_value!r}, observed {observed_value!r}"
            )

    total_reads = summary.get("total_reads")
    qc_counts = summary.get("qc_bucket_counts", {})
    if isinstance(total_reads, int) and isinstance(qc_counts, dict):
        if sum(qc_counts.values()) != total_reads:
            errors.append("qc_bucket_counts do not sum to total_reads")
    barcode_counts = summary.get("barcode_qc", {}).get("status_counts", {})
    if isinstance(total_reads, int) and isinstance(barcode_counts, dict):
        if sum(barcode_counts.values()) != total_reads:
            errors.append("barcode status counts do not sum to total_reads")

    for name in expected["required_outputs"]:
        path = run_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"missing or empty output: {name}")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate the reproducible public 10x ONT AnchorScope example."
    )
    parser.add_argument("--run-dir", type=Path, default=HERE / "anchorscope_report")
    parser.add_argument(
        "--fastq",
        type=Path,
        default=HERE / "SC5pv2_GEX_Human_Lung_Carcinoma_DTC.first_10000.fastq",
    )
    parser.add_argument("--provenance", type=Path, default=HERE / "provenance.json")
    parser.add_argument("--expected", type=Path, default=HERE / "expected_results.json")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()

    errors = validate(
        args.run_dir,
        args.fastq,
        args.provenance,
        args.expected,
        tolerance=args.tolerance,
    )
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        raise SystemExit(1)
    print("PASS: public 10x ONT example matches all expected results")


if __name__ == "__main__":
    main()
