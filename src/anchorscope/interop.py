from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Any

from .anchors import levenshtein_distance
from .multiqc import write_tag_audit_multiqc


_BARCODE_SUFFIX = re.compile(r"-\d+$")


def _normalized_barcode(value: str) -> str:
    return _BARCODE_SUFFIX.sub("", value.upper())


@dataclass
class TagAuditAccumulator:
    total_reads: int = 0
    tag_counts: Counter[str] = field(default_factory=Counter)
    barcode_edit_distances: Counter[int] = field(default_factory=Counter)
    umi_edit_distances: Counter[int] = field(default_factory=Counter)
    poly_a_lengths: Counter[int] = field(default_factory=Counter)
    ambiguous_barcodes: int = 0

    def update(self, tags: Mapping[str, Any]) -> None:
        self.total_reads += 1
        for tag in ("CR", "CB", "CY", "UR", "UB", "UY", "pt"):
            if tag in tags:
                self.tag_counts[tag] += 1
        raw_barcode = tags.get("CR")
        corrected_barcode = tags.get("CB")
        if raw_barcode and corrected_barcode:
            self.barcode_edit_distances[
                levenshtein_distance(
                    _normalized_barcode(str(raw_barcode)),
                    _normalized_barcode(str(corrected_barcode)),
                )
            ] += 1
        raw_umi = tags.get("UR")
        corrected_umi = tags.get("UB")
        if raw_umi and corrected_umi:
            self.umi_edit_distances[
                levenshtein_distance(str(raw_umi).upper(), str(corrected_umi).upper())
            ] += 1
        if tags.get("XA") == "ambiguous" or tags.get("barcode_status") == "ambiguous":
            self.ambiguous_barcodes += 1
        if "pt" in tags:
            try:
                length = int(tags["pt"])
            except (TypeError, ValueError):
                pass
            else:
                if length >= 0:
                    self.poly_a_lengths[length] += 1

    @staticmethod
    def _median(counter: Counter[int]) -> float | None:
        total = sum(counter.values())
        if not total:
            return None
        target = (total - 1) // 2
        cumulative = 0
        for value in sorted(counter):
            cumulative += counter[value]
            if cumulative > target:
                return float(value)
        return None

    def summary(self, sample_name: str, input_path: str) -> dict[str, Any]:
        total = self.total_reads
        return {
            "sample_name": sample_name,
            "input_path": input_path,
            "total_reads": total,
            "tag_counts": dict(self.tag_counts),
            "raw_barcode_ratio": self.tag_counts.get("CR", 0) / total if total else 0.0,
            "corrected_barcode_ratio": self.tag_counts.get("CB", 0) / total if total else 0.0,
            "corrected_umi_ratio": self.tag_counts.get("UB", 0) / total if total else 0.0,
            "ambiguous_barcode_ratio": self.ambiguous_barcodes / total if total else 0.0,
            "poly_a_tag_ratio": self.tag_counts.get("pt", 0) / total if total else 0.0,
            "barcode_edit_distance_counts": dict(sorted(self.barcode_edit_distances.items())),
            "umi_edit_distance_counts": dict(sorted(self.umi_edit_distances.items())),
            "median_poly_a_length": self._median(self.poly_a_lengths),
            "poly_a_length_counts": dict(sorted(self.poly_a_lengths.items())),
        }


def summarize_tag_records(
    records: Iterable[Mapping[str, Any]], sample_name: str, input_path: str = ""
) -> dict[str, Any]:
    accumulator = TagAuditAccumulator()
    for tags in records:
        accumulator.update(tags)
    return accumulator.summary(sample_name, input_path)


def _parse_sam_tags(fields: list[str]) -> dict[str, Any]:
    tags: dict[str, Any] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) != 3:
            continue
        tag, value_type, value = parts
        if value_type in {"i", "I"}:
            try:
                tags[tag] = int(value)
            except ValueError:
                continue
        elif value_type == "f":
            try:
                tags[tag] = float(value)
            except ValueError:
                continue
        else:
            tags[tag] = value
    return tags


def _write_audit_outputs(
    accumulator: TagAuditAccumulator,
    input_path: Path,
    output_path: Path,
    sample_name: str | None,
) -> dict[str, Any]:
    summary = accumulator.summary(sample_name or input_path.stem, str(input_path))
    (output_path / "tag_audit_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_tag_audit_multiqc(summary, output_path / "anchorscope_tag_audit_mqc.json")
    return summary


def audit_bam(path: str, outdir: str, sample_name: str | None = None) -> dict[str, Any]:
    input_path = Path(path)
    output_path = Path(outdir)
    output_path.mkdir(parents=True, exist_ok=True)
    accumulator = TagAuditAccumulator()
    if input_path.suffix.lower() == ".sam":
        with input_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line or line.startswith("@"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) >= 11:
                    accumulator.update(_parse_sam_tags(fields[11:]))
        return _write_audit_outputs(
            accumulator, input_path, output_path, sample_name
        )
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Binary BAM auditing requires the optional dependency: pip install 'AnchorScope[bam]'. "
            "Text SAM input is supported without optional dependencies."
        ) from exc
    with pysam.AlignmentFile(str(input_path), "rb", check_sq=False) as handle:
        for record in handle.fetch(until_eof=True):
            accumulator.update(dict(record.get_tags()))
    return _write_audit_outputs(accumulator, input_path, output_path, sample_name)
