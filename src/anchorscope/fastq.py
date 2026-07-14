from __future__ import annotations

import gzip
from dataclasses import dataclass
from math import log10
from pathlib import Path
from typing import Iterator

VALID_BASES = frozenset("ACGTN")


@dataclass
class FastqRead:
    name: str
    sequence: str
    quality: str
    qscore_method: str = "conservative"  # "conservative" or "arithmetic_mean"

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def phred_scores(self) -> list[int]:
        return [ord(ch) - 33 for ch in self.quality]

    @property
    def normalized_sequence(self) -> str:
        return self.sequence.upper()

    @property
    def n_count(self) -> int:
        return self.normalized_sequence.count("N")

    @property
    def invalid_base_count(self) -> int:
        return sum(1 for base in self.normalized_sequence if base not in VALID_BASES)

    @property
    def read_qscore(self) -> float:
        if not self.quality:
            return 0.0

        if self.qscore_method == "arithmetic_mean":
            return sum(self.phred_scores) / len(self.phred_scores)

        mean_error_rate = sum(10 ** (-q / 10) for q in self.phred_scores) / len(self.quality)
        if mean_error_rate <= 0:
            return 0.0
        return -10 * log10(mean_error_rate)

    @property
    def mean_q(self) -> float:
        if not self.quality:
            return 0.0
        return sum(self.phred_scores) / len(self.phred_scores)

    @property
    def n_fraction(self) -> float:
        if not self.length:
            return 0.0
        return self.n_count / self.length

    @property
    def invalid_base_fraction(self) -> float:
        if not self.length:
            return 0.0
        return self.invalid_base_count / self.length


def open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def open_output_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("w", encoding="utf-8")


def fastq_output_path(outdir: str | Path, stem: str, gzip_output: bool) -> Path:
    suffix = ".fastq.gz" if gzip_output else ".fastq"
    return Path(outdir) / f"{stem}{suffix}"


def read_fastq(path: str | Path, qscore_method: str = "conservative") -> Iterator[FastqRead]:
    with open_text(path) as handle:
        record_index = 0
        while True:
            header = handle.readline()
            if not header:
                break
            record_index += 1
            if not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ header at record #{record_index}: {header.strip()!r}")
            sequence_line = handle.readline()
            plus = handle.readline()
            quality_line = handle.readline()
            if not sequence_line or not plus or not quality_line:
                raise ValueError(f"Malformed FASTQ record #{record_index}: incomplete 4-line record for {header.strip()}")
            sequence = sequence_line.rstrip("\r\n")
            quality = quality_line.rstrip("\r\n")
            if not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record for {header.strip()}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"Malformed FASTQ record #{record_index} for {header.strip()}: sequence and quality lengths differ"
                )
            yield FastqRead(name=header.strip()[1:], sequence=sequence, quality=quality, qscore_method=qscore_method)
