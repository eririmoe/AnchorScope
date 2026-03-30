from __future__ import annotations

import gzip
from dataclasses import dataclass
from math import log10
from pathlib import Path
from typing import Iterator


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
    def read_qscore(self) -> float:
        if not self.quality:
            return 0.0
        
        if self.qscore_method == "arithmetic_mean":
            # Traditional method: arithmetic mean of Phred scores
            return sum(self.phred_scores) / len(self.phred_scores)
        else:
            # Conservative method: average error rate then convert back to Qscore
            mean_error_rate = sum(10 ** (-q / 10) for q in self.phred_scores) / len(self.quality)
            if mean_error_rate <= 0:
                return 0.0
            return -10 * log10(mean_error_rate)

    @property
    def mean_q(self) -> float:
        # Preserve the historical API semantics: mean_q is the arithmetic
        # mean of per-base Phred scores regardless of read-level Qscore mode.
        if not self.quality:
            return 0.0
        return sum(self.phred_scores) / len(self.phred_scores)


def open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fastq(path: str | Path, qscore_method: str = "conservative") -> Iterator[FastqRead]:
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record for {header.strip()}")
            yield FastqRead(name=header.strip()[1:], sequence=sequence, quality=quality, qscore_method=qscore_method)
