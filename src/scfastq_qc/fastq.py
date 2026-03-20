from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass
class FastqRead:
    name: str
    sequence: str
    quality: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def mean_q(self) -> float:
        if not self.quality:
            return 0.0
        return sum(ord(ch) - 33 for ch in self.quality) / len(self.quality)


def open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fastq(path: str | Path) -> Iterator[FastqRead]:
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\n")
            if not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record for {header.strip()}")
            yield FastqRead(name=header.strip()[1:], sequence=sequence, quality=quality)
