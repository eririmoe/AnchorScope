from __future__ import annotations

from dataclasses import dataclass

from .anchors import AnchorHit
from .classify import find_structure_cycles


@dataclass(frozen=True)
class SplitMolecule:
    parent_read_id: str
    cycle_index: int
    sequence: str
    quality: str
    start: int
    end: int

    @property
    def read_id(self) -> str:
        return f"{self.parent_read_id}|cycle={self.cycle_index}"


def split_complete_cycles(
    read_id: str,
    oriented_sequence: str,
    oriented_quality: str,
    hits: list[AnchorHit],
    expected_order: list[str],
) -> list[SplitMolecule]:
    molecules: list[SplitMolecule] = []
    for cycle_index, cycle in enumerate(
        find_structure_cycles(hits, expected_order), start=1
    ):
        start = cycle[0].start
        end = cycle[-1].end
        if end <= start:
            continue
        molecules.append(
            SplitMolecule(
                parent_read_id=read_id,
                cycle_index=cycle_index,
                sequence=oriented_sequence[start:end],
                quality=oriented_quality[start:end],
                start=start,
                end=end,
            )
        )
    return molecules
