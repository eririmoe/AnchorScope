from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from anchorscope.anchors import AnchorHit
from anchorscope.barcode import (
    BarcodeAccumulator,
    BarcodeMatcher,
    evaluate_barcode_recoverability,
)
from anchorscope.config import BarcodeConfig


def hit(name: str, start: int, end: int) -> AnchorHit:
    return AnchorHit(name, start, end, 0, "A" * (end - start))


class BarcodeRecoverabilityTests(unittest.TestCase):
    def setUp(self):
        self.config = BarcodeConfig(
            enabled=True,
            left_anchor="left",
            right_anchor="right",
            barcode_length=4,
            umi_length=2,
            min_base_quality=15,
            max_edit_distance=1,
            min_distance_margin=1,
        )
        self.hits = {"left": hit("left", 0, 4), "right": hit("right", 10, 14)}

    def _matcher(self, values: list[str]) -> BarcodeMatcher:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "whitelist.txt"
            path.write_text("\n".join(values), encoding="utf-8")
            return BarcodeMatcher.from_path(str(path))

    def test_exact_barcode_and_umi_are_extracted(self):
        observation = evaluate_barcode_recoverability(
            "AAAAACGTGGCCCC",
            "I" * 14,
            self.hits,
            self.config,
            self._matcher(["ACGT", "TTTT"]),
        )
        self.assertEqual(observation.status, "exact")
        self.assertEqual(observation.raw_barcode, "ACGT")
        self.assertEqual(observation.raw_umi, "GG")

    def test_unique_and_ambiguous_candidates_are_distinguished(self):
        unique = evaluate_barcode_recoverability(
            "AAAAACGAGGCCCC",
            "I" * 14,
            self.hits,
            self.config,
            self._matcher(["ACGT", "TTTT"]),
        )
        self.assertEqual(unique.status, "uniquely_recoverable")
        ambiguous = evaluate_barcode_recoverability(
            "AAAAAGGTGGCCCC",
            "I" * 14,
            self.hits,
            self.config,
            self._matcher(["ACGT", "ATGT"]),
        )
        self.assertEqual(ambiguous.status, "ambiguous")

    def test_low_quality_is_not_called_recoverable(self):
        observation = evaluate_barcode_recoverability(
            "AAAAACGTGGCCCC",
            "!" * 14,
            self.hits,
            self.config,
            self._matcher(["ACGT"]),
        )
        self.assertEqual(observation.status, "low_quality")
        accumulator = BarcodeAccumulator()
        accumulator.update(observation)
        self.assertEqual(accumulator.summary()["low_quality_ratio"], 1.0)


if __name__ == "__main__":
    unittest.main()
