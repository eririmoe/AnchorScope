from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from anchorscope.anchors import find_anchor_hits, prepare_anchors
from anchorscope.config import AnchorConfig


class EditDistanceAnchorTests(unittest.TestCase):
    def _hits(self, sequence: str, motif: str = "ACGT"):
        prepared = prepare_anchors(
            [AnchorConfig(name="adapter", type="fixed", sequence=motif, max_edits=1)]
        )
        return find_anchor_hits(sequence, prepared)

    def test_substitution_is_reported(self):
        hit = self._hits("TTACCTAA")[0]
        self.assertEqual((hit.start, hit.end), (2, 6))
        self.assertEqual(hit.mismatches, 1)
        self.assertEqual(hit.substitutions, 1)
        self.assertEqual(hit.cigar, "2M1X1M")

    def test_insertion_is_reported(self):
        hit = self._hits("TTACAGTAA")[0]
        self.assertEqual(hit.matched_sequence, "ACAGT")
        self.assertEqual(hit.insertions, 1)
        self.assertEqual(hit.mismatches, 1)

    def test_deletion_is_reported(self):
        hit = self._hits("TTACTAA")[0]
        self.assertEqual(hit.matched_sequence, "ACT")
        self.assertEqual(hit.deletions, 1)
        self.assertEqual(hit.mismatches, 1)

    def test_search_window_limits_hits(self):
        prepared = prepare_anchors(
            [
                AnchorConfig(
                    name="adapter",
                    type="fixed",
                    sequence="ACGT",
                    max_edits=0,
                    search_region="5p",
                    search_window_bp=6,
                )
            ]
        )
        self.assertEqual(len(find_anchor_hits("ACGTTTACGT", prepared)), 1)
        self.assertEqual(find_anchor_hits("TTTTTTACGT", prepared), [])


if __name__ == "__main__":
    unittest.main()
