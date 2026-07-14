from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from anchorscope.anchors import AnchorHit
from anchorscope.classify import classify_best_orientation, classify_structure, find_structure_cycles
from anchorscope.config import AnchorRule, SegmentConfig, StructureConfig
from anchorscope.concatemers import split_complete_cycles
from anchorscope.segments import SegmentAccumulator, observe_segments


def hit(name: str, start: int, end: int) -> AnchorHit:
    return AnchorHit(name, start, end, 0, "A" * (end - start))


class StructureGrammarTests(unittest.TestCase):
    def test_terminal_and_distance_rules(self):
        config = StructureConfig(
            expected_order=["left", "right"],
            anchor_rules={
                "left": AnchorRule(terminal="5p", max_terminal_offset=0.1),
                "right": AnchorRule(
                    terminal="3p",
                    max_terminal_offset=0.1,
                    min_distance_from_previous=5,
                    max_distance_from_previous=20,
                ),
            },
        )
        valid = classify_structure([hit("left", 1, 5), hit("right", 15, 19)], config, 20)
        self.assertEqual(valid.label, "full_structure")
        invalid = classify_structure([hit("left", 5, 9), hit("right", 10, 14)], config, 20)
        self.assertEqual(invalid.label, "structure_rule_violation")
        self.assertIn("left:not_5p_terminal", invalid.violations)
        self.assertIn("right:distance_below_min", invalid.violations)

    def test_two_complete_cycles_are_concatemer(self):
        hits = [
            hit("left", 0, 4),
            hit("right", 10, 14),
            hit("left", 20, 24),
            hit("right", 30, 34),
        ]
        config = StructureConfig(expected_order=["left", "right"])
        structure = classify_structure(hits, config, 40)
        self.assertEqual(structure.label, "concatemer_candidate")
        self.assertEqual(structure.cycle_count, 2)
        self.assertEqual(len(find_structure_cycles(hits, config.expected_order)), 2)
        molecules = split_complete_cycles("read1", "A" * 40, "I" * 40, hits, config.expected_order)
        self.assertEqual([m.read_id for m in molecules], ["read1|cycle=1", "read1|cycle=2"])
        self.assertEqual([(m.start, m.end) for m in molecules], [(0, 14), (20, 34)])

    def test_orientation_is_only_failed_when_protocol_specifies_it(self):
        forward: list[AnchorHit] = []
        reverse = [hit("left", 0, 4), hit("right", 10, 14)]
        either, _ = classify_best_orientation(
            forward, reverse, StructureConfig(expected_order=["left", "right"]), 14
        )
        self.assertEqual(either.label, "full_structure")
        forward_only, _ = classify_best_orientation(
            forward,
            reverse,
            StructureConfig(expected_order=["left", "right"], expected_orientation="forward"),
            14,
        )
        self.assertEqual(forward_only.label, "orientation_unexpected")


class SegmentQCTests(unittest.TestCase):
    def test_segment_metrics_and_aggregation(self):
        segments = [
            SegmentConfig(
                name="insert",
                start_anchor="left",
                end_anchor="right",
                min_length=4,
                max_length=8,
            )
        ]
        observations = observe_segments(
            "AAAACCGTGGGG",
            "I" * 12,
            {"left": hit("left", 0, 4), "right": hit("right", 8, 12)},
            segments,
        )
        observation = observations["insert"]
        self.assertEqual(observation.length, 4)
        self.assertEqual(observation.status, "pass")
        self.assertAlmostEqual(observation.gc_fraction or 0, 0.75)
        accumulator = SegmentAccumulator()
        accumulator.update(observation)
        summary = accumulator.summary()
        self.assertEqual(summary["observed_ratio"], 1.0)
        self.assertEqual(summary["median_length"], 4.0)


if __name__ == "__main__":
    unittest.main()
