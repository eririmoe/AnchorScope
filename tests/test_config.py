from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from anchorscope.config import load_config


class ConfigValidationTests(unittest.TestCase):
    def _load(self, payload: dict):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            return load_config(path)

    def test_nested_structure_rules_and_segments_load(self):
        config = self._load(
            {
                "anchors": [
                    {"name": "a", "type": "fixed", "sequence": "AAAA", "max_edits": 1},
                    {"name": "b", "type": "fixed", "sequence": "CCCC"},
                ],
                "structure": {
                    "expected_order": ["a", "b"],
                    "anchor_rules": {"a": {"terminal": "5p", "max_terminal_offset": 0.1}},
                    "segments": [
                        {"name": "insert", "start_anchor": "a", "end_anchor": "b"}
                    ],
                },
            }
        )
        self.assertEqual(config.structure.anchor_rules["a"].terminal, "5p")
        self.assertEqual(config.structure.segments[0].name, "insert")

    def test_unknown_anchor_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown anchor"):
            self._load(
                {
                    "anchors": [{"name": "a", "type": "fixed", "sequence": "AAAA"}],
                    "structure": {"expected_order": ["missing"]},
                }
            )

    def test_invalid_threshold_order_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "fail threshold"):
            self._load(
                {
                    "thresholds": {
                        "warn_long_high_quality_ratio": 0.5,
                        "fail_long_high_quality_ratio": 0.7,
                    }
                }
            )


if __name__ == "__main__":
    unittest.main()
