from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from benchmarks.validate_comprehensive_synthetic import run_validation


class ComprehensiveSyntheticValidationTests(unittest.TestCase):
    def test_truth_suite_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_validation(Path(tmp), replicates=1, seed=20260713, threads=1)
            self.assertTrue(result["overall_passed"])
            self.assertEqual(result["per_read_accuracy"], 1.0)
            self.assertGreaterEqual(result["scenario_types"], 30)


if __name__ == "__main__":
    unittest.main()
