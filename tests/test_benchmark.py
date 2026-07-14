from __future__ import annotations

import unittest

from benchmarks.benchmark_synthetic import run_benchmark


class SyntheticBenchmarkTests(unittest.TestCase):
    def test_indel_aware_truth_set_is_reproducible(self):
        first = run_benchmark(reads_per_class=3, seed=17)
        second = run_benchmark(reads_per_class=3, seed=17)
        self.assertEqual(
            first["indel_aware"]["confusion_matrix"],
            second["indel_aware"]["confusion_matrix"],
        )
        self.assertEqual(first["indel_aware"]["accuracy"], 1.0)
        self.assertGreater(
            first["indel_aware"]["accuracy"], first["legacy_hamming"]["accuracy"]
        )


if __name__ == "__main__":
    unittest.main()
