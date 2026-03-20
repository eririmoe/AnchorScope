from pathlib import Path
import json
import shutil
import tempfile
import unittest

from scfastq_qc.anchors import find_anchor_hits, get_rust_anchor_engine, prepare_anchors
from scfastq_qc.config import AnchorConfig, load_config
from scfastq_qc.report import run_qc


class SmokeTest(unittest.TestCase):
    def test_run_qc_generates_report(self):
        repo = Path(__file__).resolve().parents[1]
        config = load_config(repo / 'examples' / 'config.json')
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = run_qc(str(repo / 'examples' / 'example.fastq'), config, tmpdir)
            self.assertEqual(summary['total_reads'], 4)
            self.assertTrue((Path(tmpdir) / 'report.html').exists())
            self.assertTrue((Path(tmpdir) / 'summary.json').exists())
            data = json.loads((Path(tmpdir) / 'summary.json').read_text())
            self.assertIn('anchor_detection_ratio', data)

    def test_prepare_anchors_exact_and_approximate_paths(self):
        anchors = prepare_anchors(
            [
                AnchorConfig(name='exact', type='fixed', sequence='ACGT', max_mismatches=0),
                AnchorConfig(name='approx', type='fixed', sequence='AACCGG', max_mismatches=1),
                AnchorConfig(name='polyT', type='regex', pattern='T{4,}'),
            ]
        )
        hits = find_anchor_hits('TTACGTGGAACCAGTTTT', anchors)
        self.assertEqual([hit.anchor_name for hit in hits], ['exact', 'approx', 'polyT'])
        self.assertEqual(hits[0].mismatches, 0)
        self.assertEqual(hits[1].mismatches, 1)

    @unittest.skipUnless(shutil.which('cargo'), 'cargo not available')
    def test_rust_engine_available_after_auto_build(self):
        self.assertIsNotNone(get_rust_anchor_engine())


if __name__ == '__main__':
    unittest.main()
