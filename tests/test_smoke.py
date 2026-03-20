from pathlib import Path
import json
import shutil
import tempfile
import unittest
from unittest import mock

from scfastq_qc import anchors as anchors_module
from scfastq_qc.anchors import find_anchor_hits, get_rust_anchor_engine, prepare_anchors
from scfastq_qc.config import AnchorConfig, load_config
from scfastq_qc.fastq import read_fastq
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
        engine = get_rust_anchor_engine()
        if engine is None:
            self.skipTest('Rust accelerator build/load did not succeed; optional dependency remains unavailable')
        self.assertIsNotNone(engine)

    def test_prepare_anchors_rejects_negative_max_mismatches(self):
        with self.assertRaisesRegex(ValueError, 'max_mismatches must be non-negative'):
            prepare_anchors([AnchorConfig(name='bad', type='fixed', sequence='ACGT', max_mismatches=-1)])

    def test_prepare_anchors_requires_sequence_and_pattern(self):
        with self.assertRaisesRegex(ValueError, "requires a non-empty sequence"):
            prepare_anchors([AnchorConfig(name='fixed_missing', type='fixed')])
        with self.assertRaisesRegex(ValueError, "requires a non-empty pattern"):
            prepare_anchors([AnchorConfig(name='regex_missing', type='regex')])

    def test_read_fastq_handles_crlf_sequences(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = Path(tmpdir) / 'crlf.fastq'
            fastq_path.write_bytes(b'@read1\r\nACGT\r\n+\r\n!!!!\r\n')
            reads = list(read_fastq(fastq_path))
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].sequence, 'ACGT')
        self.assertEqual(reads[0].quality, '!!!!')

    def test_get_rust_anchor_engine_caches_loaded_engine(self):
        fake_engine = object()
        with mock.patch.object(anchors_module, '_RUST_ENGINE', None), \
             mock.patch.object(anchors_module, '_RUST_BUILD_ATTEMPTED', True), \
             mock.patch.object(anchors_module, '_find_rust_library', return_value=Path('/tmp/libscfastq_qc_rs.so')), \
             mock.patch.object(anchors_module, 'RustAnchorEngine', return_value=fake_engine) as engine_cls:
            first = anchors_module.get_rust_anchor_engine()
            second = anchors_module.get_rust_anchor_engine()
        self.assertIs(first, fake_engine)
        self.assertIs(second, fake_engine)
        self.assertEqual(engine_cls.call_count, 1)


if __name__ == '__main__':
    unittest.main()
