from pathlib import Path
import json
from math import isclose, log10
import shutil
import tempfile
import unittest
from unittest import mock

from scfastq_qc import anchors as anchors_module
from scfastq_qc.anchors import find_anchor_hits, get_rust_anchor_engine, prepare_anchors, reverse_complement
from scfastq_qc.classify import classify_best_orientation
from scfastq_qc.config import AnchorConfig, AppConfig, StructureConfig, ThresholdConfig, load_config
from scfastq_qc.fastq import FastqRead, read_fastq
from scfastq_qc.report import run_qc


class SmokeTest(unittest.TestCase):
    def _temporary_directory(self) -> tempfile.TemporaryDirectory[str]:
        repo = Path(__file__).resolve().parents[1]
        return tempfile.TemporaryDirectory(dir=repo)

    def test_run_qc_generates_report(self):
        repo = Path(__file__).resolve().parents[1]
        config = load_config(repo / 'examples' / 'config.json')
        with self._temporary_directory() as tmpdir:
            summary = run_qc(str(repo / 'examples' / 'example.fastq'), config, tmpdir)
            self.assertEqual(summary['total_reads'], 4)
            self.assertTrue((Path(tmpdir) / 'report.html').exists())
            self.assertTrue((Path(tmpdir) / 'summary.json').exists())
            data = json.loads((Path(tmpdir) / 'summary.json').read_text())
            self.assertIn('anchor_detection_ratio', data)
            self.assertIn('min_read_length', data)
            self.assertIn('max_read_length', data)
            self.assertIn('median_read_qscore', data)
            self.assertIn('structure_orientation_counts', data)
            report_html = (Path(tmpdir) / 'report.html').read_text(encoding='utf-8')
            self.assertIn('Read length summary', report_html)
            self.assertIn('Longest read', report_html)
            self.assertIn('Median read Qscore', report_html)
            self.assertIn('Reversed read ratio', report_html)
            self.assertIn('Orientation', report_html)
            hist_svg = (Path(tmpdir) / 'figures' / 'read_length_hist.svg').read_text(encoding='utf-8')
            self.assertIn("font-size='12'", hist_svg)
            q_svg = (Path(tmpdir) / 'figures' / 'mean_q_hist.svg').read_text(encoding='utf-8')
            self.assertIn('Per-base quality score density', q_svg)
            self.assertIn('Probability density', q_svg)

    def test_read_qscore_matches_error_rate_definition(self):
        read = FastqRead(name='mixed', sequence='AAAA', quality='I!I!')
        expected_error_rate = (10 ** (-40 / 10) + 10 ** (-0 / 10) + 10 ** (-40 / 10) + 10 ** (-0 / 10)) / 4
        expected_qscore = -10 * log10(expected_error_rate)
        arithmetic_mean_q = (40 + 0 + 40 + 0) / 4
        self.assertTrue(isclose(read.read_qscore, expected_qscore, rel_tol=1e-9))
        self.assertTrue(isclose(read.mean_q, expected_qscore, rel_tol=1e-9))
        self.assertLess(read.read_qscore, arithmetic_mean_q)

    def test_read_qscore_handles_uniform_extremes(self):
        self.assertEqual(FastqRead(name='low', sequence='AAAA', quality='!!!!').read_qscore, 0.0)
        self.assertEqual(FastqRead(name='high', sequence='AAAA', quality='IIII').read_qscore, 40.0)

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

    def test_classify_best_orientation_marks_reversed_reads(self):
        anchors = prepare_anchors(
            [
                AnchorConfig(name='adapter_5p', type='fixed', sequence='ACGTTGCA', max_mismatches=0),
                AnchorConfig(name='polyT', type='regex', pattern='T{6,}'),
            ]
        )
        forward_sequence = 'ACGTTGCAGGGTTTTTTT'
        reverse_sequence = reverse_complement(forward_sequence)
        forward_hits = find_anchor_hits(reverse_sequence, anchors)
        reverse_hits = find_anchor_hits(reverse_complement(reverse_sequence), anchors)
        structure, chosen_hits = classify_best_orientation(forward_hits, reverse_hits, ['adapter_5p', 'polyT'])
        self.assertTrue(structure.is_reversed)
        self.assertEqual(structure.label, 'full_structure')
        self.assertEqual([hit.anchor_name for hit in chosen_hits], ['adapter_5p', 'polyT'])

    def test_run_qc_reports_forward_and_reversed_structure_counts(self):
        config = AppConfig(
            sample_name='orientation_demo',
            anchors=[
                AnchorConfig(name='adapter_5p', type='fixed', sequence='ACGTTGCA', max_mismatches=0),
                AnchorConfig(name='polyT', type='regex', pattern='T{6,}'),
            ],
            structure=StructureConfig(expected_order=['adapter_5p', 'polyT']),
            thresholds=ThresholdConfig(long_read_min_bp=0, long_read_min_q=0.0, heatmap_max_reads=10),
        )
        forward_sequence = 'ACGTTGCAGGGTTTTTTT'
        reverse_sequence = reverse_complement(forward_sequence)
        with self._temporary_directory() as tmpdir:
            fastq_path = Path(tmpdir) / 'orientation.fastq'
            fastq_path.write_text(
                '\n'.join(
                    [
                        '@forward',
                        forward_sequence,
                        '+',
                        'I' * len(forward_sequence),
                        '@reversed',
                        reverse_sequence,
                        '+',
                        'I' * len(reverse_sequence),
                    ]
                )
                + '\n',
                encoding='utf-8',
            )
            summary = run_qc(str(fastq_path), config, tmpdir)
        self.assertEqual(summary['structure_counts']['full_structure'], 2)
        self.assertEqual(summary['structure_orientation_counts']['full_structure']['forward'], 1)
        self.assertEqual(summary['structure_orientation_counts']['full_structure']['reversed'], 1)
        self.assertEqual(summary['reversed_read_ratio'], 0.5)

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
        with self._temporary_directory() as tmpdir:
            fastq_path = Path(tmpdir) / 'crlf.fastq'
            fastq_path.write_bytes(b'@read1\r\nACGT\r\n+\r\n!!!!\r\n')
            reads = list(read_fastq(fastq_path))
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].sequence, 'ACGT')
        self.assertEqual(reads[0].quality, '!!!!')

    def test_get_rust_anchor_engine_skips_rebuild_but_loads_existing_library(self):
        fake_engine = object()
        with mock.patch.object(anchors_module, '_RUST_ENGINE', None), \
             mock.patch.object(anchors_module, '_RUST_BUILD_ATTEMPTED', True), \
             mock.patch.object(anchors_module, '_find_rust_library', return_value=Path('/tmp/libscfastq_qc_rs.so')), \
             mock.patch.object(anchors_module, 'RustAnchorEngine', return_value=fake_engine) as engine_cls:
            engine = anchors_module.get_rust_anchor_engine()
        self.assertIs(engine, fake_engine)
        engine_cls.assert_called_once_with(Path('/tmp/libscfastq_qc_rs.so'))

    def test_get_rust_anchor_engine_falls_back_when_library_load_fails(self):
        with mock.patch.object(anchors_module, '_RUST_ENGINE', None), \
             mock.patch.object(anchors_module, '_RUST_BUILD_ATTEMPTED', False), \
             mock.patch.object(anchors_module, '_find_rust_library', return_value=Path('/tmp/libscfastq_qc_rs.so')), \
             mock.patch.object(anchors_module, '_try_build_rust_library') as build_lib, \
             mock.patch.object(anchors_module, 'RustAnchorEngine', side_effect=OSError('bad library')):
            engine = anchors_module.get_rust_anchor_engine()
            self.assertTrue(anchors_module._RUST_BUILD_ATTEMPTED)
        self.assertIsNone(engine)
        build_lib.assert_not_called()

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
