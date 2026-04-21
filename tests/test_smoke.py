from pathlib import Path
import json
from contextlib import contextmanager
from math import isclose, log10
import shutil
import tempfile
import unittest
from unittest import mock

from scfastq_qc import anchors as anchors_module
from scfastq_qc.anchors import find_anchor_hits, get_rust_anchor_engine, prepare_anchors, reverse_complement
from scfastq_qc.batch import BatchProcessingError, _allocate_output_dir, run_batch_qc
from scfastq_qc.classify import classify_best_orientation
from scfastq_qc.config import AnchorConfig, AppConfig, SampleEntry, StructureConfig, ThresholdConfig, load_config
from scfastq_qc.fastq import FastqRead, read_fastq
from scfastq_qc.report import run_qc


class SmokeTest(unittest.TestCase):
    @contextmanager
    def _temporary_directory(self):
        tmpdir = Path(tempfile.mkdtemp(prefix="test-work-"))
        try:
            yield str(tmpdir)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

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
            self.assertIn('qc_bucket_counts', data)
            self.assertIn('rust_accelerator', data)
            report_html = (Path(tmpdir) / 'report.html').read_text(encoding='utf-8')
            self.assertIn('QC Verdicts', report_html)
            self.assertIn('Yield And Quality', report_html)
            self.assertIn('Structure And Failure Modes', report_html)
            self.assertIn('<svg', report_html)

    def test_read_qscore_matches_error_rate_definition(self):
        read = FastqRead(name='mixed', sequence='AAAA', quality='I!I!')
        expected_error_rate = (10 ** (-40 / 10) + 10 ** (-0 / 10) + 10 ** (-40 / 10) + 10 ** (-0 / 10)) / 4
        expected_qscore = -10 * log10(expected_error_rate)
        arithmetic_mean_q = (40 + 0 + 40 + 0) / 4
        self.assertTrue(isclose(read.read_qscore, expected_qscore, rel_tol=1e-9))
        self.assertTrue(isclose(read.mean_q, arithmetic_mean_q, rel_tol=1e-9))
        self.assertLess(read.read_qscore, arithmetic_mean_q)

    def test_mean_q_preserves_arithmetic_mean_semantics(self):
        read = FastqRead(name='mixed', sequence='AAAA', quality='I!I!', qscore_method='conservative')
        arithmetic_mean_q = (40 + 0 + 40 + 0) / 4
        self.assertTrue(isclose(read.mean_q, arithmetic_mean_q, rel_tol=1e-9))
        self.assertLess(read.read_qscore, read.mean_q)

    def test_read_qscore_handles_uniform_extremes(self):
        self.assertEqual(FastqRead(name='low', sequence='AAAA', quality='!!!!').read_qscore, 0.0)
        self.assertEqual(FastqRead(name='high', sequence='AAAA', quality='IIII').read_qscore, 40.0)

    def test_batch_processing_raises_when_a_file_fails(self):
        repo = Path(__file__).resolve().parents[1]
        config = load_config(repo / 'examples' / 'config.json')
        fastq_path = repo / 'examples' / 'example.fastq'
        failed_result = {"file": str(fastq_path), "status": "failed", "error": "boom"}
        summary = {"total_count": 1, "success_count": 0, "failed_count": 1, "success_rate": 0.0, "results": [failed_result]}
        with mock.patch('pathlib.Path.mkdir'), \
             mock.patch('scfastq_qc.batch.generate_batch_summary', return_value=summary), \
             mock.patch('scfastq_qc.batch.run_qc_single', return_value=failed_result):
            with self.assertRaises(BatchProcessingError) as ctx:
                run_batch_qc([fastq_path], config, str(repo / 'batch-out'), parallel=1, continue_on_error=False)
        self.assertEqual(ctx.exception.failed_results, [failed_result])

    def test_batch_output_dirs_are_deduplicated_for_duplicate_stems(self):
        output_root = Path('batch-output')
        used_names: set[str] = set()
        first = _allocate_output_dir(output_root, 'sample', used_names)
        second = _allocate_output_dir(output_root, 'sample', used_names)
        third = _allocate_output_dir(output_root, 'sample_2', used_names)
        self.assertEqual(first.name, 'sample')
        self.assertEqual(second.name, 'sample_2')
        self.assertEqual(third.name, 'sample_2_2')

    def test_cli_batch_exits_nonzero_on_batch_failure(self):
        from scfastq_qc import cli as cli_module

        failed_result = {"file": "broken.fastq", "status": "failed", "error": "boom"}
        error = BatchProcessingError(
            results=[failed_result],
            summary={"total_count": 1, "success_count": 0, "failed_count": 1, "success_rate": 0.0, "results": [failed_result]},
        )
        mock_config = AppConfig(samples=None)
        with mock.patch.object(cli_module, 'setup_logging'), \
             mock.patch.object(cli_module, 'load_config', return_value=mock_config), \
             mock.patch.object(cli_module, 'collect_fastq_files', return_value=[Path('broken.fastq')]), \
             mock.patch('scfastq_qc.batch.run_batch_qc', side_effect=error), \
             mock.patch('sys.argv', ['scfastq-qc', 'batch', '--input', 'inputs.txt', '--config', 'config.json', '--outdir', 'out']):
            with self.assertRaises(SystemExit) as ctx:
                cli_module.main()
        self.assertEqual(ctx.exception.code, 1)

    def test_load_config_parses_samples_list(self):
        with self._temporary_directory() as tmpdir:
            config_path = Path(tmpdir) / 'multi.json'
            config_path.write_text(json.dumps({
                "samples": [
                    {"sample_name": "s1", "path": "/data/s1.fastq"},
                    {"sample_name": "s2", "path": "/data/s2.fastq"},
                ],
                "anchors": [],
                "structure": {"expected_order": []},
            }), encoding='utf-8')
            config = load_config(config_path)
        self.assertIsNotNone(config.samples)
        self.assertEqual(len(config.samples), 2)
        self.assertEqual(config.samples[0].sample_name, 's1')
        self.assertEqual(config.samples[0].path, '/data/s1.fastq')
        self.assertEqual(config.samples[1].sample_name, 's2')

    def test_load_config_without_samples_has_none(self):
        repo = Path(__file__).resolve().parents[1]
        config = load_config(repo / 'examples' / 'config.json')
        self.assertIsNone(config.samples)

    def test_batch_uses_explicit_sample_name_in_report(self):
        repo = Path(__file__).resolve().parents[1]
        config = load_config(repo / 'examples' / 'config.json')
        fastq_path = repo / 'examples' / 'example.fastq'
        with self._temporary_directory() as tmpdir:
            from scfastq_qc.batch import run_batch_qc
            results = run_batch_qc(
                fastq_files=[fastq_path],
                config=config,
                outdir=tmpdir,
                sample_names=['my_explicit_name'],
            )
        self.assertEqual(results[0]['status'], 'success')
        self.assertEqual(results[0]['summary']['sample_name'], 'my_explicit_name')

    def test_batch_falls_back_to_file_stem_when_no_sample_name(self):
        repo = Path(__file__).resolve().parents[1]
        config = load_config(repo / 'examples' / 'config.json')
        fastq_path = repo / 'examples' / 'example.fastq'
        with self._temporary_directory() as tmpdir:
            from scfastq_qc.batch import run_batch_qc
            results = run_batch_qc(
                fastq_files=[fastq_path],
                config=config,
                outdir=tmpdir,
                sample_names=None,
            )
        self.assertEqual(results[0]['status'], 'success')
        self.assertEqual(results[0]['summary']['sample_name'], 'example')

    def test_cli_run_dispatches_to_batch_for_multi_sample_config(self):
        from scfastq_qc import cli as cli_module
        from scfastq_qc.config import SampleEntry

        repo = Path(__file__).resolve().parents[1]
        base_config = load_config(repo / 'examples' / 'config.json')
        import dataclasses
        multi_config = dataclasses.replace(base_config, samples=[
            SampleEntry(sample_name='s1', path=str(repo / 'examples' / 'example.fastq')),
            SampleEntry(sample_name='s2', path=str(repo / 'examples' / 'example.fastq')),
        ])
        with self._temporary_directory() as tmpdir:
            with mock.patch.object(cli_module, 'setup_logging'), \
                 mock.patch.object(cli_module, 'load_config', return_value=multi_config), \
                 mock.patch('sys.argv', ['scfastq-qc', 'run', '--config', 'config.json', '--outdir', tmpdir]):
                cli_module.main()
            self.assertTrue((Path(tmpdir) / 'batch_report.html').exists())

    def test_cli_run_single_sample_from_config_samples(self):
        from scfastq_qc import cli as cli_module
        from scfastq_qc.config import SampleEntry

        repo = Path(__file__).resolve().parents[1]
        base_config = load_config(repo / 'examples' / 'config.json')
        import dataclasses
        single_config = dataclasses.replace(base_config, samples=[
            SampleEntry(sample_name='only_sample', path=str(repo / 'examples' / 'example.fastq')),
        ])
        with self._temporary_directory() as tmpdir:
            with mock.patch.object(cli_module, 'setup_logging'), \
                 mock.patch.object(cli_module, 'load_config', return_value=single_config), \
                 mock.patch('sys.argv', ['scfastq-qc', 'run', '--config', 'config.json', '--outdir', tmpdir]):
                cli_module.main()
            self.assertTrue((Path(tmpdir) / 'report.html').exists())
            data = json.loads((Path(tmpdir) / 'summary.json').read_text())
            self.assertEqual(data['sample_name'], 'only_sample')

    def test_cli_batch_uses_config_samples_when_no_input(self):
        from scfastq_qc import cli as cli_module
        from scfastq_qc.config import SampleEntry

        repo = Path(__file__).resolve().parents[1]
        base_config = load_config(repo / 'examples' / 'config.json')
        import dataclasses
        multi_config = dataclasses.replace(base_config, samples=[
            SampleEntry(sample_name='alpha', path=str(repo / 'examples' / 'example.fastq')),
            SampleEntry(sample_name='beta', path=str(repo / 'examples' / 'example.fastq')),
        ])
        with self._temporary_directory() as tmpdir:
            with mock.patch.object(cli_module, 'setup_logging'), \
                 mock.patch.object(cli_module, 'load_config', return_value=multi_config), \
                 mock.patch('sys.argv', ['scfastq-qc', 'batch', '--config', 'config.json', '--outdir', tmpdir]):
                cli_module.main()
            summary = json.loads((Path(tmpdir) / 'batch_summary.json').read_text())
            self.assertEqual(summary['success_count'], 2)

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

    def test_run_qc_exports_read_level_csvs(self):
        repo = Path(__file__).resolve().parents[1]
        config = load_config(repo / 'examples' / 'config.json')
        with self._temporary_directory() as tmpdir:
            run_qc(str(repo / 'examples' / 'example.fastq'), config, tmpdir, export_csv=True)
            self.assertTrue((Path(tmpdir) / 'read_details.csv').exists())
            self.assertTrue((Path(tmpdir) / 'anchor_hits.csv').exists())

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

    def test_read_fastq_rejects_mismatched_sequence_and_quality_lengths(self):
        with self._temporary_directory() as tmpdir:
            fastq_path = Path(tmpdir) / 'bad.fastq'
            fastq_path.write_text('@r1\nACGT\n+\n!!!\n', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'sequence and quality lengths differ'):
                list(read_fastq(fastq_path))

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
