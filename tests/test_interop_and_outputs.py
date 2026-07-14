from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from anchorscope.fastq import fastq_output_path, open_output_text
from anchorscope.interop import audit_bam, summarize_tag_records
from anchorscope.multiqc import write_multiqc_summary


class TagAuditTests(unittest.TestCase):
    def test_external_tags_and_poly_a_are_summarized(self):
        summary = summarize_tag_records(
            [
                {"CR": "ACGT", "CB": "ACGT-1", "UR": "AAAA", "UB": "AAAT", "pt": 35},
                {"CR": "ACGA", "CB": "ACGT-1", "UR": "CCCC", "UB": "CCCC", "pt": 25},
                {"CR": "TTTT"},
            ],
            "sample",
        )
        self.assertEqual(summary["total_reads"], 3)
        self.assertAlmostEqual(summary["corrected_barcode_ratio"], 2 / 3)
        self.assertEqual(summary["barcode_edit_distance_counts"], {0: 1, 1: 1})
        self.assertEqual(summary["umi_edit_distance_counts"], {0: 1, 1: 1})
        self.assertEqual(summary["median_poly_a_length"], 25.0)

    def test_sam_file_audit_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sam = root / "tagged.sam"
            sam.write_text(
                "@HD\tVN:1.6\tSO:unknown\n"
                "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tCR:Z:ACGA\tCB:Z:ACGT-1\tUR:Z:AAAA\tUB:Z:AAAT\tpt:i:35\n",
                encoding="utf-8",
            )
            outdir = root / "audit"
            summary = audit_bam(str(sam), str(outdir), sample_name="sample")
            self.assertEqual(summary["barcode_edit_distance_counts"], {1: 1})
            self.assertEqual(summary["median_poly_a_length"], 35.0)
            self.assertTrue((outdir / "tag_audit_summary.json").is_file())
            self.assertTrue((outdir / "anchorscope_tag_audit_mqc.json").is_file())


class WorkflowOutputTests(unittest.TestCase):
    def test_gzip_fastq_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = fastq_output_path(tmp, "passed", True)
            with open_output_text(path) as handle:
                handle.write("@r1\nACGT\n+\nIIII\n")
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "@r1\nACGT\n+\nIIII\n")

    def test_multiqc_custom_content_is_emitted(self):
        summary = {
            "sample_name": "sample",
            "total_reads": 10,
            "median_read_length": 100.0,
            "median_read_qscore": 20.0,
            "correct_anchor_order_ratio": 0.8,
            "qc_bucket_ratios": {"pass": 0.7, "no_anchor": 0.1},
            "barcode_qc": None,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "anchorscope_mqc.json"
            write_multiqc_summary(summary, path)
            payload = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(payload["plot_type"], "generalstats")
            self.assertEqual(payload["data"]["sample"]["total_reads"], 10)


if __name__ == "__main__":
    unittest.main()
