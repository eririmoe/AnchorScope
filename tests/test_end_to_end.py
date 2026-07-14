from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from anchorscope.config import (
    AnchorConfig,
    AnchorRule,
    AppConfig,
    BarcodeConfig,
    SegmentConfig,
    StructureConfig,
    ThresholdConfig,
)
from anchorscope.report import run_qc


class EndToEndTests(unittest.TestCase):
    def test_complete_pipeline_and_parallel_consistency(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            whitelist = root / "whitelist.txt"
            whitelist.write_text("ACGT\nTTTT\n", encoding="utf-8")
            full = "AAAAACGTGGCCCC"
            concatemer = full + full
            fastq = root / "reads.fastq"
            fastq.write_text(
                "".join(
                    [
                        f"@full\n{full}\n+\n{'I' * len(full)}\n",
                        f"@concatemer\n{concatemer}\n+\n{'I' * len(concatemer)}\n",
                    ]
                ),
                encoding="utf-8",
            )
            config = AppConfig(
                sample_name="synthetic",
                anchors=[
                    AnchorConfig("left", "fixed", sequence="AAAA", max_edits=0),
                    AnchorConfig("right", "fixed", sequence="CCCC", max_edits=0),
                ],
                structure=StructureConfig(
                    expected_order=["left", "right"],
                    anchor_rules={
                        "left": AnchorRule(max_count=1),
                        "right": AnchorRule(max_count=1),
                    },
                    segments=[SegmentConfig("payload", "left", "right", min_length=6, max_length=6)],
                    split_concatemers=True,
                ),
                barcode=BarcodeConfig(
                    enabled=True,
                    left_anchor="left",
                    right_anchor="right",
                    barcode_length=4,
                    umi_length=2,
                    whitelist_path=str(whitelist),
                ),
                thresholds=ThresholdConfig(long_read_min_bp=1, long_read_min_q=0),
            )
            serial_dir = root / "serial"
            parallel_dir = root / "parallel"
            serial = run_qc(
                str(fastq),
                config,
                str(serial_dir),
                output_passed_fastq=True,
                output_failed_fastq=True,
                export_csv=True,
                gzip_output=True,
            )
            parallel = run_qc(
                str(fastq), config, str(parallel_dir), threads=2, export_csv=True
            )

            for key in (
                "total_reads",
                "structure_counts",
                "qc_bucket_counts",
                "segment_qc",
                "barcode_qc",
                "correct_anchor_order_ratio",
            ):
                self.assertEqual(serial[key], parallel[key])
            self.assertEqual(serial["structure_counts"]["concatemer_candidate"], 1)
            self.assertEqual(serial["matcher_backend"]["mode"], "python-edit")
            self.assertEqual(serial["split_molecules_exported"], 2)
            self.assertEqual(serial["barcode_qc"]["status_counts"]["exact"], 2)
            self.assertTrue((serial_dir / "report.html").is_file())
            self.assertTrue((serial_dir / "anchorscope_mqc.json").is_file())
            self.assertTrue((serial_dir / "split_concatemers.fastq.gz").is_file())
            with gzip.open(serial_dir / "split_concatemers.fastq.gz", "rt", encoding="utf-8") as handle:
                self.assertIn("@concatemer|cycle=2", handle.read())
            json.loads((serial_dir / "summary.json").read_text(encoding="utf-8"))
            for csv_name in ("read_details.csv", "anchor_hits.csv"):
                self.assertEqual(
                    (serial_dir / csv_name).read_text(encoding="utf-8"),
                    (parallel_dir / csv_name).read_text(encoding="utf-8"),
                )
            report_text = (serial_dir / "report.html").read_text(encoding="utf-8")
            self.assertTrue(report_text.startswith("<!doctype html>"))
            self.assertIn("Anchor-relative segment QC", report_text)
            self.assertIn("Barcode and UMI recoverability", report_text)
            self.assertIn("Structure and failure modes", report_text)
            self.assertGreaterEqual(report_text.count("<svg"), 5)
            self.assertTrue(report_text.rstrip().endswith("</html>"))


if __name__ == "__main__":
    unittest.main()
