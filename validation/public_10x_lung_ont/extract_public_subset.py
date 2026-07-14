from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path


DATASET_ID = "SC5pv2_GEX_Human_Lung_Carcinoma_DTC"
SOURCE_URL = (
    "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/7.0.1/"
    "SC5pv2_GEX_Human_Lung_Carcinoma_DTC/"
    "SC5pv2_GEX_Human_Lung_Carcinoma_DTC_ONT.fastq.gz"
)
DATASET_PAGE = (
    "https://www.10xgenomics.com/datasets/"
    "3k-human-squamous-cell-lung-carcinoma-dtcs-chromium-x-2-standard"
)


def extract_subset(archive_prefix: Path, output_dir: Path, reads: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fastq_path = output_dir / f"{DATASET_ID}.first_{reads}.fastq"
    provenance_path = output_dir / "provenance.json"

    digest = hashlib.sha256()
    written = 0
    with archive_prefix.open("rb") as raw:
        with gzip.GzipFile(fileobj=raw) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                with fastq_path.open("w", encoding="utf-8", newline="\n") as output:
                    for _ in range(reads):
                        record = [text.readline() for _ in range(4)]
                        if not record[0]:
                            break
                        if any(line == "" for line in record):
                            raise EOFError("The downloaded gzip prefix ended within a FASTQ record")
                        if not record[0].startswith("@") or not record[2].startswith("+"):
                            raise ValueError(f"Malformed FASTQ record at index {written + 1}")
                        sequence = record[1].rstrip("\r\n")
                        quality = record[3].rstrip("\r\n")
                        if len(sequence) != len(quality):
                            raise ValueError(
                                f"Sequence/quality length mismatch at record {written + 1}"
                            )
                        normalized = "".join(line.rstrip("\r\n") + "\n" for line in record)
                        output.write(normalized)
                        digest.update(normalized.encode("utf-8"))
                        written += 1

    if written != reads:
        raise EOFError(f"Requested {reads} reads but extracted only {written}")

    provenance = {
        "dataset_id": DATASET_ID,
        "dataset_page": DATASET_PAGE,
        "source_url": SOURCE_URL,
        "source_byte_range": "0-67108863",
        "source_prefix_file": archive_prefix.name,
        "selection": "first complete FASTQ records in archive order",
        "requested_reads": reads,
        "written_reads": written,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "subset_sha256": digest.hexdigest(),
        "license": "CC BY 4.0",
        "notes": (
            "The subset is used to demonstrate AnchorScope report behavior. "
            "It does not provide curated read-level structural truth."
        ),
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(f"fastq={fastq_path}")
    print(f"reads={written}")
    print(f"sha256={digest.hexdigest()}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--archive-prefix",
        type=Path,
        default=Path(__file__).resolve().parent / f"{DATASET_ID}.prefix_64MiB.fastq.gz",
    )
    parser.add_argument("--reads", type=int, default=10_000)
    args = parser.parse_args()
    if args.reads <= 0:
        parser.error("--reads must be positive")
    extract_subset(args.archive_prefix, args.output_dir, args.reads)


if __name__ == "__main__":
    main()
