from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .config import load_config
from .logging_config import setup_logging
from .report import run_qc

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Structure-aware FASTQ QC for single-cell long-read libraries"
    )

    parser.add_argument("--log-file", help="Log file path (optional)")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Log level (default: INFO)",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run QC on a single FASTQ file")
    run_parser.add_argument("--fastq", required=True, help="Input FASTQ/FASTQ.GZ file")
    run_parser.add_argument("--config", required=True, help="JSON config file")
    run_parser.add_argument("--outdir", required=True, help="Output directory")
    run_parser.add_argument(
        "--export-csv",
        action="store_true",
        help="Export results to CSV files",
    )

    batch_parser = subparsers.add_parser("batch", help="Run QC on multiple FASTQ files")
    batch_parser.add_argument(
        "--input",
        required=True,
        help="Input directory or file list (one path per line)",
    )
    batch_parser.add_argument("--config", required=True, help="JSON config file")
    batch_parser.add_argument("--outdir", required=True, help="Output directory")
    batch_parser.add_argument(
        "--pattern",
        default="*.fastq*",
        help="File pattern for directory input (default: *.fastq*)",
    )
    batch_parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel processes (default: 1)",
    )
    batch_parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue processing if one file fails",
    )
    batch_parser.add_argument(
        "--export-csv",
        action="store_true",
        help="Export results to CSV files",
    )

    return parser


def collect_fastq_files(input_path: str, pattern: str) -> list[Path]:
    source = Path(input_path)

    if source.is_file():
        files: list[Path] = []
        with source.open("r", encoding="utf-8") as handle:
            for line in handle:
                entry = line.strip()
                if not entry:
                    continue
                path = Path(entry)
                if not path.is_absolute():
                    path = source.parent / path
                files.append(path)
        return files

    if source.is_dir():
        return sorted(source.glob(pattern))

    raise ValueError(f"Input path does not exist: {source}")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    log_file = Path(args.log_file) if args.log_file else None
    setup_logging(log_file=log_file, level=args.log_level)

    if args.command == "run":
        config = load_config(args.config)
        run_qc(args.fastq, config, args.outdir, export_csv=args.export_csv)
        return

    from .batch import BatchProcessingError, run_batch_qc

    config = load_config(args.config)
    fastq_files = collect_fastq_files(args.input, args.pattern)
    try:
        run_batch_qc(
            fastq_files=fastq_files,
            config=config,
            outdir=args.outdir,
            parallel=args.parallel,
            continue_on_error=args.continue_on_error,
            export_csv=args.export_csv,
        )
    except BatchProcessingError as exc:
        logger.error("%s", exc)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
