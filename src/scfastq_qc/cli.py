from __future__ import annotations

import argparse
import dataclasses
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
    run_parser.add_argument(
        "--fastq",
        default=None,
        help="Input FASTQ/FASTQ.GZ file (required when 'samples' is not provided in config)",
    )
    run_parser.add_argument("--config", required=True, help="JSON config file")
    run_parser.add_argument("--outdir", required=True, help="Output directory")
    run_parser.add_argument(
        "--export-csv",
        action="store_true",
        help="Export results to CSV files",
    )
    run_parser.add_argument(
        "--output-passed-fastq",
        action="store_true",
        help="Export reads with QC bucket 'pass' to passed.fastq",
    )
    run_parser.add_argument(
        "--output-failed-fastq",
        action="store_true",
        help="Export reads with QC bucket != 'pass' to failed.fastq",
    )
    run_parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for single-sample processing (default: 1)",
    )

    batch_parser = subparsers.add_parser("batch", help="Run QC on multiple FASTQ files")
    batch_parser.add_argument(
        "--input",
        default=None,
        help="Input directory or file list (one path per line); required when 'samples' is not provided in config",
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
    batch_parser.add_argument(
        "--output-passed-fastq",
        action="store_true",
        help="Export reads with QC bucket 'pass' to passed.fastq",
    )
    batch_parser.add_argument(
        "--output-failed-fastq",
        action="store_true",
        help="Export reads with QC bucket != 'pass' to failed.fastq",
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
        if config.samples:
            from .batch import BatchProcessingError, run_batch_qc
            if len(config.samples) == 1:
                entry = config.samples[0]
                effective_config = dataclasses.replace(config, sample_name=entry.sample_name)
                run_qc(
                    entry.path, effective_config, args.outdir,
                    export_csv=args.export_csv,
                    output_passed_fastq=args.output_passed_fastq,
                    output_failed_fastq=args.output_failed_fastq,
                    threads=args.threads,
                )
            else:
                fastq_files = [Path(e.path) for e in config.samples]
                sample_names: list[str] = [e.sample_name for e in config.samples]
                try:
                    run_batch_qc(
                        fastq_files=fastq_files,
                        config=config,
                        outdir=args.outdir,
                        parallel=1,
                        continue_on_error=False,
                        export_csv=args.export_csv,
                        output_passed_fastq=args.output_passed_fastq,
                        output_failed_fastq=args.output_failed_fastq,
                        sample_names=sample_names,
                    )
                except BatchProcessingError as exc:
                    logger.error("%s", exc)
                    raise SystemExit(1) from None
        else:
            if not args.fastq:
                parser.error("run: --fastq is required when 'samples' is not provided in config")
            run_qc(
                args.fastq, config, args.outdir,
                export_csv=args.export_csv,
                output_passed_fastq=args.output_passed_fastq,
                output_failed_fastq=args.output_failed_fastq,
                threads=args.threads,
            )
        return

    from .batch import BatchProcessingError, run_batch_qc

    config = load_config(args.config)

    if config.samples and not args.input:
        fastq_files = [Path(e.path) for e in config.samples]
        batch_sample_names: list[str] = [e.sample_name for e in config.samples]
    elif args.input:
        fastq_files = collect_fastq_files(args.input, args.pattern)
        batch_sample_names = []
    else:
        parser.error("batch: --input is required when 'samples' is not provided in config")

    try:
        run_batch_qc(
            fastq_files=fastq_files,
            config=config,
            outdir=args.outdir,
            parallel=args.parallel,
            continue_on_error=args.continue_on_error,
            export_csv=args.export_csv,
            output_passed_fastq=args.output_passed_fastq,
            output_failed_fastq=args.output_failed_fastq,
            sample_names=batch_sample_names if batch_sample_names else None,
        )
    except BatchProcessingError as exc:
        logger.error("%s", exc)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
