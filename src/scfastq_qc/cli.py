from __future__ import annotations

import argparse

from .config import load_config
from .report import run_qc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Structure-aware FASTQ QC for single-cell long-read libraries")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run QC and generate HTML report")
    run_parser.add_argument("--fastq", required=True, help="Input FASTQ/FASTQ.GZ file")
    run_parser.add_argument("--config", required=True, help="JSON config file")
    run_parser.add_argument("--outdir", required=True, help="Output directory")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "run":
        config = load_config(args.config)
        run_qc(args.fastq, config, args.outdir)


if __name__ == "__main__":
    main()
