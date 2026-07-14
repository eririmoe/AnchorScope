from __future__ import annotations

import argparse
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

from extract_public_subset import DATASET_ID, SOURCE_URL, extract_subset
from validate_public_results import validate


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
PREFIX_BYTES = 64 * 1024 * 1024


def download_prefix(url: str, destination: Path, byte_count: int, force: bool = False) -> None:
    if destination.is_file() and destination.stat().st_size == byte_count and not force:
        print(f"reuse_prefix={destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(
        url,
        headers={"Range": f"bytes=0-{byte_count - 1}", "User-Agent": "AnchorScope-validation/0.2"},
    )
    written = 0
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with temporary.open("wb") as output:
                while written < byte_count:
                    chunk = response.read(min(1024 * 1024, byte_count - written))
                    if not chunk:
                        break
                    output.write(chunk)
                    written += len(chunk)
        if written != byte_count:
            raise EOFError(f"Expected {byte_count} source bytes, downloaded {written}")
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    print(f"downloaded_prefix={destination}")


def run_anchorscope(fastq: Path, config: Path, run_dir: Path, threads: int) -> None:
    environment = os.environ.copy()
    source_path = str(REPOSITORY / "src")
    environment["PYTHONPATH"] = source_path + os.pathsep + environment.get("PYTHONPATH", "")
    command = [
        sys.executable,
        "-B",
        "-m",
        "anchorscope.cli",
        "run",
        "--fastq",
        str(fastq),
        "--config",
        str(config),
        "--outdir",
        str(run_dir),
        "--threads",
        str(threads),
        "--export-csv",
    ]
    print("command=" + " ".join(command))
    subprocess.run(command, cwd=REPOSITORY, env=environment, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download, run and validate the public 10x ONT AnchorScope example."
    )
    parser.add_argument("--work-dir", type=Path, default=HERE)
    parser.add_argument("--reads", type=int, default=10_000)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--skip-run", action="store_true", help="Only validate existing outputs")
    args = parser.parse_args()
    if args.reads != 10_000:
        parser.error("The published expected-results contract currently requires --reads 10000")
    if args.threads < 1:
        parser.error("--threads must be positive")

    work_dir = args.work_dir.resolve()
    prefix = work_dir / f"{DATASET_ID}.prefix_64MiB.fastq.gz"
    fastq = work_dir / f"{DATASET_ID}.first_{args.reads}.fastq"
    provenance = work_dir / "provenance.json"
    run_dir = work_dir / "anchorscope_report"

    if not args.skip_run:
        download_prefix(SOURCE_URL, prefix, PREFIX_BYTES, force=args.force_download)
        extract_subset(prefix, work_dir, args.reads)
        run_anchorscope(fastq, HERE / "config.json", run_dir, args.threads)

    errors = validate(run_dir, fastq, provenance, HERE / "expected_results.json")
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        raise SystemExit(1)
    print(f"report={run_dir / 'report.html'}")
    print("PASS: download, analysis and result contract are reproducible")


if __name__ == "__main__":
    main()
