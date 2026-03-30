from __future__ import annotations

import json
import logging
from multiprocessing import Pool
from pathlib import Path
from typing import Any

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback for minimally provisioned environments
    def tqdm(iterable, **_: Any):
        return iterable

from .config import AppConfig
from .report import run_qc

logger = logging.getLogger(__name__)


class BatchProcessingError(RuntimeError):
    def __init__(self, results: list[dict[str, Any]], summary: dict[str, Any]):
        self.results = results
        self.summary = summary
        self.failed_results = [result for result in results if result.get("status") == "failed"]
        preview = ", ".join(result["file"] for result in self.failed_results[:3])
        suffix = "..." if len(self.failed_results) > 3 else ""
        detail = f": {preview}{suffix}" if preview else ""
        super().__init__(f"Batch processing failed for {len(self.failed_results)} file(s){detail}")


def run_qc_single(args: tuple[Path, AppConfig, Path, bool]) -> dict[str, Any]:
    """Process a single FASTQ file for batch execution."""
    fastq_path, config, outdir, export_csv = args
    try:
        logger.info("Processing: %s", fastq_path)
        summary = run_qc(str(fastq_path), config, str(outdir), export_csv=export_csv)
        logger.info("Completed: %s", fastq_path)
        return {"file": str(fastq_path), "status": "success", "summary": summary}
    except Exception as exc:
        logger.error("Failed to process %s: %s", fastq_path, exc)
        return {"file": str(fastq_path), "status": "failed", "error": str(exc)}


def _allocate_output_dir(output_root: Path, fastq_path: Path, used_names: set[str]) -> Path:
    base_name = fastq_path.stem or fastq_path.name
    candidate = base_name
    suffix = 2
    while candidate in used_names:
        candidate = f"{base_name}_{suffix}"
        suffix += 1
    used_names.add(candidate)
    return output_root / candidate


def _build_batch_tasks(
    fastq_files: list[Path],
    config: AppConfig,
    output_root: Path,
    export_csv: bool,
    continue_on_error: bool,
) -> tuple[list[tuple[Path, AppConfig, Path, bool]], list[dict[str, Any]]]:
    tasks: list[tuple[Path, AppConfig, Path, bool]] = []
    failures: list[dict[str, Any]] = []
    used_output_names: set[str] = set()

    for fastq_path in fastq_files:
        if not fastq_path.exists():
            error = f"File not found: {fastq_path}"
            logger.warning(error)
            failures.append({"file": str(fastq_path), "status": "failed", "error": error})
            if not continue_on_error:
                break
            continue

        file_outdir = _allocate_output_dir(output_root, fastq_path, used_output_names)
        tasks.append((fastq_path, config, file_outdir, export_csv))

    return tasks, failures


def run_batch_qc(
    fastq_files: list[Path],
    config: AppConfig,
    outdir: str,
    parallel: int = 1,
    continue_on_error: bool = False,
    export_csv: bool = False,
) -> list[dict[str, Any]]:
    """Run QC across multiple FASTQ files."""
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    tasks, results = _build_batch_tasks(fastq_files, config, out_path, export_csv, continue_on_error)
    logger.info("Starting batch processing: %s files", len(tasks))

    if tasks and (continue_on_error or not results):
        if parallel > 1:
            pool = Pool(processes=parallel)
            try:
                iterator = pool.imap(run_qc_single, tasks, chunksize=1)
                for result in tqdm(iterator, total=len(tasks), desc="Processing files"):
                    results.append(result)
                    if result["status"] == "failed" and not continue_on_error:
                        logger.error("Stopping batch processing due to error in %s", result["file"])
                        pool.terminate()
                        break
                else:
                    pool.close()
            except Exception:
                pool.terminate()
                raise
            finally:
                pool.join()
        else:
            for task in tqdm(tasks, desc="Processing files"):
                result = run_qc_single(task)
                results.append(result)
                if result["status"] == "failed" and not continue_on_error:
                    logger.error("Stopping batch processing due to error in %s", result["file"])
                    break

    summary = generate_batch_summary(results, out_path)
    logger.info(
        "Batch processing completed: %s/%s files succeeded",
        summary["success_count"],
        summary["total_count"],
    )

    failed_results = [result for result in results if result["status"] == "failed"]
    if failed_results and not continue_on_error:
        raise BatchProcessingError(results, summary)

    return results


def generate_batch_summary(results: list[dict[str, Any]], outdir: Path) -> dict[str, Any]:
    """Generate a batch-level JSON summary."""
    success_count = sum(1 for result in results if result["status"] == "success")
    failed_count = sum(1 for result in results if result["status"] == "failed")

    summary = {
        "total_count": len(results),
        "success_count": success_count,
        "failed_count": failed_count,
        "success_rate": success_count / len(results) if results else 0.0,
        "results": results,
    }

    summary_path = outdir / "batch_summary.json"
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    logger.info("Batch summary saved to: %s", summary_path)
    return summary
