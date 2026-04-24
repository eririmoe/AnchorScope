from __future__ import annotations

import dataclasses
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


def _write_batch_html_report(summary: dict[str, Any], outdir: Path) -> None:
    successful = [item for item in summary["results"] if item.get("status") == "success" and item.get("summary")]
    rows = "".join(
        f"<tr><td>{Path(item['file']).name}</td><td>{item['summary']['total_reads']:,}</td><td>{item['summary']['long_high_quality_ratio']:.1%}</td><td>{item['summary']['correct_anchor_order_ratio']:.1%}</td><td>{item['summary']['qc_verdicts']['overall']['status']}</td></tr>"
        for item in successful
    )
    html = f"""<!DOCTYPE html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>scfastq-qc batch summary</title>
  <style>
    body {{ font-family:'IBM Plex Sans','Segoe UI',sans-serif; margin:0; background:#f7f5ef; color:#12232f; }}
    .page {{ max-width:1100px; margin:0 auto; padding:28px 18px 40px; }}
    .hero, .panel {{ background:rgba(255,255,255,.86); border:1px solid rgba(18,35,47,.1); border-radius:22px; box-shadow:0 18px 48px rgba(21,40,54,.08); }}
    .hero {{ padding:24px; }}
    .panel {{ margin-top:20px; padding:18px; }}
    .cards {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:14px; margin-top:18px; }}
    .card {{ background:rgba(255,255,255,.78); border:1px solid rgba(18,35,47,.08); border-radius:18px; padding:14px; }}
    .label {{ color:#5b6b78; font-size:.82rem; text-transform:uppercase; letter-spacing:.05em; }}
    .value {{ font-size:1.8rem; font-weight:700; margin-top:8px; }}
    table {{ width:100%; border-collapse:collapse; }}
    th, td {{ padding:12px 14px; text-align:left; border-bottom:1px solid rgba(18,35,47,.08); }}
    th {{ color:#5b6b78; font-size:.78rem; text-transform:uppercase; letter-spacing:.06em; }}
    tr:last-child td {{ border-bottom:none; }}
    @media (max-width:900px) {{ .cards {{ grid-template-columns:1fr 1fr; }} }}
  </style>
</head>
<body>
  <div class='page'>
    <div class='hero'>
      <h1>Batch QC Summary</h1>
      <p>Cross-sample comparison for the basic QC metrics produced by scfastq-qc.</p>
      <div class='cards'>
        <div class='card'><div class='label'>Samples</div><div class='value'>{summary['total_count']}</div></div>
        <div class='card'><div class='label'>Successful</div><div class='value'>{summary['success_count']}</div></div>
        <div class='card'><div class='label'>Failed</div><div class='value'>{summary['failed_count']}</div></div>
        <div class='card'><div class='label'>Success rate</div><div class='value'>{summary['success_rate']:.1%}</div></div>
      </div>
    </div>
    <div class='panel'>
      <h2>Per-sample comparison</h2>
      <table>
        <thead><tr><th>Sample</th><th>Total reads</th><th>Long high-quality</th><th>Anchor order</th><th>Overall QC</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </div>
  </div>
</body>
</html>"""
    (outdir / "batch_report.html").write_text(html, encoding="utf-8")


class BatchProcessingError(RuntimeError):
    def __init__(self, results: list[dict[str, Any]], summary: dict[str, Any]):
        self.results = results
        self.summary = summary
        self.failed_results = [result for result in results if result.get("status") == "failed"]
        preview = ", ".join(result["file"] for result in self.failed_results[:3])
        suffix = "..." if len(self.failed_results) > 3 else ""
        detail = f": {preview}{suffix}" if preview else ""
        super().__init__(f"Batch processing failed for {len(self.failed_results)} file(s){detail}")


def _write_batch_report(summary: dict[str, Any], outdir: Path) -> None:
    _write_batch_html_report(summary, outdir)


def run_qc_single(args: tuple[Path, AppConfig, Path, bool, bool, bool, str | None]) -> dict[str, Any]:
    """Process a single FASTQ file for batch execution."""
    fastq_path, config, outdir, export_csv, output_passed_fastq, output_failed_fastq, sample_name = args
    effective_name = sample_name if sample_name is not None else fastq_path.stem
    effective_config = dataclasses.replace(config, sample_name=effective_name)
    try:
        logger.info("Processing: %s", fastq_path)
        summary = run_qc(
            str(fastq_path), effective_config, str(outdir),
            export_csv=export_csv,
            output_passed_fastq=output_passed_fastq,
            output_failed_fastq=output_failed_fastq,
        )
        logger.info("Completed: %s", fastq_path)
        return {"file": str(fastq_path), "status": "success", "summary": summary}
    except Exception as exc:
        logger.error("Failed to process %s: %s", fastq_path, exc)
        return {"file": str(fastq_path), "status": "failed", "error": str(exc)}


def _allocate_output_dir(output_root: Path, base_name: str, used_names: set[str]) -> Path:
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
    output_passed_fastq: bool = False,
    output_failed_fastq: bool = False,
    sample_names: list[str | None] | None = None,
) -> tuple[list[tuple[Path, AppConfig, Path, bool, bool, bool, str | None]], list[dict[str, Any]]]:
    if sample_names is not None and len(sample_names) != len(fastq_files):
        raise ValueError(
            "sample_names length must match fastq_files length: "
            f"{len(sample_names)} != {len(fastq_files)}"
        )

    tasks: list[tuple[Path, AppConfig, Path, bool, bool, bool, str | None]] = []
    failures: list[dict[str, Any]] = []
    used_output_names: set[str] = set()

    for i, fastq_path in enumerate(fastq_files):
        sample_name = sample_names[i] if sample_names is not None else None
        base_name = sample_name if sample_name is not None else (fastq_path.stem or fastq_path.name)

        if not fastq_path.exists():
            error = f"File not found: {fastq_path}"
            logger.warning(error)
            failures.append({"file": str(fastq_path), "status": "failed", "error": error})
            if not continue_on_error:
                break
            continue

        file_outdir = _allocate_output_dir(output_root, base_name, used_output_names)
        tasks.append((fastq_path, config, file_outdir, export_csv, output_passed_fastq, output_failed_fastq, sample_name))

    return tasks, failures


def run_batch_qc(
    fastq_files: list[Path],
    config: AppConfig,
    outdir: str,
    parallel: int = 1,
    continue_on_error: bool = False,
    export_csv: bool = False,
    output_passed_fastq: bool = False,
    output_failed_fastq: bool = False,
    sample_names: list[str | None] | None = None,
) -> list[dict[str, Any]]:
    """Run QC across multiple FASTQ files."""
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    tasks, results = _build_batch_tasks(
        fastq_files, config, out_path, export_csv, continue_on_error,
        output_passed_fastq=output_passed_fastq,
        output_failed_fastq=output_failed_fastq,
        sample_names=sample_names,
    )
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
    _write_batch_report(summary, outdir)

    logger.info("Batch summary saved to: %s", summary_path)
    return summary
