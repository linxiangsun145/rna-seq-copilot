"""
R runner service — calls DESeq2 via Rscript subprocess.
"""
from __future__ import annotations

import json
import logging
import subprocess
import sys
from pathlib import Path

from config import settings
from models.schemas import AnalysisSummary

logger = logging.getLogger(__name__)


def run_deseq2(
    job_dir: Path,
    formula: str,
    contrast: list[str],
) -> AnalysisSummary:
    """
    Invoke the R DESeq2 script, collect summary JSON, return AnalysisSummary.
    Raises RuntimeError on failure.
    """
    r_script = Path(settings.r_scripts_dir) / "run_deseq2.R"
    if not r_script.exists():
        raise FileNotFoundError(f"DESeq2 R script not found: {r_script}")

    (job_dir / "plots").mkdir(exist_ok=True)
    (job_dir / "results").mkdir(exist_ok=True)

    cmd = [
        settings.rscript_path,
        str(r_script),
        "--counts", str(job_dir / "counts.csv"),
        "--metadata", str(job_dir / "metadata.csv"),
        "--formula", formula,
        "--contrast", ",".join(contrast),
        "--outdir", str(job_dir),
    ]

    logger.info("Running R command: %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=600,  # 10-minute ceiling
    )

    stderr_text = result.stderr or ""

    # Primary success indicator: summary.json exists (R may exit 1 due to package warnings)
    summary_path = job_dir / "results" / "summary.json"
    if not summary_path.exists():
        logger.error("R stderr:\n%s", stderr_text[-3000:])
        raise RuntimeError(
            f"DESeq2 script failed (exit {result.returncode}).\n"
            f"Last R stderr:\n{stderr_text[-2000:]}"
        )

    if result.returncode != 0:
        logger.warning("R exited with code %d (non-fatal, summary.json exists)", result.returncode)

    logger.debug("R stdout:\n%s", (result.stdout or "")[-2000:])

    # summary.json already checked above
    with open(summary_path) as f:
        raw = json.load(f)

    return AnalysisSummary(**raw)
