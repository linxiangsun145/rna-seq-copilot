"""
Analysis router — POST /validate/{job_id}  &  POST /run-analysis/{job_id}
"""
import json
import logging
from pathlib import Path

from fastapi import APIRouter, BackgroundTasks, HTTPException

from config import settings
from db.database import get_job, update_job_status
from models.schemas import JobStatus, RunAnalysisRequest, ValidationReport
from services.validator import validate_inputs
from services.r_runner import run_deseq2
from services.llm_client import generate_interpretation
from services.realism_validator import validate_realism
from services.report_builder import build_report

logger = logging.getLogger(__name__)
router = APIRouter()


def _job_dir(job_id: str) -> Path:
    p = Path(settings.jobs_dir) / job_id
    if not p.exists():
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return p


# ─── Validate ────────────────────────────────────────────────────────────────

@router.post("/validate/{job_id}", response_model=ValidationReport)
def validate_job(job_id: str):
    """Validate uploaded files and return a structured report."""
    job_dir = _job_dir(job_id)
    counts_path = job_dir / "counts.csv"
    meta_path = job_dir / "metadata.csv"

    if not counts_path.exists() or not meta_path.exists():
        raise HTTPException(status_code=400, detail="Upload files not found. Did you upload them?")

    report = validate_inputs(counts_path, meta_path)
    return report


# ─── Run analysis ─────────────────────────────────────────────────────────────

@router.post("/run-analysis/{job_id}")
def run_analysis(
    job_id: str,
    request: RunAnalysisRequest,
    background_tasks: BackgroundTasks,
):
    """Kick off DESeq2 analysis as a background task."""
    row = get_job(job_id)
    if row is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    if row["status"] == JobStatus.running.value:
        raise HTTPException(status_code=409, detail="Analysis already running")

    from datetime import datetime, timezone
    now = datetime.now(timezone.utc).isoformat()
    contrast_str = f"{request.contrast[0]}_{request.contrast[1]}_vs_{request.contrast[2]}"

    update_job_status(
        job_id,
        JobStatus.running,
        now,
        formula=request.formula,
        contrast=contrast_str,
    )

    background_tasks.add_task(
        _run_pipeline,
        job_id=job_id,
        formula=request.formula,
        contrast=list(request.contrast),
    )

    return {
        "job_id": job_id,
        "status": JobStatus.running.value,
        "created_at": row["created_at"],
        "updated_at": now,
        "formula": request.formula,
        "contrast": contrast_str,
    }


def _run_pipeline(job_id: str, formula: str, contrast: list[str]) -> None:
    """Background pipeline: DESeq2 → summary → LLM → report."""
    from datetime import datetime, timezone

    job_dir = Path(settings.jobs_dir) / job_id

    def _now():
        return datetime.now(timezone.utc).isoformat()

    try:
        logger.info("[%s] Running DESeq2…", job_id)
        summary = run_deseq2(job_dir, formula, contrast)

        logger.info("[%s] Running strict realism validation…", job_id)
        realism = validate_realism(job_dir, summary.model_dump())
        summary.realism_validation = realism

        qc_report = None
        qc_path = job_dir / "results" / "qc_report.json"
        if qc_path.exists():
            try:
                qc_report = json.loads(qc_path.read_text(encoding="utf-8"))
            except Exception:
                qc_report = None

        # Keep filesystem summary.json in sync with backend-enriched summary.
        summary_path = job_dir / "results" / "summary.json"
        summary_path.write_text(
            json.dumps(summary.model_dump(), indent=2),
            encoding="utf-8",
        )

        logger.info("[%s] Generating LLM interpretation…", job_id)
        llm = None
        try:
            llm = generate_interpretation(summary, qc_report=qc_report)
        except Exception as exc:
            logger.warning("[%s] LLM interpretation failed: %s", job_id, exc)

        logger.info("[%s] Building HTML report…", job_id)
        build_report(job_dir, summary, llm, formula=formula, contrast=contrast)

        update_job_status(
            job_id,
            JobStatus.done,
            _now(),
            summary=summary.model_dump(),
            llm_json=llm.model_dump() if llm else None,
        )
        logger.info("[%s] Pipeline complete.", job_id)

    except Exception as exc:
        logger.exception("[%s] Pipeline failed: %s", job_id, exc)
        update_job_status(job_id, JobStatus.failed, _now(), error=str(exc))
