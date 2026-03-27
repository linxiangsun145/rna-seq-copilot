"""
Results router
GET /results/{job_id}           — polling endpoint (returns ResultsPayload)
GET /results/{job_id}/deg.csv   — DEG table download
GET /results/{job_id}/plots/{plot_name} — plot image serving
"""
import json
import logging
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from config import settings
from db.database import get_job
from models.schemas import (
    AnalysisSummary,
    JobStatus,
    LLMInterpretation,
    QCReport,
    ResultsPayload,
)

logger = logging.getLogger(__name__)
router = APIRouter()

BASE_URL = ""  # resolved at runtime


def _job_dir(job_id: str) -> Path:
    p = Path(settings.jobs_dir) / job_id
    if not p.exists():
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return p


@router.get("/results/{job_id}", response_model=ResultsPayload)
def get_results(job_id: str):
    row = get_job(job_id)
    if row is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    job_dir = _job_dir(job_id)
    status = JobStatus(row["status"])

    summary: AnalysisSummary | None = None
    if row["summary"]:
        summary = AnalysisSummary(**json.loads(row["summary"]))

    llm: LLMInterpretation | None = None
    if row["llm_json"]:
        llm = LLMInterpretation(**json.loads(row["llm_json"]))

    qc_report: QCReport | None = None
    qc_path = job_dir / "results" / "qc_report.json"
    if qc_path.exists():
        qc_report = QCReport(**json.loads(qc_path.read_text(encoding="utf-8")))

    # Build plot presence map
    plot_dir = job_dir / "plots"
    plots: dict = {}
    for name in [
        "pca",
        "volcano",
        "ma",
        "heatmap",
        "sample_distance_heatmap",
        "sample_correlation_heatmap",
        "library_size",
        "count_distribution",
        "zero_fraction",
    ]:
        p = plot_dir / f"{name}.png"
        plots[name] = f"/results/{job_id}/plots/{name}.png" if p.exists() else None

    deg_path = job_dir / "results" / "deg_results.csv"
    deg_url = f"/results/{job_id}/deg.csv" if deg_path.exists() else None

    report_path = job_dir / "report.html"
    report_url = f"/report/{job_id}" if report_path.exists() else None

    return ResultsPayload(
        job_id=job_id,
        status=status,
        summary=summary,
        qc_report=qc_report,
        deg_table_url=deg_url,
        plots=plots,
        llm_interpretation=llm,
        report_url=report_url,
        error=row["error"],
    )


@router.get("/results/{job_id}/deg.csv")
def download_deg(job_id: str):
    job_dir = _job_dir(job_id)
    deg_path = job_dir / "results" / "deg_results.csv"
    if not deg_path.exists():
        raise HTTPException(status_code=404, detail="DEG results not yet available")
    return FileResponse(
        str(deg_path),
        media_type="text/csv",
        filename=f"deg_results_{job_id}.csv",
    )


@router.get("/results/{job_id}/plots/{plot_name}")
def get_plot(job_id: str, plot_name: str):
    # Sanitize to prevent path traversal
    safe_name = Path(plot_name).name
    job_dir = _job_dir(job_id)
    plot_path = job_dir / "plots" / safe_name
    if not plot_path.exists() or plot_path.suffix != ".png":
        raise HTTPException(status_code=404, detail="Plot not found")
    return FileResponse(str(plot_path), media_type="image/png")
