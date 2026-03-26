"""
Report router — GET /report/{job_id}
Serves the generated HTML report.
"""
import logging
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from config import settings

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/report/{job_id}")
def get_report(job_id: str):
    report_path = Path(settings.jobs_dir) / job_id / "report.html"
    if not report_path.exists():
        raise HTTPException(status_code=404, detail="Report not yet available")
    return FileResponse(str(report_path), media_type="text/html")
