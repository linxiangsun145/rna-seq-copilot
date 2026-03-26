"""
File upload router — POST /upload
Saves uploaded files and creates a job record.
"""
import logging
import uuid
from datetime import datetime, timezone
from pathlib import Path

from fastapi import APIRouter, File, HTTPException, UploadFile

from config import settings
from db.database import create_job

logger = logging.getLogger(__name__)
router = APIRouter()

ALLOWED_EXTENSIONS = {".csv", ".tsv", ".txt"}
MAX_FILE_SIZE = 100 * 1024 * 1024  # 100 MB


def _validate_extension(filename: str) -> None:
    ext = Path(filename).suffix.lower()
    if ext not in ALLOWED_EXTENSIONS:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported file type '{ext}'. Allowed: {ALLOWED_EXTENSIONS}",
        )


@router.post("/upload")
async def upload_files(
    counts_file: UploadFile = File(...),
    metadata_file: UploadFile = File(...),
):
    """Receive counts + metadata files, persist to disk, return job_id."""
    _validate_extension(counts_file.filename or "file.csv")
    _validate_extension(metadata_file.filename or "file.csv")

    job_id = str(uuid.uuid4())
    job_dir = Path(settings.jobs_dir) / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # Save counts
    counts_path = job_dir / "counts.csv"
    counts_content = await counts_file.read()
    if len(counts_content) > MAX_FILE_SIZE:
        raise HTTPException(status_code=413, detail="counts_file exceeds 100 MB limit")
    counts_path.write_bytes(counts_content)

    # Save metadata
    meta_path = job_dir / "metadata.csv"
    meta_content = await metadata_file.read()
    if len(meta_content) > MAX_FILE_SIZE:
        raise HTTPException(status_code=413, detail="metadata_file exceeds 100 MB limit")
    meta_path.write_bytes(meta_content)

    # Create job record
    now = datetime.now(timezone.utc).isoformat()
    create_job(job_id, now)

    logger.info("Job %s created — counts=%d bytes, meta=%d bytes", job_id, len(counts_content), len(meta_content))
    return {"job_id": job_id}
