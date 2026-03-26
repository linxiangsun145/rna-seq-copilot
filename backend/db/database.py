"""
SQLite-backed job store.
Stores job metadata; analysis artefacts are written to the filesystem.
"""
import json
import logging
import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Optional

from config import settings
from models.schemas import JobStatus

logger = logging.getLogger(__name__)

DB_PATH = Path(settings.jobs_dir) / "jobs.db"


def get_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH), check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


@contextmanager
def db_cursor():
    conn = get_connection()
    try:
        cur = conn.cursor()
        yield cur
        conn.commit()
    finally:
        conn.close()


def init_db() -> None:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    with db_cursor() as cur:
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                job_id      TEXT PRIMARY KEY,
                status      TEXT NOT NULL DEFAULT 'pending',
                created_at  TEXT NOT NULL,
                updated_at  TEXT NOT NULL,
                formula     TEXT,
                contrast    TEXT,
                error       TEXT,
                summary     TEXT,
                llm_json    TEXT
            )
            """
        )
    logger.info("Database initialised at %s", DB_PATH)


def create_job(job_id: str, created_at: str) -> None:
    with db_cursor() as cur:
        cur.execute(
            "INSERT INTO jobs (job_id, status, created_at, updated_at) VALUES (?, ?, ?, ?)",
            (job_id, JobStatus.pending.value, created_at, created_at),
        )


def update_job_status(
    job_id: str,
    status: JobStatus,
    updated_at: str,
    error: Optional[str] = None,
    formula: Optional[str] = None,
    contrast: Optional[str] = None,
    summary: Optional[dict] = None,
    llm_json: Optional[dict] = None,
) -> None:
    fields = ["status = ?", "updated_at = ?"]
    values: list = [status.value, updated_at]

    if error is not None:
        fields.append("error = ?")
        values.append(error)
    if formula is not None:
        fields.append("formula = ?")
        values.append(formula)
    if contrast is not None:
        fields.append("contrast = ?")
        values.append(contrast)
    if summary is not None:
        fields.append("summary = ?")
        values.append(json.dumps(summary))
    if llm_json is not None:
        fields.append("llm_json = ?")
        values.append(json.dumps(llm_json))

    values.append(job_id)
    with db_cursor() as cur:
        cur.execute(
            f"UPDATE jobs SET {', '.join(fields)} WHERE job_id = ?",
            values,
        )


def get_job(job_id: str) -> Optional[sqlite3.Row]:
    with db_cursor() as cur:
        cur.execute("SELECT * FROM jobs WHERE job_id = ?", (job_id,))
        return cur.fetchone()
