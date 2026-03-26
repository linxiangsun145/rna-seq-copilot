"""
RNA-seq Copilot — FastAPI Backend
Entry point
"""
import logging
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from routers import upload, analysis, results, report
from db.database import init_db
from config import settings

logging.basicConfig(
    level=settings.log_level,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup / shutdown lifecycle."""
    Path(settings.jobs_dir).mkdir(parents=True, exist_ok=True)
    init_db()
    logger.info("RNA-seq Copilot backend started — jobs dir: %s", settings.jobs_dir)
    yield
    logger.info("Shutting down.")


app = FastAPI(
    title="RNA-seq Copilot API",
    version="0.1.0",
    description="Intelligent RNA-seq count matrix analysis backend",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Tighten in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files (plots, reports)
jobs_path = Path(settings.jobs_dir)
jobs_path.mkdir(parents=True, exist_ok=True)
app.mount("/static/jobs", StaticFiles(directory=str(jobs_path)), name="jobs")

# Routers
app.include_router(upload.router, tags=["upload"])
app.include_router(analysis.router, tags=["analysis"])
app.include_router(results.router, tags=["results"])
app.include_router(report.router, tags=["report"])


@app.get("/healthz")
def health():
    return {"status": "ok", "version": "0.1.0"}
