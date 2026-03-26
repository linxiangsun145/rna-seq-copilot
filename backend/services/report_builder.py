"""
HTML report builder — renders a Jinja2 template with analysis artefacts.
"""
from __future__ import annotations

import base64
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from models.schemas import AnalysisSummary, LLMInterpretation

logger = logging.getLogger(__name__)

TEMPLATES_DIR = Path(__file__).parent.parent / "templates"


def _img_to_base64(path: Path) -> Optional[str]:
    if path.exists():
        return base64.b64encode(path.read_bytes()).decode()
    return None


def build_report(
    job_dir: Path,
    summary: AnalysisSummary,
    llm: Optional[LLMInterpretation],
) -> None:
    """Render HTML report and write to job_dir/report.html."""
    env = Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        autoescape=select_autoescape(["html"]),
    )

    template = env.get_template("report.html.j2")

    plots_dir = job_dir / "plots"
    plots = {
        name: _img_to_base64(plots_dir / f"{name}.png")
        for name in ["pca", "volcano", "ma", "heatmap"]
    }

    html = template.render(
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        summary=summary.model_dump(),
        llm=llm.model_dump() if llm else None,
        plots=plots,
        job_id=job_dir.name,
    )

    out = job_dir / "report.html"
    out.write_text(html, encoding="utf-8")
    logger.info("Report written to %s", out)
