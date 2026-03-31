"""
Analysis router — POST /validate/{job_id}  &  POST /run-analysis/{job_id}
"""
import json
import logging
import re
from pathlib import Path

from fastapi import APIRouter, BackgroundTasks, HTTPException
import pandas as pd

from config import settings
from db.database import get_job, update_job_status
from models.schemas import JobStatus, RunAnalysisRequest, ValidationReport
from services.validator import validate_inputs
from services.r_runner import run_deseq2
from services.llm_client import generate_interpretation
from services.realism_validator import validate_realism
from services.ncrna_analysis import analyze_ncrna_assessment
from services.confidence_score import compute_confidence_from_pipeline
from services.statistical_guardrails import (
    adjust_confidence_for_statistical_validity,
    build_analysis_status,
    compute_realism_stability,
    compute_sparse_ncrna_metrics,
    detect_deg_status,
    run_exploratory_deg_fallback,
    validate_expression_matrix_orientation,
)
from services.report_builder import build_report

logger = logging.getLogger(__name__)
router = APIRouter()


def _read_tabular(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, index_col=0)


def _extract_formula_vars(formula: str) -> list[str]:
    # Keep parsing simple and robust for common formulas (~ a, ~ a + b, ~ a:b).
    tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", formula)
    return [t for t in tokens if t not in {"C"}]


def _resolve_formula_and_contrast(job_dir: Path, formula: str, contrast: list[str]) -> tuple[str, list[str]]:
    meta_path = job_dir / "metadata.csv"
    if not meta_path.exists():
        return formula, contrast

    meta = _read_tabular(meta_path)
    columns = set(meta.columns.astype(str))
    resolved_formula = formula.strip()
    resolved_contrast = list(contrast)

    # If formula variables are absent, try a conservative condition<->group fallback.
    vars_in_formula = _extract_formula_vars(resolved_formula)
    if vars_in_formula and any(v not in columns for v in vars_in_formula):
        if "condition" in vars_in_formula and "condition" not in columns and "group" in columns:
            resolved_formula = resolved_formula.replace("condition", "group")
        elif "group" in vars_in_formula and "group" not in columns and "condition" in columns:
            resolved_formula = resolved_formula.replace("group", "condition")

    # Contrast factor fallback uses the same mapping.
    if len(resolved_contrast) == 3:
        factor, numerator, denominator = resolved_contrast
        if factor not in columns:
            if factor == "condition" and "group" in columns:
                factor = "group"
            elif factor == "group" and "condition" in columns:
                factor = "condition"
            resolved_contrast = [factor, numerator, denominator]

    return resolved_formula, resolved_contrast


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
    job_dir = _job_dir(job_id)
    resolved_formula, resolved_contrast = _resolve_formula_and_contrast(
        job_dir=job_dir,
        formula=request.formula,
        contrast=list(request.contrast),
    )
    contrast_str = f"{resolved_contrast[0]}_{resolved_contrast[1]}_vs_{resolved_contrast[2]}"

    update_job_status(
        job_id,
        JobStatus.running,
        now,
        formula=resolved_formula,
        contrast=contrast_str,
    )

    background_tasks.add_task(
        _run_pipeline,
        job_id=job_id,
        formula=resolved_formula,
        contrast=resolved_contrast,
    )

    return {
        "job_id": job_id,
        "status": JobStatus.running.value,
        "created_at": row["created_at"],
        "updated_at": now,
        "formula": resolved_formula,
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

        logger.info("[%s] Running ncRNA-aware annotation and diagnostics…", job_id)
        ncrna_assessment = analyze_ncrna_assessment(job_dir)
        summary.ncrna_assessment = ncrna_assessment

        fallback = run_exploratory_deg_fallback(job_dir)
        strongest_candidate_fc = float(fallback.get("strongest_abs_log2fc", 0.0) or 0.0)

        provisional_qc_status = "critical" if ((qc_report or {}).get("qc_critical", []) or []) else (
            "warning" if ((qc_report or {}).get("qc_warnings", []) or []) else "pass"
        )
        deg_status = detect_deg_status(
            primary_deg=int(fallback.get("primary_deg", 0) or 0),
            exploratory_deg=int(fallback.get("exploratory_deg", 0) or 0),
            has_any_padj_lt_005=bool(fallback.get("has_any_padj_lt_005", False)),
            qc_context={"qc_status": provisional_qc_status},
            pca_separation=str(summary.pca_separation or ""),
            strongest_candidate_log2fc=strongest_candidate_fc,
        )

        realism_stability = compute_realism_stability(
            n_deg=int(fallback.get("primary_deg", 0) or 0),
            realism_metrics=(realism.model_dump().get("metrics") if realism else {}) or {},
        )
        if realism:
            realism.realism_status = str(realism_stability.get("realism_status", "unstable"))
            realism.gated_metrics = list(realism_stability.get("gated_metrics", []))
            realism.gating_note = str(realism_stability.get("gating_note", ""))

        sparse_ncrna = compute_sparse_ncrna_metrics(qc_report or {}, ncrna_assessment)
        # Merge guardrail warnings into ncRNA section while preserving deterministic order.
        ncrna_warnings = [str(x) for x in (ncrna_assessment.get("warnings", []) or []) if str(x).strip()]
        for extra in sparse_ncrna.get("warnings", []):
            t = str(extra).strip()
            if t and t not in ncrna_warnings:
                ncrna_warnings.append(t)
        ncrna_assessment["warnings"] = ncrna_warnings
        summary.ncrna_assessment = ncrna_assessment

        analysis_status = build_analysis_status(
            qc_report=qc_report or {},
            deg_status=deg_status,
            realism_status=str(realism_stability.get("realism_status", "unstable")),
        )

        summary.deg_status = deg_status
        summary.near_sig_genes = int(fallback.get("near_sig_genes", 0) or 0)
        summary.exploratory_deg_candidates = fallback.get("candidates", [])
        summary.sparse_ncrna_metrics = sparse_ncrna
        summary.analysis_status = analysis_status

        statistical_validation = validate_expression_matrix_orientation(qc_report or {})
        if realism_stability.get("gating_note"):
            statistical_validation.append(str(realism_stability.get("gating_note")))
        summary.statistical_validation = list(dict.fromkeys([str(x) for x in statistical_validation if str(x).strip()]))

        confidence = compute_confidence_from_pipeline(
            summary_dict=summary.model_dump(),
            qc_report=qc_report,
            realism_dict=realism.model_dump() if realism else None,
            ncrna_assessment=ncrna_assessment,
        )
        confidence = adjust_confidence_for_statistical_validity(
            confidence=confidence,
            analysis_status=analysis_status,
            realism_stability=realism_stability,
        )
        summary.confidence_assessment = confidence

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
