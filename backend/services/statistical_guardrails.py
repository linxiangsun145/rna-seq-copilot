"""Statistical guardrails for lab-grade RNA-seq reliability assessment.

This module is deterministic and rule-based. It extends the existing pipeline
without replacing DESeq2 or introducing black-box models.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PRIMARY_PADJ = 0.05
PRIMARY_ABS_LOG2FC = 1.0
FALLBACK_PADJ = 0.10
FALLBACK_PVALUE = 0.05
FALLBACK_ABS_LOG2FC = 1.0


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in (None, "", "NA", "NaN"):
            return default
        return float(value)
    except Exception:
        return default


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        if value in (None, "", "NA", "NaN"):
            return default
        return int(value)
    except Exception:
        return default


def _load_deg_table(job_dir: Path) -> pd.DataFrame:
    deg_path = job_dir / "results" / "deg_results.csv"
    if not deg_path.exists():
        return pd.DataFrame(columns=["gene_id", "log2FoldChange", "pvalue", "padj"])
    df = pd.read_csv(deg_path)
    for col in ["gene_id", "log2FoldChange", "pvalue", "padj"]:
        if col not in df.columns:
            df[col] = np.nan
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    return df


def validate_expression_matrix_orientation(qc_report: dict[str, Any]) -> list[str]:
    """Validate sample-correlation related outputs from QC layer."""
    issues: list[str] = []
    qcr = qc_report or {}

    group_qc = qcr.get("group_qc") if isinstance(qcr.get("group_qc"), dict) else {}
    means: list[float] = []
    for obj in group_qc.values():
        if not isinstance(obj, dict):
            continue
        if obj.get("mean_correlation") is None:
            continue
        m = _safe_float(obj.get("mean_correlation"), float("nan"))
        if math.isfinite(m):
            means.append(m)

    if means and float(np.mean(means)) < 0:
        issues.append("Correlation computation or data structure may be problematic.")

    qc_warnings = [str(x).lower() for x in (qcr.get("qc_warnings", []) or []) if str(x).strip()]
    if any("implausible control correlation" in x for x in qc_warnings):
        issues.append("Correlation computation or data structure may be problematic.")

    return list(dict.fromkeys(issues))


def run_exploratory_deg_fallback(job_dir: Path, limit: int = 30) -> dict[str, Any]:
    """Return exploratory DEG candidates when primary DEG calls are sparse."""
    df = _load_deg_table(job_dir)
    if df.empty:
        return {
            "primary_deg": 0,
            "exploratory_deg": 0,
            "near_sig_genes": 0,
            "has_any_padj_lt_005": False,
            "strongest_abs_log2fc": 0.0,
            "candidates": [],
        }

    abs_fc = df["log2FoldChange"].abs()
    strongest_abs_log2fc = float(abs_fc.max(skipna=True)) if not abs_fc.empty else 0.0
    primary_mask = (df["padj"] < PRIMARY_PADJ) & (abs_fc > PRIMARY_ABS_LOG2FC)
    near_sig_mask = (df["padj"] < PRIMARY_PADJ) & (abs_fc > 0.5)
    has_any_padj_lt_005 = bool((df["padj"] < PRIMARY_PADJ).fillna(False).any())
    exploratory_mask = ((df["padj"] < FALLBACK_PADJ) | (df["pvalue"] < FALLBACK_PVALUE)) & (abs_fc > FALLBACK_ABS_LOG2FC)

    primary_deg = int(primary_mask.fillna(False).sum())
    near_sig_genes = int(near_sig_mask.fillna(False).sum())
    exploratory = df[exploratory_mask.fillna(False)].copy()
    exploratory_deg = int(exploratory.shape[0])

    if exploratory.empty:
        return {
            "primary_deg": primary_deg,
            "exploratory_deg": 0,
            "near_sig_genes": near_sig_genes,
            "has_any_padj_lt_005": has_any_padj_lt_005,
            "strongest_abs_log2fc": strongest_abs_log2fc,
            "candidates": [],
        }

    exploratory["rank_score"] = exploratory.apply(
        lambda r: float(abs(_safe_float(r.get("log2FoldChange"), 0.0)) * max(-math.log10(max(_safe_float(r.get("pvalue"), 1.0), 1e-300)), 0.0)),
        axis=1,
    )
    exploratory = exploratory.sort_values(["rank_score", "padj", "pvalue"], ascending=[False, True, True])

    candidates: list[dict[str, Any]] = []
    for _, row in exploratory.head(limit).iterrows():
        candidates.append(
            {
                "gene": str(row.get("gene_id", "")),
                "log2fc": _safe_float(row.get("log2FoldChange"), 0.0),
                "padj": _safe_float(row.get("padj"), 1.0),
                "pvalue": _safe_float(row.get("pvalue"), 1.0),
                "rank_score": _safe_float(row.get("rank_score"), 0.0),
                "label": "low-confidence DEG candidate",
            }
        )

    return {
        "primary_deg": primary_deg,
        "exploratory_deg": exploratory_deg,
        "near_sig_genes": near_sig_genes,
        "has_any_padj_lt_005": has_any_padj_lt_005,
        "strongest_abs_log2fc": strongest_abs_log2fc,
        "candidates": candidates,
    }


def detect_deg_status(
    primary_deg: int,
    exploratory_deg: int,
    has_any_padj_lt_005: bool,
    qc_context: dict[str, Any],
    pca_separation: str,
    strongest_candidate_log2fc: float,
) -> str:
    """Classify DEG evidence quality without conflating biology and reliability."""
    if primary_deg > 0:
        return "confirmed"

    clear_pca = str(pca_separation or "").lower() == "clear"
    if bool(has_any_padj_lt_005):
        return "low_effect_size"
    if clear_pca:
        return "low_power"

    return "no_signal"


def compute_realism_stability(n_deg: int, realism_metrics: dict[str, Any]) -> dict[str, Any]:
    """Gate unstable realism-fraction metrics for low-DEG scenarios."""
    gated_metrics = [
        "canonical_fraction_top20",
        "housekeeping_genes_in_top20",
        "top_gene_dominance",
    ]
    metrics = realism_metrics or {}
    if n_deg < 10:
        return {
            "realism_status": "not_applicable",
            "gating_note": "Realism fraction metrics were not evaluated due to insufficient DEG count.",
            "gated_metrics": gated_metrics,
            "is_stable": False,
        }

    canonical = metrics.get("canonical_fraction_top20")
    if canonical is None:
        return {
            "realism_status": "unstable",
            "gating_note": "realism fractions are incomplete or unstable for this run",
            "gated_metrics": ["canonical_fraction_top20"],
            "is_stable": False,
        }

    return {
        "realism_status": "stable",
        "gating_note": "",
        "gated_metrics": [],
        "is_stable": True,
    }


def compute_sparse_ncrna_metrics(
    qc_report: dict[str, Any],
    ncrna_assessment: dict[str, Any],
) -> dict[str, Any]:
    """Compute sparse/ncRNA robustness metrics and conservative warnings."""
    qcr = qc_report or {}
    ncrna = ncrna_assessment or {}
    qcm = qcr.get("qc_metrics") if isinstance(qcr.get("qc_metrics"), dict) else {}
    zero_obj = qcm.get("zero_fraction") if isinstance(qcm.get("zero_fraction"), dict) else {}
    nqc = ncrna.get("qc_metrics") if isinstance(ncrna.get("qc_metrics"), dict) else {}

    zero_fraction_global = _safe_float(zero_obj.get("mean"), 0.0)
    zero_fraction_ncrna = _safe_float(nqc.get("ncrna_zero_fraction"), 0.0)
    low_count_fraction_ncrna = _safe_float(nqc.get("ncrna_low_count_fraction"), 0.0)
    detected_miRNA_count = _safe_int(nqc.get("miRNA_detected_count"), 0)
    detected_lncRNA_count = _safe_int(nqc.get("lncRNA_detected_count"), 0)
    n_ncrna_features = _safe_int(nqc.get("n_ncrna_features"), 0)

    warnings: list[str] = []
    dataset_label = "ncRNA-rich dataset"
    if zero_fraction_ncrna > 0.7 or low_count_fraction_ncrna > 0.7:
        warnings.append("High ncRNA sparsity may limit interpretability.")
        dataset_label = "ncRNA-sparse unreliable dataset"
    if low_count_fraction_ncrna > 0.6:
        warnings.append("lncRNA findings are dominated by low-abundance features.")
    if detected_miRNA_count < 5:
        warnings.append("miRNA signal is weak across most samples.")
    if n_ncrna_features >= 50 and zero_fraction_ncrna < 0.5 and low_count_fraction_ncrna < 0.5:
        dataset_label = "ncRNA-rich dataset"

    return {
        "zero_fraction_global": zero_fraction_global,
        "zero_fraction_ncrna": zero_fraction_ncrna,
        "low_count_fraction_ncrna": low_count_fraction_ncrna,
        "detected_miRNA_count": detected_miRNA_count,
        "detected_lncRNA_count": detected_lncRNA_count,
        "dataset_label": dataset_label,
        "warnings": list(dict.fromkeys(warnings)),
    }


def build_analysis_status(
    qc_report: dict[str, Any],
    deg_status: str,
    realism_status: str,
) -> dict[str, str]:
    qcr = qc_report or {}
    qc_critical = qcr.get("qc_critical", []) if isinstance(qcr.get("qc_critical", []), list) else []
    qc_warnings = qcr.get("qc_warnings", []) if isinstance(qcr.get("qc_warnings", []), list) else []

    if qc_critical:
        qc_status = "critical"
    elif qc_warnings:
        qc_status = "warning"
    else:
        qc_status = "pass"

    rs = str(realism_status or "").lower()
    if rs == "not_applicable":
        normalized_realism_status = "not_applicable"
    elif rs in {"high", "low"}:
        normalized_realism_status = rs
    elif rs == "stable":
        normalized_realism_status = "low"
    else:
        normalized_realism_status = "high"

    return {
        "qc_status": qc_status,
        "deg_status": str(deg_status),
        "realism_status": normalized_realism_status,
    }


def adjust_confidence_for_statistical_validity(
    confidence: dict[str, Any],
    analysis_status: dict[str, str],
    realism_stability: dict[str, Any],
) -> dict[str, Any]:
    """Adjust confidence conservatively when downstream metrics are unstable."""
    out = dict(confidence or {})
    breakdown = out.get("score_breakdown") if isinstance(out.get("score_breakdown"), dict) else {}

    qc_penalty = _safe_float(breakdown.get("qc_penalty"), 0.0)
    design_penalty = _safe_float(breakdown.get("design_penalty"), 0.0)
    realism_penalty = _safe_float(breakdown.get("realism_penalty"), 0.0)

    deg_status = str((analysis_status or {}).get("deg_status", "no_signal"))
    realism_status = str((analysis_status or {}).get("realism_status", "unstable"))
    qc_status = str((analysis_status or {}).get("qc_status", "warning"))

    validation = [str(x) for x in (out.get("confidence_validation", []) or []) if str(x).strip()]
    explanations = [str(x) for x in (out.get("explanations", []) or []) if str(x).strip()]

    if realism_status == "not_applicable":
        if realism_penalty != 0:
            validation.append("Realism penalty forced to 0 because realism_status is not_applicable.")
        realism_penalty = 0.0
        note = "Realism fraction metrics were not evaluated due to insufficient DEG count."
        if note not in explanations:
            explanations.append(note)
    elif realism_status == "low":
        if realism_penalty > 8:
            validation.append("Realism penalty reduced to small range because realism_status is low.")
            realism_penalty = 8.0
    elif realism_status == "high":
        if realism_penalty < 20:
            validation.append("Realism penalty increased to high range because realism_status is high.")
            realism_penalty = 20.0

    if deg_status == "low_power":
        note = "No confirmed DEGs under primary thresholds; interpretation is limited by statistical power."
        if note not in explanations:
            explanations.append(note)

    if deg_status == "failed_detection" or qc_status == "critical":
        note = "Technical reliability concerns are the primary driver of reduced confidence."
        if note not in explanations:
            explanations.append(note)
        # Conservative floor for technically unreliable runs.
        score = max(0, int(round(100 - qc_penalty - design_penalty - realism_penalty)))
        if score > 45:
            qc_penalty += float(score - 45)

    score = max(0, int(round(100 - qc_penalty - design_penalty - realism_penalty)))
    level = "HIGH" if score >= 80 else ("MEDIUM" if score >= 50 else "LOW")

    out["score_breakdown"] = {
        "qc_penalty": float(qc_penalty),
        "design_penalty": float(design_penalty),
        "realism_penalty": float(realism_penalty),
    }
    out["confidence_score"] = score
    out["confidence_level"] = level
    out["explanations"] = explanations
    out["confidence_breakdown_text"] = [
        f"Data Quality Penalty: -{int(round(qc_penalty))}",
        f"Experimental Design Penalty: -{int(round(design_penalty))}",
        f"Data Realism Penalty: -{int(round(realism_penalty))}",
    ]
    if explanations:
        out["confidence_explanation"] = explanations[-1] if str(explanations[-1]).endswith(".") else f"{explanations[-1]}."
    out["confidence_validation"] = list(dict.fromkeys(validation))
    return out


def render_analysis_status_html(analysis_status: dict[str, str]) -> dict[str, Any]:
    status = analysis_status or {}
    return {
        "qc_status": str(status.get("qc_status", "warning")),
        "deg_status": str(status.get("deg_status", "no_signal")),
        "realism_status": str(status.get("realism_status", "high")),
    }
