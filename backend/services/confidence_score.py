"""Deterministic confidence scoring for RNA-seq report reliability.

This module is rule-based and fully explainable. No machine learning is used.
"""
from __future__ import annotations

import re
from typing import Any


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


def _group_ratio_from_sizes(group_sizes: dict[str, int]) -> float:
    vals = [int(v) for v in (group_sizes or {}).values() if _safe_int(v, -1) > 0]
    if len(vals) < 2:
        return 1.0
    low = min(vals)
    high = max(vals)
    if low <= 0:
        return 1.0
    return float(high) / float(low)


def _extract_top_gene_fraction(realism_metrics: dict[str, Any], realism_flags: list[str]) -> float:
    if "top_gene_fraction" in realism_metrics:
        return _safe_float(realism_metrics.get("top_gene_fraction"), 0.0)

    for f in realism_flags or []:
        m = re.search(r"top\s*5[^\d]*([0-9]*\.?[0-9]+)", str(f or ""), flags=re.IGNORECASE)
        if m:
            return _safe_float(m.group(1), 0.0)
        m = re.search(r"([0-9]*\.?[0-9]+)\s+of\s+total\s+-log10", str(f or ""), flags=re.IGNORECASE)
        if m:
            return _safe_float(m.group(1), 0.0)

    return 0.0


def compute_qc_penalty(qc_metrics: dict[str, Any], qc_flags: list[str]) -> tuple[float, list[str]]:
    penalty = 0.0
    explanations: list[str] = []

    mean_corr = _safe_float(qc_metrics.get("mean_correlation"), 1.0)
    if mean_corr < 0.5:
        penalty += 25
        explanations.append("Low sample correlation reduces confidence")
    elif mean_corr < 0.75:
        penalty += 15
        explanations.append("Moderate sample correlation reduces confidence")

    lib_ratio = _safe_float(qc_metrics.get("min_library_size_ratio"), 1.0)
    if lib_ratio < 0.3:
        penalty += 15
        explanations.append("Severe library-size imbalance reduces confidence")
    elif lib_ratio < 0.5:
        penalty += 10
        explanations.append("Library-size imbalance reduces confidence")

    zero_fraction = _safe_float(qc_metrics.get("zero_fraction"), 0.0)
    if zero_fraction > 0.6:
        penalty += 10
        explanations.append("High zero-count fraction reduces confidence")

    flag_penalty = min(len(qc_flags or []) * 5, 20)
    if flag_penalty > 0:
        penalty += flag_penalty
        explanations.append("Accumulated QC flags reduce confidence")

    return penalty, explanations


def compute_design_penalty(n_samples: int, groups: dict[str, int]) -> tuple[float, list[str]]:
    penalty = 0.0
    explanations: list[str] = []

    if n_samples < 4:
        penalty += 25
        explanations.append("Very small sample size reduces confidence")
    elif n_samples < 6:
        penalty += 15
        explanations.append("Limited sample size reduces confidence")

    ratio = _group_ratio_from_sizes(groups)
    if ratio > 2.0:
        penalty += 15
        explanations.append("Strong group imbalance weakens statistical reliability")
    elif ratio > 1.5:
        penalty += 10
        explanations.append("Group imbalance weakens statistical reliability")

    if any(int(v) < 3 for v in (groups or {}).values() if _safe_int(v, -1) >= 0):
        penalty += 10
        explanations.append("Insufficient replicates per group reduce confidence")

    return penalty, explanations


def compute_realism_penalty(realism_metrics: dict[str, Any], realism_flags: list[str]) -> tuple[float, list[str]]:
    penalty = 0.0
    explanations: list[str] = []

    canonical_fraction = _safe_float(realism_metrics.get("canonical_fraction"), 0.0)
    if canonical_fraction > 0.5:
        penalty += 25
        explanations.append("High canonical gene fraction suggests potential bias")
    elif canonical_fraction > 0.3:
        penalty += 15
        explanations.append("Elevated canonical gene fraction lowers confidence")

    housekeeping_count = _safe_int(realism_metrics.get("housekeeping_count"), 0)
    if housekeeping_count >= 2:
        penalty += 10
        explanations.append("Housekeeping-gene enrichment lowers confidence")

    top_gene_fraction = _safe_float(realism_metrics.get("top_gene_fraction"), 0.0)
    if top_gene_fraction > 0.6:
        penalty += 10
        explanations.append("Top-gene dominance lowers confidence")

    flag_penalty = min(len(realism_flags or []) * 5, 20)
    if flag_penalty > 0:
        penalty += flag_penalty
        explanations.append("Accumulated realism flags reduce confidence")

    return penalty, explanations


def generate_confidence_explanations(
    qc_explanations: list[str],
    design_explanations: list[str],
    realism_explanations: list[str],
) -> list[str]:
    out: list[str] = []
    for text in [*(qc_explanations or []), *(design_explanations or []), *(realism_explanations or [])]:
        t = str(text or "").strip()
        if t and t not in out:
            out.append(t)
    return out


def compute_confidence_score(
    n_samples: int,
    groups: dict[str, int],
    qc_metrics: dict[str, Any],
    qc_flags: list[str],
    realism_metrics: dict[str, Any],
    realism_flags: list[str],
) -> dict[str, Any]:
    qc_penalty, qc_exp = compute_qc_penalty(qc_metrics=qc_metrics, qc_flags=qc_flags)
    design_penalty, design_exp = compute_design_penalty(n_samples=n_samples, groups=groups)
    realism_penalty, realism_exp = compute_realism_penalty(realism_metrics=realism_metrics, realism_flags=realism_flags)

    score = max(0, int(round(100 - qc_penalty - design_penalty - realism_penalty)))
    if score >= 80:
        level = "HIGH"
    elif score >= 50:
        level = "MEDIUM"
    else:
        level = "LOW"

    return {
        "confidence_score": score,
        "confidence_level": level,
        "score_breakdown": {
            "qc_penalty": float(qc_penalty),
            "realism_penalty": float(realism_penalty),
            "design_penalty": float(design_penalty),
        },
        "explanations": generate_confidence_explanations(qc_exp, design_exp, realism_exp),
    }


def compute_confidence_from_pipeline(
    summary_dict: dict[str, Any],
    qc_report: dict[str, Any] | None,
    realism_dict: dict[str, Any] | None,
) -> dict[str, Any]:
    """Adapter that maps existing pipeline outputs into confidence-score input schema."""
    qcr = qc_report or {}
    realism = realism_dict or {}
    realism_metrics_obj = realism.get("metrics") if isinstance(realism.get("metrics"), dict) else {}

    group_qc = qcr.get("group_qc") if isinstance(qcr.get("group_qc"), dict) else {}
    mean_corr = None
    for obj in group_qc.values():
        if isinstance(obj, dict) and obj.get("mean_correlation") is not None:
            mean_corr = _safe_float(obj.get("mean_correlation"), None)
            break

    qc_metrics_obj = qcr.get("qc_metrics") if isinstance(qcr.get("qc_metrics"), dict) else {}
    library_obj = qc_metrics_obj.get("library_size") if isinstance(qc_metrics_obj.get("library_size"), dict) else {}
    zero_obj = qc_metrics_obj.get("zero_fraction") if isinstance(qc_metrics_obj.get("zero_fraction"), dict) else {}

    qc_metrics = {
        "mean_correlation": 1.0 if mean_corr is None else mean_corr,
        "min_library_size_ratio": _safe_float(library_obj.get("min_median_ratio"), 1.0),
        "zero_fraction": _safe_float(zero_obj.get("mean"), 0.0),
    }

    reps = qcr.get("replicates_per_group") if isinstance(qcr.get("replicates_per_group"), dict) else {}
    if not reps:
        groups = {str(g): 0 for g in (summary_dict.get("groups", []) or [])}
    else:
        groups = {str(k): _safe_int(v, 0) for k, v in reps.items()}

    qc_flags = [str(x) for x in (qcr.get("qc_critical", []) or []) if str(x).strip()]
    qc_flags.extend([str(x) for x in (qcr.get("qc_warnings", []) or []) if str(x).strip()])

    realism_flags = [str(x) for x in (realism.get("critical", []) or []) if str(x).strip()]
    realism_flags.extend([str(x) for x in (realism.get("warnings", []) or []) if str(x).strip()])

    realism_metrics = {
        "canonical_fraction": _safe_float(realism_metrics_obj.get("canonical_fraction_top20"), 0.0),
        "housekeeping_count": _safe_int(realism_metrics_obj.get("housekeeping_genes_in_top20"), 0),
        "top_gene_fraction": _extract_top_gene_fraction(realism_metrics_obj, realism_flags),
    }

    return compute_confidence_score(
        n_samples=_safe_int(summary_dict.get("n_samples"), 0),
        groups=groups,
        qc_metrics=qc_metrics,
        qc_flags=qc_flags,
        realism_metrics=realism_metrics,
        realism_flags=realism_flags,
    )
