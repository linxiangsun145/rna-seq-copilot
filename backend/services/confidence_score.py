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


def format_confidence_breakdown(score_breakdown: dict[str, Any]) -> list[str]:
    """Render deterministic confidence-penalty breakdown lines."""
    qc = int(round(_safe_float(score_breakdown.get("qc_penalty"), 0.0)))
    design = int(round(_safe_float(score_breakdown.get("design_penalty"), 0.0)))
    realism = int(round(_safe_float(score_breakdown.get("realism_penalty"), 0.0)))
    return [
        f"Data Quality Penalty: -{qc}",
        f"Experimental Design Penalty: -{design}",
        f"Data Realism Penalty: -{realism}",
    ]


def generate_penalty_explanations(
    qc_metrics: dict[str, Any],
    realism_metrics: dict[str, Any],
    flags: dict[str, list[str]],
    score_breakdown: dict[str, Any],
) -> dict[str, list[str]]:
    """Map each non-zero penalty category to 1-2 deterministic causal explanations."""
    out: dict[str, list[str]] = {"qc": [], "design": [], "realism": []}

    qc_penalty = _safe_float(score_breakdown.get("qc_penalty"), 0.0)
    design_penalty = _safe_float(score_breakdown.get("design_penalty"), 0.0)
    realism_penalty = _safe_float(score_breakdown.get("realism_penalty"), 0.0)

    qc_flags = [str(x) for x in (flags.get("qc_flags", []) or []) if str(x).strip()]
    realism_flags = [str(x) for x in (flags.get("realism_flags", []) or []) if str(x).strip()]
    group_sizes = flags.get("group_sizes", {}) if isinstance(flags.get("group_sizes", {}), dict) else {}
    n_samples = _safe_int(flags.get("n_samples"), 0)

    if qc_penalty > 0:
        mean_corr = _safe_float(qc_metrics.get("mean_correlation"), 1.0)
        lib_ratio = _safe_float(qc_metrics.get("min_library_size_ratio"), 1.0)
        zero_fraction = _safe_float(qc_metrics.get("zero_fraction"), 0.0)

        if mean_corr < 0.75:
            out["qc"].append("Low sample correlation reduces data reliability")
        if lib_ratio < 0.5 and len(out["qc"]) < 2:
            out["qc"].append("Library-size imbalance weakens comparability across samples")
        if zero_fraction > 0.6 and len(out["qc"]) < 2:
            out["qc"].append("High zero fraction indicates sparse expression profiles")
        if qc_flags and len(out["qc"]) < 2:
            out["qc"].append("Multiple QC flags indicate elevated technical risk")

    if design_penalty > 0:
        ratio = _group_ratio_from_sizes({str(k): _safe_int(v, 0) for k, v in group_sizes.items()})
        if n_samples < 6:
            out["design"].append("Limited sample size reduces robustness")
        if ratio > 1.5 and len(out["design"]) < 2:
            out["design"].append("Group imbalance weakens statistical power")
        if any(_safe_int(v, 0) < 3 for v in group_sizes.values()) and len(out["design"]) < 2:
            out["design"].append("Insufficient replicates per group reduce robustness")

    if realism_penalty > 0:
        canonical_fraction = _safe_float(realism_metrics.get("canonical_fraction"), 0.0)
        housekeeping_count = _safe_int(realism_metrics.get("housekeeping_count"), 0)
        top_gene_fraction = _safe_float(realism_metrics.get("top_gene_fraction"), 0.0)

        if canonical_fraction > 0.3:
            out["realism"].append("High canonical gene fraction suggests potential bias")
        if top_gene_fraction > 0.6 and len(out["realism"]) < 2:
            out["realism"].append("Dominance of top-ranked genes indicates skewed signal")
        if housekeeping_count >= 2 and len(out["realism"]) < 2:
            out["realism"].append("Housekeeping-gene enrichment indicates non-specific signal")
        if realism_flags and len(out["realism"]) < 2:
            out["realism"].append("Multiple realism flags indicate atypical signal structure")

    return out


def generate_confidence_explanation(
    score_breakdown: dict[str, Any],
    explanations: dict[str, list[str]],
) -> str:
    """Generate one deterministic causal sentence based on largest penalty contributors."""
    def _to_driver_phrase(text: str) -> str:
        t = str(text or "").strip().rstrip(".")
        if not t:
            return ""
        for marker in [" reduces ", " weakens ", " indicates ", " suggests ", " lowers "]:
            idx = t.lower().find(marker)
            if idx > 0:
                return t[:idx].strip().lower()
        return t[0].lower() + t[1:] if len(t) > 1 else t.lower()

    categories = [
        ("qc", _safe_float(score_breakdown.get("qc_penalty"), 0.0), "data-quality issues"),
        ("realism", _safe_float(score_breakdown.get("realism_penalty"), 0.0), "realism concerns"),
        ("design", _safe_float(score_breakdown.get("design_penalty"), 0.0), "design limitations"),
    ]
    ordered = sorted(categories, key=lambda x: (-x[1], x[0]))

    selected: list[str] = []
    for key, val, fallback in ordered:
        if val <= 0:
            continue
        phrase = fallback
        if explanations.get(key):
            first = _to_driver_phrase(explanations[key][0])
            if first:
                phrase = first
        selected.append(phrase)
        if len(selected) >= 3:
            break

    if not selected:
        return "The confidence score is supported by the absence of major QC, design, and realism penalties."

    if len(selected) == 1:
        return f"The confidence score is primarily driven by {selected[0]}."
    if len(selected) == 2:
        return f"The confidence score is primarily driven by {selected[0]} and {selected[1]}."
    return f"The confidence score is primarily driven by {selected[0]}, {selected[1]}, and {selected[2]}."


def inject_confidence_into_summary(summary_text: str, confidence_data: dict[str, Any]) -> str:
    """Append confidence score and causal explanation to summary deterministically."""
    base = str(summary_text or "").strip()
    score = _safe_int(confidence_data.get("confidence_score"), 0)
    level = str(confidence_data.get("confidence_level", "LOW") or "LOW").upper()
    expl = str(confidence_data.get("confidence_explanation", "") or "").strip()

    score_sentence = f"Global confidence score = {score}/100 ({level})."
    parts = [base]
    if score_sentence not in base:
        parts.append(score_sentence)
    if expl and expl not in base:
        parts.append(expl if expl.endswith(".") else f"{expl}.")

    return " ".join([p for p in parts if p]).strip()


def render_confidence_section_html(confidence_data: dict[str, Any]) -> dict[str, Any]:
    """Return deterministic UI-friendly structure for confidence subsection rendering."""
    breakdown = confidence_data.get("score_breakdown", {}) if isinstance(confidence_data.get("score_breakdown", {}), dict) else {}
    ordered_categories = sorted(
        [
            ("qc", _safe_float(breakdown.get("qc_penalty"), 0.0)),
            ("design", _safe_float(breakdown.get("design_penalty"), 0.0)),
            ("realism", _safe_float(breakdown.get("realism_penalty"), 0.0)),
        ],
        key=lambda x: (-x[1], x[0]),
    )

    return {
        "score": _safe_int(confidence_data.get("confidence_score"), 0),
        "level": str(confidence_data.get("confidence_level", "LOW") or "LOW").upper(),
        "breakdown": [str(x) for x in (confidence_data.get("confidence_breakdown_text", []) or []) if str(x).strip()],
        "explanation": str(confidence_data.get("confidence_explanation", "") or "").strip(),
        "penalty_explanations": confidence_data.get("confidence_penalty_explanations", {}),
        "ordered_categories": [k for k, v in ordered_categories if v > 0],
    }


def _validate_breakdown_consistency(confidence_score: int, score_breakdown: dict[str, Any]) -> list[str]:
    logs: list[str] = []
    qc = _safe_float(score_breakdown.get("qc_penalty"), 0.0)
    design = _safe_float(score_breakdown.get("design_penalty"), 0.0)
    realism = _safe_float(score_breakdown.get("realism_penalty"), 0.0)
    expected = max(0, int(round(100 - qc - design - realism)))
    if expected != _safe_int(confidence_score, 0):
        logs.append(
            f"Confidence breakdown mismatch: score={confidence_score}, expected={expected} from penalties (qc={qc}, design={design}, realism={realism})."
        )
    return logs


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

    score_breakdown = {
        "qc_penalty": float(qc_penalty),
        "realism_penalty": float(realism_penalty),
        "design_penalty": float(design_penalty),
    }

    penalty_explanations = generate_penalty_explanations(
        qc_metrics=qc_metrics,
        realism_metrics=realism_metrics,
        flags={
            "qc_flags": qc_flags,
            "realism_flags": realism_flags,
            "group_sizes": groups,
            "n_samples": n_samples,
        },
        score_breakdown=score_breakdown,
    )
    confidence_explanation = generate_confidence_explanation(score_breakdown, penalty_explanations)
    confidence_breakdown_text = format_confidence_breakdown(score_breakdown)
    consistency_logs = _validate_breakdown_consistency(score, score_breakdown)

    return {
        "confidence_score": score,
        "confidence_level": level,
        "score_breakdown": score_breakdown,
        "explanations": generate_confidence_explanations(qc_exp, design_exp, realism_exp),
        "confidence_breakdown_text": confidence_breakdown_text,
        "confidence_penalty_explanations": penalty_explanations,
        "confidence_explanation": confidence_explanation,
        "confidence_validation": consistency_logs,
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
