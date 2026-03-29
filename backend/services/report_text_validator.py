"""Post-generation validator and rewriter for scientific report text.

This module runs after text generation and before HTML rendering.
It enforces deterministic reporting rules for executive summary,
assessment basis bullets, and AI interpretation blocks.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Any


FORBIDDEN_VAGUE_PATTERNS = [
    r"\bdata quality concerns were identified\b",
    r"\blow correlation detected\b",
    r"\bmay indicate\b",
    r"\bsuggests\b",
    r"\bpossible issue\b",
    r"\btechnical warning\b",
    r"\brealism-related warning\b",
]

SUPPRESSION_PATTERNS = [
    r"sample-level warnings suppressed",
    r"suppressed because group-level inconsistency is primary",
]


def _clean_text(text: Any) -> str:
    t = str(text or "").strip()
    if not t:
        return ""
    t = re.sub(r"\s+", " ", t)
    t = re.sub(r"\s+([,.;:])", r"\1", t)
    return t.strip()


def _contains_number(text: str) -> bool:
    return bool(re.search(r"[-+]?\d*\.?\d+", str(text or "")))


def _contains_threshold(text: str) -> bool:
    return bool(re.search(r"\bthreshold\b|<=|>=|<|>", str(text or ""), flags=re.IGNORECASE))


def _is_quantified_statement(text: str) -> bool:
    t = _clean_text(text)
    return bool(t and "=" in t and _contains_number(t) and _contains_threshold(t))


def _is_qc_statement(text: str) -> bool:
    t = str(text or "").lower()
    return any(
        k in t
        for k in [
            "qc",
            "correlation",
            "library",
            "group imbalance",
            "zero",
            "outlier",
            "sample",
        ]
    )


def _is_realism_statement(text: str) -> bool:
    t = str(text or "").lower()
    return any(k in t for k in ["realism", "canonical", "housekeeping", "p-value", "pvalue", "top5", "top-gene"])


def _group_inconsistency_exists(analysis_json: dict[str, Any]) -> bool:
    qc_report = analysis_json.get("qc_report") if isinstance(analysis_json.get("qc_report"), dict) else {}
    group_qc = qc_report.get("group_qc") if isinstance(qc_report.get("group_qc"), dict) else {}
    for g_obj in group_qc.values():
        if isinstance(g_obj, dict) and str(g_obj.get("flag", "")).strip() == "group_inconsistency":
            return True
    return False


def _major_qc_issue_text(analysis_json: dict[str, Any]) -> str | None:
    qc_report = analysis_json.get("qc_report") if isinstance(analysis_json.get("qc_report"), dict) else {}
    group_qc = qc_report.get("group_qc") if isinstance(qc_report.get("group_qc"), dict) else {}
    for g_name, g_obj in sorted(group_qc.items(), key=lambda x: str(x[0]).lower()):
        if not isinstance(g_obj, dict) or str(g_obj.get("flag", "")) != "group_inconsistency":
            continue
        val = g_obj.get("mean_correlation")
        try:
            fval = float(val)
            return f"{str(g_name).title()} group mean correlation = {fval:.3f} (threshold < 0.75)."
        except Exception:
            return f"{str(g_name).title()} group mean correlation = 0.000 (threshold < 0.75)."

    qc_metrics = qc_report.get("qc_metrics") if isinstance(qc_report.get("qc_metrics"), dict) else {}
    lib = qc_metrics.get("library_size") if isinstance(qc_metrics.get("library_size"), dict) else {}
    ratio = lib.get("min_median_ratio")
    try:
        fr = float(ratio)
        return f"library-size ratio = {fr:.3f} (threshold < 0.50)."
    except Exception:
        pass

    zero = qc_metrics.get("zero_fraction") if isinstance(qc_metrics.get("zero_fraction"), dict) else {}
    zmean = zero.get("mean")
    try:
        fz = float(zmean)
        return f"zero fraction = {fz:.3f} (threshold > 0.40)."
    except Exception:
        return None


def _major_realism_issue_text(analysis_json: dict[str, Any]) -> str | None:
    metrics = analysis_json.get("realism_metrics") if isinstance(analysis_json.get("realism_metrics"), dict) else {}
    canonical_fraction = metrics.get("canonical_fraction")
    try:
        cf = float(canonical_fraction)
        return f"canonical fraction = {cf:.3f} (threshold > 0.30)."
    except Exception:
        pass

    extreme = metrics.get("extreme_pvalue_fraction")
    try:
        ep = float(extreme)
        return f"extreme p-value fraction = {ep:.3f} (threshold > 0.40)."
    except Exception:
        return None


def rewrite_qc_statement(statement: str, analysis_json: dict[str, Any]) -> str | None:
    t = _clean_text(statement).lower()

    if "imbalance" in t:
        qc_report = analysis_json.get("qc_report") if isinstance(analysis_json.get("qc_report"), dict) else {}
        reps = qc_report.get("replicates_per_group") if isinstance(qc_report.get("replicates_per_group"), dict) else {}
        if len(reps) >= 2:
            sorted_groups = sorted(reps.items(), key=lambda x: (int(x[1]), str(x[0])))
            low_g, low_n = sorted_groups[0]
            high_g, high_n = sorted_groups[-1]
            if int(low_n) > 0:
                ratio = float(high_n) / float(low_n)
                return f"group-size ratio = {ratio:.2f} (threshold > 1.50); {low_n} {low_g} vs {high_n} {high_g}."

    if "correlation" in t or "qc" in t or "library" in t or "zero" in t:
        return _major_qc_issue_text(analysis_json)

    return None


def rewrite_realism_statement(statement: str, analysis_json: dict[str, Any]) -> str | None:
    t = _clean_text(statement).lower()
    metrics = analysis_json.get("realism_metrics") if isinstance(analysis_json.get("realism_metrics"), dict) else {}

    if "housekeeping" in t:
        hk = metrics.get("housekeeping_genes")
        if isinstance(hk, list):
            count = len([x for x in hk if str(x).strip()])
            return f"housekeeping genes = {count} (threshold >= 2)."
        return None

    if "canonical" in t or "realism" in t or "p-value" in t or "pvalue" in t or "top" in t:
        return _major_realism_issue_text(analysis_json)

    return None


def merge_duplicate_statements(items: list[str]) -> list[str]:
    merged: list[str] = []
    seen_text: set[str] = set()
    seen_metric: set[str] = set()

    for raw in items or []:
        t = _clean_text(raw)
        if not t:
            continue
        key = re.sub(r"\s+", " ", t.lower()).strip(" .")
        if key in seen_text:
            continue

        metric_match = re.match(r"^\s*([^=]+)=", t)
        if metric_match:
            metric_key = _clean_text(metric_match.group(1)).lower()
            if metric_key in seen_metric:
                continue
            seen_metric.add(metric_key)

        seen_text.add(key)
        merged.append(t if t.endswith(".") else f"{t}.")

    return merged


def reorder_assessment_basis(items: list[str]) -> list[str]:
    buckets = {
        "group_qc": [],
        "support_qc": [],
        "diagnostic": [],
        "realism": [],
    }

    for item in items or []:
        t = _clean_text(item)
        if not t:
            continue
        tl = t.lower()
        if "mean correlation" in tl and ("group" in tl or "control" in tl or "treated" in tl):
            buckets["group_qc"].append(t)
        elif any(k in tl for k in ["library-size", "group-size ratio", "zero fraction", "outlier"]):
            buckets["support_qc"].append(t)
        elif any(k in tl for k in ["diagnostic", "pca", "batch", "system-level"]):
            buckets["diagnostic"].append(t)
        elif any(k in tl for k in ["canonical", "housekeeping", "p-value", "realism", "top5"]):
            buckets["realism"].append(t)
        else:
            buckets["diagnostic"].append(t)

    ordered: list[str] = []
    ordered.extend(buckets["group_qc"])
    ordered.extend(buckets["support_qc"])
    ordered.extend(buckets["diagnostic"])
    ordered.extend(buckets["realism"])
    return ordered


def validate_assessment_basis(items: list[str], analysis_json: dict[str, Any]) -> tuple[list[str], list[str]]:
    logs: list[str] = []
    out: list[str] = []
    group_inconsistent = _group_inconsistency_exists(analysis_json)

    for raw in items or []:
        text = _clean_text(raw)
        if not text:
            continue
        lower = text.lower()

        if group_inconsistent:
            if any(re.search(p, lower) for p in SUPPRESSION_PATTERNS):
                logs.append("Deleted suppression wording tied to sample-level QC statement")
                continue
            if ("sample" in lower or re.search(r"\bctrl[_-]?\d+\b|\btreated[_-]?\d+\b", lower)) and (
                "correlation" in lower or "outlier" in lower
            ):
                logs.append("Deleted sample-level QC warning because group-level inconsistency exists")
                continue

        if _is_qc_statement(text) and not _is_quantified_statement(text):
            rewritten = rewrite_qc_statement(text, analysis_json)
            if rewritten:
                out.append(rewritten)
                logs.append("Rewrote vague QC statement into quantified form")
            else:
                logs.append("Deleted non-quantified QC statement due to missing metric values")
            continue

        if _is_realism_statement(text) and not _is_quantified_statement(text):
            rewritten = rewrite_realism_statement(text, analysis_json)
            if rewritten:
                out.append(rewritten)
                logs.append("Rewrote vague realism statement into quantified form")
            else:
                logs.append("Deleted non-quantified realism statement due to missing metric values")
            continue

        cleaned = text
        for p in FORBIDDEN_VAGUE_PATTERNS:
            cleaned = re.sub(p, "", cleaned, flags=re.IGNORECASE)
        cleaned = _clean_text(cleaned)
        if not cleaned:
            logs.append("Deleted statement containing forbidden vague language")
            continue

        out.append(cleaned)

    deduped = merge_duplicate_statements(out)
    if len(deduped) < len(out):
        logs.append("Merged duplicated assessment statements")

    ordered = reorder_assessment_basis(deduped)
    if ordered != deduped:
        logs.append("Reordered Assessment Basis into QC to realism sequence")

    return ordered, logs


def _sample_size_groups_from_analysis(analysis_json: dict[str, Any]) -> tuple[int, list[str]]:
    n_samples = int(analysis_json.get("n_samples", 0) or 0)
    groups = [str(g) for g in (analysis_json.get("groups", []) or []) if str(g).strip()]
    return n_samples, groups


def _deg_from_analysis(analysis_json: dict[str, Any]) -> tuple[int, int, int]:
    deg_up = int(analysis_json.get("deg_up", 0) or 0)
    deg_down = int(analysis_json.get("deg_down", 0) or 0)
    return deg_up, deg_down, deg_up + deg_down


def validate_executive_summary(text: str, analysis_json: dict[str, Any]) -> tuple[str, list[str]]:
    logs: list[str] = []
    _ = text
    n_samples, groups = _sample_size_groups_from_analysis(analysis_json)
    deg_up, deg_down, total_deg = _deg_from_analysis(analysis_json)
    pca = str(analysis_json.get("pca_separation", "unknown")).strip() or "unknown"
    groups_text = ", ".join(groups) if groups else "unspecified groups"

    qc_issue = _major_qc_issue_text(analysis_json)
    realism_issue = _major_realism_issue_text(analysis_json)

    if not qc_issue:
        qc_issue = "sample size = 0 (threshold >= 6)."
        logs.append("Fallback QC sentence used because no quantitative QC metric was available")
    if not realism_issue:
        realism_issue = "canonical fraction = 0.000 (threshold > 0.30)."
        logs.append("Fallback realism sentence used because no quantitative realism metric was available")

    caution = "Interpretation should be treated as hypothesis-generating due to QC/realism threshold breaches."
    summary = " ".join(
        [
            f"The analysis included {n_samples} samples across {groups_text}, with {total_deg} DEGs ({deg_up} upregulated, {deg_down} downregulated).",
            f"PCA separation was classified as {pca}.",
            f"QC issue: {qc_issue}",
            f"Realism issue: {realism_issue}",
            caution,
        ]
    )

    if any(re.search(p, str(text or ""), flags=re.IGNORECASE) for p in FORBIDDEN_VAGUE_PATTERNS):
        logs.append("Rewrote vague Executive Summary wording into quantified form")
    else:
        logs.append("Rebuilt Executive Summary to enforce compact quantified structure")

    return summary, logs


def _replace_forbidden_language(text: str) -> tuple[str, bool]:
    original = _clean_text(text)
    cleaned = original
    replacements = {
        r"\bmay indicate\b": "is consistent with",
        r"\bsuggests\b": "is consistent with",
        r"\bpossible issue\b": "measured issue",
        r"\btechnical warning\b": "quantified technical anomaly",
        r"\brealism-related warning\b": "quantified realism anomaly",
        r"\bdata quality concerns were identified\b": "quantified QC anomalies were identified",
        r"\blow correlation detected\b": "mean correlation below threshold",
    }
    for old, new in replacements.items():
        cleaned = re.sub(old, new, cleaned, flags=re.IGNORECASE)
    cleaned = _clean_text(cleaned)
    return cleaned, cleaned != original


def _contains_unsupported_biology_claim(text: str, analysis_json: dict[str, Any]) -> bool:
    t = str(text or "").lower()
    if not t:
        return False

    has_mechanism_terms = any(k in t for k in ["pathway", "mechanism", "signaling", "immune", "metabolic", "causal"])
    if not has_mechanism_terms:
        return False

    top_genes = [str(g).lower() for g in (analysis_json.get("top_genes", []) or []) if str(g).strip()]
    if not top_genes:
        return True

    return not any(g in t for g in top_genes)


def validate_ai_interpretation(section: dict[str, str], analysis_json: dict[str, Any]) -> tuple[dict[str, str], list[str]]:
    logs: list[str] = []
    out = dict(section or {})

    for key in ["pca_interpretation", "deg_summary", "biological_insight", "limitations", "recommendations"]:
        cleaned, changed = _replace_forbidden_language(out.get(key, ""))
        out[key] = cleaned
        if changed:
            logs.append(f"Rewrote forbidden vague language in AI interpretation field: {key}")

    risk_high = False
    qc_report = analysis_json.get("qc_report") if isinstance(analysis_json.get("qc_report"), dict) else {}
    if len(qc_report.get("qc_critical", []) or []) > 0:
        risk_high = True
    realism_level = str(analysis_json.get("realism_level", "LOW")).upper()
    if realism_level in {"HIGH", "MEDIUM"}:
        risk_high = True

    if _contains_unsupported_biology_claim(out.get("biological_insight", ""), analysis_json):
        out["biological_insight"] = "Observed expression patterns are reported at the gene-level only; pathway or mechanism claims are not supported by the provided inputs."
        logs.append("Removed unsupported biological pathway/mechanism claim")

    major_qc = _major_qc_issue_text(analysis_json)
    if major_qc:
        out["limitations"] = f"{major_qc} Interpretation is conservative and intended for hypothesis generation."
        logs.append("Updated AI limitations to explicitly include quantified major QC issue")

    if risk_high:
        prefix = "Conservative interpretation: "
        if out.get("deg_summary") and not out["deg_summary"].lower().startswith("conservative interpretation"):
            out["deg_summary"] = prefix + out["deg_summary"]
            logs.append("Applied conservative language to AI DEG summary due to elevated QC/realism risk")

    return out, logs


def validate_report_text(report_text: dict[str, Any], analysis_json: dict[str, Any]) -> dict[str, Any]:
    """Validate and rewrite generated report text deterministically.

    Input schema example:
    {
      "executive_summary": str,
      "assessment_basis": [str, ...],
      "ai_interpretation": {
        "pca_interpretation": str,
        "deg_summary": str,
        "biological_insight": str,
        "limitations": str,
        "recommendations": str
      }
    }
    """
    validation_log: list[str] = []
    cleaned = dict(report_text or {})

    clean_summary, summary_log = validate_executive_summary(cleaned.get("executive_summary", ""), analysis_json)
    cleaned["executive_summary"] = clean_summary
    validation_log.extend(summary_log)

    clean_basis, basis_log = validate_assessment_basis(cleaned.get("assessment_basis", []) or [], analysis_json)
    cleaned["assessment_basis"] = clean_basis
    validation_log.extend(basis_log)

    ai_section = cleaned.get("ai_interpretation", {})
    if not isinstance(ai_section, dict):
        ai_section = {}
    clean_ai, ai_log = validate_ai_interpretation(ai_section, analysis_json)
    cleaned["ai_interpretation"] = clean_ai
    validation_log.extend(ai_log)

    cleaned["validation_log"] = merge_duplicate_statements(validation_log)
    return cleaned


def build_analysis_snapshot(summary_data: dict[str, Any], qc_report: dict[str, Any], realism_metrics: dict[str, Any], realism_level: str) -> dict[str, Any]:
    """Create a deterministic analysis snapshot consumed by the text validator."""
    return {
        "n_samples": int(summary_data.get("n_samples", 0) or 0),
        "groups": [str(g) for g in (summary_data.get("groups", []) or []) if str(g).strip()],
        "deg_up": int(summary_data.get("deg_up", 0) or 0),
        "deg_down": int(summary_data.get("deg_down", 0) or 0),
        "pca_separation": str(summary_data.get("pca_separation", "unknown")),
        "top_genes": [str(g) for g in (summary_data.get("top_genes", []) or []) if str(g).strip()],
        "qc_report": qc_report or {},
        "realism_metrics": realism_metrics or {},
        "realism_level": str(realism_level or "LOW").upper(),
    }
