"""Post-generation validator and final publication-style text upgrade layer.

This module runs after deterministic report text generation and before HTML rendering.
It upgrades phrasing style without changing analytical numbers.
"""
from __future__ import annotations

import re
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

ASSESSMENT_GROUP_ORDER = [
    "QC - Group-level",
    "QC - Design",
    "QC - Diagnostic",
    "Realism",
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
            "group-size",
            "zero",
            "outlier",
            "sample",
        ]
    )


def _is_realism_statement(text: str) -> bool:
    t = str(text or "").lower()
    return any(k in t for k in ["realism", "canonical", "housekeeping", "p-value", "pvalue", "top5", "top-gene", "fraction"])


def _format_symbol(op: str) -> str:
    return op.replace(">=", "≥").replace("<=", "≤")


def _replace_metric_label(metric: str, text: str) -> str:
    m = _clean_text(metric).lower()
    t = _clean_text(text).lower()

    mapping = {
        "fraction": "canonical gene fraction",
        "ratio": "group size ratio",
        "top5 contribution": "top 5 gene contribution",
        "top-5 contribution": "top 5 gene contribution",
        "group-size ratio": "group size ratio",
        "library-size ratio": "library size ratio",
        "canonical fraction": "canonical gene fraction",
        "housekeeping genes": "housekeeping genes detected",
        "mean correlation": "mean correlation",
        "zero fraction": "zero-count fraction",
        "extreme p-value fraction": "extreme p-value fraction",
    }

    if m in mapping:
        return mapping[m]

    # Deterministic fallback for ambiguous labels.
    if m == "fraction":
        return "extreme p-value fraction" if "p-value" in t or "pvalue" in t else "canonical gene fraction"
    if m == "ratio":
        return "group size ratio"

    return metric


def upgrade_metric_name(text: str) -> str:
    """Replace generic metric labels with self-explanatory domain labels."""
    t = _clean_text(text)
    if not t or "=" not in t:
        return t

    m = re.match(r"^\s*([^=]+?)\s*=\s*(.+)$", t)
    if not m:
        return t

    raw_metric = _clean_text(m.group(1))
    rest = m.group(2)
    upgraded_metric = _replace_metric_label(raw_metric, t)
    return f"{upgraded_metric} = {rest}"


def _metric_key_from_text(text: str) -> str:
    t = str(text or "").lower()
    if "mean correlation" in t:
        return "mean_correlation"
    if "library size ratio" in t:
        return "library_size_ratio"
    if "group size ratio" in t:
        return "group_size_ratio"
    if "zero-count fraction" in t or "zero fraction" in t:
        return "zero_fraction"
    if "canonical gene fraction" in t or "canonical fraction" in t:
        return "canonical_fraction"
    if "housekeeping genes" in t:
        return "housekeeping_genes"
    if "top 5 gene contribution" in t or "top5 contribution" in t:
        return "top5_contribution"
    if "extreme p-value fraction" in t:
        return "extreme_pvalue_fraction"
    return ""


def rewrite_threshold_phrasing(text: str, metric_code: str | None = None) -> str:
    """Convert threshold wording into publication-style expectation phrasing."""
    t = _clean_text(text)
    if not t:
        return t

    metric_key = metric_code or _metric_key_from_text(t)

    def replacement(match: re.Match[str]) -> str:
        op = match.group(1)
        value = match.group(2)

        # Metric-aware phrasing first.
        if metric_key in {"group_size_ratio", "housekeeping_genes", "top5_contribution"}:
            return f"(warning threshold {_format_symbol(op)} {value})"
        if metric_key in {"mean_correlation", "library_size_ratio"}:
            return f"(expected ≥ {value})"
        if metric_key in {"canonical_fraction", "zero_fraction", "extreme_pvalue_fraction"}:
            return f"(expected ≤ {value})"

        # Generic deterministic inversion fallback.
        invert = {"<": "≥", "<=": ">", ">": "≤", ">=": "<"}
        return f"(expected {invert.get(op, _format_symbol(op))} {value})"

    t = re.sub(r"\(\s*threshold\s*(<=|>=|<|>)\s*([-+]?\d*\.?\d+)\s*\)", replacement, t, flags=re.IGNORECASE)
    return _clean_text(t)


def polish_scientific_phrasing(text: str) -> str:
    """Normalize tone into publication-style scientific English."""
    original = _clean_text(text)
    if not original:
        return ""

    cleaned = original
    replacements = {
        r"\bQC issue\s*:\s*": "Data quality assessment identified ",
        r"\bRealism issue\s*:\s*": "Realism evaluation identified ",
        r"\bissue\s*:\s*": "identified ",
        r"\bwarning\s*:\s*": "observed ",
        r"\bdetected\b": "identified",
        r"\bpossible issue\b": "measured deviation",
        r"\btechnical warning\b": "quantified technical anomaly",
        r"\brealism-related warning\b": "quantified realism anomaly",
        r"\bmay indicate\b": "is consistent with",
        r"\bsuggests\b": "is consistent with",
    }
    for old, new in replacements.items():
        cleaned = re.sub(old, new, cleaned, flags=re.IGNORECASE)

    cleaned = _clean_text(cleaned)
    return cleaned


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
    groups = group_assessment_basis(items)
    ordered: list[str] = []
    for section in ASSESSMENT_GROUP_ORDER:
        ordered.extend(groups.get(section, []))
    return ordered


def group_assessment_basis(items: list[str]) -> dict[str, list[str]]:
    """Deterministically group Assessment Basis items into semantic sections."""
    grouped = {k: [] for k in ASSESSMENT_GROUP_ORDER}

    for raw in items or []:
        t = _clean_text(raw)
        if not t:
            continue
        tl = t.lower()

        if "mean correlation" in tl and ("group" in tl or "control" in tl or "treated" in tl):
            grouped["QC - Group-level"].append(t)
        elif any(k in tl for k in ["group size ratio", "library size ratio"]):
            grouped["QC - Design"].append(t)
        elif any(k in tl for k in ["zero-count fraction", "batch", "pca", "diagnostic", "outlier"]):
            grouped["QC - Diagnostic"].append(t)
        elif any(k in tl for k in ["canonical", "housekeeping", "top 5 gene contribution", "extreme p-value"]):
            grouped["Realism"].append(t)
        elif _is_realism_statement(t):
            grouped["Realism"].append(t)
        else:
            grouped["QC - Diagnostic"].append(t)

    return {k: v for k, v in grouped.items() if v}


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
        except Exception:
            fval = 0.0
        return f"{str(g_name).lower()} group mean correlation = {fval:.3f} (threshold < 0.75)."

    qc_metrics = qc_report.get("qc_metrics") if isinstance(qc_report.get("qc_metrics"), dict) else {}
    lib = qc_metrics.get("library_size") if isinstance(qc_metrics.get("library_size"), dict) else {}
    ratio = lib.get("min_median_ratio")
    try:
        fr = float(ratio)
        return f"library size ratio = {fr:.3f} (threshold < 0.50)."
    except Exception:
        pass

    zero = qc_metrics.get("zero_fraction") if isinstance(qc_metrics.get("zero_fraction"), dict) else {}
    zmean = zero.get("mean")
    try:
        fz = float(zmean)
        return f"zero-count fraction = {fz:.3f} (threshold > 0.40)."
    except Exception:
        return None


def _major_realism_issue_text(analysis_json: dict[str, Any]) -> str | None:
    metrics = analysis_json.get("realism_metrics") if isinstance(analysis_json.get("realism_metrics"), dict) else {}

    canonical_fraction = metrics.get("canonical_fraction")
    try:
        cf = float(canonical_fraction)
        return f"canonical gene fraction = {cf:.3f} (threshold > 0.30)."
    except Exception:
        pass

    hk = metrics.get("housekeeping_genes")
    if isinstance(hk, list):
        hk_count = len([x for x in hk if str(x).strip()])
        return f"housekeeping genes detected = {hk_count} (threshold >= 2)."

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
                return f"group size ratio = {ratio:.2f} (threshold > 1.50); {low_n} {low_g} vs {high_n} {high_g}."

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
            return f"housekeeping genes detected = {count} (threshold >= 2)."
        return None

    if "canonical" in t or "realism" in t or "p-value" in t or "pvalue" in t or "top" in t or "fraction" in t:
        return _major_realism_issue_text(analysis_json)

    return None


def _apply_metric_threshold_upgrade(text: str) -> str:
    t = upgrade_metric_name(text)
    metric_key = _metric_key_from_text(t)
    t = rewrite_threshold_phrasing(t, metric_key)
    t = polish_scientific_phrasing(t)
    return t


def validate_assessment_basis(items: list[str], analysis_json: dict[str, Any]) -> tuple[list[str], dict[str, list[str]], list[str]]:
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

    upgraded = [_apply_metric_threshold_upgrade(x) for x in out]

    deduped = merge_duplicate_statements(upgraded)
    if len(deduped) < len(upgraded):
        logs.append("Merged duplicated assessment statements")

    ordered = reorder_assessment_basis(deduped)
    grouped = group_assessment_basis(ordered)
    if grouped:
        logs.append("Grouped Assessment Basis into semantic sections")

    if ordered != deduped:
        logs.append("Reordered Assessment Basis into QC to realism sequence")

    return ordered, grouped, logs


def _sample_size_groups_from_analysis(analysis_json: dict[str, Any]) -> tuple[int, list[str]]:
    n_samples = int(analysis_json.get("n_samples", 0) or 0)
    groups = [str(g) for g in (analysis_json.get("groups", []) or []) if str(g).strip()]
    return n_samples, groups


def _deg_from_analysis(analysis_json: dict[str, Any]) -> tuple[int, int, int]:
    deg_up = int(analysis_json.get("deg_up", 0) or 0)
    deg_down = int(analysis_json.get("deg_down", 0) or 0)
    return deg_up, deg_down, deg_up + deg_down


def rewrite_executive_summary(text: str, analysis_json: dict[str, Any]) -> str:
    """Rewrite summary into compact publication-style scientific prose."""
    _ = text
    n_samples, groups = _sample_size_groups_from_analysis(analysis_json)
    deg_up, deg_down, total_deg = _deg_from_analysis(analysis_json)
    pca = str(analysis_json.get("pca_separation", "unknown")).strip() or "unknown"

    if len(groups) > 1:
        group_text = " and ".join(groups)
    elif groups:
        group_text = groups[0]
    else:
        group_text = "unspecified"

    qc_finding = _apply_metric_threshold_upgrade(_major_qc_issue_text(analysis_json) or "sample size = 0 (threshold >= 6).")
    realism_finding = _apply_metric_threshold_upgrade(_major_realism_issue_text(analysis_json) or "canonical gene fraction = 0.000 (threshold > 0.30).")

    summary = (
        f"The analysis included {n_samples} samples across {group_text} groups, with {total_deg} differentially expressed genes identified "
        f"({deg_up} upregulated and {deg_down} downregulated). "
        f"PCA separation was classified as {pca}. "
        f"Data quality assessment identified {qc_finding.rstrip('.')}. "
        f"Realism evaluation identified {realism_finding.rstrip('.')}. "
        "These findings are consistent with a hypothesis-generating interpretation rather than a confirmatory conclusion."
    )
    return _clean_text(summary)


def validate_executive_summary(text: str, analysis_json: dict[str, Any]) -> tuple[str, list[str]]:
    logs: list[str] = []
    summary = rewrite_executive_summary(text, analysis_json)

    if any(re.search(p, str(text or ""), flags=re.IGNORECASE) for p in FORBIDDEN_VAGUE_PATTERNS):
        logs.append("Rewrote vague Executive Summary wording into publication-style quantified prose")
    else:
        logs.append("Rebuilt Executive Summary into publication-style quantified prose")

    return summary, logs


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
        out[key] = polish_scientific_phrasing(out.get(key, ""))

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
        out["limitations"] = (
            f"{_apply_metric_threshold_upgrade(major_qc).rstrip('.')} "
            "Interpretation is conservative and intended for hypothesis generation."
        )
        logs.append("Updated AI limitations to explicitly include quantified major QC issue")

    if risk_high:
        prefix = "Conservative interpretation: "
        if out.get("deg_summary") and not out["deg_summary"].lower().startswith("conservative interpretation"):
            out["deg_summary"] = prefix + out["deg_summary"]
            logs.append("Applied conservative language to AI DEG summary due to elevated QC/realism risk")

    for k, v in list(out.items()):
        out[k] = _clean_text(v)

    return out, logs


def validate_report_text(report_text: dict[str, Any], analysis_json: dict[str, Any]) -> dict[str, Any]:
    """Validate and rewrite generated report text deterministically."""
    validation_log: list[str] = []
    cleaned = dict(report_text or {})

    clean_summary, summary_log = validate_executive_summary(cleaned.get("executive_summary", ""), analysis_json)
    cleaned["executive_summary"] = clean_summary
    validation_log.extend(summary_log)

    clean_basis, grouped_basis, basis_log = validate_assessment_basis(cleaned.get("assessment_basis", []) or [], analysis_json)
    cleaned["assessment_basis"] = clean_basis
    cleaned["assessment_basis_grouped"] = grouped_basis
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
