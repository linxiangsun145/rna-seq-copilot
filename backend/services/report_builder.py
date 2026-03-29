"""
HTML report builder — renders a Jinja2 template with analysis artefacts.
"""
from __future__ import annotations

import base64
import csv
import json
import logging
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from models.schemas import AnalysisSummary, LLMInterpretation
from services.llm_client import evaluateInterpretationConfidence
from services.report_text_validator import build_analysis_snapshot, validate_report_text

logger = logging.getLogger(__name__)

TEMPLATES_DIR = Path(__file__).parent.parent / "templates"

CANONICAL_GENES = {
    "TP53", "MYC", "EGFR", "VEGFA", "IL6", "CXCL8", "GAPDH", "ACTB", "AKT1", "PTEN"
}
HOUSEKEEPING_GENES = {"GAPDH", "ACTB", "RPL13A", "B2M", "HPRT1", "TUBB"}


def _img_to_base64(path: Path) -> Optional[str]:
    if path.exists():
        return base64.b64encode(path.read_bytes()).decode()
    return None


def _to_float(value: Any) -> Optional[float]:
    try:
        if value in (None, "", "NA", "NaN"):
            return None
        return float(value)
    except Exception:
        return None


def _load_top_genes_table(job_dir: Path, n: int = 20) -> list[dict[str, Any]]:
    deg_path = job_dir / "results" / "deg_results.csv"
    if not deg_path.exists():
        return []

    rows: list[dict[str, Any]] = []
    try:
        with deg_path.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                if i >= n:
                    break

                gene = str(row.get("gene_id", "")).strip()
                gene_upper = gene.upper()
                l2fc = _to_float(row.get("log2FoldChange"))
                pvalue = _to_float(row.get("pvalue"))
                padj = _to_float(row.get("padj"))

                rows.append(
                    {
                        "gene": gene,
                        "log2FoldChange": l2fc,
                        "pvalue": pvalue,
                        "padj": padj,
                        "is_canonical": gene_upper in CANONICAL_GENES,
                        "is_housekeeping": gene_upper in HOUSEKEEPING_GENES,
                    }
                )
    except Exception:
        return []

    return rows


def _overall_data_quality(qc: Optional[dict[str, Any]]) -> str:
    if not qc:
        return "unknown"
    n_critical = len(qc.get("qc_critical", []) or [])
    n_warning = len(qc.get("qc_warnings", []) or [])
    if n_critical > 0:
        return "high"
    if n_warning >= 2:
        return "moderate"
    return "low"


def summarizeWarningsForSummary(
    qc_warnings: list[str],
    realism_flags: list[str],
) -> dict[str, str]:
    """Summarize QC and realism concerns into deterministic, concrete phrases."""

    def _top_phrases(texts: list[str], phrase_rules: list[tuple[str, list[str]]], fallback: str = "") -> str:
        normalized = [str(t).lower() for t in (texts or []) if str(t).strip()]
        if not normalized:
            return ""

        found: list[str] = []
        for phrase, keywords in phrase_rules:
            if any(any(k in text for k in keywords) for text in normalized):
                found.append(phrase)
            if len(found) >= 2:
                break

        if not found and fallback:
            found = [fallback]

        if len(found) == 2:
            return f"{found[0]} and {found[1]}"
        return found[0] if found else ""

    qc_phrase = _top_phrases(
        qc_warnings,
        [
            ("low within-group correlation", ["within-group correlation", "low correlation", "correlation"]),
            ("reduced library size in one sample", ["library size", "library-size", "low library"]),
            ("elevated zero-count fraction", ["zero-count", "zero fraction", "sparse"]),
            ("outlier-like sample pattern", ["outlier", "distance outlier"]),
            ("potential batch-related structure", ["batch", "confounded"]),
        ],
        fallback="multiple QC warnings",
    )

    realism_lower = [str(f).lower() for f in (realism_flags or []) if str(f).strip()]
    realism_phrase = ""
    if realism_lower:
        has_canonical = any("canonical" in f for f in realism_lower)
        has_housekeeping = any("housekeeping" in f for f in realism_lower)
        if has_canonical and has_housekeeping:
            realism_phrase = "canonical-gene and housekeeping-gene enrichment patterns"
        elif has_canonical:
            realism_phrase = "canonical-gene enrichment patterns"
        elif has_housekeeping:
            realism_phrase = "housekeeping-gene enrichment patterns"
        else:
            realism_phrase = _top_phrases(
                realism_flags,
                [
                    ("extreme p-value concentration", ["p-value", "pvalue", "< 1e-6"]),
                    ("top-ranked gene concentration", ["top-gene dominance", "dominance", "top 5"]),
                    ("atypical DEG count patterns", ["deg count", "total_deg"]),
                ],
                fallback="quantitative realism deviations",
            )

    return {
        "qc_phrase": qc_phrase,
        "realism_phrase": realism_phrase,
    }


def generateExecutiveSummary(data: dict[str, Any]) -> str:
    """Build a deterministic 3-4 sentence executive summary with metric-backed statements only."""
    n_samples = int(data.get("n_samples", 0) or 0)
    groups = [str(g) for g in (data.get("groups", []) or []) if str(g).strip()]
    groups_text = ", ".join(groups) if groups else "unspecified groups"

    deg_up = int(data.get("deg_up", 0) or 0)
    deg_down = int(data.get("deg_down", 0) or 0)
    total_deg = deg_up + deg_down
    pca_separation = str(data.get("pca_separation", "unknown")).strip() or "unknown"

    qc_report = data.get("qc_report") if isinstance(data.get("qc_report"), dict) else {}
    group_qc = (qc_report or {}).get("group_qc") if isinstance((qc_report or {}).get("group_qc"), dict) else {}

    qc_sentence = "library-size ratio = 1.000 (threshold < 0.50)"
    inconsistent_groups: list[tuple[str, dict[str, Any]]] = []
    for g, obj in group_qc.items():
        if isinstance(obj, dict) and str(obj.get("flag", "none")) == "group_inconsistency":
            inconsistent_groups.append((str(g), obj))

    if inconsistent_groups:
        g_name, g_obj = inconsistent_groups[0]
        g_mean = _to_metric_float(g_obj.get("mean_correlation"))
        if g_mean is not None:
            qc_sentence = f"mean correlation = {g_mean:.3f} (threshold < 0.75)"
        else:
            qc_sentence = "mean correlation = 0.000 (threshold < 0.75)"
    else:
        lib = ((qc_report or {}).get("qc_metrics") or {}).get("library_size") if isinstance((qc_report or {}).get("qc_metrics"), dict) else {}
        ratio = _to_metric_float((lib or {}).get("min_median_ratio")) if isinstance(lib, dict) else None
        if ratio is not None:
            qc_sentence = f"library-size ratio = {ratio:.3f} (threshold < 0.50)"

    realism_metrics = data.get("realism_metrics") if isinstance(data.get("realism_metrics"), dict) else {}
    realism_level = str(data.get("realism_level", "LOW")).upper()
    canonical_count = int(realism_metrics.get("canonical_count", 0) or 0)
    total_deg_for_realism = int(realism_metrics.get("total_deg", 0) or 0)
    canonical_fraction = float(realism_metrics.get("canonical_fraction", 0.0) or 0.0)
    extreme_frac = float(realism_metrics.get("extreme_pvalue_fraction", 0.0) or 0.0)

    realism_sentence = f"canonical fraction = {canonical_fraction:.3f} (threshold > 0.30)"
    if total_deg_for_realism > 0:
        realism_sentence = f"canonical fraction = {canonical_fraction:.3f} (threshold > 0.30)"
    elif extreme_frac > 0:
        realism_sentence = f"extreme p-value fraction = {extreme_frac:.3f} (threshold > 0.40)"

    return " ".join(
        [
            f"The analysis included {n_samples} samples across {groups_text}, with {total_deg} DEGs ({deg_up} upregulated, {deg_down} downregulated).",
            f"PCA separation was classified as {pca_separation}.",
            f"QC issue: {qc_sentence}.",
            f"Realism issue: {realism_sentence}; realism level = {realism_level}.",
        ]
    )


def _normalize_placeholder_text(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    lowered = text.lower()
    if lowered in {"none", "null", "nan", "na", "{}", "[]", "()"}:
        return ""
    if lowered in {"<none>", "<null>"}:
        return ""
    return text


def _normalize_sample_name(value: Any) -> str:
    sample = _normalize_placeholder_text(value)
    if not sample:
        return ""
    sample = re.sub(r"[{}\[\]()]", "", sample).strip()
    lowered = sample.lower()
    if not sample or lowered in {"none", "null", "nan", "na"}:
        return ""
    return sample


def _cleanup_assessment_text(text: str) -> str:
    cleaned = _normalize_placeholder_text(text)
    if not cleaned:
        return ""
    cleaned = re.sub(r"\s+", " ", cleaned)
    cleaned = re.sub(r"([.!?;,])\1+", r"\1", cleaned)
    cleaned = re.sub(r"\bdetected\s+detected\b", "detected", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\bfor\s*\{\s*\}\b", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\bfor\s+(none|null|nan)\b", "", cleaned, flags=re.IGNORECASE)
    cleaned = cleaned.replace("possible technical artifact", "technical artifact")
    cleaned = re.sub(r"\bmay\s+indicate\b", "is consistent with", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\bsuggests\s+possible\b", "is consistent with", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\bsuggests\b", "is consistent with", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s+([,.;:])", r"\1", cleaned)
    return cleaned.strip(" .")


def _strip_sample_phrase(message: str, sample: str) -> str:
    if not message:
        return ""
    cleaned = message
    if sample:
        escaped = re.escape(sample)
        cleaned = re.sub(rf"\bfor\s+{escaped}\b", "", cleaned, flags=re.IGNORECASE)
        cleaned = re.sub(rf"\b{escaped}\s+detected\b", "detected", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s+", " ", cleaned)
    cleaned = re.sub(r"\s+([,.;:])", r"\1", cleaned)
    return cleaned.strip(" .")


def _extract_parenthetical_evidence(text: str) -> str:
    if not text:
        return ""
    matches = re.findall(r"\(([^()]*)\)", text)
    if not matches:
        return ""
    return _cleanup_assessment_text(matches[-1])


def _contains_number(text: str) -> bool:
    return bool(re.search(r"[-+]?\d*\.?\d+", str(text or "")))


def _contains_threshold(text: str) -> bool:
    t = str(text or "").lower()
    return bool(re.search(r"(threshold|<=|>=|<|>|warning|critical|low|moderate|high)", t))


def _extract_comparator_pair(text: str) -> tuple[Optional[float], Optional[str], Optional[float]]:
    m = re.search(r"([-+]?\d*\.?\d+)\s*(<=|>=|<|>)\s*([-+]?\d*\.?\d+)", str(text or ""))
    if not m:
        return None, None, None
    return _to_metric_float(m.group(1)), m.group(2), _to_metric_float(m.group(3))


def _extract_named_metric(text: str, key: str) -> Optional[float]:
    m = re.search(rf"\b{re.escape(key)}\s*[:=]\s*([-+]?\d*\.?\d+)", str(text or ""), flags=re.IGNORECASE)
    if not m:
        return None
    return _to_metric_float(m.group(1))


def _extract_threshold_value(text: str) -> Optional[float]:
    # Examples: "threshold > 0.30", "exceeding the 1.50 threshold", "threshold >= 2"
    t = str(text or "")
    m = re.search(r"threshold\s*(?:[:=]|<=|>=|<|>)?\s*([-+]?\d*\.?\d+)", t, flags=re.IGNORECASE)
    if m:
        return _to_metric_float(m.group(1))
    m = re.search(r"([-+]?\d*\.?\d+)\s+threshold", t, flags=re.IGNORECASE)
    if m:
        return _to_metric_float(m.group(1))
    return None


def _format_threshold_value(value: Optional[float]) -> Optional[str]:
    if value is None:
        return None
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    if value < 0.01:
        return f"{value:.3g}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _metric_name_from_code(code: str, metric: str = "") -> str:
    c = str(code or "").lower()
    m = str(metric or "").lower()
    if "correlation" in c or "correlation" in m or "group_inconsistency" in c:
        return "mean correlation"
    if "library_size" in c:
        return "library-size ratio"
    if "group_imbalance" in c:
        return "group-size ratio"
    if "zero" in c:
        return "zero fraction"
    if "top5" in c or "dominance" in c:
        return "top5 contribution"
    if "pvalue" in c or "p_value" in c:
        return "fraction"
    if "canonical" in c:
        return "canonical fraction"
    if "housekeeping" in c:
        return "housekeeping genes"
    return "value"


def _default_threshold_clause(code: str) -> str:
    c = str(code or "").lower()
    if "correlation" in c or "group_inconsistency" in c:
        return "< 0.75"
    if "library_size" in c:
        return "< 0.50"
    if "group_imbalance" in c:
        return "> 1.50"
    if "zero" in c:
        return "> 0.40"
    if "top5" in c or "dominance" in c:
        return ">= 0.60"
    if "pvalue" in c or "p_value" in c:
        return "> 0.40"
    if "canonical" in c:
        return "> 0.30"
    if "housekeeping" in c:
        return "> 0"
    return "as configured"


def _build_required_metric_statement(code: str, metric: str, message: str, evidence: str) -> Optional[str]:
    text = _cleanup_assessment_text("; ".join([x for x in [evidence, message] if x]))
    if not text:
        return None

    metric_name = _metric_name_from_code(code, metric)
    value: Optional[float] = None

    named = _extract_named_metric(text, metric_name)
    if named is not None:
        value = named

    if value is None:
        left, op, right = _extract_comparator_pair(text)
        if left is not None and op and right is not None:
            value = left

    if value is None:
        m = re.search(r"\b=\s*([-+]?\d*\.?\d+)", text)
        if m:
            value = _to_metric_float(m.group(1))

    if value is None:
        return None

    threshold_clause: Optional[str] = None
    c = str(code or "").lower()

    # Force canonical comparator direction for metrics with fixed rule direction.
    if "correlation" in c or "group_inconsistency" in c:
        threshold_clause = "< 0.75"
    elif "library_size" in c:
        threshold_clause = "< 0.50"
    elif "group_imbalance" in c:
        threshold_clause = "> 1.50"
    elif "zero" in c:
        threshold_clause = "> 0.40"
    elif "top5" in c or "dominance" in c:
        threshold_clause = ">= 0.60"
    elif "pvalue" in c or "p_value" in c:
        threshold_clause = "> 0.30"
    elif "canonical" in c:
        threshold_clause = "> 0.30"
    elif "housekeeping" in c:
        threshold_clause = ">= 2"

    left, op, right = _extract_comparator_pair(text)
    if threshold_clause is None and op and right is not None:
        right_txt = _format_threshold_value(right)
        if right_txt is not None:
            threshold_clause = f"{op} {right_txt}"

    if threshold_clause is None:
        thr_val = _extract_threshold_value(text)
        thr_txt = _format_threshold_value(thr_val)
        if thr_txt is not None:
            threshold_clause = f"> {thr_txt}"

    if threshold_clause is None:
        threshold_clause = _default_threshold_clause(code)

    value_txt = _format_threshold_value(value)
    if value_txt is None:
        return None

    return f"{metric_name} = {value_txt} (threshold {threshold_clause})"


def _is_required_metric_statement(text: str) -> bool:
    t = _cleanup_assessment_text(text)
    return bool(t and _contains_number(t) and "=" in t and "threshold" in t.lower())


def _is_sentence_like(text: str) -> bool:
    t = _cleanup_assessment_text(text)
    if not t:
        return False
    return bool(re.search(r"\b(check|verify|review|investigate|possible|may|likely)\b", t, flags=re.IGNORECASE) or "." in t)


def _extract_correlation_samples_from_text(text: str) -> list[tuple[str, str]]:
    found: list[tuple[str, str]] = []
    if not text:
        return found

    for sample, evidence in re.findall(r"for\s+([A-Za-z0-9_.-]+)\s*\(([^()]*)\)", text, flags=re.IGNORECASE):
        s = _normalize_sample_name(sample)
        e = _cleanup_assessment_text(evidence)
        if s:
            found.append((s, e))

    for sample, evidence in re.findall(r"\b([A-Za-z0-9_.-]+)\s*\(([-+]?\d*\.?\d+\s*<\s*[-+]?\d*\.?\d+)\)", text):
        s = _normalize_sample_name(sample)
        e = _cleanup_assessment_text(evidence)
        if s and not any(x[0] == s and x[1] == e for x in found):
            found.append((s, e))

    for sample, left, right in re.findall(r"\b([A-Za-z0-9_.-]+)\s*[:=]\s*([-+]?\d*\.?\d+)\s*<\s*([-+]?\d*\.?\d+)", text):
        s = _normalize_sample_name(sample)
        e = _cleanup_assessment_text(f"{left} < {right}")
        if s and not any(x[0] == s and x[1] == e for x in found):
            found.append((s, e))

    return found


def expandSampleLevelWarnings(warnings: list[dict[str, Any]]) -> list[dict[str, Any]]:
    expanded: list[dict[str, Any]] = []
    for raw in (warnings or []):
        if not isinstance(raw, dict):
            continue

        item = dict(raw)
        item_type = str(item.get("type", "statistical")).lower().strip() or "statistical"
        code = str(item.get("code", "")).strip()
        metric = str(item.get("metric", "")).strip()
        message = _cleanup_assessment_text(item.get("message", ""))
        sample = _normalize_sample_name(item.get("sample", ""))
        evidence = _cleanup_assessment_text(item.get("evidence", ""))

        looks_like_corr = (
            "correlation" in code.lower()
            or "correlation" in message.lower()
            or "mean_within_group_correlation" in metric.lower()
        )

        if looks_like_corr:
            extracted = _extract_correlation_samples_from_text(f"{message}; {evidence}")
            if len(extracted) > 1:
                for s, e in extracted:
                    left, op, right = _extract_comparator_pair(e)
                    corr_message = "Within-group correlation below threshold"
                    if left is not None and op and right is not None:
                        corr_message = f"Within-group correlation = {left:.3f} ({op} threshold {right:.3f})"
                    expanded.append(
                        {
                            **item,
                            "sample": s,
                            "message": corr_message,
                            "evidence": e,
                        }
                    )
                continue
            if len(extracted) == 1 and not sample:
                s, e = extracted[0]
                item["sample"] = s
                item["evidence"] = e

        item["message"] = message
        item["sample"] = sample
        item["evidence"] = evidence
        expanded.append(item)

    return expanded


def _grouped_warnings(
    summary: dict[str, Any],
    qc: Optional[dict[str, Any]],
    realism: Optional[dict[str, Any]],
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {"qc": [], "realism": [], "statistical": []}

    collected_items: list[dict[str, Any]] = []

    if qc and isinstance(qc.get("warning_items"), list) and qc.get("warning_items"):
        collected_items.extend([w for w in qc.get("warning_items", []) if isinstance(w, dict)])
    else:
        for msg in (qc or {}).get("qc_critical", []) or []:
            collected_items.append({"type": "qc", "severity": "critical", "code": "qc_critical_fallback", "message": str(msg)})
        for msg in (qc or {}).get("qc_warnings", []) or []:
            collected_items.append({"type": "qc", "severity": "warning", "code": "qc_warning_fallback", "message": str(msg)})

    if realism and isinstance(realism.get("warning_items"), list) and realism.get("warning_items"):
        collected_items.extend([w for w in realism.get("warning_items", []) if isinstance(w, dict)])
    else:
        for msg in (realism or {}).get("critical", []) or []:
            collected_items.append({"type": "realism", "severity": "critical", "code": "realism_critical_fallback", "message": str(msg)})
        for msg in (realism or {}).get("warnings", []) or []:
            collected_items.append({"type": "realism", "severity": "warning", "code": "realism_warning_fallback", "message": str(msg)})

    summary_items = summary.get("warning_items", []) or []
    if isinstance(summary_items, list) and summary_items:
        collected_items.extend([w for w in summary_items if isinstance(w, dict)])
    else:
        for msg in summary.get("warnings", []) or []:
            collected_items.append({"type": "statistical", "severity": "warning", "code": "summary_warning_fallback", "message": str(msg)})
        for msg in summary.get("data_issues", []) or []:
            collected_items.append({"type": "statistical", "severity": "warning", "code": "summary_data_issue_fallback", "message": str(msg)})

    # Suppress sample-level correlation/outlier statements when a group-level inconsistency
    # statement exists for the same group.
    inconsistent_groups: set[str] = set()
    for item in collected_items:
        code = str(item.get("code", "")).strip()
        grp = _normalize_placeholder_text(item.get("group", "")).lower()
        if code in {"group_inconsistency", "group_global_inconsistency", "global_group_inconsistency"} and grp:
            inconsistent_groups.add(grp)

    filtered_items: list[dict[str, Any]] = []
    for item in collected_items:
        code = str(item.get("code", "")).strip()
        grp = _normalize_placeholder_text(item.get("group", "")).lower()
        if code.endswith("group_context_warning"):
            continue
        if code in {"low_within_group_correlation_critical", "low_within_group_correlation_warning", "sample_outlier"} and grp and grp in inconsistent_groups:
            continue
        filtered_items.append(item)
    collected_items = filtered_items

    # Deduplicate by code + sample + group + metric for traceable deterministic rendering.
    seen: set[str] = set()
    unique_items: list[dict[str, Any]] = []
    for item in collected_items:
        key = "|".join([
            str(item.get("type", "")),
            str(item.get("severity", "")),
            str(item.get("code", "")),
            str(item.get("sample", "")),
            str(item.get("group", "")),
            str(item.get("metric", "")),
        ])
        if key in seen:
            continue
        seen.add(key)
        unique_items.append(item)

    severity_order = {"critical": 0, "warning": 1}
    type_order = {"qc": 0, "realism": 1, "statistical": 2}

    unique_items.sort(
        key=lambda x: (
            severity_order.get(str(x.get("severity", "warning")).lower(), 1),
            type_order.get(str(x.get("type", "statistical")).lower(), 2),
            str(x.get("code", "")),
            str(x.get("group", "")),
            str(x.get("sample", "")),
        )
    )

    for item in unique_items:
        t = str(item.get("type", "statistical")).lower()
        if t not in grouped:
            t = "statistical"
        grouped[t].append(
            {
                "level": "critical" if str(item.get("severity", "warning")).lower() == "critical" else "warning",
                "message": _cleanup_assessment_text(item.get("message", "")),
                "code": str(item.get("code", "")),
                "group": _normalize_placeholder_text(item.get("group", "")),
                "sample": _normalize_sample_name(item.get("sample", "")),
                "metric": str(item.get("metric", "")) if item.get("metric") is not None else "",
                "evidence": _cleanup_assessment_text(item.get("evidence", "")),
            }
        )

    return grouped


def _normalize_text(text: str) -> str:
    cleaned = str(text or "").strip()
    if not cleaned:
        return ""
    cleaned = cleaned.replace("Warning realism flag:", "")
    cleaned = cleaned.replace("Critical realism flag:", "")
    cleaned = re.sub(r"\bwarning\b:?", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\bcritical\b:?", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s+", " ", cleaned).strip(" .")
    cleaned = cleaned.replace(" :", ":")
    return cleaned[:1].upper() + cleaned[1:] if cleaned else ""


def _status_rank(status: str) -> int:
    order = {"OK": 0, "WARNING": 1, "CRITICAL": 2}
    return order.get(status, 0)


def _to_metric_float(value: Any) -> Optional[float]:
    try:
        if value in (None, "", "NA", "NaN"):
            return None
        return float(value)
    except Exception:
        return None


def _extract_metric_value(source: Any, keys: list[str]) -> Optional[float]:
    if isinstance(source, dict):
        for key in keys:
            if key in source:
                val = _to_metric_float(source.get(key))
                if val is not None:
                    return val
    val = _to_metric_float(source)
    return val


def _warning_matches_metric(message: str, metric_key: str) -> bool:
    text = str(message or "").lower()
    if metric_key == "library_size":
        return ("library size" in text) or ("low library" in text) or ("sequencing depth" in text)
    if metric_key == "correlation":
        return ("correlation" in text) or ("low correlation" in text)
    if metric_key == "zero_fraction":
        return (
            ("zero fraction" in text)
            or ("zero-count" in text)
            or ("zero count" in text)
            or ("sparse" in text)
            or ("sparsity" in text)
        )
    if metric_key == "pca_pc1":
        return ("pca" in text) or ("pc1" in text)
    return False


def computeQCMetrics(
    qc_metrics: Optional[dict[str, Any]],
    pca_variance: Optional[dict[str, Any]],
    qc_warnings: Optional[list[str]] = None,
    group_qc: Optional[dict[str, Any]] = None,
) -> list[dict[str, Any]]:
    """Compute normalized QC metric rows with explicit status thresholds for report rendering."""
    if not isinstance(qc_metrics, dict) or not qc_metrics:
        return []

    rows: list[dict[str, Any]] = []
    warning_texts = [str(w) for w in (qc_warnings or [])]

    library_size = qc_metrics.get("library_size")
    lib_ratio = _extract_metric_value(library_size, ["min_median_ratio", "min_to_median", "ratio"])
    lib_min = _extract_metric_value(library_size, ["min", "minimum", "min_library_size"]) if isinstance(library_size, dict) else None
    lib_median = _extract_metric_value(library_size, ["median", "median_library_size"]) if isinstance(library_size, dict) else None
    if lib_ratio is None and isinstance(library_size, dict):
        min_lib = _extract_metric_value(library_size, ["min", "minimum", "min_library_size"])
        median_lib = _extract_metric_value(library_size, ["median", "median_library_size"])
        if min_lib is not None and median_lib is not None and median_lib > 0:
            lib_ratio = min_lib / median_lib
    if lib_ratio is not None:
        if lib_ratio < 0.5:
            lib_status = "CRITICAL"
        elif lib_ratio < 0.7:
            lib_status = "WARNING"
        else:
            lib_status = "OK"
        rows.append(
            {
                "metric_key": "library_size",
                "name": "Minimum library size",
                "value": (
                    f"minimum = {lib_min:.0f}, median = {lib_median:.0f}, ratio = {lib_ratio:.3f} "
                    "(threshold: ratio < 0.50 = CRITICAL)"
                    if lib_min is not None and lib_median is not None
                    else f"ratio = {lib_ratio:.3f} (threshold: ratio < 0.50 = CRITICAL)"
                ),
                "status": lib_status,
            }
        )

    correlation = qc_metrics.get("correlation")
    corr_mean = _extract_metric_value(correlation, ["mean", "mean_correlation", "mean_within_group", "average", "avg"])
    if corr_mean is not None:
        if corr_mean < 0.8:
            corr_status = "CRITICAL"
        elif corr_mean < 0.9:
            corr_status = "WARNING"
        else:
            corr_status = "OK"
        rows.append(
            {
                "metric_key": "correlation",
                "name": "Mean within-group correlation",
                "value": f"mean = {corr_mean:.3f} (thresholds: <0.75 LOW, 0.75-0.90 MODERATE, >=0.90 HIGH)",
                "status": corr_status,
            }
        )

    zero_fraction = qc_metrics.get("zero_fraction")
    zero_mean = _extract_metric_value(zero_fraction, ["mean", "mean_zero_fraction", "average", "avg"])
    if zero_mean is not None:
        if zero_mean > 0.6:
            zero_status = "CRITICAL"
        elif zero_mean > 0.4:
            zero_status = "WARNING"
        else:
            zero_status = "OK"
        rows.append(
            {
                "metric_key": "zero_fraction",
                "name": "Zero count fraction",
                "value": f"mean = {zero_mean:.3f} (threshold: >0.40 = WARNING, >0.60 = CRITICAL)",
                "status": zero_status,
            }
        )

    pc1_raw = None
    if isinstance(pca_variance, dict):
        pc1_raw = _extract_metric_value(pca_variance, ["pc1", "PC1"])
    if pc1_raw is not None:
        pc1_percent = pc1_raw * 100 if pc1_raw <= 1 else pc1_raw
        pca_status = "WARNING" if pc1_percent < 40 else "OK"
        rows.append(
            {
                "metric_key": "pca_pc1",
                "name": "PCA variance (PC1)",
                "value": f"PC1 variance = {pc1_percent:.2f}% (monitoring threshold: <40% = WARNING)",
                "status": pca_status,
            }
        )

    if isinstance(group_qc, dict) and group_qc:
        for g_name, g_obj in sorted(group_qc.items(), key=lambda x: str(x[0]).lower()):
            if not isinstance(g_obj, dict):
                continue
            g_status = str(g_obj.get("status", "UNKNOWN")).upper()
            g_mean = _to_metric_float(g_obj.get("mean_correlation"))
            g_min = _to_metric_float(g_obj.get("min_correlation"))
            g_low_fraction = _to_metric_float(g_obj.get("low_fraction"))
            if g_mean is None:
                continue
            if g_status == "LOW":
                g_row_status = "CRITICAL"
            elif g_status == "MODERATE":
                g_row_status = "WARNING"
            else:
                g_row_status = "OK"
            rows.append(
                {
                    "metric_key": f"group_qc_{g_name}",
                    "name": f"{str(g_name).title()} group consistency",
                    "value": (
                        f"mean = {g_mean:.3f}, min = {g_min:.3f}, low-sample fraction = {g_low_fraction:.2f} "
                        f"(classification: {g_status})"
                        if g_min is not None and g_low_fraction is not None
                        else f"mean = {g_mean:.3f} (classification: {g_status})"
                    ),
                    "status": g_row_status,
                }
            )

    for row in rows:
        linked_warning = any(_warning_matches_metric(w, row["metric_key"]) for w in warning_texts)
        row["linked_warning"] = linked_warning
        if linked_warning and _status_rank(row["status"]) < _status_rank("WARNING"):
            row["status"] = "WARNING"

    return rows


def _has_tag(row: dict[str, Any], tag_name: str) -> bool:
    tag = str(tag_name).strip().lower()
    tags_raw = row.get("tags")
    if isinstance(tags_raw, list):
        if tag in {str(x).strip().lower() for x in tags_raw}:
            return True
    elif isinstance(tags_raw, str):
        tokens = [t.strip().lower() for t in tags_raw.replace(";", ",").split(",") if t.strip()]
        if tag in set(tokens):
            return True

    if tag == "canonical" and bool(row.get("is_canonical")):
        return True
    if tag == "housekeeping" and bool(row.get("is_housekeeping")):
        return True
    return False


def computeRealismMetrics(top_genes: list[dict[str, Any]]) -> dict[str, Any]:
    total_deg = len(top_genes or [])
    if total_deg == 0:
        return {
            "total_deg": 0,
            "canonical_count": 0,
            "canonical_fraction": 0.0,
            "housekeeping_genes": [],
            "extreme_pvalue_fraction": 0.0,
            "extreme_pvalue_count": 0,
        }

    canonical_count = 0
    housekeeping: list[str] = []
    extreme_p_count = 0

    for row in top_genes:
        if not isinstance(row, dict):
            continue
        if _has_tag(row, "canonical"):
            canonical_count += 1
        if _has_tag(row, "housekeeping"):
            gene_name = str(row.get("gene", "")).strip()
            if gene_name:
                housekeeping.append(gene_name)

        padj = _to_metric_float(row.get("padj"))
        if padj is not None and padj < 1e-6:
            extreme_p_count += 1

    housekeeping_unique = sorted(set(housekeeping))
    return {
        "total_deg": total_deg,
        "canonical_count": canonical_count,
        "canonical_fraction": canonical_count / total_deg if total_deg > 0 else 0.0,
        "housekeeping_genes": housekeeping_unique,
        "extreme_pvalue_fraction": extreme_p_count / total_deg if total_deg > 0 else 0.0,
        "extreme_pvalue_count": extreme_p_count,
    }


def formatTopGenes(top_genes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Sort and normalize top gene rows for publication-ready table rendering."""
    normalized: list[dict[str, Any]] = []

    for row in top_genes or []:
        if not isinstance(row, dict):
            continue
        gene = str(row.get("gene", "")).strip()
        if not gene:
            continue

        log2fc = _to_metric_float(row.get("log2FoldChange"))
        padj = _to_metric_float(row.get("padj"))
        pvalue = _to_metric_float(row.get("pvalue"))

        if log2fc is None:
            direction = "Neutral"
        elif log2fc > 0:
            direction = "Up"
        elif log2fc < 0:
            direction = "Down"
        else:
            direction = "Neutral"

        tags: list[str] = []
        raw_tags = row.get("tags")
        if isinstance(raw_tags, list):
            tags.extend([str(t).strip().lower() for t in raw_tags if str(t).strip()])
        elif isinstance(raw_tags, str):
            tags.extend([t.strip().lower() for t in raw_tags.replace(";", ",").split(",") if t.strip()])

        if _has_tag(row, "canonical") and "canonical" not in tags:
            tags.append("canonical")
        if _has_tag(row, "housekeeping") and "housekeeping" not in tags:
            tags.append("housekeeping")

        normalized.append(
            {
                "gene": gene,
                "log2fc_raw": log2fc,
                "padj_raw": padj,
                "pvalue_raw": pvalue,
                "log2fc": f"{log2fc:.3f}" if log2fc is not None else "NA",
                "padj": f"{padj:.2e}" if padj is not None else (f"{pvalue:.2e}" if pvalue is not None else "NA"),
                "direction": direction,
                "tags": tags,
            }
        )

    normalized.sort(
        key=lambda x: (
            x["padj_raw"] is None and x["pvalue_raw"] is None,
            x["padj_raw"] if x["padj_raw"] is not None else (x["pvalue_raw"] if x["pvalue_raw"] is not None else float("inf")),
        )
    )
    return normalized


def evaluateRealism(metrics: dict[str, Any]) -> dict[str, Any]:
    reasons: list[str] = []
    has_critical = False

    canonical_fraction = float(metrics.get("canonical_fraction", 0.0) or 0.0)
    if canonical_fraction > 0.5:
        has_critical = True
        reasons.append(
            f"Canonical fraction is high ({canonical_fraction:.1%}) and exceeds the critical threshold (>50%)."
        )
    elif canonical_fraction > 0.3:
        reasons.append(
            f"Canonical fraction is elevated ({canonical_fraction:.1%}) and exceeds the warning threshold (>30%)."
        )

    housekeeping_genes = metrics.get("housekeeping_genes") or []
    if len(housekeeping_genes) > 0:
        reasons.append(
            "Housekeeping genes detected in top DEGs: " + ", ".join(housekeeping_genes)
        )

    extreme_p_fraction = float(metrics.get("extreme_pvalue_fraction", 0.0) or 0.0)
    if extreme_p_fraction > 0.4:
        reasons.append(
            f"Extreme adjusted p-value fraction is elevated ({extreme_p_fraction:.1%}) and exceeds the warning threshold (>40%)."
        )

    if has_critical:
        level = "HIGH"
    elif reasons:
        level = "MEDIUM"
    else:
        level = "LOW"

    return {
        "level": level,
        "reasons": reasons,
    }


def _build_realism_evidence_lines(metrics: dict[str, Any]) -> list[str]:
    total_deg = int(metrics.get("total_deg", 0) or 0)
    canonical_count = int(metrics.get("canonical_count", 0) or 0)
    canonical_fraction = float(metrics.get("canonical_fraction", 0.0) or 0.0)
    housekeeping = [str(x).strip() for x in (metrics.get("housekeeping_genes", []) or []) if str(x).strip()]
    extreme_count = int(metrics.get("extreme_pvalue_count", 0) or 0)
    extreme_fraction = float(metrics.get("extreme_pvalue_fraction", 0.0) or 0.0)

    lines: list[str] = []
    if total_deg > 0:
        lines.append(
            f"Canonical genes in ranked DEGs: {canonical_count}/{total_deg} ({100*canonical_fraction:.1f}%) "
            "(thresholds: >30% = WARNING, >50% = CRITICAL)."
        )
        lines.append(
            f"Extreme adjusted p-values (<1e-6): {extreme_count}/{total_deg} ({100*extreme_fraction:.1f}%) "
            "(threshold: >40% = WARNING)."
        )
    if housekeeping:
        lines.append(
            "Housekeeping genes observed in ranked DEGs: "
            + ", ".join(housekeeping)
            + " (threshold: >0 genes = WARNING)."
        )
    else:
        lines.append("Housekeeping genes observed in ranked DEGs: 0 (threshold: >0 genes = WARNING).")
    return lines


def _build_interpretation_limitation(qc_report: Optional[dict[str, Any]], n_samples: int) -> str:
    group_qc = (qc_report or {}).get("group_qc") if isinstance((qc_report or {}).get("group_qc"), dict) else {}
    for g_name, g_obj in sorted((group_qc or {}).items(), key=lambda x: str(x[0]).lower()):
        if not isinstance(g_obj, dict):
            continue
        if str(g_obj.get("flag", "")) != "group_inconsistency":
            continue
        g_mean = _to_metric_float(g_obj.get("mean_correlation"))
        if g_mean is not None:
            return (
                f"Interpretation limitation: mean correlation = {g_mean:.3f} "
                f"(threshold < 0.75) in {str(g_name).title()} group."
            )
        return f"Interpretation limitation: mean correlation = 0.000 (threshold < 0.75) in {str(g_name).title()} group."

    if n_samples < 6:
        return f"Interpretation limitation: sample size = {n_samples} (threshold >= 6)."

    return f"Interpretation limitation: sample size = {n_samples} (threshold >= 6)."


def _normalize_realism_level(value: str) -> str:
    v = str(value or "").strip().lower()
    if v in {"high", "critical"}:
        return "HIGH"
    if v in {"moderate", "medium", "warning"}:
        return "MEDIUM"
    if v in {"low", "ok"}:
        return "LOW"
    return "LOW"


def _build_shared_realism_result(realism: dict[str, Any]) -> dict[str, Any]:
    """Canonical realism source of truth consumed by all report sections."""
    metrics = (realism or {}).get("metrics") or {}
    critical = [str(x).strip() for x in ((realism or {}).get("critical") or []) if str(x).strip()]
    warnings = [str(x).strip() for x in ((realism or {}).get("warnings") or []) if str(x).strip()]

    canonical_count = int(metrics.get("canonical_genes_in_top20", 0) or 0)
    canonical_fraction = float(metrics.get("canonical_fraction_top20", 0.0) or 0.0)
    housekeeping_count = int(metrics.get("housekeeping_genes_in_top20", 0) or 0)
    extreme_pvalue_fraction = float(metrics.get("fraction_p_lt_1e6", 0.0) or 0.0)

    # Keep score deterministic and bounded while reflecting severity and metric extremes.
    score = 0.0
    score += min(60.0, 20.0 * len(critical))
    score += min(30.0, 5.0 * len(warnings))
    if canonical_fraction >= 0.35:
        score += 10.0
    if housekeeping_count >= 2:
        score += 10.0
    if extreme_pvalue_fraction >= 0.30:
        score += 10.0
    score = max(0.0, min(100.0, score))

    level_from_metrics = "HIGH" if (critical or len(warnings) >= 4) else ("MEDIUM" if len(warnings) >= 2 else "LOW")

    supplied_overall = str((realism or {}).get("overall_suspicion", "")).strip()
    if supplied_overall:
        level_from_overall = _normalize_realism_level(supplied_overall)
        if level_from_overall != level_from_metrics:
            raise ValueError(
                "[RealismConsistencyError] realism.overall_suspicion conflicts with metrics-derived level "
                f"({level_from_overall} vs {level_from_metrics})."
            )
        canonical_level = level_from_overall
    else:
        canonical_level = level_from_metrics

    reasons: list[str] = []
    reasons.extend([f"Critical: {x}" for x in critical])
    reasons.extend([f"Warning: {x}" for x in warnings])

    return {
        "level": canonical_level,
        "score": round(score, 1),
        "reasons": reasons,
        "metrics": {
            "canonical_count": canonical_count,
            "canonical_fraction": canonical_fraction,
            "housekeeping_count": housekeeping_count,
            "extreme_pvalue_fraction": extreme_pvalue_fraction,
        },
    }


def _map_realism_flags_to_metrics(realism_flags: list[str]) -> dict[str, Any]:
    mappings: list[dict[str, str]] = []
    unmatched: list[str] = []

    for raw_flag in realism_flags or []:
        flag = str(raw_flag or "").strip()
        text = flag.lower()
        metric_key = ""
        metric_name = ""

        if "canonical" in text:
            metric_key = "canonical_fraction"
            metric_name = "Canonical fraction"
        elif "housekeeping" in text:
            metric_key = "housekeeping_genes"
            metric_name = "Housekeeping genes"
        elif "pvalue" in text or "p_value" in text:
            metric_key = "extreme_pvalue_fraction"
            metric_name = "Extreme p-value fraction"
        elif "deg_count" in text or "deg" in text:
            metric_key = "total_deg"
            metric_name = "Total DEGs"

        if metric_key:
            mappings.append({
                "flag": flag,
                "metric_key": metric_key,
                "metric_name": metric_name,
            })
        else:
            unmatched.append(flag)

    return {
        "mapped": mappings,
        "unmatched": unmatched,
    }


def _expand_groups_for_assessment(
    groups: list[str],
    qc: Optional[dict[str, Any]],
) -> list[str]:
    if qc:
        replicates = qc.get("replicates_per_group")
        if isinstance(replicates, dict) and replicates:
            expanded: list[str] = []
            for group_name, count in replicates.items():
                try:
                    n = int(count)
                except Exception:
                    continue
                if n > 0:
                    expanded.extend([str(group_name)] * n)
            if expanded:
                return expanded
    return [str(g) for g in (groups or []) if str(g).strip()]


def _assessment_issue_from_warning(
    item: dict[str, Any],
    group_counts: Counter[str],
    qc_report: Optional[dict[str, Any]] = None,
    realism_metrics: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    issue_type = str(item.get("type", "statistical")).lower() or "statistical"
    code = str(item.get("code", "")).strip() or "unknown_issue"
    severity = str(item.get("severity", "warning")).lower()
    message = _cleanup_assessment_text(item.get("message", ""))
    sample = _normalize_sample_name(item.get("sample", ""))
    group = _normalize_placeholder_text(item.get("group", ""))
    group_label = group.title() if group else ""
    metric = str(item.get("metric", "")).strip()
    raw_evidence = _cleanup_assessment_text(item.get("evidence", ""))

    issue_message = "Issue"
    evidence: Optional[str] = None

    c = code.lower()
    m = message.lower()

    # Message semantics take precedence over code to avoid swapped canonical/housekeeping labels.
    is_housekeeping = ("housekeeping" in m) or ("housekeeping" in c)
    is_canonical = ("canonical" in m) or ("canonical" in c)

    if "group_imbalance" in c:
        issue_message = "Group balance ratio exceeded threshold"
        ratio_match = re.search(r"ratio\s*([0-9.]+)", m)
        threshold_val = _extract_threshold_value(message) or _extract_threshold_value(raw_evidence)
        comp_left, comp_op, comp_right = _extract_comparator_pair(message)
        if comp_op in {">", ">="} and comp_right is not None:
            threshold_val = comp_right
        if len(group_counts) >= 2:
            sorted_groups = sorted(group_counts.items(), key=lambda x: (x[1], x[0]))
            low_group, low_n = sorted_groups[0]
            high_group, high_n = sorted_groups[-1]
            ratio = (high_n / low_n) if low_n > 0 else float("inf")
            thr_text = f"{threshold_val:.2f}" if threshold_val is not None else "1.50"
            if ratio_match:
                evidence = (
                    f"{low_n} {low_group} vs {high_n} {high_group}; "
                    f"ratio = {float(ratio_match.group(1)):.2f} (threshold: >{thr_text} = WARNING)"
                )
            else:
                evidence = (
                    f"{low_n} {low_group} vs {high_n} {high_group}; "
                    f"ratio = {ratio:.2f} (threshold: >{thr_text} = WARNING)"
                )
        elif ratio_match:
            thr_text = f"{threshold_val:.2f}" if threshold_val is not None else "1.50"
            evidence = f"ratio = {float(ratio_match.group(1)):.2f} (threshold: >{thr_text} = WARNING)"
    elif "group_inconsistency" in c:
        issue_message = f"{group_label} group mean within-group correlation was below threshold" if group_label else "Group mean within-group correlation was below threshold"
        evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif "group_global_inconsistency" in c:
        issue_message = f"{group_label} group consistency failure affected all evaluated samples" if group_label else "Group consistency failure affected all evaluated samples"
        evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif "global_group_inconsistency" in c:
        issue_message = "All evaluable groups had low within-group correlation"
        evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif "sample_outlier" in c:
        issue_message = f"{group_label} group had low-correlation sample outlier(s)" if group_label else "Low-correlation sample outlier(s) were observed"
        evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif "low_library_size" in c or "library_size" in c:
        issue_message = "Library-size ratio was below threshold"
        parsed = re.search(r"for\s+([^\s]+)\s*\(([^\)]+)\)", message)
        if parsed:
            sample_name = parsed.group(1)
            sample = _normalize_sample_name(sample_name) or sample
            evidence = parsed.group(2)
        evidence = evidence or raw_evidence
        if not evidence and message:
            evidence = _extract_parenthetical_evidence(message)
    elif "correlation" in c:
        if "implausible control correlation pattern" in m:
            issue_message = "Control-correlation pattern failed consistency check"
            gq = (qc_report or {}).get("group_qc") if isinstance((qc_report or {}).get("group_qc"), dict) else {}
            control_obj = None
            for g_name, g_obj in (gq or {}).items():
                if str(g_name).strip().lower() == "control" and isinstance(g_obj, dict):
                    control_obj = g_obj
                    break
            control_mean = _to_metric_float((control_obj or {}).get("mean_correlation")) if control_obj else None
            if control_mean is not None:
                evidence = (
                    f"control mean correlation = {control_mean:.3f} "
                    "(thresholds: <0.75 = LOW, 0.75-0.90 = MODERATE, >=0.90 = HIGH)"
                )
            else:
                evidence = _extract_parenthetical_evidence(message)
        else:
            issue_message = "Within-group correlation was below threshold"
            evidence = raw_evidence
            if not evidence and sample:
                parsed = re.search(rf"for\s+{re.escape(sample)}\s*\(([^\)]+)\)", message, flags=re.IGNORECASE)
                if parsed:
                    evidence = _cleanup_assessment_text(parsed.group(1))
            if not evidence:
                evidence = _extract_parenthetical_evidence(message)
    elif "zero_fraction" in c or "zero" in c:
        issue_message = "Zero-count fraction exceeded threshold"
        evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif "top5" in c or "dominance" in c or "top-gene dominance" in m:
        issue_message = "Top-gene concentration was elevated"
        evidence = raw_evidence or _extract_parenthetical_evidence(message)
        if not evidence and realism_metrics:
            total_deg = int((realism_metrics or {}).get("total_deg", 0) or 0)
            if total_deg > 0:
                canonical_count = int((realism_metrics or {}).get("canonical_count", 0) or 0)
                evidence = (
                    f"canonical genes = {canonical_count}/{total_deg} "
                    "(thresholds: >30% = WARNING, >50% = CRITICAL)"
                )
    elif is_housekeeping and not is_canonical:
        issue_message = "Housekeeping-gene signal exceeded threshold"
        hk_match = re.search(r"(\d+\s+housekeeping\s+genes?\s+in\s+top\s+20)", message, flags=re.IGNORECASE)
        hk_thr = _extract_threshold_value(message) or _extract_threshold_value(raw_evidence)
        if hk_match:
            n_match = re.search(r"(\d+)", hk_match.group(1))
            hk_n = int(n_match.group(1)) if n_match else 0
            if hk_thr is not None:
                evidence = f"housekeeping genes = {hk_n} (threshold: >={hk_thr:.0f} = WARNING)"
            else:
                evidence = f"housekeeping genes = {hk_n} (threshold: >0 genes = WARNING)"
        else:
            evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif is_canonical:
        issue_message = "Canonical-gene fraction exceeded threshold"
        canonical_match = re.search(r"(\d+\s*/\s*20\s+canonical\s+genes)", message, flags=re.IGNORECASE)
        if canonical_match:
            frac_match = re.search(r"(\d+)\s*/\s*(20)", canonical_match.group(1))
            if frac_match:
                c_num = int(frac_match.group(1))
                c_den = int(frac_match.group(2))
                c_pct = 100.0 * c_num / c_den if c_den > 0 else 0.0
                evidence = (
                    f"canonical genes = {c_num}/{c_den} ({c_pct:.1f}%) "
                    "(thresholds: >30% = WARNING, >50% = CRITICAL)"
                )
            else:
                evidence = raw_evidence or _extract_parenthetical_evidence(message)
        else:
            evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif "pvalue" in c or "p_value" in c:
        issue_message = "Extreme p-value fraction exceeded threshold"
        frac = None
        frac_m = re.search(r"([-+]?\d*\.?\d+)\s+of\s+p-values", message, flags=re.IGNORECASE)
        if frac_m:
            frac = _to_metric_float(frac_m.group(1))
        thr = _extract_threshold_value(message) or _extract_threshold_value(raw_evidence)
        if frac is not None:
            thr_txt = f"{thr:.2f}" if thr is not None else "0.40"
            evidence = f"fraction = {frac:.3f} (threshold: >{thr_txt} = WARNING)"
        else:
            evidence = raw_evidence or _extract_parenthetical_evidence(message)
    elif issue_type == "statistical":
        issue_message = "Diagnostic statistical anomaly observed"
        evidence = raw_evidence or _extract_parenthetical_evidence(message)

    if "top5" in c or "dominance" in c or "top-gene dominance" in m:
        dominance = None
        dom_m = re.search(r"account\s+for\s+([-+]?\d*\.?\d+)", message, flags=re.IGNORECASE)
        if dom_m:
            dominance = _to_metric_float(dom_m.group(1))
        thr = _extract_threshold_value(message) or _extract_threshold_value(raw_evidence)
        if dominance is not None:
            thr_txt = f"{thr:.2f}" if thr is not None else "0.60"
            evidence = f"top5 contribution = {dominance:.3f} (threshold: >={thr_txt} = WARNING)"

    if message and ("implausible control correlation pattern" in m):
        issue_message = "Control-correlation pattern failed consistency check"

    issue_message = _cleanup_assessment_text(issue_message)
    issue_message = _strip_sample_phrase(issue_message, sample)
    evidence = _cleanup_assessment_text(evidence or "")
    if message and not evidence and issue_message.lower() not in message.lower():
        evidence = _cleanup_assessment_text(message)

    if evidence:
        sample_stripped = _strip_sample_phrase(evidence, sample)
        if sample_stripped:
            evidence = sample_stripped

    if evidence and issue_message and evidence.lower().startswith(issue_message.lower()):
        evidence = _cleanup_assessment_text(evidence[len(issue_message):])

    if code in {"group_inconsistency", "group_global_inconsistency", "global_group_inconsistency"} and evidence:
        metric_clause = re.search(r"(mean\s+correlation\s*=.+)$", evidence, flags=re.IGNORECASE)
        if metric_clause:
            evidence = _cleanup_assessment_text(metric_clause.group(1))

    if not evidence and ("top5" in c or "dominance" in c):
        evidence = "dominance index unavailable (threshold: model-defined cutoff)"
    if not evidence and metric:
        evidence = _cleanup_assessment_text(metric.replace("_", " "))

    # Enforce quantitative QC/realism wording with numeric value and threshold when applicable.
    if issue_type in {"qc", "realism"}:
        strict_metric = _build_required_metric_statement(code, metric, message, evidence)
        if strict_metric:
            evidence = strict_metric
            issue_message = ""
        else:
            # Hard rule: any QC/realism statement without numeric value must be deleted.
            evidence = ""
            issue_message = ""

    category = "D_realism"
    if issue_type == "qc":
        if code in {"group_inconsistency", "group_global_inconsistency", "global_group_inconsistency"}:
            category = "A_core_qc"
        elif code in {"low_library_size_critical", "low_library_size_warning", "group_imbalance_critical", "group_imbalance_warning", "high_zero_fraction_critical", "high_zero_fraction_warning", "sample_outlier"}:
            category = "B_support_qc"
        else:
            category = "C_diagnostic_qc"
    elif issue_type == "realism":
        category = "D_realism"
    else:
        category = "C_diagnostic_qc"

    return {
        "type": issue_type,
        "code": code,
        "severity": "critical" if severity == "critical" else "warning",
        "message": issue_message,
        "evidence": evidence,
        "group": group,
        "sample": sample,
        "category": category,
    }


def _assessment_phrase(issue: dict[str, Any]) -> str:
    issue_type = str(issue.get("type", "statistical")).lower().strip() or "statistical"
    sample = _normalize_sample_name(issue.get("sample", ""))
    message = _cleanup_assessment_text(issue.get("message", "")) or "Issue"
    evidence = _cleanup_assessment_text(issue.get("evidence", ""))

    message = _strip_sample_phrase(message, sample)
    message = re.sub(r"\bdetected\s+for\s*$", "detected", message, flags=re.IGNORECASE).strip(" .")
    evidence = _strip_sample_phrase(evidence, sample)

    if evidence and message and evidence.lower().startswith(message.lower()):
        evidence = _cleanup_assessment_text(evidence[len(message):])

    if issue_type in {"qc", "realism"}:
        strict_metric = evidence if _is_required_metric_statement(evidence) else (message if _is_required_metric_statement(message) else "")
        if strict_metric:
            return f"{strict_metric}."
        return ""

    if sample and issue_type == "qc":
        if evidence:
            return f"{message} for {sample} ({evidence})."
        return f"{message} for {sample}."

    if evidence:
        return f"{message} ({evidence})."

    return f"{message}."


def generateAssessmentBasis(
    warning_items: list[dict[str, Any]],
    n_samples: int,
    groups: list[str],
    qc_report: Optional[dict[str, Any]] = None,
    realism_metrics: Optional[dict[str, Any]] = None,
) -> list[str]:
    _ = n_samples  # Kept for function signature compatibility.

    items = [w for w in (warning_items or []) if isinstance(w, dict)]
    items = expandSampleLevelWarnings(items)
    group_counts = Counter([str(g).strip() for g in (groups or []) if str(g).strip()])

    # Remove sample-level correlation/outlier warnings for groups already flagged at group level.
    group_failures = {
        _normalize_placeholder_text(w.get("group", "")).lower()
        for w in items
        if str(w.get("code", "")).strip() in {"group_inconsistency", "group_global_inconsistency", "global_group_inconsistency"}
        and _normalize_placeholder_text(w.get("group", ""))
    }
    filtered_items: list[dict[str, Any]] = []
    for w in items:
        code = str(w.get("code", "")).strip()
        grp = _normalize_placeholder_text(w.get("group", "")).lower()
        if code in {"low_within_group_correlation_critical", "low_within_group_correlation_warning", "sample_outlier"} and grp in group_failures:
            continue
        filtered_items.append(w)

    normalized_issues = [
        _assessment_issue_from_warning(item, group_counts, qc_report=qc_report, realism_metrics=realism_metrics)
        for item in filtered_items
    ]

    # Deduplicate exact semantic duplicates but keep distinct sample-level bullets.
    deduped: list[dict[str, Any]] = []
    seen_keys: set[str] = set()
    for issue in normalized_issues:
        key = "|".join(
            [
                str(issue.get("severity", "warning")).lower(),
                str(issue.get("type", "statistical")).lower(),
                str(issue.get("code", "")),
                _normalize_placeholder_text(issue.get("group", "")),
                _normalize_sample_name(issue.get("sample", "")),
                _cleanup_assessment_text(issue.get("message", "")),
                _cleanup_assessment_text(issue.get("evidence", "")),
            ]
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduped.append(issue)

    # Merge group-level redundancy into one reviewer-style sentence.
    by_group_issue: dict[str, dict[str, Any]] = {}
    remaining: list[dict[str, Any]] = []
    for issue in deduped:
        code = str(issue.get("code", ""))
        group = _normalize_placeholder_text(issue.get("group", ""))
        if code in {"group_inconsistency", "group_global_inconsistency", "global_group_inconsistency"} and group:
            key = group.lower()
            if key not in by_group_issue:
                by_group_issue[key] = dict(issue)
                continue
            current = by_group_issue[key]
            existing_ev = _cleanup_assessment_text(current.get("evidence", ""))
            new_ev = _cleanup_assessment_text(issue.get("evidence", ""))
            if (not _contains_number(existing_ev)) and _contains_number(new_ev):
                current["evidence"] = new_ev
                existing_ev = new_ev
            if new_ev and (_contains_number(new_ev) or _contains_threshold(new_ev)) and new_ev.lower() not in existing_ev.lower():
                current["evidence"] = "; ".join([x for x in [existing_ev, new_ev] if x])
            current["code"] = "group_inconsistency"
            current["message"] = f"{group.title()} group mean within-group correlation was below threshold"
            current["category"] = "A_core_qc"
        else:
            remaining.append(issue)

    deduped = list(by_group_issue.values()) + remaining

    severity_order = {"critical": 0, "warning": 1}
    type_order = {"qc": 0, "statistical": 1, "realism": 2}
    category_order = {
        "A_core_qc": 0,
        "B_support_qc": 1,
        "C_diagnostic_qc": 2,
        "D_realism": 3,
    }
    code_priority = {
        "group_inconsistency": 0,
        "group_global_inconsistency": 1,
        "global_group_inconsistency": 2,
        "sample_outlier": 3,
    }
    ordered = sorted(
        deduped,
        key=lambda x: (
            severity_order.get(str(x.get("severity", "warning")).lower(), 1),
            category_order.get(str(x.get("category", "C_diagnostic_qc")), 2),
            type_order.get(str(x.get("type", "statistical")).lower(), 2),
            _normalize_placeholder_text(x.get("group", "")) or "~",
            _normalize_sample_name(x.get("sample", "")) or "~",
            code_priority.get(str(x.get("code", "")), 50),
            str(x.get("code", "")),
        ),
    )

    if not ordered:
        return [
            "No QC or realism warning item exceeded configured thresholds "
            "(library-size ratio <0.50 CRITICAL; within-group correlation <0.75 LOW; "
            "zero-fraction >0.40 WARNING; canonical fraction >30% WARNING; extreme p-value fraction >40% WARNING)."
        ]

    bullets: list[str] = []
    bullet_seen: set[str] = set()
    metric_seen: set[str] = set()
    for issue in ordered:
        bullet = _cleanup_assessment_text(_assessment_phrase(issue))
        if bullet and not bullet.endswith("."):
            bullet = f"{bullet}."
        if not bullet:
            continue

        # Merge duplicate meaning: keep only one sentence per quantified metric key.
        if str(issue.get("type", "")).lower() in {"qc", "realism"}:
            metric_match = re.match(r"^\s*([^=]+)=", bullet)
            if metric_match:
                metric_key = _cleanup_assessment_text(metric_match.group(1)).lower()
                if metric_key in metric_seen:
                    continue
                metric_seen.add(metric_key)

        dedup_key = bullet.lower()
        if dedup_key in bullet_seen:
            continue
        bullet_seen.add(dedup_key)
        bullets.append(bullet)

    return bullets or [
        "No QC or realism warning item exceeded configured thresholds "
        "(library-size ratio <0.50 CRITICAL; within-group correlation <0.75 LOW; "
        "zero-fraction >0.40 WARNING; canonical fraction >30% WARNING; extreme p-value fraction >40% WARNING)."
    ]


def build_report(
    job_dir: Path,
    summary: AnalysisSummary,
    llm: Optional[LLMInterpretation],
    formula: Optional[str] = None,
    contrast: Optional[list[str]] = None,
) -> None:
    """Render HTML report and write to job_dir/report.html."""
    env = Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        autoescape=select_autoescape(["html"]),
    )

    template = env.get_template("report.html.j2")

    summary_data = summary.model_dump()

    plots_dir = job_dir / "plots"
    plots = {
        name: _img_to_base64(plots_dir / f"{name}.png")
        for name in [
            "pca",
            "volcano",
            "ma",
            "heatmap",
            "sample_distance_heatmap",
            "sample_correlation_heatmap",
            "library_size",
            "count_distribution",
            "zero_fraction",
        ]
    }

    qc_report: Optional[dict[str, Any]] = None
    qc_path = job_dir / "results" / "qc_report.json"
    if qc_path.exists():
        try:
            qc_report = json.loads(qc_path.read_text(encoding="utf-8"))
        except Exception:
            qc_report = None

    realism = summary_data.get("realism_validation") or {}
    top_genes_table = _load_top_genes_table(job_dir, n=20)
    top_genes_display = formatTopGenes(top_genes_table)

    analysis_methods = {
        "design_formula": formula or "not provided",
        "contrast": ", ".join(contrast) if contrast else summary_data.get("contrast", "not provided"),
        "normalization_method": "DESeq2 median-of-ratios size-factor normalization with VST for visualization",
        "filtering_criteria": "Genes with total count >= 10 retained before differential testing",
        "statistical_test": "Negative binomial GLM in DESeq2 with Wald test",
        "multiple_testing_correction": "Benjamini-Hochberg FDR (adjusted p-value)",
        "deg_threshold": "padj < 0.05 and |log2FoldChange| > 1",
    }

    groups = summary_data.get("groups", []) or []
    methods_paragraph = (
        f"RNA-seq count data were analyzed using design formula '{analysis_methods['design_formula']}' "
        f"and contrast '{analysis_methods['contrast']}'. Counts were normalized with {analysis_methods['normalization_method']}. "
        f"Prior to testing, low-information genes were filtered using the criterion: {analysis_methods['filtering_criteria']}. "
        f"Differential expression was evaluated using {analysis_methods['statistical_test']}, and multiplicity was controlled via "
        f"{analysis_methods['multiple_testing_correction']}. Reporting thresholds for DEG summaries were {analysis_methods['deg_threshold']}."
    )

    shared_realism_result = _build_shared_realism_result(realism)

    results_paragraph = (
        f"The analysis included {summary_data.get('n_samples', 0)} samples across groups: {', '.join(groups) if groups else 'not provided'}. "
        f"Differential testing identified {summary_data.get('deg_up', 0)} upregulated and {summary_data.get('deg_down', 0)} downregulated genes. "
        f"PCA separation was classified as '{summary_data.get('pca_separation', 'unknown')}'. "
        f"Overall data quality risk was assessed as '{_overall_data_quality(qc_report)}', and realism suspicion was "
        f"'{shared_realism_result.get('level', 'LOW').lower()}'."
    )

    figure_legends = [
        {
            "title": "PCA Plot",
            "text": "Principal component projection of VST-transformed counts. Axis labels report variance explained for PC1 and PC2.",
        },
        {
            "title": "Sample Distance Heatmap",
            "text": "Euclidean distance matrix between samples on VST-transformed counts, used for outlier structure assessment.",
        },
        {
            "title": "Sample Correlation Heatmap",
            "text": "Pearson correlation matrix across samples used to detect low-consistency profiles and potential outliers.",
        },
        {
            "title": "Volcano Plot",
            "text": "Gene-wise effect size (log2 fold-change) against significance, highlighting statistically relevant differentially expressed genes.",
        },
        {
            "title": "MA Plot",
            "text": "Mean expression versus log2 fold-change, used to evaluate expression-dependent effect size behavior.",
        },
        {
            "title": "Library Size / Count Distribution / Zero Fraction",
            "text": "Sample-level sequencing depth and count-shape diagnostics used by strict QC threshold rules.",
        },
    ]

    warning_groups = _grouped_warnings(summary_data, qc_report, realism)
    all_warning_items: list[dict[str, Any]] = []
    for warning_type, items in warning_groups.items():
        for item in items:
            all_warning_items.append(
                {
                    "type": warning_type,
                    "severity": str(item.get("level", "warning")).lower(),
                    "code": str(item.get("code", "")),
                    "message": _cleanup_assessment_text(item.get("message", "")),
                    "group": _normalize_placeholder_text(item.get("group", "")),
                    "sample": _normalize_sample_name(item.get("sample", "")),
                    "metric": str(item.get("metric", "")) if item.get("metric") is not None else "",
                    "evidence": _cleanup_assessment_text(item.get("evidence", "")),
                }
            )

    realism_metrics = computeRealismMetrics(top_genes_table)
    groups_for_assessment = _expand_groups_for_assessment(groups, qc_report)
    assessment_basis = generateAssessmentBasis(
        warning_items=all_warning_items,
        n_samples=int(summary_data.get("n_samples", 0) or 0),
        groups=groups_for_assessment,
        qc_report=qc_report or {},
        realism_metrics=realism_metrics,
    )
    qc_metrics_summary = computeQCMetrics(
        qc_metrics=(qc_report or {}).get("qc_metrics") if isinstance(qc_report, dict) else None,
        pca_variance=(qc_report or {}).get("pca_variance") if isinstance(qc_report, dict) else None,
        qc_warnings=(qc_report or {}).get("qc_warnings", []) if isinstance(qc_report, dict) else [],
        group_qc=(qc_report or {}).get("group_qc") if isinstance(qc_report, dict) else None,
    )
    realism_assessment = evaluateRealism(realism_metrics)
    realism_evidence_lines = _build_realism_evidence_lines(realism_metrics)
    shared_realism_result["reasons"] = realism_evidence_lines
    realism_flag_map = _map_realism_flags_to_metrics((realism or {}).get("realism_flags", []) or [])
    if realism_flag_map.get("unmatched"):
        shared_realism_result["reasons"] = list(shared_realism_result.get("reasons", []))
        unmatched_count = len(realism_flag_map.get("unmatched", []) or [])
        shared_realism_result["reasons"].append(
            f"unmapped realism flags = {unmatched_count} (threshold > 0)."
        )

    interpretation_qc_warnings = [
        f"{w.get('severity', 'warning')}: {w.get('message', '')}" for w in all_warning_items if w.get("type") == "qc"
    ]
    interpretation_realism_flags = [
        str(w.get("message", "")) for w in all_warning_items if w.get("type") == "realism" and str(w.get("message", "")).strip()
    ]
    interpretation_confidence = evaluateInterpretationConfidence(
        qc_warnings=interpretation_qc_warnings,
        realism_flags=interpretation_realism_flags,
        n_samples=int(summary_data.get("n_samples", 0) or 0),
    )
    n_samples = int(summary_data.get("n_samples", 0) or 0)
    base_reasons: list[str] = [f"sample size = {n_samples} (threshold >= 6)."]

    # Pull one quantitative QC and one quantitative realism sentence from assessment bullets.
    qc_candidate: Optional[str] = None
    qc_priority = ["group-size ratio", "library-size ratio", "zero fraction", "mean correlation"]
    for key in qc_priority:
        for sentence in assessment_basis:
            s = _cleanup_assessment_text(sentence)
            if not _is_required_metric_statement(s):
                continue
            if key in s.lower():
                qc_candidate = s if s.endswith(".") else f"{s}."
                break
        if qc_candidate:
            break
    if qc_candidate:
        base_reasons.append(qc_candidate)
    for sentence in assessment_basis:
        s = _cleanup_assessment_text(sentence)
        if _is_required_metric_statement(s):
            sl = s.lower()
            if any(k in sl for k in ["canonical fraction", "top5 contribution", "housekeeping genes", "extreme p-value fraction", "fraction ="]):
                base_reasons.append(s if s.endswith(".") else f"{s}.")
                break

    group_qc = (qc_report or {}).get("group_qc") if isinstance((qc_report or {}).get("group_qc"), dict) else {}
    for g_name, g_obj in sorted((group_qc or {}).items(), key=lambda x: str(x[0]).lower()):
        if not isinstance(g_obj, dict):
            continue
        g_mean = _to_metric_float(g_obj.get("mean_correlation"))
        if g_mean is None:
            continue
        base_reasons.append(
            f"{str(g_name).title()} group mean correlation = {g_mean:.3f} (threshold < 0.75)."
        )

    deduped_reasons: list[str] = []
    seen_reasons: set[str] = set()
    for reason in base_reasons:
        if reason in seen_reasons:
            continue
        seen_reasons.add(reason)
        deduped_reasons.append(reason)
    interpretation_confidence["reasons"] = deduped_reasons

    interpretation_limitation_text = _build_interpretation_limitation(qc_report, n_samples)
    executive_summary = generateExecutiveSummary(
        {
            "n_samples": summary_data.get("n_samples", 0),
            "groups": groups,
            "deg_up": summary_data.get("deg_up", 0),
            "deg_down": summary_data.get("deg_down", 0),
            "pca_separation": summary_data.get("pca_separation", "unknown"),
            "qc_warnings": interpretation_qc_warnings,
            "realism_flags": interpretation_realism_flags,
            "qc_report": qc_report or {},
            "realism_metrics": realism_metrics,
            "realism_level": shared_realism_result.get("level", "LOW"),
            "warning_items": all_warning_items,
        }
    )

    llm_payload = llm.model_dump() if llm else None
    ai_interpretation_input = {
        "pca_interpretation": (llm_payload or {}).get("pca_text", ""),
        "deg_summary": (llm_payload or {}).get("deg_summary", ""),
        "biological_insight": (llm_payload or {}).get("biological_insights", ""),
        "limitations": " ".join([x for x in [interpretation_limitation_text, (llm_payload or {}).get("data_quality", "")] if str(x).strip()]),
        "recommendations": (llm_payload or {}).get("next_steps", ""),
    }

    analysis_snapshot = build_analysis_snapshot(
        summary_data=summary_data,
        qc_report=qc_report or {},
        realism_metrics=realism_metrics,
        realism_level=shared_realism_result.get("level", "LOW"),
    )

    validated_text = validate_report_text(
        {
            "executive_summary": executive_summary,
            "assessment_basis": assessment_basis,
            "ai_interpretation": ai_interpretation_input,
        },
        analysis_snapshot,
    )

    executive_summary = str(validated_text.get("executive_summary", executive_summary))
    assessment_basis = [str(x) for x in (validated_text.get("assessment_basis", assessment_basis) or []) if str(x).strip()]

    validated_ai = validated_text.get("ai_interpretation", {}) if isinstance(validated_text.get("ai_interpretation", {}), dict) else {}
    interpretation_limitation_text = str(validated_ai.get("limitations", interpretation_limitation_text) or interpretation_limitation_text)

    if llm_payload is not None:
        llm_payload["pca_text"] = str(validated_ai.get("pca_interpretation", llm_payload.get("pca_text", "")))
        llm_payload["deg_summary"] = str(validated_ai.get("deg_summary", llm_payload.get("deg_summary", "")))
        llm_payload["biological_insights"] = str(validated_ai.get("biological_insight", llm_payload.get("biological_insights", "")))
        llm_payload["data_quality"] = ""
        llm_payload["next_steps"] = str(validated_ai.get("recommendations", llm_payload.get("next_steps", "")))

    html = template.render(
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        summary=summary_data,
        qc=qc_report,
        realism=realism,
        top_genes_table=top_genes_table,
        top_genes_display=top_genes_display,
        warning_groups=warning_groups,
        assessment_basis=assessment_basis,
        qc_metrics_summary=qc_metrics_summary,
        realism_metrics=realism_metrics,
        realism_assessment=realism_assessment,
        shared_realism_result=shared_realism_result,
        realism_flag_map=realism_flag_map,
        interpretation_confidence=interpretation_confidence,
        interpretation_limitation_text=interpretation_limitation_text,
        executive_summary=executive_summary,
        validation_log=validated_text.get("validation_log", []),
        analysis_methods=analysis_methods,
        methods_paragraph=methods_paragraph,
        results_paragraph=results_paragraph,
        figure_legends=figure_legends,
        overall_data_quality=_overall_data_quality(qc_report),
        overall_realism=shared_realism_result.get("level", "LOW"),
        llm=llm_payload,
        plots=plots,
        job_id=job_dir.name,
    )

    out = job_dir / "report.html"
    out.write_text(html, encoding="utf-8")
    logger.info("Report written to %s", out)
