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
from services.llm_client import evaluateInterpretationConfidence, generateLimitationText

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


def summarizeWarningsForSummary(warning_items: list[dict[str, Any]]) -> str:
    """Compress warning items into 1-2 concrete phrases for executive summary."""
    warnings = [w for w in (warning_items or []) if str(w.get("type", "")).lower() in {"qc", "realism"}]
    if not warnings:
        return ""

    phrases: list[str] = []
    seen: set[str] = set()

    for item in warnings:
        code = str(item.get("code", "")).lower()
        text = str(item.get("message", "")).lower()
        if "correlation" in code or "correlation" in text:
            phrase = "low sample correlation"
        elif "library" in code or "library" in text:
            phrase = "low library-size pattern"
        elif "zero" in code or "zero fraction" in text or "sparse" in text:
            phrase = "high zero-count fraction"
        elif "outlier" in code or "outlier" in text:
            phrase = "sample outlier pattern"
        elif "batch" in code or "batch" in text:
            phrase = "possible batch structure"
        elif "canonical" in code or "housekeeping" in code:
            phrase = "canonical-gene and housekeeping-gene enrichment patterns"
        elif "pvalue" in code:
            phrase = "extreme p-value concentration"
        else:
            phrase = "additional technical warning"

        if phrase not in seen:
            seen.add(phrase)
            phrases.append(phrase)
        if len(phrases) >= 2:
            break

    return " and ".join(phrases)


def generateExecutiveSummary(data: dict[str, Any]) -> str:
    """Build a deterministic, publication-style executive summary from structured inputs."""
    n_samples = int(data.get("n_samples", 0) or 0)
    groups = [str(g) for g in (data.get("groups", []) or []) if str(g).strip()]
    groups_text = ", ".join(groups) if groups else "unspecified groups"

    deg_up = int(data.get("deg_up", 0) or 0)
    deg_down = int(data.get("deg_down", 0) or 0)
    total_deg = deg_up + deg_down

    pca_separation = str(data.get("pca_separation", "unknown")).strip() or "unknown"
    qc_warnings = [str(w).strip() for w in (data.get("qc_warnings", []) or []) if str(w).strip()]
    realism_flags = [str(f).strip() for f in (data.get("realism_flags", []) or []) if str(f).strip()]
    warning_items = [w for w in (data.get("warning_items", []) or []) if isinstance(w, dict)]

    sentences: list[str] = []

    sentences.append(f"This analysis includes {n_samples} samples across {groups_text}.")

    if total_deg > 0:
        sentences.append(
            f"A total of {total_deg} differentially expressed genes were identified ({deg_up} upregulated, {deg_down} downregulated)."
        )
    else:
        sentences.append("No significant differentially expressed genes were detected.")

    sentences.append(f"PCA indicates {pca_separation} separation between groups.")

    if qc_warnings:
        summarized = summarizeWarningsForSummary(warning_items)
        if summarized:
            sentences.append(f"Data quality concerns were identified, including {summarized}.")

    if realism_flags:
        sentences.append("The DEG profile shows patterns that may indicate potential data realism issues.")

    if n_samples < 4:
        sentences.append("The small sample size limits statistical power.")

    if qc_warnings or realism_flags:
        sentences.append("Therefore, results should be interpreted with caution.")
    else:
        sentences.append("No major technical concerns were detected.")

    return " ".join(sentences[:5])


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

    # Deduplicate by code + sample + metric for traceable deterministic rendering.
    seen: set[str] = set()
    unique_items: list[dict[str, Any]] = []
    for item in collected_items:
        key = "|".join([
            str(item.get("type", "")),
            str(item.get("severity", "")),
            str(item.get("code", "")),
            str(item.get("sample", "")),
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
                "message": str(item.get("message", "")),
                "code": str(item.get("code", "")),
                "sample": str(item.get("sample", "")) if item.get("sample") is not None else "",
                "metric": str(item.get("metric", "")) if item.get("metric") is not None else "",
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
) -> list[dict[str, Any]]:
    """Compute normalized QC metric rows with explicit status thresholds for report rendering."""
    if not isinstance(qc_metrics, dict) or not qc_metrics:
        return []

    rows: list[dict[str, Any]] = []
    warning_texts = [str(w) for w in (qc_warnings or [])]

    library_size = qc_metrics.get("library_size")
    lib_ratio = _extract_metric_value(library_size, ["min_median_ratio", "min_to_median", "ratio"])
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
                "value": f"{lib_ratio:.3f} (min/median)",
                "status": lib_status,
            }
        )

    correlation = qc_metrics.get("correlation")
    corr_mean = _extract_metric_value(correlation, ["mean", "mean_correlation", "average", "avg"])
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
                "name": "Sample correlation",
                "value": f"{corr_mean:.3f}",
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
                "value": f"{zero_mean:.3f}",
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
                "value": f"{pc1_percent:.2f}%",
                "status": pca_status,
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
) -> dict[str, Any]:
    issue_type = str(item.get("type", "statistical")).lower() or "statistical"
    code = str(item.get("code", "")).strip() or "unknown_issue"
    severity = str(item.get("severity", "warning")).lower()
    message = str(item.get("message", "")).strip()
    sample = str(item.get("sample", "")).strip()
    metric = str(item.get("metric", "")).strip()

    issue_name = "Data quality issue"
    evidence: Optional[str] = None

    c = code.lower()
    m = message.lower()

    # Message semantics take precedence over code to avoid swapped canonical/housekeeping labels.
    is_housekeeping = ("housekeeping" in m) or ("housekeeping" in c)
    is_canonical = ("canonical" in m) or ("canonical" in c)

    if "group_imbalance" in c:
        issue_name = "Group imbalance"
        ratio_match = re.search(r"ratio\s*([0-9.]+)", m)
        if len(group_counts) >= 2:
            sorted_groups = sorted(group_counts.items(), key=lambda x: (x[1], x[0]))
            low_group, low_n = sorted_groups[0]
            high_group, high_n = sorted_groups[-1]
            ratio = (high_n / low_n) if low_n > 0 else float("inf")
            if ratio_match:
                evidence = f"{low_n} {low_group} vs {high_n} {high_group}; ratio {float(ratio_match.group(1)):.2f}"
            else:
                evidence = f"{low_n} {low_group} vs {high_n} {high_group}; ratio {ratio:.2f}"
        elif ratio_match:
            evidence = f"ratio {float(ratio_match.group(1)):.2f}"
    elif "low_library_size" in c or "library_size" in c:
        issue_name = f"Low library size{' for ' + sample if sample else ''}"
        parsed = re.search(r"for\s+([^\s]+)\s*\(([^\)]+)\)", message)
        if parsed:
            sample_name = parsed.group(1)
            issue_name = f"Low library size for {sample_name}"
            evidence = parsed.group(2)
        elif message:
            evidence = _normalize_text(message)
    elif "correlation" in c:
        issue_name = f"Low within-group correlation{' for ' + sample if sample else ''}"
        if message:
            evidence = _normalize_text(message)
    elif "zero_fraction" in c or "zero" in c:
        issue_name = f"High zero-count fraction{' for ' + sample if sample else ''}"
        if message:
            evidence = _normalize_text(message)
    elif "top5" in c or "dominance" in c or "top-gene dominance" in m:
        issue_name = "Top-ranked gene concentration"
        if message:
            evidence = _normalize_text(message)
        if not evidence:
            evidence = "dominance among leading genes"
    elif is_housekeeping and not is_canonical:
        issue_name = "Housekeeping-gene enrichment"
        hk_match = re.search(r"(\d+\s+housekeeping\s+genes?\s+in\s+top\s+20)", message, flags=re.IGNORECASE)
        if hk_match:
            evidence = f"{hk_match.group(1)} DEGs"
        elif message:
            evidence = _normalize_text(message)
    elif is_canonical:
        issue_name = "Canonical-gene enrichment"
        canonical_match = re.search(r"(\d+\s*/\s*20\s+canonical\s+genes)", message, flags=re.IGNORECASE)
        if canonical_match:
            evidence = f"{canonical_match.group(1)} among top-ranked DEGs"
        elif message:
            evidence = _normalize_text(message)
    elif "pvalue" in c or "p_value" in c:
        issue_name = "P-value distribution anomaly"
        if message:
            evidence = _normalize_text(message)
    elif issue_type == "statistical":
        issue_name = "Statistical consistency issue"
        if message:
            evidence = _normalize_text(message)

    if not evidence and message:
        evidence = _normalize_text(message)
    if not evidence and ("top5" in c or "dominance" in c):
        evidence = "dominance among leading genes"
    if not evidence and metric:
        evidence = _normalize_text(metric.replace("_", " "))

    return {
        "type": issue_type,
        "code": code,
        "severity": "critical" if severity == "critical" else "warning",
        "message": message,
        "evidence": evidence,
        "issue": issue_name,
    }


def _assessment_phrase(issue: dict[str, Any]) -> str:
    issue_name = str(issue.get("issue", "Data quality issue")).strip() or "Data quality issue"
    evidence = str(issue.get("evidence", "")).strip()
    if evidence:
        return f"{issue_name} detected ({evidence})."
    if "top-ranked gene concentration" in issue_name.lower():
        return "Top-ranked gene concentration detected (dominance among leading genes)."
    return f"{issue_name} detected."


def generateAssessmentBasis(
    warning_items: list[dict[str, Any]],
    n_samples: int,
    groups: list[str],
) -> list[str]:
    _ = n_samples  # Kept for function signature compatibility.

    items = [w for w in (warning_items or []) if isinstance(w, dict)]
    group_counts = Counter([str(g).strip() for g in (groups or []) if str(g).strip()])

    normalized_issues = [_assessment_issue_from_warning(item, group_counts) for item in items]

    # Deduplicate by code while preserving code-level evidence aggregation.
    by_code: dict[str, dict[str, Any]] = {}
    for issue in normalized_issues:
        code = str(issue.get("code", "")).strip() or "unknown_issue"
        if code not in by_code:
            by_code[code] = issue
        else:
            current = by_code[code]
            current_evidence = str(current.get("evidence", "")).strip()
            new_evidence = str(issue.get("evidence", "")).strip()
            if new_evidence and new_evidence not in current_evidence:
                merged = "; ".join([x for x in [current_evidence, new_evidence] if x])
                current["evidence"] = merged
            if str(current.get("severity", "warning")).lower() != "critical" and str(issue.get("severity", "warning")).lower() == "critical":
                current["severity"] = "critical"

    severity_order = {"critical": 0, "warning": 1}
    type_order = {"qc": 0, "realism": 1, "statistical": 2}
    ordered = sorted(
        by_code.values(),
        key=lambda x: (
            severity_order.get(str(x.get("severity", "warning")).lower(), 1),
            type_order.get(str(x.get("type", "statistical")).lower(), 2),
            str(x.get("code", "")),
        ),
    )

    if not ordered:
        return ["No major quality or realism concerns detected."]

    return [_assessment_phrase(issue) for issue in ordered]


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

    results_paragraph = (
        f"The analysis included {summary_data.get('n_samples', 0)} samples across groups: {', '.join(groups) if groups else 'not provided'}. "
        f"Differential testing identified {summary_data.get('deg_up', 0)} upregulated and {summary_data.get('deg_down', 0)} downregulated genes. "
        f"PCA separation was classified as '{summary_data.get('pca_separation', 'unknown')}'. "
        f"Overall data quality risk was assessed as '{_overall_data_quality(qc_report)}', and realism suspicion was "
        f"'{realism.get('overall_suspicion', 'unknown')}'."
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
                    "message": str(item.get("message", "")),
                    "sample": str(item.get("sample", "")) if item.get("sample") is not None else "",
                    "metric": str(item.get("metric", "")) if item.get("metric") is not None else "",
                }
            )

    groups_for_assessment = _expand_groups_for_assessment(groups, qc_report)
    assessment_basis = generateAssessmentBasis(
        warning_items=all_warning_items,
        n_samples=int(summary_data.get("n_samples", 0) or 0),
        groups=groups_for_assessment,
    )
    qc_metrics_summary = computeQCMetrics(
        qc_metrics=(qc_report or {}).get("qc_metrics") if isinstance(qc_report, dict) else None,
        pca_variance=(qc_report or {}).get("pca_variance") if isinstance(qc_report, dict) else None,
        qc_warnings=(qc_report or {}).get("qc_warnings", []) if isinstance(qc_report, dict) else [],
    )
    realism_metrics = computeRealismMetrics(top_genes_table)
    realism_assessment = evaluateRealism(realism_metrics)
    realism_flag_map = _map_realism_flags_to_metrics((realism or {}).get("realism_flags", []) or [])
    if realism_flag_map.get("unmatched"):
        realism_assessment["reasons"] = list(realism_assessment.get("reasons", []))
        realism_assessment["reasons"].append(
            "Some realism flags do not map to quantitative metrics; metric-based assessment is prioritized."
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
    interpretation_limitation_text = generateLimitationText(
        interpretation_confidence.get("level", "MEDIUM"),
        interpretation_confidence.get("reasons", []),
    )
    executive_summary = generateExecutiveSummary(
        {
            "n_samples": summary_data.get("n_samples", 0),
            "groups": groups,
            "deg_up": summary_data.get("deg_up", 0),
            "deg_down": summary_data.get("deg_down", 0),
            "pca_separation": summary_data.get("pca_separation", "unknown"),
            "qc_warnings": interpretation_qc_warnings,
            "realism_flags": interpretation_realism_flags,
            "warning_items": all_warning_items,
        }
    )

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
        realism_flag_map=realism_flag_map,
        interpretation_confidence=interpretation_confidence,
        interpretation_limitation_text=interpretation_limitation_text,
        executive_summary=executive_summary,
        analysis_methods=analysis_methods,
        methods_paragraph=methods_paragraph,
        results_paragraph=results_paragraph,
        figure_legends=figure_legends,
        overall_data_quality=_overall_data_quality(qc_report),
        overall_realism=realism.get("overall_suspicion", "unknown"),
        llm=llm.model_dump() if llm else None,
        plots=plots,
        job_id=job_dir.name,
    )

    out = job_dir / "report.html"
    out.write_text(html, encoding="utf-8")
    logger.info("Report written to %s", out)
