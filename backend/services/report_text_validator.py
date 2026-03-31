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


# Centralized display-precision policy for all rendered report text.
METRIC_PRECISION: dict[str, int] = {
    "canonical_fraction": 3,
    "mean_correlation": 3,
    "extreme_pvalue_fraction": 3,
    "zero_fraction": 3,
    "library_size_ratio": 3,
    "top5_contribution": 3,
    "group_size_ratio": 2,
    "housekeeping_genes": 0,
    "deg_up": 0,
    "deg_down": 0,
    "deg_total": 0,
}

METRIC_LABEL_TO_CODE: dict[str, str] = {
    "canonical gene fraction": "canonical_fraction",
    "canonical fraction": "canonical_fraction",
    "control group mean correlation": "mean_correlation",
    "treated group mean correlation": "mean_correlation",
    "mean correlation": "mean_correlation",
    "group mean correlation": "mean_correlation",
    "top 5 gene contribution": "top5_contribution",
    "top5 contribution": "top5_contribution",
    "group size ratio": "group_size_ratio",
    "library size ratio": "library_size_ratio",
    "zero-count fraction": "zero_fraction",
    "zero fraction": "zero_fraction",
    "extreme p-value fraction": "extreme_pvalue_fraction",
    "housekeeping genes detected": "housekeeping_genes",
}


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


def get_metric_precision(metric_code: str) -> int:
    """Return the central display precision for a metric code."""
    return int(METRIC_PRECISION.get(str(metric_code or "").strip(), 3))


def format_metric_value(metric_code: str, value: Any) -> str:
    """Format numeric value using the central metric precision policy."""
    code = str(metric_code or "").strip()
    precision = get_metric_precision(code)

    try:
        num = float(value)
    except Exception:
        return _clean_text(value)

    if precision <= 0:
        return str(int(round(num)))
    return f"{num:.{precision}f}"


def _metric_code_from_label(label: str) -> str:
    key = _clean_text(label).lower()
    if key in METRIC_LABEL_TO_CODE:
        return METRIC_LABEL_TO_CODE[key]
    return _metric_key_from_text(key)


def _safe_float(value: Any) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def _extract_group_size_ratio_from_qc(qc_report: dict[str, Any]) -> float | None:
    reps = qc_report.get("replicates_per_group") if isinstance(qc_report.get("replicates_per_group"), dict) else {}
    if len(reps) < 2:
        return None
    try:
        nums = [int(v) for v in reps.values()]
        low = min(nums)
        high = max(nums)
        if low <= 0:
            return None
        return float(high) / float(low)
    except Exception:
        return None


def _build_metric_value_registry(analysis_json: dict[str, Any]) -> dict[str, Any]:
    """Build a canonical metric value registry from structured analysis outputs."""
    registry: dict[str, Any] = {}

    realism = analysis_json.get("realism_metrics") if isinstance(analysis_json.get("realism_metrics"), dict) else {}
    qc_report = analysis_json.get("qc_report") if isinstance(analysis_json.get("qc_report"), dict) else {}

    if "canonical_fraction" in realism:
        registry["canonical_fraction"] = realism.get("canonical_fraction")
    if "top5_contribution" in realism:
        registry["top5_contribution"] = realism.get("top5_contribution")
    if "extreme_pvalue_fraction" in realism:
        registry["extreme_pvalue_fraction"] = realism.get("extreme_pvalue_fraction")

    hk = realism.get("housekeeping_genes")
    if isinstance(hk, list):
        registry["housekeeping_genes"] = len([x for x in hk if str(x).strip()])

    group_qc = qc_report.get("group_qc") if isinstance(qc_report.get("group_qc"), dict) else {}
    for group_obj in group_qc.values():
        if isinstance(group_obj, dict) and "mean_correlation" in group_obj:
            registry["mean_correlation"] = group_obj.get("mean_correlation")
            break

    qc_metrics = qc_report.get("qc_metrics") if isinstance(qc_report.get("qc_metrics"), dict) else {}
    zero_obj = qc_metrics.get("zero_fraction") if isinstance(qc_metrics.get("zero_fraction"), dict) else {}
    if "mean" in zero_obj:
        registry["zero_fraction"] = zero_obj.get("mean")

    library_obj = qc_metrics.get("library_size") if isinstance(qc_metrics.get("library_size"), dict) else {}
    if "min_median_ratio" in library_obj:
        registry["library_size_ratio"] = library_obj.get("min_median_ratio")

    ratio = _extract_group_size_ratio_from_qc(qc_report)
    if ratio is not None:
        registry["group_size_ratio"] = ratio

    registry["deg_up"] = int(analysis_json.get("deg_up", 0) or 0)
    registry["deg_down"] = int(analysis_json.get("deg_down", 0) or 0)
    registry["deg_total"] = registry["deg_up"] + registry["deg_down"]

    return registry


def normalize_metric_display(text: str, metric_registry: dict[str, Any]) -> tuple[str, list[tuple[str, str, str]]]:
    """Normalize metric value display using central policy and canonical metric registry."""
    t = str(text or "")
    changes: list[tuple[str, str, str]] = []
    if not t:
        return t, changes

    def _repl(match: re.Match[str]) -> str:
        raw_label = str(match.group(1) or "")
        leading_ws_match = re.match(r"^\s*", raw_label)
        leading_ws = leading_ws_match.group(0) if leading_ws_match else ""
        label_body = raw_label[len(leading_ws):]

        conjunction = ""
        for c in ["and ", "or "]:
            if label_body.lower().startswith(c):
                conjunction = label_body[: len(c)]
                label_body = label_body[len(c):]
                break

        metric_label = _clean_text(label_body)
        shown_value = _clean_text(match.group(2))
        code = _metric_code_from_label(metric_label)
        if not code:
            return match.group(0)

        canonical_raw = metric_registry.get(code, shown_value)
        canonical_fmt = format_metric_value(code, canonical_raw)
        if shown_value != canonical_fmt:
            changes.append((code, shown_value, canonical_fmt))
        return f"{leading_ws}{conjunction}{metric_label} = {canonical_fmt}"

    t = re.sub(r"([A-Za-z0-9\- ]+?)\s*=\s*([-+]?\d*\.?\d+)", _repl, t)
    return _clean_text(t), changes


def validate_metric_display_consistency(report_text: dict[str, Any], metric_registry: dict[str, Any]) -> list[str]:
    """Validate rendered metric precision consistency across sections."""
    logs: list[str] = []
    observed: dict[str, set[str]] = {}

    def _collect_from_text(text: str, section: str) -> None:
        for m in re.finditer(r"([A-Za-z0-9\- ]+?)\s*=\s*([-+]?\d*\.?\d+)", str(text or "")):
            code = _metric_code_from_label(m.group(1))
            if not code:
                continue
            shown = _clean_text(m.group(2))
            observed.setdefault(code, set()).add(shown)
            expected = format_metric_value(code, metric_registry.get(code, shown))
            if shown != expected:
                logs.append(f"Metric display mismatch for {code} in {section}: {shown} vs expected {expected}")

    _collect_from_text(str(report_text.get("executive_summary", "")), "Executive Summary")
    for item in report_text.get("assessment_basis", []) or []:
        _collect_from_text(str(item), "Assessment Basis")

    ai = report_text.get("ai_interpretation", {})
    if isinstance(ai, dict):
        for key, val in ai.items():
            _collect_from_text(str(val), f"AI Interpretation/{key}")

    for code, values in observed.items():
        if len(values) > 1:
            logs.append(f"Inconsistent display precision for {code}: {sorted(values)}")

    return logs


def _extract_metric_value_threshold(text: str, metric_code: str | None = None) -> tuple[str, str, str, str] | None:
    """Extract metric key, value, threshold operator, and threshold value from one statement."""
    t = _clean_text(text)
    if not t:
        return None

    metric_key = metric_code or _metric_key_from_text(t)
    m_value = re.search(r"=\s*([-+]?\d*\.?\d+)", t)
    m_threshold = re.search(
        r"\(\s*(?:warning\s+)?threshold\s*(<=|>=|<|>)\s*([-+]?\d*\.?\d+)\s*\)",
        t,
        flags=re.IGNORECASE,
    )
    if not (m_value and m_threshold):
        return None
    return metric_key, m_value.group(1), m_threshold.group(1), m_threshold.group(2)


def normalize_threshold_phrasing(text: str, metric_code: str | None = None) -> str:
    """Normalize threshold wording into consistent expected/warning prose."""
    t = _clean_text(text)
    if not t:
        return t

    parsed = _extract_metric_value_threshold(t, metric_code)
    if not parsed:
        return t

    metric_key, metric_value_str, op, threshold_value_str = parsed

    if metric_key in {"mean_correlation", "library_size_ratio"}:
        replacement = f"(expected ≥ {threshold_value_str})"
    elif metric_key in {"canonical_fraction", "zero_fraction", "extreme_pvalue_fraction"}:
        replacement = f"(expected ≤ {threshold_value_str})"
    elif metric_key in {"group_size_ratio", "housekeeping_genes", "top5_contribution"}:
        relation = "exceeds"
        try:
            mv = float(metric_value_str)
            tv = float(threshold_value_str)
            if abs(mv - tv) < 1e-12:
                relation = "meets"
            elif op in {">", ">="} and mv < tv:
                relation = "below"
            elif op in {"<", "<="} and mv > tv:
                relation = "below"
        except Exception:
            relation = "exceeds"

        if relation == "below":
            replacement = f"(below warning threshold of {threshold_value_str})"
        else:
            replacement = f"({relation} warning threshold of {threshold_value_str})"
    else:
        invert = {"<": "≥", "<=": ">", ">": "≤", ">=": "<"}
        replacement = f"(expected {invert.get(op, _format_symbol(op))} {threshold_value_str})"

    t = re.sub(
        r"\(\s*(?:warning\s+)?threshold\s*(<=|>=|<|>)\s*([-+]?\d*\.?\d+)\s*\)",
        replacement,
        t,
        flags=re.IGNORECASE,
    )
    return _clean_text(t)


def rewrite_threshold_phrasing(text: str, metric_code: str | None = None) -> str:
    """Backward-compatible wrapper for threshold phrasing normalization."""
    return normalize_threshold_phrasing(text, metric_code)


def infer_meaning_phrase(metric_code: str, value: float | None, threshold: float | None) -> str:
    """Infer deterministic scientific meaning phrase from metric semantics."""
    code = str(metric_code or "").strip()

    if code == "mean_correlation":
        if value is not None and threshold is not None and value < threshold:
            return "Low control-group consistency was observed"
        return "Control-group consistency was within expected range"
    if code == "library_size_ratio":
        if value is not None and threshold is not None and value < threshold:
            return "Library-size imbalance was observed"
        return "Library-size balance was within expected range"
    if code == "group_size_ratio":
        if value is not None and threshold is not None and value > threshold:
            return "Group imbalance was observed"
        if value is not None and threshold is not None and abs(value - threshold) < 1e-12:
            return "Group imbalance met the warning threshold"
        return "Group-size distribution was within expected range"
    if code == "canonical_fraction":
        if value is not None and threshold is not None and value > threshold:
            return "Borderline canonical-gene enrichment was observed"
        if value is not None and threshold is not None and abs(value - threshold) < 1e-12:
            return "Borderline canonical-gene enrichment was observed"
        return "Canonical-gene enrichment was within expected range"
    if code == "housekeeping_genes":
        if value is not None and threshold is not None and value > threshold:
            return "Elevated housekeeping-gene signal was observed"
        if value is not None and threshold is not None and abs(value - threshold) < 1e-12:
            return "Housekeeping-gene signal met the warning threshold"
        return "Housekeeping-gene signal was limited"
    if code == "top5_contribution":
        if value is not None and threshold is not None and value > threshold:
            return "Expression dominance by top-ranked genes was observed"
        return "Top-gene contribution was within expected range"
    if code == "zero_fraction":
        if value is not None and threshold is not None and value > threshold:
            return "High zero-count burden was observed"
        return "Zero-count burden was within expected range"
    if code == "extreme_pvalue_fraction":
        if value is not None and threshold is not None and value > threshold:
            return "P-value distribution anomaly was observed"
        return "P-value distribution was within expected range"

    return "A quantified finding was observed"


def polish_metric_sentence(text: str, metric_code: str | None = None) -> str:
    """Convert metric-log style statements into scientific Results-style prose."""
    t = _clean_text(text)
    if not t:
        return t

    normalized = normalize_threshold_phrasing(upgrade_metric_name(t), metric_code)
    metric_key = metric_code or _metric_key_from_text(normalized)

    match = re.search(
        r"([^=()]+?)\s*=\s*([-+]?\d*\.?\d+)\s*\(([^)]*)\)",
        normalized,
        flags=re.IGNORECASE,
    )
    if not match:
        return polish_scientific_phrasing(normalized if normalized.endswith(".") else f"{normalized}.")

    metric_label = _clean_text(match.group(1))
    value_str = match.group(2)
    threshold_text = _clean_text(match.group(3))

    value_num = None
    threshold_num = None
    try:
        value_num = float(value_str)
    except Exception:
        value_num = None

    tnum = re.search(r"([-+]?\d*\.?\d+)", threshold_text)
    if tnum:
        try:
            threshold_num = float(tnum.group(1))
        except Exception:
            threshold_num = None

    meaning = infer_meaning_phrase(metric_key, value_num, threshold_num)
    sentence = f"{meaning} ({metric_label} = {value_str}; {threshold_text})."
    return polish_scientific_phrasing(sentence)


def rewrite_summary_sentence(text: str, analysis_json: dict[str, Any]) -> str:
    """Rewrite one summary finding sentence into meaning-first scientific prose."""
    _ = analysis_json
    t = _clean_text(text)
    if not t:
        return ""

    t = re.sub(r"^(Data quality assessment (identified|indicated)\s*)", "", t, flags=re.IGNORECASE)
    t = re.sub(r"^(Realism evaluation (identified|indicated)\s*)", "", t, flags=re.IGNORECASE)
    return polish_metric_sentence(t)


def fix_double_verbs(text: str) -> str:
    """Remove redundant double-verb constructions in scientific sentences."""
    t = _clean_text(text)
    if not t:
        return t

    t = re.sub(
        r"\bData quality assessment identified that\s+(.+?)\s+was observed\s*\(",
        r"Data quality assessment indicated \1 (",
        t,
        flags=re.IGNORECASE,
    )
    t = re.sub(
        r"\bRealism evaluation indicated that\s+(.+?)\s+was observed\s*\(",
        r"Realism evaluation indicated \1 (",
        t,
        flags=re.IGNORECASE,
    )
    t = re.sub(
        r"\b(identified|indicated) that\s+(.+?)\s+was observed\s*\(",
        r"indicated \2 (",
        t,
        flags=re.IGNORECASE,
    )
    t = re.sub(r"\bindicated that\s+", "indicated ", t, flags=re.IGNORECASE)
    return _clean_text(t)


def normalize_scientific_phrasing(text: str, metric_code: str | None = None) -> str:
    """Apply deterministic, concise scientific phrasing normalization."""
    _ = metric_code
    t = _clean_text(text)
    if not t:
        return t

    t = re.sub(r"\bmet the warning threshold\b", "reached the warning threshold", t, flags=re.IGNORECASE)
    t = re.sub(r"\bmeets warning threshold of\b", "reached warning threshold of", t, flags=re.IGNORECASE)
    t = re.sub(r"\bis at the warning threshold of\b", "reached warning threshold of", t, flags=re.IGNORECASE)
    t = re.sub(r"\bwas observed\s*\(", "was observed (", t)
    t = re.sub(r"\bindicated\s*\(", "indicated (", t)
    return _clean_text(t)


def rewrite_pca_sentence(text: str) -> str:
    """Rewrite PCA sentence into concise Results-style phrasing without changing meaning."""
    t = _clean_text(text)
    if not t:
        return t

    def _repl(match: re.Match[str]) -> str:
        level = _clean_text(match.group(1)).lower()
        return f"PCA showed {level} separation between groups."

    t = re.sub(
        r"\bPCA separation was classified as\s+([a-zA-Z-]+)\s*\.",
        _repl,
        t,
        flags=re.IGNORECASE,
    )
    return _clean_text(t)


def unify_metric_precision(text: str, metric_registry: dict[str, Any]) -> tuple[str, list[tuple[str, str, str]]]:
    """Unify metric display using centralized precision policy."""
    return normalize_metric_display(text, metric_registry)


def reduce_sentence_redundancy(lines: list[str]) -> list[str]:
    """Reduce repeated sentence patterns deterministically across consecutive lines."""
    out: list[str] = []
    observed_streak = 0

    for raw in lines or []:
        line = _clean_text(raw)
        if not line:
            continue

        if re.search(r"\bwas observed\s*\(", line, flags=re.IGNORECASE):
            observed_streak += 1
        else:
            observed_streak = 0

        if observed_streak >= 3:
            if observed_streak % 2 == 1:
                line = re.sub(r"\bwas observed\s*\(", "was detected (", line, flags=re.IGNORECASE)
            else:
                line = re.sub(r"\bwas observed\s*\(", "was noted (", line, flags=re.IGNORECASE)

        out.append(_clean_text(line if line.endswith(".") else f"{line}."))

    return out


def apply_final_micro_polish(report_text: dict[str, Any], analysis_json: dict[str, Any]) -> dict[str, Any]:
    """Apply final deterministic micro-polish without changing structure or numeric values."""
    out = dict(report_text or {})

    metric_registry = _build_metric_value_registry(analysis_json if isinstance(analysis_json, dict) else {})
    normalization_log: list[str] = []

    summary = _clean_text(out.get("executive_summary", ""))
    summary = rewrite_pca_sentence(summary)
    summary = fix_double_verbs(summary)
    summary = normalize_scientific_phrasing(summary)
    summary, summary_changes = unify_metric_precision(summary, metric_registry)
    for code, old, new in summary_changes:
        normalization_log.append(f"Normalized {code} display from {old} to {new} in Executive Summary")
    out["executive_summary"] = summary

    basis_lines = [str(x) for x in (out.get("assessment_basis", []) or [])]
    polished_basis: list[str] = []
    for line in basis_lines:
        v = _clean_text(line)
        v = fix_double_verbs(v)
        v = normalize_scientific_phrasing(v, _metric_key_from_text(v))
        v, basis_changes = unify_metric_precision(v, metric_registry)
        for code, old, new in basis_changes:
            normalization_log.append(f"Normalized {code} display from {old} to {new} in Assessment Basis")
        polished_basis.append(v)
    out["assessment_basis"] = reduce_sentence_redundancy(polished_basis)

    grouped = out.get("assessment_basis_grouped", {})
    if isinstance(grouped, dict):
        polished_grouped: dict[str, list[str]] = {}
        for section, lines in grouped.items():
            section_lines = [str(x) for x in (lines or [])]
            local: list[str] = []
            for line in section_lines:
                v = _clean_text(line)
                v = fix_double_verbs(v)
                v = normalize_scientific_phrasing(v, _metric_key_from_text(v))
                v, grouped_changes = unify_metric_precision(v, metric_registry)
                for code, old, new in grouped_changes:
                    normalization_log.append(f"Normalized {code} display from {old} to {new} in Assessment Basis Grouped")
                local.append(v)
            polished_grouped[str(section)] = reduce_sentence_redundancy(local)
        out["assessment_basis_grouped"] = polished_grouped

    ai = out.get("ai_interpretation", {})
    if isinstance(ai, dict):
        polished_ai: dict[str, str] = {}
        for key, value in ai.items():
            s = normalize_scientific_phrasing(_clean_text(value))
            s, ai_changes = unify_metric_precision(s, metric_registry)
            for code, old, new in ai_changes:
                normalization_log.append(f"Normalized {code} display from {old} to {new} in AI Interpretation/{key}")
            polished_ai[str(key)] = s
        out["ai_interpretation"] = polished_ai

    out["precision_normalization_log"] = merge_duplicate_statements(normalization_log)

    return out


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


def sort_realism_items(items: list[str]) -> list[str]:
    """Sort realism statements as canonical -> housekeeping -> distribution/dominance -> others."""
    def sort_key(statement: str) -> tuple[int, str]:
        t = str(statement or "").lower()
        if "canonical" in t:
            return (0, t)
        if "housekeeping" in t:
            return (1, t)
        if any(k in t for k in ["distribution", "dominance", "top 5 gene contribution", "top5 contribution", "p-value", "pvalue", "extreme p-value"]):
            return (2, t)
        return (3, t)

    return sorted([_clean_text(x) for x in (items or []) if _clean_text(x)], key=sort_key)


def relabel_assessment_groups(groups: dict[str, list[str]]) -> dict[str, list[str]]:
    """Ensure correlation/consistency findings are always under QC - Group-level."""
    relabeled = {k: list(v) for k, v in (groups or {}).items()}
    relabeled.setdefault("QC - Group-level", [])

    correlation_markers = [
        "mean correlation",
        "group consistency",
        "internal consistency",
        "control-group consistency",
        "treated-group consistency",
    ]

    for source_key in list(relabeled.keys()):
        if source_key == "QC - Group-level":
            continue

        kept: list[str] = []
        for item in relabeled.get(source_key, []):
            tl = str(item or "").lower()
            if any(marker in tl for marker in correlation_markers):
                relabeled["QC - Group-level"].append(_clean_text(item))
            else:
                kept.append(_clean_text(item))
        relabeled[source_key] = kept

    ordered: dict[str, list[str]] = {}
    for section in ASSESSMENT_GROUP_ORDER:
        vals = merge_duplicate_statements(relabeled.get(section, []))
        if vals:
            ordered[section] = vals
    for section, vals in relabeled.items():
        if section not in ordered:
            dedup = merge_duplicate_statements(vals)
            if dedup:
                ordered[section] = dedup
    return ordered


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

    compact = {k: v for k, v in grouped.items() if v}
    compact = relabel_assessment_groups(compact)
    if compact.get("Realism"):
        compact["Realism"] = sort_realism_items(compact["Realism"])
    return compact


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
    t = polish_metric_sentence(text)
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

    qc_finding = rewrite_summary_sentence(
        _major_qc_issue_text(analysis_json) or "sample size = 0 (threshold >= 6).",
        analysis_json,
    )
    realism_finding = rewrite_summary_sentence(
        _major_realism_issue_text(analysis_json) or "canonical gene fraction = 0.000 (threshold > 0.30).",
        analysis_json,
    )

    confidence = analysis_json.get("confidence_assessment") if isinstance(analysis_json.get("confidence_assessment"), dict) else {}
    confidence_score = int(confidence.get("confidence_score", 0) or 0)
    confidence_level = str(confidence.get("confidence_level", "LOW") or "LOW").upper()
    confidence_explanation = str(confidence.get("confidence_explanation", "") or "").strip()

    summary = (
        f"The analysis included {n_samples} samples across {group_text} groups, with {total_deg} differentially expressed genes identified "
        f"({deg_up} upregulated and {deg_down} downregulated). "
        f"PCA separation was classified as {pca}. "
        f"Data quality assessment identified that {qc_finding.rstrip('.')[:1].lower() + qc_finding.rstrip('.')[1:]}. "
        f"Realism evaluation indicated that {realism_finding.rstrip('.')[:1].lower() + realism_finding.rstrip('.')[1:]}. "
        f"Global confidence score = {confidence_score}/100 ({confidence_level}). "
        f"{confidence_explanation if confidence_explanation else ''} "
        "These findings are consistent with a hypothesis-generating interpretation rather than a confirmatory conclusion."
    )
    return _clean_text(summary)


def apply_final_final_patch(report_text: dict[str, Any], analysis_json: dict[str, Any]) -> dict[str, Any]:
    """Apply final deterministic micro-polish for wording and grouping semantics."""
    out = dict(report_text or {})

    out["executive_summary"] = rewrite_executive_summary(out.get("executive_summary", ""), analysis_json)

    basis_items = [polish_metric_sentence(x) for x in (out.get("assessment_basis", []) or []) if _clean_text(x)]
    basis_items = merge_duplicate_statements(basis_items)
    grouped = group_assessment_basis(basis_items)
    grouped = relabel_assessment_groups(grouped)

    flattened: list[str] = []
    for section in ASSESSMENT_GROUP_ORDER:
        flattened.extend(grouped.get(section, []))

    out["assessment_basis"] = flattened
    out["assessment_basis_grouped"] = grouped

    ai = out.get("ai_interpretation", {})
    if isinstance(ai, dict):
        out["ai_interpretation"] = {k: polish_scientific_phrasing(v) for k, v in ai.items()}

    return out


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

    cleaned = apply_final_final_patch(cleaned, analysis_json)
    validation_log.append("Applied FINAL FINAL wording/grouping micro-polish layer")

    cleaned = apply_final_micro_polish(cleaned, analysis_json)
    validation_log.append("Applied FINAL MICRO phrasing consistency polish")

    for item in cleaned.get("precision_normalization_log", []) or []:
        validation_log.append(str(item))
    cleaned.pop("precision_normalization_log", None)

    metric_registry = _build_metric_value_registry(analysis_json)
    precision_check_logs = validate_metric_display_consistency(cleaned, metric_registry)
    validation_log.extend(precision_check_logs)
    if precision_check_logs:
        validation_log.append("Metric display precision consistency check reported mismatches")
    else:
        validation_log.append("Metric display precision consistency check passed")

    cleaned["validation_log"] = merge_duplicate_statements(validation_log)
    return cleaned


def build_analysis_snapshot(
    summary_data: dict[str, Any],
    qc_report: dict[str, Any],
    realism_metrics: dict[str, Any],
    realism_level: str,
    confidence_assessment: dict[str, Any] | None = None,
    ncrna_assessment: dict[str, Any] | None = None,
    analysis_status: dict[str, Any] | None = None,
    deg_status: str | None = None,
) -> dict[str, Any]:
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
        "analysis_status": analysis_status or {},
        "deg_status": deg_status or "",
        "confidence_assessment": confidence_assessment or {},
        "ncrna_assessment": ncrna_assessment or {},
    }
