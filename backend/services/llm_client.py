"""
LLM interpretation service.
Accepts a structured AnalysisSummary JSON — never raw counts.
"""
from __future__ import annotations

import json
import logging
import re
from typing import Any, Optional

from openai import OpenAI

from config import settings
from models.schemas import AnalysisSummary, LLMInterpretation

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are a senior bioinformatics analyst writing a Results-style interpretation from structured RNA-seq summary JSON only.

Mission:
Produce scientifically rigorous, evidence-based interpretation with explicit uncertainty control based on data quality and realism risk.

Hard constraints:
1) Use only the provided JSON fields.
2) Do not invent pathways, mechanisms, diseases, or external biology.
3) Do not mention genes not present in top_genes.
4) If evidence is insufficient, explicitly state insufficient evidence.
5) If qc_critical is non-empty or overall_qc_status is high risk, downgrade confidence and state unreliability.
6) If overall_realism is high suspicion or realism_flags is non-empty, explicitly state possible artificial or biased signal.
7) Every conclusion must be traceable to specific input fields.
8) You MUST explicitly reference sample size (n_samples) in the interpretation.
9) You MUST mention QC limitations whenever confidence is not HIGH.
10) Avoid strong causal claims (e.g., "causes", "proves", "demonstrates").
11) Do not infer pathways unless pathway information is explicitly provided in input.

Style:
- Formal scientific English.
- Concise and specific.
- No fluff.

Return JSON with exactly these keys:
{
  "pca_text": "<Section 1: PCA Interpretation>",
  "deg_summary": "<Section 2: Differential Expression Summary>",
  "biological_insights": "<Section 3: Biological Insight>",
  "data_quality": "<Section 4: Limitations>",
  "next_steps": "<Section 5: Recommendations>",
  "methods_paragraph": "<Publication text: Methods paragraph>",
  "results_paragraph": "<Publication text: Results paragraph>",
  "figure_legend": "<Publication text: PCA/Volcano figure legend>"
}

Section requirements:
- Keep each section concise and evidence-bound.
- If confidence level is LOW, use extra-cautious language and short statements.
- Limitations section must align with provided confidence reasons.

Critical quality controls:
- If qc_critical exists: use language such as "may", "suggests", "uncertain".
- If overall_qc_status is high risk: explicitly say conclusions are unreliable.
- If overall_realism is high suspicion: explicitly say results may reflect artificial or biased data.
- Never make strong biological claims under elevated QC/realism risk.
"""


def _qc_status(qc_report: Optional[dict[str, Any]]) -> str:
    if not qc_report:
        return "low risk"
    n_critical = len(qc_report.get("qc_critical", []) or [])
    n_warning = len(qc_report.get("qc_warnings", []) or [])
    if n_critical > 0:
        return "high risk"
    if n_warning > 0:
        return "moderate risk"
    return "low risk"


def evaluateInterpretationConfidence(
    qc_warnings: list[str],
    realism_flags: list[str],
    n_samples: int,
) -> dict[str, Any]:
    warnings = [str(w).strip() for w in (qc_warnings or []) if str(w).strip()]
    flags = [str(f).strip() for f in (realism_flags or []) if str(f).strip()]

    reasons: list[str] = [f"Sample size considered in interpretation (n={int(n_samples)})."]
    has_critical_qc = any(
        any(token in w.lower() for token in ["critical", "high risk", "severe"]) for w in warnings
    )

    if has_critical_qc or len(flags) > 0:
        level = "LOW"
        if has_critical_qc:
            reasons.append("Critical QC issue detected; interpretation reliability is reduced.")
        if flags:
            reasons.append("Realism flags detected; potential artificial or biased signal cannot be excluded.")
    elif warnings or int(n_samples) < 6:
        level = "MEDIUM"
        if warnings:
            reasons.append("QC warnings present without critical findings.")
        if int(n_samples) < 6:
            reasons.append(f"Small sample size (n={int(n_samples)}) limits statistical stability.")
    else:
        level = "HIGH"
        reasons.append("No major technical limitations detected.")

    return {
        "level": level,
        "reasons": reasons,
    }


def generateLimitationText(confidence: str, reasons: list[str]) -> str:
    level = str(confidence or "").upper()
    if level == "LOW":
        base = "Interpretation is limited by data quality issues and should be treated with caution."
    elif level == "MEDIUM":
        base = "Results are generally consistent but should be interpreted with awareness of dataset limitations."
    else:
        base = "No major technical concerns detected; results are considered robust within dataset scope."

    summary = " ".join(str(r).strip() for r in (reasons or []) if str(r).strip())
    return f"{base} {summary}".strip()


def _limit_sentences(text: str, max_sentences: int = 3) -> str:
    chunks = re.split(r"(?<=[.!?])\s+", str(text or "").strip())
    chunks = [c for c in chunks if c]
    return " ".join(chunks[:max_sentences]) if chunks else str(text or "")


def _cautious_rewrite(text: str) -> str:
    replacements = {
        " proves ": " suggests ",
        " demonstrates ": " indicates ",
        " causes ": " may contribute to ",
        " confirms ": " is consistent with ",
    }
    out = f" {str(text or '')} "
    for src, dst in replacements.items():
        out = out.replace(src, dst)
    return out.strip()


def _build_llm_input(summary: AnalysisSummary, qc_report: Optional[dict[str, Any]]) -> dict[str, Any]:
    realism = summary.realism_validation
    overall_realism = "low suspicion"
    realism_flags: list[str] = []
    if realism is not None:
        overall_realism = f"{realism.overall_suspicion} suspicion"
        realism_flags = list(realism.realism_flags)

    qc_warnings = list(summary.warnings) + list(summary.data_issues)
    qc_critical = []
    if qc_report:
        qc_warnings = list(qc_report.get("qc_warnings", []) or qc_warnings)
        qc_critical = list(qc_report.get("qc_critical", []) or [])

    confidence = evaluateInterpretationConfidence(
        qc_warnings=[*qc_critical, *qc_warnings],
        realism_flags=realism_flags,
        n_samples=summary.n_samples,
    )
    limitation_text = generateLimitationText(confidence["level"], confidence["reasons"])

    structured_summary = {
        "deg_up": summary.deg_up,
        "deg_down": summary.deg_down,
        "pca_separation": summary.pca_separation,
        "top_genes": summary.top_genes,
    }

    return {
        "n_samples": summary.n_samples,
        "groups": summary.groups,
        "deg_up": summary.deg_up,
        "deg_down": summary.deg_down,
        "top_genes": summary.top_genes,
        "pca_separation": summary.pca_separation,
        "qc_warnings": qc_warnings,
        "qc_critical": qc_critical,
        "realism_flags": realism_flags,
        "overall_qc_status": _qc_status(qc_report),
        "overall_realism": overall_realism,
        "interpretation_confidence": confidence,
        "limitation_text": limitation_text,
        "structured_summary": structured_summary,
    }


def _extract_json(text: str) -> dict:
    """Extract JSON from raw LLM output, handling markdown code blocks and <think> tags."""
    text = text.strip()
    # Strip <think>...</think> reasoning blocks (some models like MiniMax emit these)
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL).strip()
    # Strip markdown code fences: ```json ... ``` or ``` ... ```
    match = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.DOTALL)
    if match:
        text = match.group(1)
    else:
        # Find last (outermost) { ... } block
        brace_match = re.search(r"\{.*\}", text, re.DOTALL)
        if brace_match:
            text = brace_match.group(0)
    return json.loads(text)


def generate_interpretation(summary: AnalysisSummary, qc_report: Optional[dict[str, Any]] = None) -> LLMInterpretation:
    """Call LLM with summary JSON and parse response."""
    if not settings.llm_api_key:
        raise ValueError("LLM_API_KEY not configured")

    client = OpenAI(
        api_key=settings.llm_api_key,
        base_url=settings.llm_base_url,
    )

    llm_input = _build_llm_input(summary, qc_report)
    user_content = (
        "Input JSON:\n"
        f"{json.dumps(llm_input, indent=2)}\n\n"
        "Generate the output JSON following all constraints in the system prompt."
    )

    # Build request kwargs; some providers don't support response_format
    request_kwargs: dict = dict(
        model=settings.llm_model,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_content},
        ],
        temperature=0.3,
        timeout=60,
    )

    # Try with json_object response format first; fall back if provider rejects it
    try:
        response = client.chat.completions.create(
            **request_kwargs,
            response_format={"type": "json_object"},
        )
    except Exception as e:
        logger.warning("response_format=json_object not supported (%s), retrying without it", e)
        response = client.chat.completions.create(**request_kwargs)

    raw = response.choices[0].message.content or "{}"

    try:
        parsed = _extract_json(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"LLM returned invalid JSON: {exc}\nRaw: {raw[:500]}") from exc

    required_keys = {
        "pca_text",
        "deg_summary",
        "biological_insights",
        "data_quality",
        "next_steps",
        "methods_paragraph",
        "results_paragraph",
        "figure_legend",
    }
    missing = required_keys - set(parsed.keys())
    if missing:
        raise ValueError(f"LLM response missing keys: {missing}")

    payload = {k: str(parsed[k]) for k in required_keys}

    confidence = llm_input.get("interpretation_confidence", {})
    limitation_text = generateLimitationText(
        confidence.get("level", "MEDIUM"),
        confidence.get("reasons", []),
    )

    # Guardrail: always inject deterministic limitation context and sample size reference.
    payload["data_quality"] = f"{limitation_text} {payload['data_quality']} Sample size: n={summary.n_samples}.".strip()

    # Guardrail: for low confidence, shorten and soften interpretive language.
    if str(confidence.get("level", "")).upper() == "LOW":
        for key in ["pca_text", "deg_summary", "biological_insights", "next_steps"]:
            payload[key] = _cautious_rewrite(_limit_sentences(payload[key], max_sentences=3))

    return LLMInterpretation(**payload)
