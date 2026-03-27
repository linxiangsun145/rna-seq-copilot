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
    return LLMInterpretation(**payload)
