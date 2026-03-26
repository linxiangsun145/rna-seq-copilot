"""
LLM interpretation service.
Accepts a structured AnalysisSummary JSON — never raw counts.
"""
from __future__ import annotations

import json
import logging
import re

from openai import OpenAI

from config import settings
from models.schemas import AnalysisSummary, LLMInterpretation

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are a bioinformatics expert assisting wet-lab researchers in interpreting RNA-seq differential expression results.

You will receive a structured JSON summary of a DESeq2 analysis. Based ONLY on the values in this summary, you must generate clear, accurate biological interpretations.

RULES:
- Do NOT fabricate statistics, p-values, or gene names not present in the summary.
- Do NOT claim specific biology unless it is directly supported by the numbers.
- Use clear, non-technical language where possible.
- Be concise — each section should be 2–4 sentences.

Return a JSON object with exactly these keys:
{
  "pca_text": "<PCA interpretation>",
  "deg_summary": "<DEG summary>",
  "biological_insights": "<biology speculations with caveats>",
  "data_quality": "<data quality assessment>",
  "next_steps": "<recommended follow-up steps>"
}"""


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


def generate_interpretation(summary: AnalysisSummary) -> LLMInterpretation:
    """Call LLM with summary JSON and parse response."""
    if not settings.llm_api_key:
        raise ValueError("LLM_API_KEY not configured")

    client = OpenAI(
        api_key=settings.llm_api_key,
        base_url=settings.llm_base_url,
    )

    user_content = f"Here is the DESeq2 analysis summary:\n\n{json.dumps(summary.model_dump(), indent=2)}"

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

    required_keys = {"pca_text", "deg_summary", "biological_insights", "data_quality", "next_steps"}
    missing = required_keys - set(parsed.keys())
    if missing:
        raise ValueError(f"LLM response missing keys: {missing}")

    return LLMInterpretation(**{k: str(parsed[k]) for k in required_keys})
