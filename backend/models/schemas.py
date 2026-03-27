"""Pydantic schemas shared across routers and services."""
from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel, field_validator


class JobStatus(str, Enum):
    pending = "pending"
    validating = "validating"
    running = "running"
    done = "done"
    failed = "failed"


# ─── Validation ──────────────────────────────────────────────────────────────

class ValidationIssue(BaseModel):
    level: str  # "error" | "warning" | "info"
    field: str
    message: str


class ValidationReport(BaseModel):
    valid: bool
    n_genes: int
    n_samples: int
    sample_names: List[str]
    groups: Dict[str, List[str]]
    issues: List[ValidationIssue]


# ─── Job ─────────────────────────────────────────────────────────────────────

class JobInfo(BaseModel):
    job_id: str
    status: JobStatus
    created_at: str
    updated_at: str
    error: Optional[str] = None
    contrast: Optional[str] = None
    formula: Optional[str] = None


# ─── Analysis request ────────────────────────────────────────────────────────

class RunAnalysisRequest(BaseModel):
    formula: str = "~ condition"
    contrast: Tuple[str, str, str]  # (factor, numerator, denominator)

    @field_validator("formula")
    @classmethod
    def formula_starts_with_tilde(cls, v: str) -> str:
        if not v.strip().startswith("~"):
            raise ValueError("Design formula must start with ~")
        return v.strip()


# ─── Analysis summary ────────────────────────────────────────────────────────

class AnalysisSummary(BaseModel):
    n_samples: int
    groups: List[str]
    contrast: str
    outliers: List[str]
    pca_separation: str  # "clear" | "weak" | "none"
    deg_up: int
    deg_down: int
    top_genes: List[str]
    warnings: List[str]
    data_issues: List[str]


class GroupBalance(BaseModel):
    status: str  # balanced | warning | critical
    ratio: float


class SampleFlag(BaseModel):
    sample: str
    level: str  # warning | critical
    value: Optional[float] = None
    threshold: Optional[float] = None
    rule: str


class CorrelationFlag(BaseModel):
    sample: str
    level: str  # warning | critical
    mean_correlation: float
    threshold: float
    rule: str


class GenericFlag(BaseModel):
    level: str  # warning | critical
    rule: str
    detail: Optional[str] = None
    sample: Optional[str] = None
    value: Optional[float] = None
    threshold: Optional[float] = None
    genes: Optional[List[str]] = None
    ks_pvalue: Optional[float] = None
    mean_pvalue: Optional[float] = None
    fraction: Optional[float] = None
    z_median: Optional[float] = None
    iqr: Optional[float] = None


class QCReport(BaseModel):
    outliers: List[str]
    low_quality_samples: List[str]
    group_balance: GroupBalance
    replicates_per_group: Dict[str, int]
    library_size_flags: List[SampleFlag]
    zero_fraction_flags: List[SampleFlag]
    correlation_flags: List[CorrelationFlag]
    batch_flags: List[GenericFlag]
    realism_flags: List[GenericFlag]
    qc_warnings: List[str]
    qc_critical: List[str]
    pca_variance: Dict[str, float]


# ─── LLM interpretation ──────────────────────────────────────────────────────

class LLMInterpretation(BaseModel):
    pca_text: str
    deg_summary: str
    biological_insights: str
    data_quality: str
    next_steps: str


# ─── Results payload ─────────────────────────────────────────────────────────

class ResultsPayload(BaseModel):
    job_id: str
    status: JobStatus
    summary: Optional[AnalysisSummary] = None
    qc_report: Optional[QCReport] = None
    deg_table_url: Optional[str] = None
    plots: Dict[str, Optional[str]] = {}
    llm_interpretation: Optional[LLMInterpretation] = None
    report_url: Optional[str] = None
    error: Optional[str] = None
