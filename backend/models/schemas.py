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
    warning_items: List["WarningItem"] = []
    realism_validation: Optional["RealismValidation"] = None
    confidence_assessment: Optional["ConfidenceAssessment"] = None
    ncrna_assessment: Optional["NcRNAAssessment"] = None
    deg_status: str = "no_signal"
    near_sig_genes: int = 0
    exploratory_deg_candidates: List["ExploratoryDEGCandidate"] = []
    sparse_ncrna_metrics: Optional["SparseNcRNAMetrics"] = None
    analysis_status: Optional["AnalysisStatus"] = None
    statistical_validation: List[str] = []
    stability_assessment: Optional["StabilityAssessment"] = None


class ConfidencePenaltyBreakdown(BaseModel):
    qc_penalty: float
    realism_penalty: float
    design_penalty: float
    stability_penalty: float = 0.0


class StabilityPenalty(BaseModel):
    stability_penalty: float = 0.0
    penalty_level: str = "none"
    penalty_reason: str = ""


class StabilityAssessment(BaseModel):
    stability_score: Optional[float] = None
    deg_stability_score: Optional[float] = None
    effect_stability_score: Optional[float] = None
    final_stability_score: Optional[float] = None
    stability_level: str = "unknown"  # high | moderate | low | low_signal | unknown
    stability_badge: str = "unknown"  # high | moderate | low | low_signal | not_applicable | unknown
    signal_state: str = "no_detectable_signal"  # no_detectable_signal | weak_signal | strong_signal | unknown
    stability_run_status: str = "limited"  # completed | limited | failed
    deg_metrics_applicable: bool = False
    stability_mode: str = "not_applicable"  # deg_based | effect_only | not_applicable
    influence_mode: str = "not_applicable"  # deg_based | effect_or_rank_based | not_applicable
    stability_headline: str = ""
    key_stability_findings: List[str] = []
    warnings: List[str] = []
    summary_text: str = ""
    directionality_text: str = ""
    top_rank_definition_text: str = ""
    pca_qc_conflict_text: str = ""
    sample_influence_text: str = ""
    reference_deg_count: Optional[int] = None
    mean_deg_recovery_rate: Optional[float] = None
    mean_top_n_overlap: Optional[float] = None
    top_rank_overlap: Optional[float] = None
    mean_log2fc_correlation: Optional[float] = None
    mean_log2fc_rmse: Optional[float] = None
    signal_collapse_fraction: Optional[float] = None
    sample_influence: Optional[Dict[str, Any]] = None
    stability_penalty: StabilityPenalty = StabilityPenalty()

    @field_validator(
        "stability_score",
        "deg_stability_score",
        "effect_stability_score",
        "final_stability_score",
        "mean_deg_recovery_rate",
        "mean_top_n_overlap",
        "top_rank_overlap",
        "mean_log2fc_correlation",
        "mean_log2fc_rmse",
        "signal_collapse_fraction",
        mode="before",
    )
    @classmethod
    def _coerce_optional_float(cls, v: Any) -> Optional[float]:
        if v in (None, "", "NA", "NaN", "nan", "N/A"):
            return None
        try:
            return float(v)
        except Exception:
            return None


class ConfidenceAssessment(BaseModel):
    confidence_score: int
    confidence_level: str  # HIGH | MEDIUM | LOW
    score_breakdown: ConfidencePenaltyBreakdown
    explanations: List[str]
    confidence_breakdown_text: List[str] = []
    confidence_explanation: str = ""
    confidence_penalty_explanations: Dict[str, List[str]] = {}
    confidence_validation: List[str] = []


class NcRNASummaryMetrics(BaseModel):
    n_total_deg: int = 0
    n_ncrna_deg: int = 0
    n_lncRNA_deg: int = 0
    n_miRNA_deg: int = 0
    n_protein_coding_deg: int = 0
    n_unknown_biotype_deg: int = 0
    n_circRNA_deg: int = 0
    n_other_ncRNA_deg: int = 0
    ncrna_fraction_deg: float = 0.0
    lncRNA_fraction_deg: float = 0.0
    miRNA_fraction_deg: float = 0.0


class NcRNAQCMetrics(BaseModel):
    ncrna_zero_fraction: float = 0.0
    ncrna_low_count_fraction: float = 0.0
    lncrna_low_count_fraction: float = 0.0
    miRNA_detected_count: int = 0
    lncRNA_detected_count: int = 0
    n_ncrna_features: int = 0


class NcRNACandidate(BaseModel):
    gene: str
    biotype: str
    priority_score: float
    padj: float
    abs_log2fc: float


class NcRNAAssessment(BaseModel):
    summary_metrics: NcRNASummaryMetrics = NcRNASummaryMetrics()
    qc_metrics: NcRNAQCMetrics = NcRNAQCMetrics()
    warnings: List[str] = []
    top_candidates: List[NcRNACandidate] = []


class ExploratoryDEGCandidate(BaseModel):
    gene: str
    log2fc: float
    padj: float
    pvalue: float
    rank_score: float
    label: str = "low-confidence DEG candidate"


class SparseNcRNAMetrics(BaseModel):
    zero_fraction_global: float = 0.0
    zero_fraction_ncrna: float = 0.0
    low_count_fraction_ncrna: float = 0.0
    detected_miRNA_count: int = 0
    detected_lncRNA_count: int = 0
    dataset_label: str = "ncRNA-rich dataset"
    warnings: List[str] = []


class AnalysisStatus(BaseModel):
    qc_status: str = "warning"  # pass | warning | critical
    deg_status: str = "no_signal"  # confirmed | low_effect_size | low_power | no_signal | failed_detection
    realism_status: str = "unstable"  # stable | unstable | not_applicable


class WarningItem(BaseModel):
    type: str  # qc | realism | statistical
    severity: str  # warning | critical
    code: str
    message: str
    level: Optional[str] = None
    group: Optional[str] = None
    sample: Optional[str] = None
    metric: Optional[str] = None
    evidence: Optional[str] = None
    metrics: Optional[Dict[str, Any]] = None


class RealismMetrics(BaseModel):
    canonical_genes_in_top20: int
    canonical_fraction_top20: float
    housekeeping_genes_in_top20: int
    total_deg: int
    fraction_p_lt_1e6: float
    fraction_p_lt_1e3: float
    fraction_p_gt_0_9: float
    fraction_deg_abs_log2fc_gt_5: float
    fraction_deg_abs_log2fc_gt_10: float


class RealismValidation(BaseModel):
    realism_flags: List[str]
    suspicious_patterns: List[str]
    warnings: List[str]
    critical: List[str]
    warning_items: List[WarningItem] = []
    metrics: RealismMetrics
    overall_suspicion: str  # low | moderate | high
    realism_status: str = "stable"
    gated_metrics: List[str] = []
    gating_note: str = ""


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
    qc_metrics: Dict[str, Any] = {}
    group_qc: Dict[str, Any] = {}
    per_sample_qc_metrics: List[Dict[str, Any]] = []
    warning_items: List[WarningItem] = []


# ─── LLM interpretation ──────────────────────────────────────────────────────

class LLMInterpretation(BaseModel):
    pca_text: str
    deg_summary: str
    biological_insights: str
    data_quality: str
    next_steps: str
    methods_paragraph: Optional[str] = None
    results_paragraph: Optional[str] = None
    figure_legend: Optional[str] = None


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
