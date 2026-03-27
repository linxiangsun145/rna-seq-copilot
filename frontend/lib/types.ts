// ─── API types shared between frontend components ───────────────────────────

export type JobStatus =
  | "pending"
  | "validating"
  | "running"
  | "done"
  | "failed";

export interface ValidationIssue {
  level: "error" | "warning" | "info";
  field: string;
  message: string;
}

export interface ValidationReport {
  valid: boolean;
  n_genes: number;
  n_samples: number;
  sample_names: string[];
  groups: Record<string, string[]>;
  issues: ValidationIssue[];
}

export interface JobInfo {
  job_id: string;
  status: JobStatus;
  created_at: string;
  updated_at: string;
  error?: string;
  contrast?: string;
  formula?: string;
}

export interface DEGRow {
  gene_id: string;
  baseMean: number;
  log2FoldChange: number;
  lfcSE: number;
  stat?: number;  // dropped by lfcShrink, may be absent
  pvalue: number;
  padj: number;
}

export interface AnalysisSummary {
  n_samples: number;
  groups: string[];
  contrast: string;
  outliers: string[];
  pca_separation: "clear" | "weak" | "none";
  deg_up: number;
  deg_down: number;
  top_genes: string[];
  warnings: string[];
  data_issues: string[];
  realism_validation?: RealismValidation;
}

export interface RealismMetrics {
  canonical_genes_in_top20: number;
  canonical_fraction_top20: number;
  housekeeping_genes_in_top20: number;
  total_deg: number;
  fraction_p_lt_1e6: number;
  fraction_p_lt_1e3: number;
  fraction_p_gt_0_9: number;
  fraction_deg_abs_log2fc_gt_5: number;
  fraction_deg_abs_log2fc_gt_10: number;
}

export interface RealismValidation {
  realism_flags: string[];
  suspicious_patterns: string[];
  warnings: string[];
  critical: string[];
  metrics: RealismMetrics;
  overall_suspicion: "low" | "moderate" | "high";
}

export interface QCGroupBalance {
  status: "balanced" | "warning" | "critical";
  ratio: number;
}

export interface QCFlagItem {
  sample?: string;
  level: "warning" | "critical";
  value?: number;
  threshold?: number;
  rule: string;
  mean_correlation?: number;
  detail?: string;
  genes?: string[];
}

export interface QCReport {
  outliers: string[];
  low_quality_samples: string[];
  group_balance: QCGroupBalance;
  replicates_per_group: Record<string, number>;
  library_size_flags: QCFlagItem[];
  zero_fraction_flags: QCFlagItem[];
  correlation_flags: QCFlagItem[];
  batch_flags: QCFlagItem[];
  realism_flags: QCFlagItem[];
  qc_warnings: string[];
  qc_critical: string[];
  pca_variance: {
    PC1: number;
    PC2: number;
  };
}

export interface ResultsPayload {
  job_id: string;
  status: JobStatus;
  error?: string;
  summary: AnalysisSummary | null;
  qc_report?: QCReport | null;
  deg_table_url: string | null;
  plots: {
    pca?: string;
    volcano?: string;
    ma?: string;
    heatmap?: string;
    sample_distance_heatmap?: string;
    sample_correlation_heatmap?: string;
    library_size?: string;
    count_distribution?: string;
    zero_fraction?: string;
  };
  llm_interpretation: LLMInterpretation | null;
  report_url: string | null;
}

export interface LLMInterpretation {
  pca_text: string;
  deg_summary: string;
  biological_insights: string;
  data_quality: string;
  next_steps: string;
  methods_paragraph?: string;
  results_paragraph?: string;
  figure_legend?: string;
}
