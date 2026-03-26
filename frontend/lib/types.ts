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
}

export interface ResultsPayload {
  job_id: string;
  status: JobStatus;
  error?: string;
  summary: AnalysisSummary | null;
  deg_table_url: string | null;
  plots: {
    pca?: string;
    volcano?: string;
    ma?: string;
    heatmap?: string;
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
}
