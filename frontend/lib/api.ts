import axios from "axios";
import type {
  ValidationReport,
  JobInfo,
  ResultsPayload,
} from "./types";

const BASE = process.env.NEXT_PUBLIC_API_BASE_URL || "http://localhost:8000";

const client = axios.create({ baseURL: BASE });

// ─── Upload files ────────────────────────────────────────────────────────────
export async function uploadFiles(
  countsFile: File,
  metadataFile: File
): Promise<{ job_id: string }> {
  const form = new FormData();
  form.append("counts_file", countsFile);
  form.append("metadata_file", metadataFile);
  const { data } = await client.post<{ job_id: string }>("/upload", form);
  return data;
}

// ─── Validate uploaded files ──────────────────────────────────────────────────
export async function validateJob(jobId: string): Promise<ValidationReport> {
  const { data } = await client.post<ValidationReport>(`/validate/${jobId}`);
  return data;
}

// ─── Run analysis ─────────────────────────────────────────────────────────────
export async function runAnalysis(
  jobId: string,
  formula: string,
  contrast: [string, string, string]
): Promise<JobInfo> {
  const { data } = await client.post<JobInfo>(`/run-analysis/${jobId}`, {
    formula,
    contrast,
  });
  return data;
}

// ─── Poll results ─────────────────────────────────────────────────────────────
export async function getResults(jobId: string): Promise<ResultsPayload> {
  const { data } = await client.get<ResultsPayload>(`/results/${jobId}`);
  return data;
}

// ─── Report URL ───────────────────────────────────────────────────────────────
export function reportUrl(jobId: string): string {
  return `${BASE}/report/${jobId}`;
}

// ─── Download DEG table URL ───────────────────────────────────────────────────
export function degDownloadUrl(jobId: string): string {
  return `${BASE}/results/${jobId}/deg.csv`;
}

// ─── Plot image URL ───────────────────────────────────────────────────────────
export function plotUrl(jobId: string, plotName: string): string {
  return `${BASE}/results/${jobId}/plots/${plotName}`;
}
