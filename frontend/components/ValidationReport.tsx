"use client";

import type { ValidationReport } from "@/lib/types";
import { cn } from "@/lib/utils";

interface Props {
  report: ValidationReport;
}

function deriveStatus(report: ValidationReport): {
  label: "OK" | "Warning" | "Critical";
  style: string;
  hint: string;
} {
  const hasError = report.issues.some((x) => x.level === "error");
  const hasWarning = report.issues.some((x) => x.level === "warning");
  if (hasError) {
    return {
      label: "Critical",
      style: "bg-red-50 border-red-200 text-red-700",
      hint: "At least one issue can invalidate differential expression (DEG) analysis.",
    };
  }
  if (hasWarning) {
    return {
      label: "Warning",
      style: "bg-yellow-50 border-yellow-200 text-yellow-700",
      hint: "Analysis can run, but result reliability or interpretation may be limited.",
    };
  }
  return {
    label: "OK",
    style: "bg-green-50 border-green-200 text-green-700",
    hint: "No blocking issues detected for DEG analysis.",
  };
}

function explainImpact(field: string, message: string): { cause: string; impact: string } {
  const text = `${field} ${message}`.toLowerCase();

  if (text.includes("integer") || text.includes("raw") || text.includes("normalized") || text.includes("tpm") || text.includes("fpkm")) {
    return {
      cause: "Input may be normalized expression values instead of raw counts.",
      impact: "This may invalidate DESeq2 differential expression (DEG) results.",
    };
  }
  if (text.includes("sample") && (text.includes("mismatch") || text.includes("missing") || text.includes("not found"))) {
    return {
      cause: "Sample identifiers are inconsistent between counts and metadata.",
      impact: "Group assignments may be incorrect, leading to unreliable DEG and QC conclusions.",
    };
  }
  if (text.includes("replicate") || text.includes("group")) {
    return {
      cause: "Insufficient or imbalanced group replication.",
      impact: "Effect-size signal and DEG-level inference can become unstable and underpowered.",
    };
  }
  if (text.includes("zero") || text.includes("sparse") || text.includes("low count")) {
    return {
      cause: "High sparsity or weak count signal in the matrix.",
      impact: "Stability / reliability may be limited, with weak reproducibility under perturbation.",
    };
  }

  return {
    cause: "Input quality or design consistency requires review.",
    impact: "This can reduce confidence in biological interpretation and reliability claims.",
  };
}

const levelColor = {
  error: "bg-red-50 border-red-200 text-red-700",
  warning: "bg-yellow-50 border-yellow-200 text-yellow-700",
  info: "bg-blue-50 border-blue-200 text-blue-700",
};

const levelIcon = { error: "✗", warning: "⚠", info: "ℹ" };

export function ValidationReport({ report }: Props) {
  const status = deriveStatus(report);

  return (
    <div className="rounded-xl border bg-white p-6 space-y-4">
      <div className="flex items-center gap-3">
        <div
          className={cn(
            "w-8 h-8 rounded-full flex items-center justify-center text-sm font-bold",
            report.valid
              ? "bg-green-100 text-green-600"
              : "bg-red-100 text-red-600"
          )}
        >
          {report.valid ? "✓" : "✗"}
        </div>
        <div>
          <h3 className="font-semibold text-gray-900">
            Validation {report.valid ? "Passed" : "Failed"}
          </h3>
          <p className="text-xs text-gray-500">
            {report.n_genes.toLocaleString()} genes · {report.n_samples} samples
          </p>
        </div>
      </div>

      {/* Validation status */}
      <div>
        <p className="text-xs font-medium text-gray-500 uppercase mb-2">Validation Status</p>
        <div className={cn("rounded-lg border px-3 py-2 text-xs", status.style)}>
          <span className="font-semibold">{status.label}</span>
          <span className="ml-2">{status.hint}</span>
        </div>
      </div>

      {/* Groups */}
      <div>
        <p className="text-xs font-medium text-gray-500 uppercase mb-2">Groups</p>
        <div className="flex flex-wrap gap-2">
          {Object.entries(report.groups).map(([group, samples]) => (
            <span
              key={group}
              className="px-2 py-1 rounded-full bg-slate-100 text-xs text-gray-700"
            >
              {group} ({samples.length})
            </span>
          ))}
        </div>
      </div>

      {/* Issues */}
      {report.issues.length > 0 && (
        <div className="space-y-2">
          <p className="text-xs font-medium text-gray-500 uppercase">Issues</p>
          {report.issues.map((issue, i) => (
            <div
              key={i}
              className={cn(
                "flex gap-2 rounded-lg border px-3 py-2 text-xs",
                levelColor[issue.level]
              )}
            >
              <span className="font-bold">{levelIcon[issue.level]}</span>
              <div>
                <div>
                  <span className="font-medium">Issue: </span>
                  <span className="font-medium">{issue.field}</span>
                  <span> — {issue.message}</span>
                </div>
              </div>
            </div>
          ))}
        </div>
      )}

      {/* Impact */}
      {report.issues.length > 0 && (
        <div className="space-y-2">
          <p className="text-xs font-medium text-gray-500 uppercase">Impact</p>
          {report.issues.map((issue, i) => {
            const impact = explainImpact(issue.field, issue.message);
            return (
              <div key={`impact-${i}`} className="rounded-lg border bg-gray-50 border-gray-200 px-3 py-2 text-xs text-gray-700">
                <div>
                  <span className="font-medium">Issue: </span>
                  {issue.field} - {issue.message}
                </div>
                <div>
                  <span className="font-medium">Possible cause: </span>
                  {impact.cause}
                </div>
                <div>
                  <span className="font-medium">Impact: </span>
                  {impact.impact}
                </div>
              </div>
            );
          })}
        </div>
      )}

      {report.issues.length === 0 && (
        <p className="text-xs text-green-600">No issues detected.</p>
      )}
    </div>
  );
}
