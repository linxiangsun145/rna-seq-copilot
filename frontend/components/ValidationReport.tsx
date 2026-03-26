"use client";

import type { ValidationReport } from "@/lib/types";
import { cn } from "@/lib/utils";

interface Props {
  report: ValidationReport;
}

const levelColor = {
  error: "bg-red-50 border-red-200 text-red-700",
  warning: "bg-yellow-50 border-yellow-200 text-yellow-700",
  info: "bg-blue-50 border-blue-200 text-blue-700",
};

const levelIcon = { error: "✗", warning: "⚠", info: "ℹ" };

export function ValidationReport({ report }: Props) {
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
                <span className="font-medium">{issue.field}: </span>
                {issue.message}
              </div>
            </div>
          ))}
        </div>
      )}

      {report.issues.length === 0 && (
        <p className="text-xs text-green-600">No issues detected.</p>
      )}
    </div>
  );
}
