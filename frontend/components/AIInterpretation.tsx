"use client";

import type { LLMInterpretation } from "@/lib/types";

interface Props {
  interpretation: LLMInterpretation;
}

const sections: { key: keyof LLMInterpretation; label: string; icon: string }[] = [
  { key: "pca_text", label: "PCA Interpretation", icon: "📍" },
  { key: "deg_summary", label: "DEG Summary", icon: "🔬" },
  { key: "biological_insights", label: "Biological Insights", icon: "💡" },
  { key: "data_quality", label: "Data Quality", icon: "✅" },
  { key: "next_steps", label: "Recommended Next Steps", icon: "→" },
];

export function AIInterpretation({ interpretation }: Props) {
  return (
    <div className="rounded-xl border bg-gradient-to-br from-violet-50 to-blue-50 p-6 space-y-5">
      <div className="flex items-center gap-2 mb-1">
        <div className="w-6 h-6 rounded-full bg-violet-600 flex items-center justify-center text-xs text-white">
          AI
        </div>
        <h2 className="font-semibold text-gray-900">AI Biological Interpretation</h2>
        <span className="ml-auto text-xs text-gray-400 bg-white px-2 py-0.5 rounded-full border">
          Based on summary JSON only
        </span>
      </div>

      <div className="space-y-4">
        {sections.map(({ key, label, icon }) => (
          <div key={key} className="bg-white rounded-lg border p-4">
            <div className="flex items-center gap-2 mb-2">
              <span>{icon}</span>
              <h3 className="text-sm font-medium text-gray-900">{label}</h3>
            </div>
            <p className="text-sm text-gray-600 leading-relaxed whitespace-pre-wrap">
              {interpretation[key]}
            </p>
          </div>
        ))}
      </div>

      <p className="text-xs text-gray-400 text-center">
        AI interpretation is generated from structured summary statistics — not raw count data.
        Always verify biological claims independently.
      </p>
    </div>
  );
}
