"use client";

import { useEffect, useState, useRef } from "react";
import { useParams } from "next/navigation";
import { toast } from "sonner";
import { getResults, degDownloadUrl, plotUrl, reportUrl } from "@/lib/api";
import type { ResultsPayload, DEGRow } from "@/lib/types";
import { PlotViewer } from "@/components/PlotViewer";
import { DEGTable } from "@/components/DEGTable";
import { AIInterpretation } from "@/components/AIInterpretation";

const POLL_INTERVAL = 3000;

async function downloadBlob(url: string, filename: string) {
  try {
    const res = await fetch(url);
    const blob = await res.blob();
    const blobUrl = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = blobUrl;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(blobUrl);
  } catch {
    window.open(url, "_blank");
  }
}

export default function ResultsPage() {
  const { jobId } = useParams<{ jobId: string }>();
  const [data, setData] = useState<ResultsPayload | null>(null);
  const [degRows, setDegRows] = useState<DEGRow[]>([]);
  const [activeTab, setActiveTab] = useState<"qc" | "plots" | "degs" | "ai">("qc");
  const intervalRef = useRef<ReturnType<typeof setInterval> | null>(null);

  useEffect(() => {
    if (!jobId) return;

    async function poll() {
      try {
        const result = await getResults(jobId);
        setData(result);
        if (result.status === "done") {
          clearInterval(intervalRef.current!);
          // Fetch DEG CSV
          const csvRes = await fetch(degDownloadUrl(jobId));
          if (csvRes.ok) {
            const text = await csvRes.text();
            setDegRows(parseDEGCsv(text));
          }
        }
        if (result.status === "failed") {
          clearInterval(intervalRef.current!);
          toast.error("Analysis failed" + (result.error ? ": " + result.error : ""));
        }
      } catch {
        // Silently retry
      }
    }

    poll();
    intervalRef.current = setInterval(poll, POLL_INTERVAL);
    return () => clearInterval(intervalRef.current!);
  }, [jobId]);

  const isDone = data?.status === "done";
  const isRunning = data?.status === "running" || data?.status === "validating" || data?.status === "pending";

  return (
    <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-10">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Analysis Results</h1>
          <p className="text-xs text-gray-400 font-mono mt-0.5">Job: {jobId}</p>
        </div>
        <div className="flex items-center gap-3">
          {data && (
            <span
              className={`px-3 py-1 rounded-full text-xs font-medium ${
                isDone
                  ? "bg-green-100 text-green-700"
                  : data.status === "failed"
                  ? "bg-red-100 text-red-700"
                  : "bg-blue-100 text-blue-700"
              }`}
            >
              {data.status.toUpperCase()}
            </span>
          )}
          {isDone && (
            <button
              onClick={() => downloadBlob(reportUrl(jobId), `report_${jobId.slice(0, 8)}.html`)}
              className="px-4 py-1.5 bg-blue-600 text-white text-xs rounded-lg hover:bg-blue-700 transition-colors"
            >
              Download Report
            </button>
          )}
        </div>
      </div>

      {/* Loading state */}
      {isRunning && (
        <div className="flex flex-col items-center justify-center py-20 space-y-4">
          <div className="w-10 h-10 border-4 border-blue-200 border-t-blue-600 rounded-full animate-spin" />
          <p className="text-gray-500 text-sm">
            {data?.status === "running"
              ? "Running DESeq2 analysis…"
              : "Preparing analysis…"}
          </p>
          <p className="text-xs text-gray-400">This may take 1–3 minutes depending on dataset size.</p>
        </div>
      )}

      {/* Results tabs */}
      {isDone && data && (
        <>
          {/* Summary cards */}
          {data.summary && (
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
              <SummaryCard label="Samples" value={data.summary.n_samples} />
              <SummaryCard
                label="Upregulated"
                value={data.summary.deg_up}
                color="text-red-600"
              />
              <SummaryCard
                label="Downregulated"
                value={data.summary.deg_down}
                color="text-blue-600"
              />
              <SummaryCard
                label="PCA Separation"
                value={data.summary.pca_separation}
                color={data.summary.pca_separation === "clear" ? "text-green-600" : "text-yellow-600"}
              />
            </div>
          )}

          {/* Tabs */}
          <div className="border-b mb-6">
            <nav className="flex gap-6 text-sm">
              {(["qc", "plots", "degs", "ai"] as const).map((t) => (
                <button
                  key={t}
                  onClick={() => setActiveTab(t)}
                  className={`pb-2 font-medium capitalize border-b-2 transition-colors ${
                    activeTab === t
                      ? "border-blue-600 text-blue-600"
                      : "border-transparent text-gray-500 hover:text-gray-700"
                  }`}
                >
                  {t === "qc" ? "QC Summary" : t === "degs" ? "DEG Table" : t === "ai" ? "AI Interpretation" : "Plots"}
                </button>
              ))}
            </nav>
          </div>

          {/* QC tab */}
          {activeTab === "qc" && data.summary && (
            <div className="space-y-4">
              {data.qc_report && (
                <div className="bg-white rounded-xl border p-5">
                  <div className="flex items-center justify-between mb-3">
                    <h3 className="font-medium text-gray-900">Strict QC Report</h3>
                    <span
                      className={`px-2.5 py-1 rounded-full text-xs font-semibold uppercase ${
                        data.qc_report.group_balance.status === "critical"
                          ? "bg-red-100 text-red-700"
                          : data.qc_report.group_balance.status === "warning"
                          ? "bg-yellow-100 text-yellow-700"
                          : "bg-green-100 text-green-700"
                      }`}
                    >
                      {data.qc_report.group_balance.status}
                    </span>
                  </div>

                  <div className="grid grid-cols-1 md:grid-cols-2 gap-3 mb-4">
                    <MiniMetric
                      label="Group balance ratio"
                      value={data.qc_report.group_balance.ratio.toFixed(3)}
                    />
                    <MiniMetric
                      label="Replicate groups"
                      value={`${Object.keys(data.qc_report.replicates_per_group || {}).length}`}
                    />
                    <MiniMetric
                      label="Outlier samples"
                      value={`${data.qc_report.outliers.length}`}
                    />
                    <MiniMetric
                      label="Low-quality samples"
                      value={`${data.qc_report.low_quality_samples.length}`}
                    />
                    <MiniMetric
                      label="PCA variance PC1"
                      value={`${(100 * (data.qc_report.pca_variance?.PC1 ?? 0)).toFixed(2)}%`}
                    />
                    <MiniMetric
                      label="PCA variance PC2"
                      value={`${(100 * (data.qc_report.pca_variance?.PC2 ?? 0)).toFixed(2)}%`}
                    />
                    <MiniMetric
                      label="Library size flags"
                      value={`${data.qc_report.library_size_flags.length}`}
                    />
                    <MiniMetric
                      label="Zero-fraction flags"
                      value={`${data.qc_report.zero_fraction_flags.length}`}
                    />
                  </div>

                  {data.qc_report.qc_critical.length > 0 && (
                    <div className="bg-red-50 border border-red-200 rounded-lg p-3 mb-3">
                      <h4 className="text-xs font-semibold text-red-800 mb-1">QC Critical</h4>
                      <ul className="space-y-1">
                        {data.qc_report.qc_critical.map((c, i) => (
                          <li key={i} className="text-xs text-red-700">✗ {c}</li>
                        ))}
                      </ul>
                    </div>
                  )}

                  {data.qc_report.qc_warnings.length > 0 && (
                    <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-3 mb-3">
                      <h4 className="text-xs font-semibold text-yellow-800 mb-1">QC Warnings</h4>
                      <ul className="space-y-1">
                        {data.qc_report.qc_warnings.map((w, i) => (
                          <li key={i} className="text-xs text-yellow-700">⚠ {w}</li>
                        ))}
                      </ul>
                    </div>
                  )}

                  {Object.keys(data.qc_report.replicates_per_group || {}).length > 0 && (
                    <div className="bg-gray-50 border border-gray-200 rounded-lg p-3">
                      <h4 className="text-xs font-semibold text-gray-700 mb-1">Replicates Per Group</h4>
                      <div className="flex flex-wrap gap-2">
                        {Object.entries(data.qc_report.replicates_per_group).map(([grp, n]) => (
                          <span key={grp} className="px-2 py-1 bg-white border border-gray-200 rounded text-xs text-gray-700">
                            {grp}: {n}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              )}

              {data.summary.realism_validation && (
                <div className="bg-white rounded-xl border p-5">
                  <div className="flex items-center justify-between mb-3">
                    <h3 className="font-medium text-gray-900">Data Realism Validation</h3>
                    <span
                      className={`px-2.5 py-1 rounded-full text-xs font-semibold uppercase ${
                        data.summary.realism_validation.overall_suspicion === "high"
                          ? "bg-red-100 text-red-700"
                          : data.summary.realism_validation.overall_suspicion === "moderate"
                          ? "bg-yellow-100 text-yellow-700"
                          : "bg-green-100 text-green-700"
                      }`}
                    >
                      {data.summary.realism_validation.overall_suspicion}
                    </span>
                  </div>

                  <div className="grid grid-cols-1 md:grid-cols-2 gap-3 mb-4">
                    <MiniMetric label="Canonical in Top20" value={`${data.summary.realism_validation.metrics.canonical_genes_in_top20}`} />
                    <MiniMetric label="Canonical Fraction Top20" value={data.summary.realism_validation.metrics.canonical_fraction_top20.toFixed(3)} />
                    <MiniMetric label="Housekeeping in Top20" value={`${data.summary.realism_validation.metrics.housekeeping_genes_in_top20}`} />
                    <MiniMetric label="Total DEG (padj < 0.05)" value={`${data.summary.realism_validation.metrics.total_deg}`} />
                    <MiniMetric label="Fraction p < 1e-6" value={data.summary.realism_validation.metrics.fraction_p_lt_1e6.toFixed(3)} />
                    <MiniMetric label="Fraction p < 1e-3" value={data.summary.realism_validation.metrics.fraction_p_lt_1e3.toFixed(3)} />
                    <MiniMetric label="Fraction p > 0.9" value={data.summary.realism_validation.metrics.fraction_p_gt_0_9.toFixed(3)} />
                    <MiniMetric label="DEG |log2FC| > 5" value={data.summary.realism_validation.metrics.fraction_deg_abs_log2fc_gt_5.toFixed(3)} />
                  </div>

                  {data.summary.realism_validation.critical.length > 0 && (
                    <div className="bg-red-50 border border-red-200 rounded-lg p-3 mb-3">
                      <h4 className="text-xs font-semibold text-red-800 mb-1">Critical</h4>
                      <ul className="space-y-1">
                        {data.summary.realism_validation.critical.map((c, i) => (
                          <li key={i} className="text-xs text-red-700">✗ {c}</li>
                        ))}
                      </ul>
                    </div>
                  )}

                  {data.summary.realism_validation.warnings.length > 0 && (
                    <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-3 mb-3">
                      <h4 className="text-xs font-semibold text-yellow-800 mb-1">Warnings</h4>
                      <ul className="space-y-1">
                        {data.summary.realism_validation.warnings.map((w, i) => (
                          <li key={i} className="text-xs text-yellow-700">⚠ {w}</li>
                        ))}
                      </ul>
                    </div>
                  )}

                  {data.summary.realism_validation.suspicious_patterns.length > 0 && (
                    <div className="bg-gray-50 border border-gray-200 rounded-lg p-3">
                      <h4 className="text-xs font-semibold text-gray-700 mb-1">Suspicious Patterns</h4>
                      <ul className="space-y-1">
                        {data.summary.realism_validation.suspicious_patterns.map((p, i) => (
                          <li key={i} className="text-xs text-gray-600">• {p}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                </div>
              )}

              {data.summary.warnings.length > 0 && (
                <div className="bg-yellow-50 border border-yellow-200 rounded-xl p-4">
                  <h3 className="text-sm font-semibold text-yellow-800 mb-2">Warnings</h3>
                  <ul className="space-y-1">
                    {data.summary.warnings.map((w, i) => (
                      <li key={i} className="text-xs text-yellow-700">⚠ {w}</li>
                    ))}
                  </ul>
                </div>
              )}
              {data.summary.data_issues.length > 0 && (
                <div className="bg-red-50 border border-red-200 rounded-xl p-4">
                  <h3 className="text-sm font-semibold text-red-800 mb-2">Data Issues</h3>
                  <ul className="space-y-1">
                    {data.summary.data_issues.map((w, i) => (
                      <li key={i} className="text-xs text-red-700">✗ {w}</li>
                    ))}
                  </ul>
                </div>
              )}
              <div className="bg-white rounded-xl border p-5">
                <h3 className="font-medium text-gray-900 mb-3">Top Differentially Expressed Genes</h3>
                <div className="flex flex-wrap gap-2">
                  {data.summary.top_genes.map((g) => (
                    <span key={g} className="px-2 py-1 bg-blue-50 text-blue-700 rounded text-xs font-mono">
                      {g}
                    </span>
                  ))}
                </div>
              </div>
            </div>
          )}

          {/* Plots tab */}
          {activeTab === "plots" && (
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              {data.plots.pca && (
                <PlotViewer
                  src={plotUrl(jobId, "pca.png")}
                  title="PCA Plot"
                  description="Principal Component Analysis of sample variation"
                />
              )}
              {data.plots.volcano && (
                <PlotViewer
                  src={plotUrl(jobId, "volcano.png")}
                  title="Volcano Plot"
                  description="log2FC vs –log10(padj), significant genes highlighted"
                />
              )}
              {data.plots.ma && (
                <PlotViewer
                  src={plotUrl(jobId, "ma.png")}
                  title="MA Plot"
                  description="Log fold change vs mean expression"
                />
              )}
              {data.plots.heatmap && (
                <PlotViewer
                  src={plotUrl(jobId, "heatmap.png")}
                  title="Sample Distance Heatmap"
                  description="Euclidean distances between samples in VST space"
                />
              )}
              {data.plots.sample_distance_heatmap && (
                <PlotViewer
                  src={plotUrl(jobId, "sample_distance_heatmap.png")}
                  title="QC: Sample Distance Heatmap"
                  description="Strict QC distance matrix using VST counts"
                />
              )}
              {data.plots.sample_correlation_heatmap && (
                <PlotViewer
                  src={plotUrl(jobId, "sample_correlation_heatmap.png")}
                  title="QC: Sample Correlation Heatmap"
                  description="Pearson correlation matrix across samples"
                />
              )}
              {data.plots.library_size && (
                <PlotViewer
                  src={plotUrl(jobId, "library_size.png")}
                  title="QC: Library Size"
                  description="Total counts per sample with strict threshold checks"
                />
              )}
              {data.plots.count_distribution && (
                <PlotViewer
                  src={plotUrl(jobId, "count_distribution.png")}
                  title="QC: Count Distribution"
                  description="Sample-level count density in log2 scale"
                />
              )}
              {data.plots.zero_fraction && (
                <PlotViewer
                  src={plotUrl(jobId, "zero_fraction.png")}
                  title="QC: Zero Fraction"
                  description="Fraction of zero-count genes per sample"
                />
              )}
            </div>
          )}

          {/* DEG table tab */}
          {activeTab === "degs" && (
            <DEGTable rows={degRows} downloadUrl={degDownloadUrl(jobId)} />
          )}

          {/* AI tab */}
          {activeTab === "ai" && data.llm_interpretation && (
            <AIInterpretation interpretation={data.llm_interpretation} />
          )}
          {activeTab === "ai" && !data.llm_interpretation && (
            <div className="text-center py-12 text-gray-400 text-sm">
              AI interpretation not available for this job.
              Check that <code>LLM_API_KEY</code> is configured in the backend.
            </div>
          )}
        </>
      )}
    </div>
  );
}

function SummaryCard({
  label,
  value,
  color = "text-gray-900",
}: {
  label: string;
  value: string | number;
  color?: string;
}) {
  return (
    <div className="bg-white rounded-xl border p-4">
      <p className="text-xs text-gray-500 mb-1">{label}</p>
      <p className={`text-2xl font-bold capitalize ${color}`}>{value}</p>
    </div>
  );
}

function MiniMetric({ label, value }: { label: string; value: string }) {
  return (
    <div className="bg-gray-50 rounded-lg border border-gray-200 p-3">
      <p className="text-[11px] text-gray-500 mb-1">{label}</p>
      <p className="text-sm font-semibold text-gray-800">{value}</p>
    </div>
  );
}

function parseDEGCsv(csv: string): DEGRow[] {
  const lines = csv.trim().split("\n");
  if (lines.length < 2) return [];
  const header = lines[0].split(",");
  const idx = (col: string) => header.indexOf(col);
  return lines.slice(1).map((line) => {
    const cols = line.split(",");
    const statIdx = idx("stat");
    return {
      gene_id: cols[idx("gene_id")] || cols[0],
      baseMean: parseFloat(cols[idx("baseMean")]),
      log2FoldChange: parseFloat(cols[idx("log2FoldChange")]),
      lfcSE: parseFloat(cols[idx("lfcSE")]),
      ...(statIdx >= 0 ? { stat: parseFloat(cols[statIdx]) } : {}),
      pvalue: parseFloat(cols[idx("pvalue")]),
      padj: parseFloat(cols[idx("padj")]),
    };
  });
}
