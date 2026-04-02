import Link from "next/link";

const features = [
  {
    icon: "🛡️",
    title: "Result Reliability",
    desc: "Quantify how trustworthy your RNA-seq results are using perturbation-based stability analysis and QC-aware interpretation.",
  },
  {
    icon: "📁",
    title: "Upload Count Matrix",
    desc: "Import genes × samples CSV/TSV count matrices alongside sample metadata.",
  },
  {
    icon: "🧪",
    title: "DESeq2 Analysis",
    desc: "Automatic differential expression using DESeq2 — normalized, tested, ranked.",
  },
  {
    icon: "📊",
    title: "Interactive Visualizations",
    desc: "PCA plot, volcano plot, MA plot, and sample-distance heatmap generated instantly.",
  },
  {
    icon: "🤖",
    title: "QC-aware Interpretation",
    desc: "Interpretation is grounded in DEG, effect-size signal, stability, and QC (data quality) diagnostics.",
  },
  {
    icon: "📄",
    title: "Publication Report",
    desc: "Auto-generated HTML report with methods text, figures, and interpretations.",
  },
  {
    icon: "✅",
    title: "Data Validation",
    desc: "Automatic checks: sample name matching, integer counts, group balance, outliers.",
  },
];

export default function HomePage() {
  return (
    <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-20">
      {/* Hero */}
      <div className="text-center mb-20">
        <div className="inline-flex items-center gap-2 px-3 py-1 rounded-full bg-blue-50 border border-blue-200 text-blue-700 text-xs font-medium mb-6">
          <span className="w-1.5 h-1.5 rounded-full bg-blue-500 animate-pulse" />
          MVP · Count Matrix Analysis Only
        </div>
        <h1 className="text-5xl font-bold text-gray-900 tracking-tight mb-6">
          RNA-seq Copilot
        </h1>
        <p className="text-xl text-gray-500 max-w-2xl mx-auto mb-10">
          Upload your count matrix and metadata. Get DESeq2 differential expression,
          publication-ready plots, and statistically grounded and QC-aware biological interpretation —{" "}
          <span className="text-gray-900 font-medium">in minutes</span>.
        </p>
        <div className="flex items-center justify-center gap-4">
          <Link
            href="/upload"
            className="px-6 py-3 bg-blue-600 text-white rounded-lg font-medium hover:bg-blue-700 transition-colors shadow-sm"
          >
            Start Analysis →
          </Link>
          <Link
            href="/upload?demo=true"
            className="px-6 py-3 bg-white text-gray-700 rounded-lg font-medium hover:bg-gray-50 transition-colors border shadow-sm"
          >
            Try Demo Dataset
          </Link>
        </div>
      </div>

      {/* Features */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 mb-20">
        {features.map((f) => (
          <div
            key={f.title}
            className="bg-white rounded-xl border p-6 hover:shadow-md transition-shadow"
          >
            <div className="text-3xl mb-3">{f.icon}</div>
            <h3 className="font-semibold text-gray-900 mb-2">{f.title}</h3>
            <p className="text-sm text-gray-500 leading-relaxed">{f.desc}</p>
          </div>
        ))}
      </div>

      {/* Workflow */}
      <div className="bg-white rounded-2xl border p-10 max-w-3xl mx-auto">
        <h2 className="text-2xl font-bold text-gray-900 text-center mb-8">How it works</h2>
        <div className="space-y-6">
          {[
            ["1", "Upload", "counts.csv + metadata.csv with group annotations"],
            ["2", "Validate", "Automatic sanity checks on sample names, count types, group balance"],
            ["3", "Analyse", "DESeq2 runs in R — log2FC, p-values, padj computed"],
            ["4", "Visualise", "PCA, volcano, MA plot, heatmap generated as PNG"],
            ["5", "Interpret", "QC-aware interpretation links differential expression (DEG), effect-size signal, and stability / reliability"],
            ["6", "Report", "Download HTML report with methods, figures, and interpretation"],
          ].map(([num, title, desc]) => (
            <div key={num} className="flex gap-4">
              <div className="w-8 h-8 rounded-full bg-blue-600 text-white flex items-center justify-center text-sm font-bold flex-shrink-0 mt-0.5">
                {num}
              </div>
              <div>
                <div className="font-semibold text-gray-900">{title}</div>
                <div className="text-sm text-gray-500">{desc}</div>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Documentation Cards */}
      <div className="max-w-5xl mx-auto mt-16">
        <div className="flex items-end justify-between mb-5">
          <div>
            <h2 className="text-2xl font-bold text-gray-900">Project Documentation</h2>
            <p className="text-sm text-gray-500 mt-1">
              Key guidance synchronized from README for quick onboarding.
            </p>
          </div>
          <div className="flex items-center gap-3 text-sm">
            <Link href="/about" className="text-blue-700 hover:text-blue-800 font-medium">About</Link>
            <Link href="/help" className="text-blue-700 hover:text-blue-800 font-medium">Help</Link>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-5">
          <div className="bg-white rounded-xl border p-6">
            <h3 className="font-semibold text-gray-900 mb-3">Scope and Architecture</h3>
            <ul className="text-sm text-gray-600 space-y-2 list-disc pl-5">
              <li>Input is processed count matrix plus sample metadata.</li>
              <li>Pipeline combines DESeq2, strict QC, realism checks, and HTML report output.</li>
              <li>FASTQ preprocessing, alignment, and transcript assembly are intentionally out of scope.</li>
            </ul>
          </div>

          <div className="bg-white rounded-xl border p-6">
            <h3 className="font-semibold text-gray-900 mb-3">Stability Interpretation Policy</h3>
            <ul className="text-sm text-gray-600 space-y-2 list-disc pl-5">
              <li>Limited status means analysis executed but DEG-level interpretation is constrained.</li>
              <li>No robust DEG signal is not treated as technical failure.</li>
              <li>
                Composite effect-size consistency is distinct from mean log2 fold-change correlation.
              </li>
            </ul>
          </div>

          <div className="bg-white rounded-xl border p-6 lg:col-span-2">
            <h3 className="font-semibold text-gray-900 mb-3">Local Validation Dataset</h3>
            <p className="text-sm text-gray-600 mb-2">
              Recommended regression dataset used in this workspace:
            </p>
            <div className="text-sm text-gray-700 bg-slate-50 border rounded-lg px-3 py-2 font-mono">
              metadata: c:/Users/Peter/Downloads/metadata_ncrna_test.csv
              <br />
              counts: c:/Users/Peter/Downloads/counts_ncrna_test.csv
            </div>
          </div>
        </div>
      </div>

      {/* Disclaimer */}
      <p className="text-center text-xs text-gray-400 mt-16">
        RNA-seq Copilot does not support FASTQ processing or read alignment.
        It operates exclusively on pre-computed count matrices.
      </p>
    </div>
  );
}
