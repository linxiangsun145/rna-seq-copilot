import Link from "next/link";

export default function HelpPage() {
  return (
    <div className="max-w-5xl mx-auto px-4 sm:px-6 lg:px-8 py-12">
      <h1 className="text-3xl font-bold text-gray-900 mb-3">Help</h1>
      <p className="text-gray-600 mb-8">
        Quick operational guidance synchronized from README for local validation and troubleshooting.
      </p>

      <div className="space-y-6">
        <section className="bg-white border rounded-xl p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-3">Quick Start Checklist</h2>
          <ol className="list-decimal pl-5 text-sm text-gray-600 space-y-2">
            <li>Start backend on http://127.0.0.1:8000 and verify /docs is reachable.</li>
            <li>Start frontend on http://127.0.0.1:3000.</li>
            <li>Upload counts and metadata, then run Validate and Analysis.</li>
            <li>Open results page and exported report for interpretation review.</li>
          </ol>
        </section>

        <section className="bg-white border rounded-xl p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-3">Default Regression Dataset</h2>
          <div className="text-sm text-gray-700 bg-slate-50 border rounded-lg px-3 py-2 font-mono">
            metadata: c:/Users/Peter/Downloads/metadata_ncrna_test.csv
            <br />
            counts: c:/Users/Peter/Downloads/counts_ncrna_test.csv
          </div>
          <p className="text-xs text-gray-500 mt-2">
            Use this pair for reproducible checks after reporting-layer changes.
          </p>
        </section>

        <section className="bg-white border rounded-xl p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-3">Interpretation Principles</h2>
          <p className="text-xs text-gray-500 mb-2">
            These rules are automatically enforced during report generation.
          </p>
          <ul className="list-disc pl-5 text-sm text-gray-600 space-y-2">
            <li>No robust DEG signal should not be described as technical failure.</li>
            <li>
              When run status is limited, prioritize QC risk and within-group inconsistency before
              discussing effect-size consistency.
            </li>
            <li>
              Composite effect-size consistency score and mean log2 fold-change correlation are related
              but different metrics.
            </li>
          </ul>
        </section>

        <section className="bg-white border rounded-xl p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-3">Common Troubleshooting</h2>
          <ul className="list-disc pl-5 text-sm text-gray-600 space-y-2">
            <li>If backend fails, verify R dependencies and RSCRIPT_PATH configuration.</li>
            <li>If report wording looks outdated, ensure latest backend commit is deployed.</li>
            <li>If push fails on github.com:443, use SSH over port 443 configuration.</li>
          </ul>
        </section>
      </div>

      <div className="mt-8 flex items-center gap-4">
        <Link href="/" className="text-blue-700 text-sm font-medium hover:text-blue-800">
          Back to Home
        </Link>
        <Link href="/upload" className="px-4 py-2 bg-blue-600 text-white rounded-lg text-sm font-medium hover:bg-blue-700">
          Go to Upload
        </Link>
      </div>
    </div>
  );
}
