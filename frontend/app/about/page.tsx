import Link from "next/link";

export default function AboutPage() {
  return (
    <div className="max-w-5xl mx-auto px-4 sm:px-6 lg:px-8 py-12">
      <h1 className="text-3xl font-bold text-gray-900 mb-3">About RNA-seq Copilot</h1>
      <p className="text-gray-600 mb-8">
        A reproducible RNA-seq analysis system that not only identifies differentially
        expressed genes, but also quantifies their reliability using perturbation-based
        stability analysis and QC-aware interpretation.
      </p>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <section className="bg-white border rounded-xl p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-3">What It Covers</h2>
          <ul className="list-disc pl-5 text-sm text-gray-600 space-y-2">
            <li>Processed count matrix and metadata validation.</li>
            <li>DESeq2 differential expression analysis.</li>
            <li>QC diagnostics, realism checks, and integrated HTML report.</li>
            <li>Optional LLM narrative layer when credentials are provided.</li>
          </ul>
        </section>

        <section className="bg-white border rounded-xl p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-3">Out of Scope</h2>
          <ul className="list-disc pl-5 text-sm text-gray-600 space-y-2">
            <li>FASTQ preprocessing and read alignment.</li>
            <li>Transcript assembly workflows.</li>
            <li>Upstream quantification pipelines.</li>
          </ul>
        </section>

        <section className="bg-white border rounded-xl p-6 lg:col-span-2">
          <h2 className="text-lg font-semibold text-gray-900 mb-3">Stability Semantics</h2>
          <ul className="list-disc pl-5 text-sm text-gray-600 space-y-2">
            <li>completed, limited, failed are runtime status states.</li>
            <li>
              limited means perturbation executed but DEG-level interpretation is constrained by
              weak or absent robust signal.
            </li>
            <li>failed is reserved for true technical failures.</li>
            <li>
              Composite effect-size consistency score is not identical to mean log2 fold-change
              correlation.
            </li>
          </ul>
        </section>
      </div>

      <div className="mt-8 flex items-center gap-4">
        <Link href="/upload" className="px-4 py-2 bg-blue-600 text-white rounded-lg text-sm font-medium hover:bg-blue-700">
          Start Analysis
        </Link>
        <Link href="/help" className="text-blue-700 text-sm font-medium hover:text-blue-800">
          Open Help
        </Link>
      </div>
    </div>
  );
}
