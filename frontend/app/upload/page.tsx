"use client";

import { useState, useEffect, Suspense } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import { toast } from "sonner";
import { FileUploader } from "@/components/FileUploader";
import { ValidationReport } from "@/components/ValidationReport";
import { uploadFiles, validateJob, runAnalysis } from "@/lib/api";
import type { ValidationReport as VR } from "@/lib/types";

const DEMO_COUNTS_URL = "/demo/counts.csv";
const DEMO_META_URL = "/demo/metadata.csv";

type Step = "upload" | "configure" | "submitting";

function UploadPageInner() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const isDemo = searchParams.get("demo") === "true";

  const [countsFile, setCountsFile] = useState<File | null>(null);
  const [metaFile, setMetaFile] = useState<File | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);
  const [validation, setValidation] = useState<VR | null>(null);
  const [step, setStep] = useState<Step>("upload");
  const [loading, setLoading] = useState(false);

  // Design formula & contrast
  const [formula, setFormula] = useState("~ condition");
  const [contrastFactor, setContrastFactor] = useState("condition");
  const [contrastNumerator, setContrastNumerator] = useState("");
  const [contrastDenominator, setContrastDenominator] = useState("");

  // Load demo files if demo=true
  useEffect(() => {
    if (!isDemo) return;
    async function loadDemo() {
      const [cRes, mRes] = await Promise.all([
        fetch(DEMO_COUNTS_URL),
        fetch(DEMO_META_URL),
      ]);
      const [cBlob, mBlob] = await Promise.all([cRes.blob(), mRes.blob()]);
      setCountsFile(new File([cBlob], "demo_counts.csv", { type: "text/csv" }));
      setMetaFile(new File([mBlob], "demo_metadata.csv", { type: "text/csv" }));
      toast.info("Demo dataset loaded");
    }
    loadDemo().catch(() => toast.error("Failed to load demo files"));
  }, [isDemo]);

  // Auto-fill contrast from groups
  useEffect(() => {
    if (!validation) return;
    const groups = Object.keys(validation.groups);
    if (groups.length >= 2) {
      setContrastNumerator(groups[1]);
      setContrastDenominator(groups[0]);
    }
  }, [validation]);

  async function handleUploadAndValidate() {
    if (!countsFile || !metaFile) {
      toast.error("Please upload both files.");
      return;
    }
    setLoading(true);
    try {
      const { job_id } = await uploadFiles(countsFile, metaFile);
      setJobId(job_id);
      toast.info("Files uploaded. Validating…");
      const report = await validateJob(job_id);
      setValidation(report);
      if (report.valid) {
        setStep("configure");
        toast.success("Validation passed!");
      } else {
        toast.error("Validation failed — review issues below.");
      }
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : "Upload failed";
      toast.error(msg);
    } finally {
      setLoading(false);
    }
  }

  async function handleRunAnalysis() {
    if (!jobId || !contrastNumerator || !contrastDenominator) {
      toast.error("Please complete the contrast configuration.");
      return;
    }
    setStep("submitting");
    setLoading(true);
    try {
      const info = await runAnalysis(jobId, formula, [
        contrastFactor,
        contrastNumerator,
        contrastDenominator,
      ]);
      toast.success("Analysis started!");
      router.push(`/results/${info.job_id}`);
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : "Failed to start analysis";
      toast.error(msg);
      setStep("configure");
    } finally {
      setLoading(false);
    }
  }

  const groups = validation ? Object.keys(validation.groups) : [];

  return (
    <div className="max-w-2xl mx-auto px-4 py-12">
      <h1 className="text-3xl font-bold text-gray-900 mb-2">Upload Data</h1>
      <p className="text-gray-500 mb-8 text-sm">
        Provide a gene count matrix and sample metadata to begin analysis.
      </p>

      {/* Step 1: Upload */}
      <div className="bg-white rounded-xl border p-6 space-y-6 mb-6">
        <h2 className="font-semibold text-gray-900 flex items-center gap-2">
          <span className="w-6 h-6 rounded-full bg-blue-600 text-white text-xs flex items-center justify-center">1</span>
          Upload Files
        </h2>
        <FileUploader
          label="Count Matrix (genes × samples)"
          hint="counts.csv — rows = genes, columns = samples, values = raw integer counts"
          onFile={setCountsFile}
          file={countsFile}
        />
        <FileUploader
          label="Sample Metadata"
          hint="metadata.csv — rows = samples, columns = variables (e.g. condition, batch)"
          onFile={setMetaFile}
          file={metaFile}
        />

        <button
          onClick={handleUploadAndValidate}
          disabled={loading || !countsFile || !metaFile}
          className="w-full py-2.5 bg-blue-600 text-white rounded-lg font-medium text-sm hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
        >
          {loading && step === "upload" ? "Uploading & Validating…" : "Upload & Validate"}
        </button>
      </div>

      {/* Validation report */}
      {validation && (
        <div className="mb-6">
          <ValidationReport report={validation} />
        </div>
      )}

      {/* Step 2: Configure */}
      {step === "configure" && validation?.valid && (
        <div className="bg-white rounded-xl border p-6 space-y-5">
          <h2 className="font-semibold text-gray-900 flex items-center gap-2">
            <span className="w-6 h-6 rounded-full bg-blue-600 text-white text-xs flex items-center justify-center">2</span>
            Configure Analysis
          </h2>

          <div className="space-y-1">
            <label className="text-xs font-medium text-gray-700 uppercase">Design Formula</label>
            <input
              type="text"
              value={formula}
              onChange={(e) => setFormula(e.target.value)}
              className="w-full border rounded-lg px-3 py-2 text-sm font-mono focus:outline-none focus:ring-2 focus:ring-blue-400"
              placeholder="~ condition"
            />
            <p className="text-xs text-gray-400">DESeq2 design formula, e.g. <code>~ condition</code></p>
          </div>

          <div className="grid grid-cols-3 gap-3">
            <div className="space-y-1">
              <label className="text-xs font-medium text-gray-700 uppercase">Factor</label>
              <input
                type="text"
                value={contrastFactor}
                onChange={(e) => setContrastFactor(e.target.value)}
                className="w-full border rounded-lg px-3 py-2 text-sm font-mono focus:outline-none focus:ring-2 focus:ring-blue-400"
              />
            </div>
            <div className="space-y-1">
              <label className="text-xs font-medium text-gray-700 uppercase">Numerator</label>
              <select
                value={contrastNumerator}
                onChange={(e) => setContrastNumerator(e.target.value)}
                className="w-full border rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-400"
              >
                {groups.map((g) => <option key={g}>{g}</option>)}
              </select>
            </div>
            <div className="space-y-1">
              <label className="text-xs font-medium text-gray-700 uppercase">Denominator</label>
              <select
                value={contrastDenominator}
                onChange={(e) => setContrastDenominator(e.target.value)}
                className="w-full border rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-400"
              >
                {groups.map((g) => <option key={g}>{g}</option>)}
              </select>
            </div>
          </div>

          {/* Contrast preview */}
          {contrastNumerator && contrastDenominator && (
            <div className="text-xs bg-slate-50 border rounded-lg px-3 py-2 font-mono text-gray-600">
              Contrast:{" "}
              <span className="text-blue-700">{contrastFactor}</span>,{" "}
              <span className="text-red-600">{contrastNumerator}</span>
              {" "}vs{" "}
              <span className="text-blue-600">{contrastDenominator}</span>
            </div>
          )}

          <button
            onClick={handleRunAnalysis}
            disabled={loading || !contrastNumerator || !contrastDenominator}
            className="w-full py-2.5 bg-green-600 text-white rounded-lg font-medium text-sm hover:bg-green-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
          >
            {loading ? "Starting Analysis…" : "Run DESeq2 Analysis →"}
          </button>
        </div>
      )}
    </div>
  );
}

export default function UploadPage() {
  return (
    <Suspense fallback={
      <div className="max-w-2xl mx-auto px-4 py-12 text-center text-gray-400">Loading…</div>
    }>
      <UploadPageInner />
    </Suspense>
  );
}
