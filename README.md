# RNA-seq Copilot

RNA-seq Copilot is a web system for count-matrix RNA-seq analysis with integrated DESeq2, QC diagnostics, realism checks, and explainable HTML reporting.

## 1. Scope

Input:
- gene-by-sample raw count matrix
- sample metadata table

Output:
- differential expression result table
- QC and realism diagnostics
- integrated summary JSON and HTML report

Out of scope:
- FASTQ preprocessing
- read alignment and quantification
- transcript assembly

## 2. Architecture

Data flow:
frontend (Next.js) -> backend (FastAPI) -> R pipeline (DESeq2, QC, stability) -> job artifacts (CSV, JSON, PNG, HTML)

Main directories:
- frontend: UI and polling workflow
- backend: API, orchestration, report rendering
- r_scripts: DESeq2 and plotting entry scripts
- R: stability and reporting modules
- backend/jobs: per-job runtime artifacts

## 3. Repository Layout

rna-seq-copilot/
  frontend/
  backend/
  r_scripts/
  R/
  sample_data/
  docker-compose.yml

## 4. Quick Start with Docker

Prerequisite:
- Docker Desktop or Docker Engine with Compose

Run:
1. cd rna-seq-copilot
2. docker compose up --build

Default endpoints:
- frontend: http://localhost:3000
- backend: http://localhost:8000
- API docs: http://localhost:8000/docs

## 5. Local Development

### 5.1 Prerequisites

- Python 3.11+
- Node.js 20+
- R 4.3+

Install required R packages:
1. Rscript -e "install.packages(c('BiocManager','optparse','jsonlite','ggplot2','ggrepel','pheatmap','RColorBrewer'))"
2. Rscript -e "BiocManager::install(c('DESeq2'))"

### 5.2 Start Backend

From rna-seq-copilot/backend:
1. pip install -r requirements.txt
2. uvicorn main:app --host 0.0.0.0 --port 8000

Alternative from repository root:
1. python dev_server.py

### 5.3 Start Frontend

From rna-seq-copilot/frontend:
1. npm install
2. create .env.local with NEXT_PUBLIC_API_BASE_URL=http://localhost:8000
3. npm run dev

## 6. Core API

- POST /upload
- POST /validate/{job_id}
- POST /run-analysis/{job_id}
- GET /results/{job_id}
- GET /results/{job_id}/deg.csv
- GET /results/{job_id}/plots/{plot_name}
- GET /report/{job_id}
- GET /healthz

Minimal run sequence:
1. upload counts and metadata
2. validate job
3. run analysis with formula and contrast
4. poll results
5. fetch report

## 7. Input Contract

counts file:
- first column is gene identifier
- other columns are sample identifiers
- values must be non-negative integer raw counts

metadata file:
- row names or first column must align to sample IDs in counts
- must include design variables referenced by formula and contrast factor
- at least two replicates per group are required for stable analysis

## 8. Stability Reporting Semantics

This project separates status logic from interpretation wording.

State machine (already enforced in code):
- stability_run_status: completed, limited, failed
- signal_state: strong_signal, weak_signal, no_detectable_signal
- stability_badge: high, moderate, low, low_signal, unknown

Important reporting behavior:
- no robust DEG signal does not imply technical failure
- limited status means analysis executed but DEG-based interpretation is constrained
- failed status is reserved for true technical failures

Metric clarity:
- public stability score may be N/A when robust DEG signal is absent
- composite effect-size consistency score is distinct from mean log2 fold-change correlation

## 9. Validation Dataset for Regression

Default local verification dataset used by this workspace:
- metadata: c:/Users/Peter/Downloads/metadata_ncrna_test.csv
- counts: c:/Users/Peter/Downloads/counts_ncrna_test.csv

Recommended quick check:
1. run analysis with the files above
2. inspect results summary JSON and report
3. verify low-signal wording is conservative and QC-first

## 10. Job Artifacts

Typical outputs in backend/jobs/{job_id}:
- results/deg_results.csv
- results/summary.json
- results/qc_report.json
- results/stability/
- plots/*.png
- report.html

## 11. Backend Configuration

Optional backend/.env keys:
- LLM_API_KEY
- LLM_BASE_URL
- LLM_MODEL
- R_SCRIPTS_DIR
- JOBS_DIR
- LOG_LEVEL
- RSCRIPT_PATH

If LLM settings are empty, the pipeline still runs and report generation remains available.

## 12. Troubleshooting

If backend starts but report generation fails:
- confirm R packages are installed
- confirm RSCRIPT_PATH is valid
- check backend logs and job-level results files

If Git push over HTTPS fails on port 443:
- switch remote to SSH and use github SSH endpoint over port 443

If stability shows limited:
- this is expected when robust DEG signal is absent or weak
- check QC metrics before biological interpretation

## 13. Production Notes

For production hardening:
- restrict CORS allow_origins in backend
- run frontend in build/start mode instead of dev mode
- persist and back up backend/jobs
- add service monitoring and log rotation

## 14. Current Limits

- no built-in multi-user auth
- no FASTQ-level upstream pipeline
- no enrichment analysis modules in current UI
