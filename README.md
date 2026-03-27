# RNA-seq Copilot

Product-grade web application for RNA-seq count-matrix analysis, differential expression, QC realism checks, and explainable HTML reporting.

## 1. Product Positioning

RNA-seq Copilot is designed for teams who need fast, reproducible analysis from processed count matrices.

- Input level: count matrix + sample metadata
- Core output: DEG table, QC diagnostics, realism signals, and integrated HTML report
- Optional capability: LLM-generated interpretation

Out of scope by design:

- FASTQ preprocessing
- read alignment / quantification
- transcript assembly

## 2. Key Features

- Guided workflow: Upload -> Validate -> Run analysis -> Review results -> Export report
- Strict validation gates before analysis starts
- DESeq2-based differential expression pipeline
- Multi-plot QC suite (PCA, volcano, MA, sample distance/correlation, library size, count distribution, zero fraction)
- Rule-based realism validation to flag suspicious statistical patterns
- Explainable "Overall Assessment" with explicit "Assessment Basis" bullets
- API-first backend suitable for web UI or automation

## 3. System Architecture

```text
frontend (Next.js) -> backend (FastAPI) -> R scripts (DESeq2/QC) -> job artifacts (CSV/JSON/PNG/HTML)
```

Main components:

- frontend/: Next.js 14 (upload, job polling, result views)
- backend/: FastAPI API, validation, orchestration, report rendering
- r_scripts/: analysis and plotting scripts
- backend/jobs/: per-job runtime artifacts and SQLite job store

## 4. Repository Layout

```text
rna-seq-copilot/
  frontend/
  backend/
  r_scripts/
  sample_data/
  docker-compose.yml
```

## 5. Quick Start (Recommended: Docker)

Prerequisites:

- Docker Engine + Docker Compose

Steps:

```bash
cd rna-seq-copilot
docker compose up --build
```

Default endpoints:

- Frontend: http://localhost:3000
- Backend API: http://localhost:8000
- OpenAPI docs: http://localhost:8000/docs

Notes:

- Backend LLM integration is optional. Service still works without LLM credentials.
- If you need LLM, create backend/.env and provide values (see section 8).

## 6. Local Development (No Docker)

### 6.1 Prerequisites

| Tool | Version |
|------|---------|
| Python | >= 3.11 |
| Node.js | >= 20 |
| R | >= 4.3 |

R packages:

```r
install.packages(c("BiocManager", "optparse", "jsonlite", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer"))
BiocManager::install(c("DESeq2"))
```

### 6.2 Backend

```bash
cd backend
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

pip install -r requirements.txt

# Create backend/.env manually if needed
uvicorn main:app --reload --port 8000
```

Alternative launcher from repo root:

```bash
python dev_server.py
```

### 6.3 Frontend

```bash
cd frontend
npm install

# Linux/macOS
echo "NEXT_PUBLIC_API_BASE_URL=http://localhost:8000" > .env.local
# Windows PowerShell
# "NEXT_PUBLIC_API_BASE_URL=http://localhost:8000" | Out-File -Encoding utf8 .env.local

npm run dev
```

## 7. API Reference

| Method | Path | Description |
|--------|------|-------------|
| POST | /upload | Upload counts_file + metadata_file |
| POST | /validate/{job_id} | Validate uploaded files |
| POST | /run-analysis/{job_id} | Start async analysis |
| GET | /results/{job_id} | Poll result payload |
| GET | /results/{job_id}/deg.csv | Download DEG CSV |
| GET | /results/{job_id}/plots/{name}.png | Get generated plot |
| GET | /report/{job_id} | Download report HTML |
| GET | /healthz | Health check |

Example run:

```bash
JOB=$(curl -s -X POST http://localhost:8000/upload \
  -F "counts_file=@sample_data/counts.csv" \
  -F "metadata_file=@sample_data/metadata.csv" | jq -r .job_id)

curl -s -X POST http://localhost:8000/validate/$JOB | jq .

curl -s -X POST http://localhost:8000/run-analysis/$JOB \
  -H "Content-Type: application/json" \
  -d '{"formula":"~ condition","contrast":["condition","treated","control"]}' | jq .

curl -s http://localhost:8000/results/$JOB | jq .status
curl -s http://localhost:8000/report/$JOB -o report.html
```

## 8. Backend Configuration

Create backend/.env when you need non-default settings.

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| LLM_API_KEY | No | empty | OpenAI-compatible key; if empty, LLM is skipped |
| LLM_BASE_URL | No | https://api.openai.com/v1 | LLM endpoint |
| LLM_MODEL | No | gpt-4o-mini | Model name |
| R_SCRIPTS_DIR | No | ../r_scripts | R scripts path |
| JOBS_DIR | No | ./jobs | Job workspace directory |
| LOG_LEVEL | No | INFO | Backend log verbosity |
| RSCRIPT_PATH | No | Rscript | Rscript executable path |

## 9. Input Contract

counts.csv requirements:

- First column must be gene_id
- Remaining columns are sample IDs
- Values must be raw integer counts (not TPM/FPKM)

metadata.csv requirements:

- Must contain sample column
- sample values must match counts.csv sample columns
- Must provide at least one grouping variable used by formula/contrast

Minimum statistical requirement:

- At least 2 replicates per group

## 10. Validation and Analysis Rules

Validation examples:

- Error: sample mismatch, non-integer/negative/missing counts, insufficient replicates
- Warning: potential normalized data signature, high sparsity, minimal replicate design

Analysis pipeline:

```text
Upload -> Validate -> DESeq2 -> Strict QC -> Realism Validation -> Summary JSON -> Optional LLM -> HTML report
```

Current DEG threshold in report summary:

- padj < 0.05 and |log2FoldChange| > 1

## 11. Job Outputs

Per job output directory:

```text
backend/jobs/<job_id>/
  results/deg_results.csv
  results/summary.json
  results/qc_report.json
  plots/*.png
  report.html
```

## 12. Production Deployment Guidance

For external sharing or team usage:

1. Deploy with Docker Compose on a Linux host.
2. Put Nginx/Caddy in front and enable HTTPS.
3. Restrict CORS allow_origins in backend for production domains.
4. Persist backend/jobs volume and implement backups.
5. Add process/service monitoring and log rotation.

Important current behavior:

- frontend Dockerfile starts Next.js in dev mode.
- For strict production mode, switch to npm run build + npm run start in frontend image.

## 13. Limitations and Roadmap

Known limitations:

- No built-in authentication/authorization
- No multi-user isolation beyond job IDs
- No enrichment modules (GO/KEGG/GSEA)
- No batch correction workflows in UI
- No FASTQ-level pipeline

Planned extension directions:

- Functional enrichment modules
- User/project-level access control
- Stronger observability and audit logs
- Optional cloud-native job queue execution
