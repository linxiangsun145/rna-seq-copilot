# RNA-seq Copilot

An intelligent web application for RNA-seq count matrix analysis, differential expression, and AI-powered biological interpretation.

> **Scope:** Count matrix level only. Does NOT support FASTQ processing, read alignment, or assembly.

---

## Project Structure

```
rna-seq-copilot/
├── frontend/               # Next.js 14 App Router + TypeScript + Tailwind + shadcn/ui
│   ├── app/
│   │   ├── layout.tsx
│   │   ├── page.tsx        # Home
│   │   ├── upload/page.tsx # Upload & config
│   │   └── results/[jobId]/page.tsx
│   ├── components/
│   │   ├── FileUploader.tsx
│   │   ├── ValidationReport.tsx
│   │   ├── PlotViewer.tsx
│   │   ├── DEGTable.tsx
│   │   └── AIInterpretation.tsx
│   ├── lib/
│   │   ├── api.ts          # Axios API client
│   │   ├── types.ts        # Shared TypeScript types
│   │   └── utils.ts
│   └── public/demo/        # Demo dataset (counts + metadata)
│
├── backend/                # Python FastAPI
│   ├── main.py             # App entry point
│   ├── config.py           # Settings (pydantic-settings)
│   ├── routers/
│   │   ├── upload.py       # POST /upload
│   │   ├── analysis.py     # POST /validate/{id}  POST /run-analysis/{id}
│   │   ├── results.py      # GET /results/{id}  GET /results/{id}/deg.csv
│   │   └── report.py       # GET /report/{id}
│   ├── services/
│   │   ├── validator.py    # Input validation logic
│   │   ├── r_runner.py     # Subprocess call to Rscript
│   │   ├── llm_client.py   # OpenAI-compatible LLM client
│   │   └── report_builder.py # Jinja2 HTML report
│   ├── models/schemas.py   # Pydantic models
│   ├── db/database.py      # SQLite job store
│   └── templates/report.html.j2
│
├── r_scripts/              # R analysis engine
│   ├── run_deseq2.R        # Main pipeline script
│   ├── plot_pca.R
│   ├── plot_volcano.R
│   ├── plot_ma.R
│   ├── plot_heatmap.R
│   └── summarize_results.R
│
├── sample_data/            # Demo input files
│   ├── counts.csv
│   └── metadata.csv
│
└── docker-compose.yml
```

---

## Quick Start

### Option A: Docker Compose (recommended)

```bash
cd rna-seq-copilot

# Copy and fill in LLM credentials
cp backend/.env.example backend/.env

# Start all services
docker-compose up --build
```

- Frontend: http://localhost:3000
- Backend API: http://localhost:8000
- API docs: http://localhost:8000/docs

---

### Option B: Manual setup

#### Prerequisites

| Tool | Version |
|------|---------|
| Python | ≥ 3.11 |
| Node.js | ≥ 20 |
| R | ≥ 4.3 |
| DESeq2 | via Bioconductor |

#### Install R dependencies

```r
install.packages(c("BiocManager", "optparse", "jsonlite", "ggplot2", "ggrepel",
                   "pheatmap", "RColorBrewer"))
BiocManager::install(c("DESeq2"))
```

#### Backend

```bash
cd backend
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
cp .env.example .env            # fill in LLM_API_KEY, LLM_MODEL etc.
uvicorn main:app --reload --port 8000
```

#### Frontend

```bash
cd frontend
npm install
echo "NEXT_PUBLIC_API_BASE_URL=http://localhost:8000" > .env.local
npm run dev
```

---

## API Reference

| Method | Path | Description |
|--------|------|-------------|
| `POST` | `/upload` | Upload `counts_file` + `metadata_file` (multipart) |
| `POST` | `/validate/{job_id}` | Validate uploaded files |
| `POST` | `/run-analysis/{job_id}` | Start DESeq2 pipeline (async) |
| `GET`  | `/results/{job_id}` | Poll for results |
| `GET`  | `/results/{job_id}/deg.csv` | Download DEG table |
| `GET`  | `/results/{job_id}/plots/{name}.png` | Get plot image |
| `GET`  | `/report/{job_id}` | Download HTML report |
| `GET`  | `/healthz` | Health check |

### Example: Run full analysis

```bash
# 1. Upload files
JOB=$(curl -s -X POST http://localhost:8000/upload \
  -F "counts_file=@sample_data/counts.csv" \
  -F "metadata_file=@sample_data/metadata.csv" | jq -r .job_id)

echo "Job ID: $JOB"

# 2. Validate
curl -s -X POST http://localhost:8000/validate/$JOB | jq .

# 3. Run analysis
curl -s -X POST http://localhost:8000/run-analysis/$JOB \
  -H "Content-Type: application/json" \
  -d '{"formula": "~ condition", "contrast": ["condition", "treated", "control"]}' | jq .

# 4. Poll results
curl -s http://localhost:8000/results/$JOB | jq .status

# 5. Download report when done
curl -s http://localhost:8000/report/$JOB -o report.html
```

---

## Environment Variables (backend)

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `LLM_API_KEY` | Yes (for AI) | — | OpenAI-compatible API key |
| `LLM_BASE_URL` | No | `https://api.openai.com/v1` | LLM endpoint |
| `LLM_MODEL` | No | `gpt-4o-mini` | Model name |
| `R_SCRIPTS_DIR` | No | `../r_scripts` | Path to R scripts |
| `JOBS_DIR` | No | `./jobs` | Job working directory |
| `LOG_LEVEL` | No | `INFO` | Logging verbosity |

> The LLM is completely optional. If `LLM_API_KEY` is not set, AI interpretation is skipped gracefully.

---

## Input File Format

### counts.csv

Rows = genes, columns = samples, values = **raw integer counts**.

```
gene_id,ctrl_1,ctrl_2,ctrl_3,treat_1,treat_2,treat_3
GAPDH,1823,1956,1788,1901,1834,1923
MYC,456,512,389,1234,1456,1189
...
```

### metadata.csv

Rows = samples, columns = experimental variables.

```
sample,condition,batch
ctrl_1,control,A
ctrl_2,control,A
treat_1,treated,A
...
```

**Requirements:**
- Sample names in `counts.csv` columns must match sample names in `metadata.csv` index
- At least 2 replicates per group
- Counts must be raw integers (not TPM/FPKM)

---

## Validation Rules

| Check | Severity |
|-------|----------|
| Sample name mismatch | Error |
| Non-integer counts | Error |
| Negative values | Error |
| Missing values | Error |
| Fewer than 2 replicates per group | Error |
| Possible TPM/FPKM detected (high mean) | Warning |
| Very sparse data (>80% zeros) | Warning |
| Only 2 replicates per group | Warning |

---

## Analysis Pipeline

```
Upload → Validate → DESeq2 → Plots → Summary JSON → LLM → HTML Report
```

1. **Filter**: genes with total raw count < 10 removed
2. **Normalise**: DESeq2 median-of-ratios normalisation
3. **Shrink**: LFC shrinkage with `ashr`
4. **Significance**: padj < 0.05, |log2FC| > 1
5. **Visualise**: PCA, volcano, MA plot, sample-distance heatmap (VST)
6. **Interpret**: LLM receives summary JSON only — never raw matrix

---

## Future Extensions (not implemented)

- GO / KEGG enrichment analysis
- GSEA
- Batch correction (ComBat, limma removeBatchEffect)
- Multi-omics integration
- Persistent user accounts
- FASTQ support (out of scope by design)
