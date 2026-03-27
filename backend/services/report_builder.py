"""
HTML report builder — renders a Jinja2 template with analysis artefacts.
"""
from __future__ import annotations

import base64
import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from models.schemas import AnalysisSummary, LLMInterpretation

logger = logging.getLogger(__name__)

TEMPLATES_DIR = Path(__file__).parent.parent / "templates"

CANONICAL_GENES = {
    "TP53", "MYC", "EGFR", "VEGFA", "IL6", "CXCL8", "GAPDH", "ACTB", "AKT1", "PTEN"
}
HOUSEKEEPING_GENES = {"GAPDH", "ACTB", "RPL13A", "B2M", "HPRT1", "TUBB"}


def _img_to_base64(path: Path) -> Optional[str]:
    if path.exists():
        return base64.b64encode(path.read_bytes()).decode()
    return None


def _to_float(value: Any) -> Optional[float]:
    try:
        if value in (None, "", "NA", "NaN"):
            return None
        return float(value)
    except Exception:
        return None


def _load_top_genes_table(job_dir: Path, n: int = 20) -> list[dict[str, Any]]:
    deg_path = job_dir / "results" / "deg_results.csv"
    if not deg_path.exists():
        return []

    rows: list[dict[str, Any]] = []
    try:
        with deg_path.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                if i >= n:
                    break

                gene = str(row.get("gene_id", "")).strip()
                gene_upper = gene.upper()
                l2fc = _to_float(row.get("log2FoldChange"))
                padj = _to_float(row.get("padj"))

                rows.append(
                    {
                        "gene": gene,
                        "log2FoldChange": l2fc,
                        "padj": padj,
                        "is_canonical": gene_upper in CANONICAL_GENES,
                        "is_housekeeping": gene_upper in HOUSEKEEPING_GENES,
                    }
                )
    except Exception:
        return []

    return rows


def _overall_data_quality(qc: Optional[dict[str, Any]]) -> str:
    if not qc:
        return "unknown"
    n_critical = len(qc.get("qc_critical", []) or [])
    n_warning = len(qc.get("qc_warnings", []) or [])
    if n_critical > 0:
        return "high"
    if n_warning >= 2:
        return "moderate"
    return "low"


def _grouped_warnings(
    summary: dict[str, Any],
    qc: Optional[dict[str, Any]],
    realism: Optional[dict[str, Any]],
) -> dict[str, list[dict[str, str]]]:
    grouped = {
        "qc": [],
        "realism": [],
        "statistical": [],
    }

    if qc:
        for msg in qc.get("qc_critical", []) or []:
            grouped["qc"].append({"level": "critical", "message": str(msg)})
        for msg in qc.get("qc_warnings", []) or []:
            grouped["qc"].append({"level": "warning", "message": str(msg)})

    if realism:
        for msg in realism.get("critical", []) or []:
            grouped["realism"].append({"level": "critical", "message": str(msg)})
        for msg in realism.get("warnings", []) or []:
            grouped["realism"].append({"level": "warning", "message": str(msg)})

    for msg in summary.get("warnings", []) or []:
        grouped["statistical"].append({"level": "warning", "message": str(msg)})
    for msg in summary.get("data_issues", []) or []:
        grouped["statistical"].append({"level": "warning", "message": str(msg)})

    return grouped


def build_report(
    job_dir: Path,
    summary: AnalysisSummary,
    llm: Optional[LLMInterpretation],
    formula: Optional[str] = None,
    contrast: Optional[list[str]] = None,
) -> None:
    """Render HTML report and write to job_dir/report.html."""
    env = Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        autoescape=select_autoescape(["html"]),
    )

    template = env.get_template("report.html.j2")

    summary_data = summary.model_dump()

    plots_dir = job_dir / "plots"
    plots = {
        name: _img_to_base64(plots_dir / f"{name}.png")
        for name in [
            "pca",
            "volcano",
            "ma",
            "heatmap",
            "sample_distance_heatmap",
            "sample_correlation_heatmap",
            "library_size",
            "count_distribution",
            "zero_fraction",
        ]
    }

    qc_report: Optional[dict[str, Any]] = None
    qc_path = job_dir / "results" / "qc_report.json"
    if qc_path.exists():
        try:
            qc_report = json.loads(qc_path.read_text(encoding="utf-8"))
        except Exception:
            qc_report = None

    realism = summary_data.get("realism_validation") or {}
    top_genes_table = _load_top_genes_table(job_dir, n=20)

    analysis_methods = {
        "design_formula": formula or "not provided",
        "contrast": ", ".join(contrast) if contrast else summary_data.get("contrast", "not provided"),
        "normalization_method": "DESeq2 median-of-ratios size-factor normalization with VST for visualization",
        "filtering_criteria": "Genes with total count >= 10 retained before differential testing",
        "statistical_test": "Negative binomial GLM in DESeq2 with Wald test",
        "multiple_testing_correction": "Benjamini-Hochberg FDR (adjusted p-value)",
        "deg_threshold": "padj < 0.05 and |log2FoldChange| > 1",
    }

    groups = summary_data.get("groups", []) or []
    methods_paragraph = (
        f"RNA-seq count data were analyzed using design formula '{analysis_methods['design_formula']}' "
        f"and contrast '{analysis_methods['contrast']}'. Counts were normalized with {analysis_methods['normalization_method']}. "
        f"Prior to testing, low-information genes were filtered using the criterion: {analysis_methods['filtering_criteria']}. "
        f"Differential expression was evaluated using {analysis_methods['statistical_test']}, and multiplicity was controlled via "
        f"{analysis_methods['multiple_testing_correction']}. Reporting thresholds for DEG summaries were {analysis_methods['deg_threshold']}."
    )

    results_paragraph = (
        f"The analysis included {summary_data.get('n_samples', 0)} samples across groups: {', '.join(groups) if groups else 'not provided'}. "
        f"Differential testing identified {summary_data.get('deg_up', 0)} upregulated and {summary_data.get('deg_down', 0)} downregulated genes. "
        f"PCA separation was classified as '{summary_data.get('pca_separation', 'unknown')}'. "
        f"Overall data quality risk was assessed as '{_overall_data_quality(qc_report)}', and realism suspicion was "
        f"'{realism.get('overall_suspicion', 'unknown')}'."
    )

    figure_legends = [
        {
            "title": "PCA Plot",
            "text": "Principal component projection of VST-transformed counts. Axis labels report variance explained for PC1 and PC2.",
        },
        {
            "title": "Sample Distance Heatmap",
            "text": "Euclidean distance matrix between samples on VST-transformed counts, used for outlier structure assessment.",
        },
        {
            "title": "Sample Correlation Heatmap",
            "text": "Pearson correlation matrix across samples used to detect low-consistency profiles and potential outliers.",
        },
        {
            "title": "Volcano Plot",
            "text": "Gene-wise effect size (log2 fold-change) against significance, highlighting statistically relevant differentially expressed genes.",
        },
        {
            "title": "MA Plot",
            "text": "Mean expression versus log2 fold-change, used to evaluate expression-dependent effect size behavior.",
        },
        {
            "title": "Library Size / Count Distribution / Zero Fraction",
            "text": "Sample-level sequencing depth and count-shape diagnostics used by strict QC threshold rules.",
        },
    ]

    warning_groups = _grouped_warnings(summary_data, qc_report, realism)

    html = template.render(
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        summary=summary_data,
        qc=qc_report,
        realism=realism,
        top_genes_table=top_genes_table,
        warning_groups=warning_groups,
        analysis_methods=analysis_methods,
        methods_paragraph=methods_paragraph,
        results_paragraph=results_paragraph,
        figure_legends=figure_legends,
        overall_data_quality=_overall_data_quality(qc_report),
        overall_realism=realism.get("overall_suspicion", "unknown"),
        llm=llm.model_dump() if llm else None,
        plots=plots,
        job_id=job_dir.name,
    )

    out = job_dir / "report.html"
    out.write_text(html, encoding="utf-8")
    logger.info("Report written to %s", out)
