"""Deterministic ncRNA-aware analysis helpers.

This module adds annotation, summary metrics, diagnostics, and candidate ranking
without changing DEG statistical inference.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd


SUPPORTED_BIOTYPES = {
    "protein_coding",
    "lncrna",
    "mirna",
    "circrna",
    "other_ncrna",
    "unknown",
}

NCRNA_SET = {"lncrna", "mirna", "circrna", "other_ncrna"}


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in (None, "", "NA", "NaN"):
            return default
        return float(value)
    except Exception:
        return default


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        if value in (None, "", "NA", "NaN"):
            return default
        return int(value)
    except Exception:
        return default


def _canonical_biotype(raw: str) -> str:
    t = str(raw or "").strip().lower().replace("-", "_")
    aliases = {
        "protein_coding": "protein_coding",
        "lncrna": "lncrna",
        "linc_rna": "lncrna",
        "lincrna": "lncrna",
        "mirna": "mirna",
        "mi_rna": "mirna",
        "micro_rna": "mirna",
        "circrna": "circrna",
        "circ_rna": "circrna",
        "other_ncrna": "other_ncrna",
        "ncrna": "other_ncrna",
        "snorna": "other_ncrna",
        "snrna": "other_ncrna",
        "rrna": "other_ncrna",
        "trna": "other_ncrna",
        "unknown": "unknown",
    }
    mapped = aliases.get(t, t)
    return mapped if mapped in SUPPORTED_BIOTYPES else "unknown"


def _load_biotype_table() -> dict[str, str]:
    """Load optional local biotype mapping table.

    Expected TSV/CSV columns (first match wins):
    - feature, gene, gene_id, transcript_id, ensembl_id
    - biotype, gene_biotype
    """
    mapping: dict[str, str] = {}
    base = Path(__file__).parent.parent / "data"
    candidates = [base / "rna_biotype_map.tsv", base / "rna_biotype_map.csv"]

    for path in candidates:
        if not path.exists():
            continue
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        try:
            df = pd.read_csv(path, sep=sep)
        except Exception:
            continue

        if df.empty:
            continue

        cols = {str(c).lower(): c for c in df.columns}
        key_col = None
        for candidate in ["feature", "gene", "gene_id", "transcript_id", "ensembl_id"]:
            if candidate in cols:
                key_col = cols[candidate]
                break
        if key_col is None:
            continue

        biotype_col = None
        for candidate in ["biotype", "gene_biotype"]:
            if candidate in cols:
                biotype_col = cols[candidate]
                break
        if biotype_col is None:
            continue

        for _, row in df.iterrows():
            key = str(row.get(key_col, "")).strip()
            if not key:
                continue
            mapping[key.upper()] = _canonical_biotype(str(row.get(biotype_col, "unknown")))

        if mapping:
            return mapping

    return mapping


def _rule_based_biotype(feature: str) -> str:
    """Conservative fallback based on common naming conventions only."""
    f = str(feature or "").strip()
    if not f:
        return "unknown"

    u = f.upper()

    if u.startswith("HSA_CIRC") or u.startswith("CIRC"):
        return "circrna"

    if u.startswith("MIR") or u.startswith("LET-") or u.startswith("MIRLET"):
        return "mirna"

    if u.startswith("LINC"):
        return "lncrna"

    if u in {"MALAT1", "NEAT1", "XIST", "H19", "MEG3"}:
        return "lncrna"

    if u.startswith("SNORD") or u.startswith("SNORA") or u.startswith("RNU") or u.startswith("RN7SK") or u.startswith("SCARNA"):
        return "other_ncrna"

    return "unknown"


def annotate_rna_biotype(features: list[str]) -> list[dict[str, str]]:
    """Annotate each feature with RNA biotype.

    The function first uses a local map, then applies conservative naming rules,
    otherwise returns unknown.
    """
    mapping = _load_biotype_table()
    out: list[dict[str, str]] = []

    for feat in features or []:
        f = str(feat or "").strip()
        if not f:
            continue
        biotype = mapping.get(f.upper())
        if not biotype:
            biotype = _rule_based_biotype(f)
        out.append({"gene": f, "biotype": _canonical_biotype(biotype)})

    return out


def attach_biotype_to_deg_results(deg_table: pd.DataFrame) -> pd.DataFrame:
    """Attach biotype column to DEG table (gene_id-driven)."""
    if deg_table is None or deg_table.empty:
        out = deg_table.copy() if isinstance(deg_table, pd.DataFrame) else pd.DataFrame()
        if isinstance(out, pd.DataFrame) and "biotype" not in out.columns:
            out["biotype"] = []
        return out

    out = deg_table.copy()
    if "gene_id" not in out.columns:
        out["biotype"] = "unknown"
        return out

    feature_list = [str(x).strip() for x in out["gene_id"].tolist()]
    ann = annotate_rna_biotype(feature_list)
    biotype_map = {a["gene"]: a["biotype"] for a in ann}
    out["biotype"] = [biotype_map.get(str(g).strip(), "unknown") for g in out["gene_id"].tolist()]
    return out


def compute_ncrna_summary_metrics(deg_table: pd.DataFrame) -> dict[str, Any]:
    """Compute ncRNA summary counts and fractions among DEGs."""
    if deg_table is None or deg_table.empty:
        total = 0
        return {
            "n_total_deg": total,
            "n_ncrna_deg": 0,
            "n_lncRNA_deg": 0,
            "n_miRNA_deg": 0,
            "n_protein_coding_deg": 0,
            "n_unknown_biotype_deg": 0,
            "n_circRNA_deg": 0,
            "n_other_ncRNA_deg": 0,
            "ncrna_fraction_deg": 0.0,
            "lncRNA_fraction_deg": 0.0,
            "miRNA_fraction_deg": 0.0,
        }

    bio = [str(x).strip().lower() for x in deg_table.get("biotype", pd.Series(["unknown"] * len(deg_table))).tolist()]
    total = int(len(bio))

    counts = {
        "protein_coding": sum(1 for b in bio if b == "protein_coding"),
        "lncrna": sum(1 for b in bio if b == "lncrna"),
        "mirna": sum(1 for b in bio if b == "mirna"),
        "circrna": sum(1 for b in bio if b == "circrna"),
        "other_ncrna": sum(1 for b in bio if b == "other_ncrna"),
        "unknown": sum(1 for b in bio if b not in SUPPORTED_BIOTYPES or b == "unknown"),
    }

    n_ncrna = counts["lncrna"] + counts["mirna"] + counts["circrna"] + counts["other_ncrna"]

    def frac(v: int) -> float:
        return float(v) / float(total) if total > 0 else 0.0

    return {
        "n_total_deg": total,
        "n_ncrna_deg": n_ncrna,
        "n_lncRNA_deg": counts["lncrna"],
        "n_miRNA_deg": counts["mirna"],
        "n_protein_coding_deg": counts["protein_coding"],
        "n_unknown_biotype_deg": counts["unknown"],
        "n_circRNA_deg": counts["circrna"],
        "n_other_ncRNA_deg": counts["other_ncrna"],
        "ncrna_fraction_deg": frac(n_ncrna),
        "lncRNA_fraction_deg": frac(counts["lncrna"]),
        "miRNA_fraction_deg": frac(counts["mirna"]),
    }


def compute_ncrna_qc_metrics(count_matrix: pd.DataFrame, biotype_map: dict[str, str]) -> dict[str, Any]:
    """Compute ncRNA-specific diagnostics from raw count matrix."""
    if count_matrix is None or count_matrix.empty:
        return {
            "ncrna_zero_fraction": 0.0,
            "ncrna_low_count_fraction": 0.0,
            "lncrna_low_count_fraction": 0.0,
            "miRNA_detected_count": 0,
            "lncRNA_detected_count": 0,
            "n_ncrna_features": 0,
        }

    if "gene_id" in count_matrix.columns:
        idx = count_matrix["gene_id"].astype(str).tolist()
        values = count_matrix.drop(columns=["gene_id"]).copy()
    else:
        idx = [str(i) for i in count_matrix.index.tolist()]
        values = count_matrix.copy()

    values = values.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    feature_bio = [str(biotype_map.get(str(g), "unknown")).lower() for g in idx]

    ncrna_rows = [i for i, b in enumerate(feature_bio) if b in NCRNA_SET]
    mirna_rows = [i for i, b in enumerate(feature_bio) if b == "mirna"]
    lncrna_rows = [i for i, b in enumerate(feature_bio) if b == "lncrna"]

    if not ncrna_rows:
        return {
            "ncrna_zero_fraction": 0.0,
            "ncrna_low_count_fraction": 0.0,
            "lncrna_low_count_fraction": 0.0,
            "miRNA_detected_count": 0,
            "lncRNA_detected_count": 0,
            "n_ncrna_features": 0,
        }

    ncrna_mat = values.iloc[ncrna_rows, :]
    ncrna_zero_fraction = float((ncrna_mat <= 0).sum().sum()) / float(ncrna_mat.size) if ncrna_mat.size > 0 else 0.0
    ncrna_low_count_fraction = float((ncrna_mat <= 5).sum().sum()) / float(ncrna_mat.size) if ncrna_mat.size > 0 else 0.0

    lncrna_low_count_fraction = 0.0
    if lncrna_rows:
        ln_mat = values.iloc[lncrna_rows, :]
        lncrna_low_count_fraction = float((ln_mat <= 5).sum().sum()) / float(ln_mat.size) if ln_mat.size > 0 else 0.0

    mirna_detected = 0
    if mirna_rows:
        mi_mat = values.iloc[mirna_rows, :]
        mirna_detected = int((mi_mat.mean(axis=1) > 1.0).sum())

    lncrna_detected = 0
    if lncrna_rows:
        ln_mat = values.iloc[lncrna_rows, :]
        lncrna_detected = int((ln_mat.mean(axis=1) > 1.0).sum())

    return {
        "ncrna_zero_fraction": round(ncrna_zero_fraction, 4),
        "ncrna_low_count_fraction": round(ncrna_low_count_fraction, 4),
        "lncrna_low_count_fraction": round(lncrna_low_count_fraction, 4),
        "miRNA_detected_count": int(mirna_detected),
        "lncRNA_detected_count": int(lncrna_detected),
        "n_ncrna_features": int(len(ncrna_rows)),
    }


def generate_ncrna_warnings(metrics: dict[str, Any]) -> list[str]:
    """Generate metric-backed ncRNA warnings (non-blocking)."""
    out: list[str] = []

    zero_frac = _safe_float(metrics.get("ncrna_zero_fraction"), 0.0)
    low_frac = _safe_float(metrics.get("ncrna_low_count_fraction"), 0.0)
    lnc_low_frac = _safe_float(metrics.get("lncrna_low_count_fraction"), 0.0)
    mirna_detected = _safe_int(metrics.get("miRNA_detected_count"), 0)

    if zero_frac > 0.8:
        out.append("High sparsity among ncRNA features may limit interpretability.")
    if mirna_detected < 10:
        out.append("miRNA signal is weak across most samples.")
    if lnc_low_frac > 0.7:
        out.append("lncRNA findings are dominated by low-abundance transcripts.")
    elif low_frac > 0.75:
        out.append("Low-count burden is high in the ncRNA subset.")

    dedup: list[str] = []
    seen: set[str] = set()
    for x in out:
        t = str(x).strip()
        if t and t not in seen:
            seen.add(t)
            dedup.append(t)
    return dedup


def rank_ncrna_candidates(deg_table: pd.DataFrame, qc_metrics: dict[str, Any]) -> list[dict[str, Any]]:
    """Rank ncRNA candidates with transparent deterministic scoring."""
    if deg_table is None or deg_table.empty:
        return []

    rows: list[dict[str, Any]] = []
    zero_frac = _safe_float(qc_metrics.get("ncrna_zero_fraction"), 0.0)
    low_frac = _safe_float(qc_metrics.get("ncrna_low_count_fraction"), 0.0)

    for _, row in deg_table.iterrows():
        biotype = str(row.get("biotype", "unknown")).strip().lower()
        if biotype not in NCRNA_SET:
            continue

        gene = str(row.get("gene_id", "")).strip()
        if not gene:
            continue

        padj = _safe_float(row.get("padj"), 1.0)
        if padj <= 0:
            padj = 1e-300
        log2fc = abs(_safe_float(row.get("log2FoldChange"), 0.0))
        base_mean = _safe_float(row.get("baseMean"), 10.0)

        sig_score = min((-1.0 * math.log10(padj)) / 10.0, 1.0)
        effect_score = min(log2fc / 4.0, 1.0)
        robust_score = min((math.log10(base_mean + 1.0)) / 3.0, 1.0)

        score = 0.5 * sig_score + 0.3 * effect_score + 0.2 * robust_score
        if zero_frac > 0.8:
            score *= 0.85
        if low_frac > 0.75:
            score *= 0.90

        rows.append(
            {
                "gene": gene,
                "biotype": biotype,
                "priority_score": round(float(score), 3),
                "padj": padj,
                "abs_log2fc": round(log2fc, 3),
            }
        )

    rows.sort(key=lambda x: (-x["priority_score"], x["padj"], x["gene"]))
    return rows[:20]


def render_ncrna_section_html(ncrna_assessment: dict[str, Any]) -> dict[str, Any]:
    """Build deterministic ncRNA section payload for HTML rendering."""
    ass = ncrna_assessment or {}
    metrics = ass.get("summary_metrics") if isinstance(ass.get("summary_metrics"), dict) else {}
    qc_metrics = ass.get("qc_metrics") if isinstance(ass.get("qc_metrics"), dict) else {}
    warnings = [str(x) for x in (ass.get("warnings", []) or []) if str(x).strip()]
    candidates = [x for x in (ass.get("top_candidates", []) or []) if isinstance(x, dict)]

    narrative = []
    n_total = _safe_int(metrics.get("n_total_deg"), 0)
    n_ncrna = _safe_int(metrics.get("n_ncrna_deg"), 0)
    lnc_frac = _safe_float(metrics.get("lncRNA_fraction_deg"), 0.0)

    narrative.append(f"Among {n_total} DEGs, {n_ncrna} were annotated as ncRNAs.")
    narrative.append(f"lncRNAs represented {lnc_frac * 100:.1f}% of the DEG set.")
    if warnings:
        narrative.append("Detected ncRNA findings should be interpreted cautiously due to ncRNA-specific QC warnings.")

    return {
        "summary_metrics": metrics,
        "qc_metrics": qc_metrics,
        "warnings": warnings,
        "top_candidates": candidates,
        "narrative": narrative,
    }


def analyze_ncrna_assessment(job_dir: Path) -> dict[str, Any]:
    """End-to-end ncRNA analysis wrapper for pipeline integration."""
    results_dir = Path(job_dir) / "results"
    deg_path = results_dir / "deg_results.csv"
    counts_path = Path(job_dir) / "counts.csv"

    if not deg_path.exists():
        return {
            "summary_metrics": compute_ncrna_summary_metrics(pd.DataFrame()),
            "qc_metrics": compute_ncrna_qc_metrics(pd.DataFrame(), {}),
            "warnings": [],
            "top_candidates": [],
        }

    deg_df = pd.read_csv(deg_path)
    deg_df = attach_biotype_to_deg_results(deg_df)
    deg_df.to_csv(deg_path, index=False)

    biotype_map = {
        str(r.get("gene_id", "")).strip(): str(r.get("biotype", "unknown")).strip().lower()
        for _, r in deg_df.iterrows()
    }

    if counts_path.exists():
        counts_df = pd.read_csv(counts_path)
    else:
        counts_df = pd.DataFrame()

    summary_metrics = compute_ncrna_summary_metrics(deg_df)
    qc_metrics = compute_ncrna_qc_metrics(counts_df, biotype_map)
    warnings = generate_ncrna_warnings(qc_metrics)
    top_candidates = rank_ncrna_candidates(deg_df, qc_metrics)

    return {
        "summary_metrics": summary_metrics,
        "qc_metrics": qc_metrics,
        "warnings": warnings,
        "top_candidates": top_candidates,
    }
