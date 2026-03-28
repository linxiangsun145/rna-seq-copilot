"""
Strict data realism validation for RNA-seq outputs.

This module is intentionally rule-based with explicit thresholds and
warning/critical severity levels. It never blocks analysis.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from models.schemas import RealismMetrics, RealismValidation

CANONICAL_GENES = {
    "TP53", "MYC", "EGFR", "VEGFA", "IL6", "CXCL8", "GAPDH", "ACTB", "AKT1", "PTEN"
}
HOUSEKEEPING_GENES = {"GAPDH", "ACTB", "RPL13A", "B2M", "HPRT1", "TUBB"}
MARKER_GENES = CANONICAL_GENES | HOUSEKEEPING_GENES


def _safe_fraction(numerator: int | float, denominator: int | float) -> float:
    if denominator is None or denominator == 0:
        return 0.0
    return float(numerator) / float(denominator)


def _empty_metrics() -> RealismMetrics:
    return RealismMetrics(
        canonical_genes_in_top20=0,
        canonical_fraction_top20=0.0,
        housekeeping_genes_in_top20=0,
        total_deg=0,
        fraction_p_lt_1e6=0.0,
        fraction_p_lt_1e3=0.0,
        fraction_p_gt_0_9=0.0,
        fraction_deg_abs_log2fc_gt_5=0.0,
        fraction_deg_abs_log2fc_gt_10=0.0,
    )


def _load_deg_table(job_dir: Path) -> pd.DataFrame:
    deg_path = job_dir / "results" / "deg_results.csv"
    if not deg_path.exists():
        return pd.DataFrame(columns=["gene_id", "log2FoldChange", "pvalue", "padj"])
    df = pd.read_csv(deg_path)
    expected = {"gene_id", "log2FoldChange", "pvalue", "padj"}
    for c in expected:
        if c not in df.columns:
            df[c] = np.nan
    return df


def _extract_qc(job_dir: Path) -> dict[str, Any] | None:
    qc_path = job_dir / "results" / "qc_report.json"
    if not qc_path.exists():
        return None
    try:
        import json

        return json.loads(qc_path.read_text(encoding="utf-8"))
    except Exception:
        return None


def validate_realism(job_dir: Path, summary_dict: dict[str, Any]) -> RealismValidation:
    warnings: list[str] = []
    critical: list[str] = []
    realism_flags: list[str] = []
    suspicious_patterns: list[str] = []

    df = _load_deg_table(job_dir)
    qc = _extract_qc(job_dir)

    if df.empty:
        metrics = _empty_metrics()
        critical.append("DEG table is missing or empty; realism checks cannot be computed from differential expression output.")
        return RealismValidation(
            realism_flags=realism_flags,
            suspicious_patterns=suspicious_patterns,
            warnings=warnings,
            critical=critical,
            metrics=metrics,
            overall_suspicion="high",
        )

    df = df.copy()
    df["gene_upper"] = df["gene_id"].astype(str).str.upper()
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")

    df_sorted = df.sort_values("padj", na_position="last")
    top20 = df_sorted.head(20)
    top5 = df_sorted.head(5)

    deg = df[df["padj"].notna() & (df["padj"] < 0.05)].copy()
    total_deg = int(deg.shape[0])

    # 1) Canonical overrepresentation
    canonical_in_top20 = int(top20["gene_upper"].isin(CANONICAL_GENES).sum())
    canonical_fraction_top20 = _safe_fraction(canonical_in_top20, 20)
    if canonical_in_top20 >= 7:
        msg = (
            f"Canonical gene overrepresentation is critical: {canonical_in_top20}/20 canonical genes in top 20 "
            f"(threshold >= 7)."
        )
        critical.append(msg)
        realism_flags.append("canonical_top20_critical")
        suspicious_patterns.append(msg)
    elif canonical_in_top20 >= 4:
        msg = (
            f"Canonical gene overrepresentation warning: {canonical_in_top20}/20 canonical genes in top 20 "
            f"(threshold >= 4)."
        )
        warnings.append(msg)
        realism_flags.append("canonical_top20_warning")
        suspicious_patterns.append(msg)

    # 2) Housekeeping genes in top20
    hk_top20 = top20[top20["gene_upper"].isin(HOUSEKEEPING_GENES)]
    hk_count = int(hk_top20.shape[0])
    if hk_count >= 2:
        msg = f"Housekeeping gene signal is critical: {hk_count} housekeeping genes in top 20 (threshold >= 2)."
        critical.append(msg)
        realism_flags.append("housekeeping_top20_critical")
        suspicious_patterns.append(msg)
    elif hk_count >= 1:
        msg = f"Housekeeping gene warning: {hk_count} housekeeping gene in top 20 (threshold >= 1)."
        warnings.append(msg)
        realism_flags.append("housekeeping_top20_warning")
        suspicious_patterns.append(msg)

    hk_strong = deg[
        deg["gene_upper"].isin(HOUSEKEEPING_GENES)
        & deg["log2FoldChange"].notna()
        & (deg["log2FoldChange"].abs() > 2)
        & (deg["padj"] < 0.01)
    ]
    if not hk_strong.empty:
        genes = ", ".join(sorted(set(hk_strong["gene_upper"].tolist())))
        msg = (
            "Housekeeping strong differential expression is critical: "
            f"abs(log2FoldChange) > 2 and padj < 0.01 detected for {genes}."
        )
        critical.append(msg)
        realism_flags.append("housekeeping_strong_effect_critical")
        suspicious_patterns.append(msg)

    # 3) DEG count sanity
    if total_deg == 0:
        msg = "DEG count sanity critical: total_deg = 0 with padj < 0.05 threshold."
        critical.append(msg)
        realism_flags.append("deg_count_zero_critical")
        suspicious_patterns.append(msg)
    elif total_deg > 10000:
        msg = f"DEG count sanity critical: total_deg = {total_deg} (> 10000)."
        critical.append(msg)
        realism_flags.append("deg_count_too_high_critical")
        suspicious_patterns.append(msg)
    elif total_deg < 10:
        msg = f"DEG count sanity warning: total_deg = {total_deg} (< 10)."
        warnings.append(msg)
        realism_flags.append("deg_count_low_warning")
        suspicious_patterns.append(msg)
    elif total_deg > 5000:
        msg = f"DEG count sanity warning: total_deg = {total_deg} (> 5000)."
        warnings.append(msg)
        realism_flags.append("deg_count_high_warning")
        suspicious_patterns.append(msg)

    # 4) P-value distribution anomalies
    valid_p = df["pvalue"].dropna()
    valid_p = valid_p[(valid_p >= 0) & (valid_p <= 1)]

    fraction_p_lt_1e6 = _safe_fraction((valid_p < 1e-6).sum(), len(valid_p))
    fraction_p_lt_1e3 = _safe_fraction((valid_p < 1e-3).sum(), len(valid_p))
    fraction_p_gt_0_9 = _safe_fraction((valid_p > 0.9).sum(), len(valid_p))

    if fraction_p_lt_1e6 > 0.5:
        msg = (
            f"P-value anomaly critical: {fraction_p_lt_1e6:.3f} of p-values are < 1e-6 "
            "(threshold > 0.50)."
        )
        critical.append(msg)
        realism_flags.append("pvalue_tiny_fraction_critical")
        suspicious_patterns.append(msg)
    elif fraction_p_lt_1e6 > 0.3:
        msg = (
            f"P-value anomaly warning: {fraction_p_lt_1e6:.3f} of p-values are < 1e-6 "
            "(threshold > 0.30)."
        )
        warnings.append(msg)
        realism_flags.append("pvalue_tiny_fraction_warning")
        suspicious_patterns.append(msg)

    # Uniformity check in non-significant region [0.1, 1.0]
    non_sig = valid_p[valid_p >= 0.1]
    if len(non_sig) >= 200:
        bins = np.linspace(0.1, 1.0, 10)
        hist, _ = np.histogram(non_sig, bins=bins)
        mean_bin = float(np.mean(hist)) if len(hist) > 0 else 0.0
        cv = float(np.std(hist) / mean_bin) if mean_bin > 0 else 0.0
        if cv < 0.2:
            msg = (
                "P-value anomaly warning: non-significant p-values are excessively uniform "
                f"across bins (coefficient of variation={cv:.3f} < 0.20)."
            )
            warnings.append(msg)
            realism_flags.append("pvalue_uniform_non_sig_warning")
            suspicious_patterns.append(msg)

    # Extreme pile-up at both ends
    if fraction_p_lt_1e3 > 0.2 and fraction_p_gt_0_9 > 0.2:
        msg = (
            "P-value anomaly warning: simultaneous pile-up near both tails detected "
            f"(fraction_p_lt_1e3={fraction_p_lt_1e3:.3f}, fraction_p_gt_0_9={fraction_p_gt_0_9:.3f})."
        )
        warnings.append(msg)
        realism_flags.append("pvalue_bimodal_tail_pileup_warning")
        suspicious_patterns.append(msg)

    # 5) Effect size plausibility (DEG only)
    if total_deg > 0:
        frac_abs_gt5 = _safe_fraction((deg["log2FoldChange"].abs() > 5).sum(), total_deg)
        frac_abs_gt10 = _safe_fraction((deg["log2FoldChange"].abs() > 10).sum(), total_deg)
    else:
        frac_abs_gt5 = 0.0
        frac_abs_gt10 = 0.0

    if frac_abs_gt10 > 0.1:
        msg = (
            f"Effect-size plausibility critical: {frac_abs_gt10:.3f} of DEGs have abs(log2FoldChange) > 10 "
            "(threshold > 0.10)."
        )
        critical.append(msg)
        realism_flags.append("effect_size_extreme_critical")
        suspicious_patterns.append(msg)

    if frac_abs_gt5 > 0.3:
        msg = (
            f"Effect-size plausibility warning: {frac_abs_gt5:.3f} of DEGs have abs(log2FoldChange) > 5 "
            "(threshold > 0.30)."
        )
        warnings.append(msg)
        realism_flags.append("effect_size_large_warning")
        suspicious_patterns.append(msg)

    # 6) Top-gene dominance
    highly_sig = df_sorted[df_sorted["padj"].notna() & (df_sorted["padj"] < 1e-6)].copy()
    top5_share = 0.0
    if not highly_sig.empty:
        clipped = highly_sig["padj"].clip(lower=1e-300)
        highly_sig["sig_weight"] = -np.log10(clipped)
        denom = float(highly_sig["sig_weight"].sum())
        numer = float(highly_sig.head(5)["sig_weight"].sum())
        top5_share = _safe_fraction(numer, denom)
        if len(highly_sig) >= 10 and top5_share >= 0.6:
            msg = (
                "Top-gene dominance warning: top 5 genes account for "
                f"{top5_share:.3f} of total -log10(padj) among highly significant genes (threshold >= 0.60)."
            )
            warnings.append(msg)
            realism_flags.append("top5_dominance_warning")
            suspicious_patterns.append(msg)

    marker_in_top5 = int(top5["gene_upper"].isin(MARKER_GENES).sum())
    if marker_in_top5 >= 3:
        msg = (
            "Top-gene dominance warning: top 5 genes are dominated by known marker/canonical genes "
            f"({marker_in_top5}/5 in marker set, threshold >= 3)."
        )
        warnings.append(msg)
        realism_flags.append("top5_marker_dominance_warning")
        suspicious_patterns.append(msg)

    # 7) Optional QC consistency cross-check
    if qc:
        qc_outliers = qc.get("outliers", []) or []
        lib_flags = qc.get("library_size_flags", []) or []
        has_low_lib = len(lib_flags) > 0
        severe_deg = total_deg >= 1000 or frac_abs_gt5 > 0.3

        if severe_deg and (len(qc_outliers) > 0 or has_low_lib):
            msg = (
                "QC consistency warning: severe DEG signal co-occurs with QC risk "
                f"(outliers={len(qc_outliers)}, low_library_flags={len(lib_flags)})."
            )
            warnings.append(msg)
            realism_flags.append("qc_consistency_deg_vs_quality_warning")
            suspicious_patterns.append(msg)

        reps = qc.get("replicates_per_group", {}) or {}
        has_low_rep_n2 = any(int(v) == 2 for v in reps.values())
        pca_var = qc.get("pca_variance", {}) or {}
        pc1 = float(pca_var.get("PC1", 0.0) or 0.0)
        strong_sep = summary_dict.get("pca_separation") == "clear" or pc1 >= 0.3

        if strong_sep and has_low_rep_n2:
            msg = (
                "QC consistency warning: strong separation signal is reported while at least one group has only n=2 replicates."
            )
            warnings.append(msg)
            realism_flags.append("qc_consistency_separation_vs_rep_warning")
            suspicious_patterns.append(msg)

        batch_flags = qc.get("batch_flags", []) or []
        confounded = any((f.get("rule") == "batch_condition_confounded") for f in batch_flags if isinstance(f, dict))
        if confounded and total_deg > 0:
            msg = (
                "QC consistency critical: condition is confounded with batch while DEG interpretation is still present."
            )
            critical.append(msg)
            realism_flags.append("qc_batch_confounded_critical")
            suspicious_patterns.append(msg)

    metrics = RealismMetrics(
        canonical_genes_in_top20=canonical_in_top20,
        canonical_fraction_top20=round(canonical_fraction_top20, 4),
        housekeeping_genes_in_top20=hk_count,
        total_deg=total_deg,
        fraction_p_lt_1e6=round(float(fraction_p_lt_1e6), 4),
        fraction_p_lt_1e3=round(float(fraction_p_lt_1e3), 4),
        fraction_p_gt_0_9=round(float(fraction_p_gt_0_9), 4),
        fraction_deg_abs_log2fc_gt_5=round(float(frac_abs_gt5), 4),
        fraction_deg_abs_log2fc_gt_10=round(float(frac_abs_gt10), 4),
    )

    # Overall suspicion scoring
    if critical or len(warnings) >= 4:
        overall = "high"
    elif len(warnings) >= 2:
        overall = "moderate"
    else:
        overall = "low"

    warning_items: list[dict[str, Any]] = []
    for idx, msg in enumerate(critical):
        code = (
            realism_flags[idx]
            if idx < len(realism_flags)
            else f"realism_critical_{idx+1}"
        )
        warning_items.append(
            {
                "type": "realism",
                "severity": "critical",
                "code": str(code),
                "message": str(msg),
                "sample": None,
                "metric": "realism",
            }
        )

    warn_offset = len(critical)
    for idx, msg in enumerate(warnings):
        flag_index = warn_offset + idx
        code = (
            realism_flags[flag_index]
            if flag_index < len(realism_flags)
            else f"realism_warning_{idx+1}"
        )
        warning_items.append(
            {
                "type": "realism",
                "severity": "warning",
                "code": str(code),
                "message": str(msg),
                "sample": None,
                "metric": "realism",
            }
        )

    return RealismValidation(
        realism_flags=realism_flags,
        suspicious_patterns=suspicious_patterns,
        warnings=warnings,
        critical=critical,
        warning_items=warning_items,
        metrics=metrics,
        overall_suspicion=overall,
    )
