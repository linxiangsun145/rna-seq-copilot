"""
Input validation service.
Validates count matrix and metadata CSV/TSV before running DESeq2.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List

import pandas as pd

from models.schemas import ValidationIssue, ValidationReport

logger = logging.getLogger(__name__)

# Heuristic thresholds
TPM_MEAN_THRESHOLD = 50.0  # If mean of all counts > this, suspect TPM/FPKM
LOW_COUNT_FRACTION = 0.8   # Warn if >80% of genes have zero counts


def _read_tabular(path: Path) -> pd.DataFrame:
    """Read CSV or TSV, auto-detect separator."""
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, index_col=0)


def validate_inputs(counts_path: Path, meta_path: Path) -> ValidationReport:
    issues: List[ValidationIssue] = []

    # ── Load files ────────────────────────────────────────────────────────────
    try:
        counts = _read_tabular(counts_path)
    except Exception as exc:
        return ValidationReport(
            valid=False,
            n_genes=0,
            n_samples=0,
            sample_names=[],
            groups={},
            issues=[ValidationIssue(level="error", field="counts_file", message=str(exc))],
        )

    try:
        meta = _read_tabular(meta_path)
    except Exception as exc:
        return ValidationReport(
            valid=False,
            n_genes=counts.shape[0],
            n_samples=0,
            sample_names=list(counts.columns),
            groups={},
            issues=[ValidationIssue(level="error", field="metadata_file", message=str(exc))],
        )

    n_genes, n_samples = counts.shape
    count_samples = list(counts.columns)
    meta_samples = list(meta.index)

    # ── Sample name matching ──────────────────────────────────────────────────
    missing_in_meta = set(count_samples) - set(meta_samples)
    missing_in_counts = set(meta_samples) - set(count_samples)
    if missing_in_meta:
        issues.append(ValidationIssue(
            level="error",
            field="sample_names",
            message=f"Samples in counts but not in metadata: {sorted(missing_in_meta)}",
        ))
    if missing_in_counts:
        issues.append(ValidationIssue(
            level="warning",
            field="sample_names",
            message=f"Samples in metadata but not in counts: {sorted(missing_in_counts)}",
        ))

    # Keep only matching samples
    common = [s for s in count_samples if s in set(meta_samples)]
    if not common:
        issues.append(ValidationIssue(
            level="error",
            field="sample_names",
            message="No common samples found between counts and metadata.",
        ))
        return ValidationReport(
            valid=False,
            n_genes=n_genes,
            n_samples=0,
            sample_names=[],
            groups={},
            issues=issues,
        )

    counts = counts[common]
    meta = meta.loc[common]

    # ── Duplicate sample check ────────────────────────────────────────────────
    dup = [s for s in common if common.count(s) > 1]
    if dup:
        issues.append(ValidationIssue(
            level="error",
            field="sample_names",
            message=f"Duplicate sample names: {list(set(dup))}",
        ))

    # ── Integer count check (vectorised, pandas 2.x compatible) ────────────
    numeric_counts = counts.apply(pd.to_numeric, errors="coerce")
    has_non_int = (numeric_counts % 1 != 0).any(axis=None)
    if has_non_int:
        issues.append(ValidationIssue(
            level="error",
            field="counts",
            message="Count matrix contains non-integer values. DESeq2 requires raw integer counts.",
        ))

    # ── Negative values ────────────────────────────────────────────────────────
    if (counts < 0).any().any():
        issues.append(ValidationIssue(
            level="error",
            field="counts",
            message="Count matrix contains negative values.",
        ))

    # ── TPM/FPKM heuristic ────────────────────────────────────────────────────
    col_means = counts.mean()
    if (col_means > TPM_MEAN_THRESHOLD).any():
        issues.append(ValidationIssue(
            level="warning",
            field="counts",
            message=(
                f"Some samples have very high mean counts (max={col_means.max():.1f}). "
                "If data is TPM/FPKM-normalised, DESeq2 results will be invalid."
            ),
        ))

    # ── Missing values ────────────────────────────────────────────────────────
    if counts.isnull().any().any():
        issues.append(ValidationIssue(
            level="error",
            field="counts",
            message="Count matrix contains missing values (NaN).",
        ))

    # ── Low count fraction ────────────────────────────────────────────────────
    zero_fraction = (counts == 0).mean(axis=1).mean()
    if zero_fraction > LOW_COUNT_FRACTION:
        issues.append(ValidationIssue(
            level="warning",
            field="counts",
            message=f"{zero_fraction:.0%} of count entries are zero — data appears very sparse.",
        ))

    # ── Group inference (first metadata column with ≤10 unique vals) ─────────
    groups: Dict[str, List[str]] = {}
    group_col = None
    for col in meta.columns:
        if meta[col].nunique() <= 10:
            group_col = col
            break

    if group_col:
        for grp in meta[group_col].unique():
            groups[str(grp)] = list(meta.index[meta[group_col] == grp])
    else:
        issues.append(ValidationIssue(
            level="warning",
            field="metadata",
            message="Could not identify a suitable grouping column (all columns have >10 unique values).",
        ))

    # ── Group balance ─────────────────────────────────────────────────────────
    if groups:
        sizes = {g: len(s) for g, s in groups.items()}
        min_size = min(sizes.values())
        if min_size < 2:
            issues.append(ValidationIssue(
                level="error",
                field="groups",
                message=f"Some groups have fewer than 2 replicates: {sizes}. DESeq2 requires ≥2 per group.",
            ))
        elif min_size < 3:
            issues.append(ValidationIssue(
                level="warning",
                field="groups",
                message=f"Some groups have only 2 replicates: {sizes}. Statistical power may be low.",
            ))

    has_error = any(i.level == "error" for i in issues)

    return ValidationReport(
        valid=not has_error,
        n_genes=n_genes,
        n_samples=len(common),
        sample_names=common,
        groups=groups,
        issues=issues,
    )
