#!/usr/bin/env python3
"""Build GSE201284 gene co-expression results for all requested metadata features.

For every valid subgroup, this script creates:

    OUTPUT/<feature>/<group>/genes/<gene>.json
    OUTPUT/<feature>/<group>/features/<gene>_feature.json

The first file contains the gene's expression values for the subgroup. The
second contains every other gene passing both |Pearson r| and p-value cutoffs.
All genes are processed; no starting-gene selection is required.

Requested metadata features:
  source_name_ch1, disease_state, preoperative_psa,
  patient_gleason_score, race, pathologic_stage, age,
  5_year_disease_status, sample_type

Age is grouped into decades. PSA is grouped into quantiles because its exact
values generally contain too few samples for stable within-group correlation.
Missing/unknown metadata and nonnumeric expression values are excluded.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


FEATURES = [
    "source_name_ch1",
    "disease_state",
    "preoperative_psa",
    "patient_gleason_score",
    "race",
    "pathologic_stage",
    "age",
    "5_year_disease_status",
    "sample_type",
]

MISSING_LABELS = {
    "", "na", "n/a", "nan", "none", "null", "unknown", "unk", "not available",
}


def safe_name(value: object) -> str:
    """Create a portable folder or filename while retaining readable labels."""
    text = str(value).strip()
    text = re.sub(r"[^A-Za-z0-9._+-]+", "_", text).strip("._")
    return text or "unnamed"


def normalize_missing(series: pd.Series) -> pd.Series:
    """Convert all configured missing labels to pandas NA."""
    result = series.astype("string").str.strip()
    missing = result.str.lower().isin(MISSING_LABELS)
    return result.mask(missing)


def make_age_groups(series: pd.Series) -> pd.Series:
    """Convert age to non-overlapping decade labels such as 60-69."""
    numeric = pd.to_numeric(normalize_missing(series), errors="coerce")
    lower = np.floor(numeric / 10.0) * 10
    labels = lower.map(lambda x: f"{int(x)}-{int(x + 9)}" if pd.notna(x) else pd.NA)
    return labels.astype("string")


def make_psa_groups(series: pd.Series, quantiles: int) -> pd.Series:
    """Create data-balanced PSA groups and preserve their numeric ranges."""
    numeric = pd.to_numeric(normalize_missing(series), errors="coerce")
    valid = numeric.dropna()
    output = pd.Series(pd.NA, index=series.index, dtype="string")
    if valid.empty:
        return output

    unique_count = valid.nunique()
    bins = min(quantiles, unique_count)
    if bins < 2:
        output.loc[valid.index] = f"PSA_{valid.iloc[0]:g}"
        return output

    # Keep identical PSA values in the same group; duplicate quantile boundaries
    # are dropped automatically if the distribution does not support all bins.
    codes = pd.qcut(valid, q=bins, labels=False, duplicates="drop")
    for code in sorted(codes.unique()):
        indices = codes.index[codes == code]
        low = numeric.loc[indices].min()
        high = numeric.loc[indices].max()
        output.loc[indices] = f"Q{int(code) + 1}_{low:g}-{high:g}"
    return output


def make_pathologic_stage_groups(series: pd.Series) -> pd.Series:
    """Collapse every T2/T3 substage into its parent stage category.

    This maps T2, T2a, T2b, T2c, T2/2a and the observed transposition Tb2
    to T2. Likewise, T3, T3a and T3b are mapped to T3.
    """
    values = normalize_missing(series).str.upper()
    compact = values.str.replace(r"[^A-Z0-9]", "", regex=True)
    output = values.copy()
    output.loc[compact.str.startswith("T2", na=False)] = "T2"
    output.loc[compact.eq("TB2")] = "T2"
    output.loc[compact.str.startswith("T3", na=False)] = "T3"
    return output


def make_gleason_groups(series: pd.Series) -> pd.Series:
    """Normalize scores and map them to standard prostate Grade Groups."""
    values = normalize_missing(series)
    base_scores = values.str.replace(
        r"(?i)\s+T\d+[A-Za-z]*\s*$", "", regex=True
    ).str.replace(" ", "", regex=False).str.strip()
    grade_group_map = {
        "3+3": "GG1",
        "3+4": "GG2",
        "4+3": "GG3",
        "4+4": "GG4",
        "3+5": "GG4",
        "5+3": "GG4",
        "4+5": "GG5",
        "5+4": "GG5",
        "5+5": "GG5",
    }
    return base_scores.map(grade_group_map).astype("string")


def build_groups(metadata: pd.DataFrame, psa_quantiles: int) -> dict[str, pd.Series]:
    missing = [feature for feature in FEATURES if feature not in metadata.columns]
    if missing:
        raise ValueError(f"Metadata is missing requested columns: {', '.join(missing)}")

    groups: dict[str, pd.Series] = {}
    for feature in FEATURES:
        if feature == "age":
            groups[feature] = make_age_groups(metadata[feature])
        elif feature == "preoperative_psa":
            groups[feature] = make_psa_groups(metadata[feature], psa_quantiles)
        elif feature == "pathologic_stage":
            groups[feature] = make_pathologic_stage_groups(metadata[feature])
        elif feature == "patient_gleason_score":
            groups[feature] = make_gleason_groups(metadata[feature])
        else:
            groups[feature] = normalize_missing(metadata[feature])
    return groups


def load_inputs(
    tpm_path: Path,
    metadata_path: Path,
    selected_genes: list[str] | None,
    max_genes: int | None,
):
    metadata = pd.read_csv(metadata_path, sep="\t", dtype="string")
    required = {"geo_accession", "expression_sample_id"}
    absent = required.difference(metadata.columns)
    if absent:
        raise ValueError(f"Metadata lacks alignment columns: {', '.join(sorted(absent))}")

    metadata["geo_accession"] = normalize_missing(metadata["geo_accession"])
    metadata["expression_sample_id"] = normalize_missing(metadata["expression_sample_id"])
    metadata = metadata.dropna(subset=["geo_accession", "expression_sample_id"]).copy()
    if metadata["expression_sample_id"].duplicated().any():
        duplicates = metadata.loc[
            metadata["expression_sample_id"].duplicated(), "expression_sample_id"
        ].tolist()
        raise ValueError(f"Duplicate expression sample IDs: {duplicates[:5]}")

    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
    tpm.index = tpm.index.astype(str).str.strip()
    if tpm.index.duplicated().any():
        raise ValueError("TPM gene identifiers must be unique")
    if max_genes is not None:
        tpm = tpm.iloc[:max_genes]

    if selected_genes:
        requested = list(dict.fromkeys(gene.strip() for gene in selected_genes))
        missing_genes = [gene for gene in requested if gene not in tpm.index]
        if missing_genes:
            raise ValueError(f"Requested genes not found in TPM: {', '.join(missing_genes)}")
    else:
        requested = tpm.index.tolist()

    sample_ids = metadata["expression_sample_id"].tolist()
    available = [sample for sample in sample_ids if sample in tpm.columns]
    if not available:
        raise ValueError("No metadata expression_sample_id values match TPM columns")

    metadata = metadata.set_index("expression_sample_id").loc[available].reset_index()
    values = tpm.loc[:, available].apply(pd.to_numeric, errors="coerce").to_numpy(
        dtype=np.float64
    )
    genes = tpm.index.to_numpy(dtype=str)
    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    anchor_indices = np.array([gene_to_index[gene] for gene in requested], dtype=int)
    return metadata, genes, values, anchor_indices


def pairwise_correlations(
    expression: np.ndarray,
    left_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate pairwise-complete Pearson r, p, and n for one gene chunk."""
    left = expression[left_indices]
    right = expression
    left_ok = np.isfinite(left).astype(np.float64)
    right_ok = np.isfinite(right).astype(np.float64)
    left_zero = np.nan_to_num(left, nan=0.0, posinf=0.0, neginf=0.0)
    right_zero = np.nan_to_num(right, nan=0.0, posinf=0.0, neginf=0.0)

    n = left_ok @ right_ok.T
    sum_x = left_zero @ right_ok.T
    sum_y = left_ok @ right_zero.T
    sum_x2 = (left_zero * left_zero) @ right_ok.T
    sum_y2 = left_ok @ (right_zero * right_zero).T
    sum_xy = left_zero @ right_zero.T

    numerator = n * sum_xy - sum_x * sum_y
    denominator = np.sqrt(
        np.maximum(n * sum_x2 - sum_x * sum_x, 0.0)
        * np.maximum(n * sum_y2 - sum_y * sum_y, 0.0)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = numerator / denominator
    corr = np.clip(corr, -1.0, 1.0)
    corr[(n < 3) | (denominator == 0)] = np.nan

    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = np.abs(corr) * np.sqrt((n - 2.0) / np.maximum(1.0 - corr * corr, 1e-15))
        pvalue = 2.0 * student_t.sf(statistic, df=n - 2.0)
    pvalue[~np.isfinite(corr)] = np.nan
    return corr, pvalue, n


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
    temporary.replace(path)


def process_group(
    *,
    feature: str,
    group: str,
    group_metadata: pd.DataFrame,
    genes: np.ndarray,
    expression: np.ndarray,
    anchor_indices: np.ndarray,
    output_root: Path,
    min_samples: int,
    correlation_cutoff: float,
    pvalue_cutoff: float,
    chunk_size: int,
    overwrite: bool,
) -> dict[str, object]:
    sample_indices = group_metadata["_sample_index"].to_numpy(dtype=int)
    sample_count = len(sample_indices)
    if sample_count < min_samples:
        return {
            "feature": feature, "group": group, "samples": sample_count,
            "status": "skipped_too_few_samples",
        }

    subgroup_expression = expression[:, sample_indices]
    gsm_ids = group_metadata["geo_accession"].astype(str).tolist()
    group_dir = output_root / safe_name(feature) / safe_name(group)
    gene_dir = group_dir / "genes"
    feature_dir = group_dir / "features"
    gene_dir.mkdir(parents=True, exist_ok=True)
    feature_dir.mkdir(parents=True, exist_ok=True)

    sample_manifest = group_metadata[
        ["expression_sample_id", "geo_accession", feature]
    ].to_dict(orient="records")
    write_json(group_dir / "samples.json", sample_manifest)

    for gene_index in anchor_indices:
        gene = genes[gene_index]
        path = gene_dir / f"{safe_name(gene)}.json"
        if path.exists() and not overwrite:
            continue
        values = subgroup_expression[gene_index]
        record: dict[str, object] = {"ENSG_ID": str(gene)}
        record.update(
            {
                gsm: None if not np.isfinite(value) else float(value)
                for gsm, value in zip(gsm_ids, values)
            }
        )
        write_json(path, record)

    comparison_gene_count = len(genes)
    anchor_count = len(anchor_indices)
    for start in range(0, anchor_count, chunk_size):
        stop = min(start + chunk_size, anchor_count)
        chunk_indices = anchor_indices[start:stop]
        corr, pvalue, paired_n = pairwise_correlations(subgroup_expression, chunk_indices)
        for local_index, gene_index in enumerate(chunk_indices):
            gene = str(genes[gene_index])
            path = feature_dir / f"{safe_name(gene)}_feature.json"
            if path.exists() and not overwrite:
                continue

            keep = (
                np.isfinite(corr[local_index])
                & (np.abs(corr[local_index]) > correlation_cutoff)
                & (pvalue[local_index] < pvalue_cutoff)
            )
            keep[gene_index] = False
            partners = np.flatnonzero(keep)
            partners = partners[np.argsort(-np.abs(corr[local_index, partners]))]
            results = [
                {
                    "Feature1": gene,
                    "Feature2": str(genes[partner]),
                    "pc": float(corr[local_index, partner]),
                    "p-value": float(pvalue[local_index, partner]),
                    "n": int(paired_n[local_index, partner]),
                }
                for partner in partners
            ]
            # Write empty lists too, making resume/completeness checks unambiguous.
            write_json(path, results)

        print(f"  {feature}/{group}: anchor genes {start + 1}-{stop} of {anchor_count}", flush=True)

    return {
        "feature": feature,
        "group": group,
        "samples": sample_count,
        "genes": anchor_count,
        "comparison_genes": comparison_gene_count,
        "status": "completed",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tpm", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("GSE201284_correlations"))
    parser.add_argument("--min-samples", type=int, default=3)
    parser.add_argument("--correlation-cutoff", type=float, default=0.3)
    parser.add_argument("--pvalue-cutoff", type=float, default=0.05)
    parser.add_argument("--psa-quantiles", type=int, default=4)
    parser.add_argument("--chunk-size", type=int, default=128)
    parser.add_argument(
        "--features", nargs="+", choices=FEATURES, default=FEATURES,
        help="Optional subset; the default processes all requested features.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--genes", nargs="+", default=None,
        help="Trial mode: process named genes such as --genes A2M. Omit for all genes.",
    )
    parser.add_argument(
        "--max-genes", type=int, default=None,
        help="Testing only: process the first N genes. Omit for all genes.",
    )
    args = parser.parse_args()

    if args.min_samples < 3:
        parser.error("--min-samples must be at least 3")
    if not 0 < args.pvalue_cutoff <= 1:
        parser.error("--pvalue-cutoff must be in (0, 1]")
    if not 0 <= args.correlation_cutoff <= 1:
        parser.error("--correlation-cutoff must be in [0, 1]")
    if args.genes and args.max_genes is not None:
        parser.error("Use either --genes or --max-genes, not both")

    metadata, genes, expression, anchor_indices = load_inputs(
        args.tpm, args.metadata, args.genes, args.max_genes
    )
    metadata["_sample_index"] = np.arange(len(metadata))
    grouped_features = build_groups(metadata, args.psa_quantiles)
    args.output_root.mkdir(parents=True, exist_ok=True)

    manifest: list[dict[str, object]] = []
    for feature in args.features:
        metadata[feature] = grouped_features[feature]
        valid = metadata.dropna(subset=[feature]).copy()
        if valid.empty:
            manifest.append({"feature": feature, "status": "skipped_no_valid_values"})
            continue

        for group, group_metadata in valid.groupby(feature, sort=True, dropna=True):
            print(f"Processing {feature}/{group} ({len(group_metadata)} samples)", flush=True)
            manifest.append(
                process_group(
                    feature=feature,
                    group=str(group),
                    group_metadata=group_metadata,
                    genes=genes,
                    expression=expression,
                    anchor_indices=anchor_indices,
                    output_root=args.output_root,
                    min_samples=args.min_samples,
                    correlation_cutoff=args.correlation_cutoff,
                    pvalue_cutoff=args.pvalue_cutoff,
                    chunk_size=args.chunk_size,
                    overwrite=args.overwrite,
                )
            )

    with (args.output_root / "run_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fieldnames = sorted({key for row in manifest for key in row})
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)

    parameters = {
        "tpm": str(args.tpm),
        "metadata": str(args.metadata),
        "features": args.features,
        "min_samples": args.min_samples,
        "correlation_cutoff": args.correlation_cutoff,
        "pvalue_cutoff": args.pvalue_cutoff,
        "psa_quantiles": args.psa_quantiles,
        "selected_genes": args.genes,
        "genes_processed": len(anchor_indices),
        "comparison_genes": len(genes),
        "samples_aligned": len(metadata),
    }
    write_json(args.output_root / "run_parameters.json", parameters)
    print(f"Finished. Results: {args.output_root.resolve()}")


if __name__ == "__main__":
    main()
