"""Input, preprocessing, and provenance helpers for CellSigN."""

from __future__ import annotations

import gzip
import hashlib
import warnings
from pathlib import Path
from typing import TextIO

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

def open_text(path: str | Path) -> TextIO:
    """Open plain-text or gzip-compressed UTF-8 input."""
    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def file_sha256(path: str | Path) -> str:
    """Calculate a file SHA-256 digest without loading the file into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_expression(
    expression_path: str | Path,
    cell_path: str | Path | None,
    cell_type_column: str,
) -> ad.AnnData:
    """Read an h5ad file or a gene-by-cell delimited expression matrix."""
    expression_path = Path(expression_path)
    if expression_path.suffix.lower() == ".h5ad":
        adata = sc.read_h5ad(expression_path)
        if cell_type_column not in adata.obs:
            raise ValueError(
                f"The h5ad file does not contain obs[{cell_type_column!r}]."
            )
    else:
        frame = pd.read_csv(expression_path, sep=None, engine="python", index_col=0)
        if frame.empty:
            raise ValueError("The expression matrix is empty.")
        frame = frame.apply(pd.to_numeric, errors="raise")
        if frame.index.has_duplicates:
            raise ValueError("Gene identifiers in the expression matrix must be unique.")
        if frame.columns.has_duplicates:
            raise ValueError("Cell identifiers in the expression matrix must be unique.")
        if cell_path is None:
            raise ValueError("--cell-path is required for a delimited expression matrix.")
        labels = _read_cell_labels(cell_path, frame.columns, cell_type_column)
        adata = ad.AnnData(
            X=csr_matrix(frame.T.to_numpy(dtype=float)),
            obs=pd.DataFrame({cell_type_column: labels}, index=frame.columns.astype(str)),
            var=pd.DataFrame(index=frame.index.astype(str)),
        )

    adata.var_names_make_unique()
    valid = adata.obs[cell_type_column].notna() & (
        adata.obs[cell_type_column].astype(str).str.strip() != ""
    )
    adata = adata[valid].copy()
    adata.obs[cell_type_column] = adata.obs[cell_type_column].astype(str).astype("category")
    if adata.n_obs == 0:
        raise ValueError("No cells with a valid cell-type label remain.")
    return adata


def _read_cell_labels(
    cell_path: str | Path,
    expected_cells: pd.Index,
    cell_type_column: str,
) -> pd.Series:
    """Read a two-column cell-label table with or without a header."""
    raw = pd.read_csv(cell_path, sep=None, engine="python", header=None, dtype=str)
    if raw.shape[1] < 2:
        raise ValueError("The cell-label file must contain at least two columns.")
    raw = raw.iloc[:, :2].copy()
    first = [str(x).strip().lower() for x in raw.iloc[0].tolist()]
    if first[0] in {"cell", "cell_id", "cellid"} or first[1] in {
        "cell_type", "celltype", cell_type_column.lower()
    }:
        raw = raw.iloc[1:]
    raw.columns = ["cell_id", cell_type_column]
    raw["cell_id"] = raw["cell_id"].astype(str)
    if raw["cell_id"].duplicated().any():
        raise ValueError("Cell identifiers in the label file must be unique.")
    mapping = raw.set_index("cell_id")[cell_type_column]
    missing = expected_cells.astype(str).difference(mapping.index)
    if len(missing):
        preview = ", ".join(missing[:5])
        raise ValueError(f"Missing labels for {len(missing)} cells, including: {preview}")
    return mapping.reindex(expected_cells.astype(str))


def preprocess(
    adata: ad.AnnData,
    min_cell_fraction: float,
    min_gene_fraction: float,
    normalize: bool,
    log1p: bool,
    target_sum: float | None,
    hvg_top_genes: int,
) -> ad.AnnData:
    """Filter, normalize, log-transform, and select highly variable genes."""
    adata = adata.copy()
    # Preserve the filtering order and integer conversion used by the original
    # CellSigN scripts: cells first, then genes using the retained cell count.
    min_genes = max(1, int(min_gene_fraction * adata.n_vars))
    sc.pp.filter_cells(adata, min_genes=min_genes)
    min_cells = max(1, int(min_cell_fraction * adata.n_obs))
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.uns["cellsign_filtered_cells"] = int(adata.n_obs)
    adata.uns["cellsign_filtered_genes"] = int(adata.n_vars)
    if normalize:
        if target_sum is None:
            sc.pp.normalize_total(adata)
        else:
            sc.pp.normalize_total(adata, target_sum=target_sum)
    if log1p:
        sc.pp.log1p(adata)
    if hvg_top_genes <= 0:
        raise ValueError("--hvg-top-genes must be a positive integer.")
    requested = min(hvg_top_genes, adata.n_vars)
    if requested == adata.n_vars:
        selected = np.ones(adata.n_vars, dtype=bool)
        adata.var["highly_variable"] = selected
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=requested)
        selected = np.asarray(adata.var["highly_variable"], dtype=bool)
    if not selected.any():
        raise ValueError("Highly variable gene selection returned no genes.")
    return adata[:, selected].copy()


def dense_matrix(adata: ad.AnnData) -> np.ndarray:
    """Return observations by variables as a finite dense float array."""
    matrix = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    matrix = np.asarray(matrix, dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError("The expression matrix contains NaN or infinite values.")
    return matrix


TOP_GENE_COLUMNS = ["Gene", "Rank", "LogFC", "P_Value"]

ONE_VS_REST_DE_COLUMNS = [
    "Gene", "CellType", "Reference", "Rank", "LogFC", "P_Value",
]


def eligible_cell_types(
    adata: ad.AnnData,
    cell_type_column: str,
    minimum_cells: int = 5,
) -> list[str]:
    """Return cell types with enough cells for one-vs-rest analysis."""
    categories = [str(x) for x in adata.obs[cell_type_column].cat.categories]
    labels = adata.obs[cell_type_column].astype(str)
    counts = labels.value_counts()
    eligible = [
        cell_type for cell_type in categories
        if int(counts.get(cell_type, 0)) >= minimum_cells
    ]
    skipped = [cell_type for cell_type in categories if cell_type not in eligible]
    if skipped:
        details = ", ".join(
            f"{cell_type} (n={int(counts.get(cell_type, 0))})"
            for cell_type in skipped
        )
        warnings.warn(
            f"Skipping cell types with fewer than {minimum_cells} cells: {details}",
            RuntimeWarning,
            stacklevel=2,
        )
    if len(eligible) < 2:
        raise ValueError(
            f"At least two cell types with {minimum_cells} or more cells are required."
        )
    return eligible


def one_vs_rest_cell_type_de_tables(
    adata: ad.AnnData,
    cell_type_column: str,
    cell_types: list[str],
) -> dict[str, pd.DataFrame]:
    """Test each requested cell type once against all remaining cells.

    Every retained HVG is ranked by a two-sided Wilcoxon rank-sum test.  Raw
    p-values are reported and used for differential-gene selection.
    """
    labels = adata.obs[cell_type_column].astype(str)
    available = set(labels)
    unknown = [cell_type for cell_type in cell_types if cell_type not in available]
    if unknown:
        raise ValueError(f"Unknown cell type(s) in one-vs-rest DE: {unknown}")
    result: dict[str, pd.DataFrame] = {}
    ranking_adata = adata.copy()
    ranking_column = "_cellsign_de_group"
    receiver_label = "selected_cell_type"
    background_label = "all_other_cells"
    for cell_type in dict.fromkeys(cell_types):
        selected = labels.to_numpy() == cell_type
        if int(selected.sum()) < 5 or int((~selected).sum()) < 5:
            raise ValueError(
                f"One-vs-rest DE for {cell_type!r} requires at least five "
                "cells in the selected type and five cells in all other types."
            )
        ranking_adata.obs[ranking_column] = pd.Categorical(
            np.where(selected, receiver_label, background_label),
            categories=[receiver_label, background_label],
        )
        sc.tl.rank_genes_groups(
            ranking_adata,
            groupby=ranking_column,
            groups=[receiver_label],
            reference=background_label,
            method="wilcoxon",
            n_genes=ranking_adata.n_vars,
            pts=False,
        )
        ranked = sc.get.rank_genes_groups_df(
            ranking_adata, group=receiver_label
        ).reset_index(drop=True)
        raw_p = pd.to_numeric(ranked["pvals"], errors="coerce").to_numpy()
        result[cell_type] = pd.DataFrame({
            "Gene": ranked["names"].astype(str),
            "CellType": cell_type,
            "Reference": "All_other_cells",
            "Rank": np.arange(1, len(ranked) + 1),
            "LogFC": pd.to_numeric(ranked["logfoldchanges"], errors="coerce"),
            "P_Value": raw_p,
        }, columns=ONE_VS_REST_DE_COLUMNS)
    return result


def select_cell_type_gene_tables(
    adata: ad.AnnData,
    cell_type_column: str,
    cell_top_genes: int,
) -> dict[str, pd.DataFrame]:
    """Return per-cell-type positive-logFC Wilcoxon rankings.

    A value of zero keeps the full HVG set and records the ranking fields as
    missing. This option is useful for small examples and diagnostic runs.
    """
    all_groups = [str(x) for x in adata.obs[cell_type_column].cat.categories]
    labels = adata.obs[cell_type_column].astype(str)
    counts = labels.value_counts()
    minimum_cells = 5
    minimum_background_cells = 4
    groups = [
        group for group in all_groups
        if int(counts.get(group, 0)) >= minimum_cells
        and int(adata.n_obs - counts.get(group, 0)) >= minimum_background_cells
    ]
    skipped = [group for group in all_groups if group not in groups]
    if skipped:
        details = ", ".join(
            f"{group} (n={int(counts.get(group, 0))}, "
            f"background={int(adata.n_obs - counts.get(group, 0))})"
            for group in skipped
        )
        warnings.warn(
            f"Skipping cell types with fewer than {minimum_cells} cells or "
            f"fewer than {minimum_background_cells} background cells: {details}",
            RuntimeWarning,
            stacklevel=2,
        )
    if not groups:
        raise ValueError(
            "No cell type has at least five receiver and four background cells "
            "after preprocessing."
        )
    if cell_top_genes < 0:
        raise ValueError("--cell-top-genes must be zero or a positive integer.")
    if cell_top_genes == 0:
        genes = list(map(str, adata.var_names))
        return {
            group: pd.DataFrame({
                "Gene": genes,
                "Rank": np.arange(1, len(genes) + 1),
                "LogFC": np.nan,
                "P_Value": np.nan,
            }, columns=TOP_GENE_COLUMNS)
            for group in groups
        }
    n_genes = min(cell_top_genes, adata.n_vars)
    # Rank each eligible cell type against a binary background separately.
    # Some Scanpy versions validate every category in ``groupby`` even when
    # ``groups`` is restricted, which makes a skipped singleton category abort
    # the entire ranking call.  A binary working label keeps skipped cells in
    # the background while ensuring that only the eligible receiver is tested.
    ranking_adata = adata.copy()
    ranking_column = "_cellsign_rank_group"
    receiver_label = "receiver"
    background_label = "background"
    result: dict[str, pd.DataFrame] = {}
    for group in groups:
        ranking_adata.obs[ranking_column] = pd.Categorical(
            np.where(labels.to_numpy() == group, receiver_label, background_label),
            categories=[receiver_label, background_label],
        )
        sc.tl.rank_genes_groups(
            ranking_adata,
            groupby=ranking_column,
            groups=[receiver_label],
            reference=background_label,
            method="wilcoxon",
            n_genes=n_genes,
            pts=False,
        )
        table = sc.get.rank_genes_groups_df(ranking_adata, group=receiver_label)
        table = table[np.isfinite(table["logfoldchanges"]) & (table["logfoldchanges"] > 0)]
        table = table.head(n_genes).reset_index(drop=True)
        result[group] = pd.DataFrame({
            "Gene": table["names"].astype(str),
            "Rank": np.arange(1, len(table) + 1),
            "LogFC": pd.to_numeric(table["logfoldchanges"], errors="coerce"),
            "P_Value": pd.to_numeric(table["pvals"], errors="coerce"),
        }, columns=TOP_GENE_COLUMNS)
    return result


def select_cell_type_genes(
    adata: ad.AnnData,
    cell_type_column: str,
    cell_top_genes: int,
) -> dict[str, list[str]]:
    """Compatibility wrapper returning gene names only."""
    tables = select_cell_type_gene_tables(adata, cell_type_column, cell_top_genes)
    return {cell_type: table["Gene"].tolist() for cell_type, table in tables.items()}
