"""Command-line interface for the staged CellSigN 0.3.5 workflow."""

from __future__ import annotations

import argparse
from importlib.metadata import version as package_version
import json
import platform
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import __version__
from .io import (
    dense_matrix,
    eligible_cell_types,
    file_sha256,
    one_vs_rest_cell_type_de_tables,
    preprocess,
    read_expression,
)
from .network import EDGE_EVIDENCE_COLUMNS, GENE_EVIDENCE_COLUMNS, build_cell_network
from .pathways import attach_sender_paths, prepare_receiver_path_cache
from .priors import load_intracellular_prior, load_ligand_receptor_prior


TEST_COUNT_COLUMNS = [
    "Direction", "Sender", "Receiver", "N_Sender_Cells", "N_Receiver_Cells",
    "N_Input_Genes", "N_Filtered_Genes", "N_HVG",
    "N_DE_Tests_Per_Cell_Type", "N_Sender_DEGs",
    "N_Receiver_DEGs", "N_Candidate_Gene_Pairs", "N_Pearson",
    "N_Nonlinear_Edges", "N_Nonlinear_Model_Fits", "N_Retained_Edges",
    "N_RTF_Pairs", "N_Path_Permutations", "DE_Threshold",
    "Pair_BH_Method", "Path_BH_Method",
    "Receptor_BH_Method",
]

EVIDENCE_PARAMETER_NAMES = [
    "cell_type_column", "sample_column", "min_cell_fraction",
    "min_gene_fraction", "target_sum", "normalize", "log1p", "hvg_top_genes",
    "de_threshold", "alpha", "min_abs_r", "min_delta_r2",
    "stability_resamples", "relation_weight", "stability_weight",
    "hub_penalty_weight", "hop_penalty", "pair", "malignant_label", "seed",
]

PATH_PARAMETER_NAMES = [
    "pair", "k_paths", "path_permutations", "path_alpha", "seed",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Infer data-driven, prior-annotated candidate "
            "L-R-M-TF-TG paths from scRNA-seq data."
        )
    )
    parser.add_argument(
        "--stage", choices=["all", "evidence", "paths"], default="all",
        help="Run both stages, evidence only, or paths from saved evidence.",
    )
    parser.add_argument("--expression", help="Gene-by-cell table or h5ad file.")
    parser.add_argument("--cell-path", help="Two-column cell ID and cell-type table.")
    parser.add_argument("--cell-type-column", default="celltype")
    parser.add_argument(
        "--sample-column", default=None,
        help="Optional h5ad obs column recorded as sample provenance.",
    )
    parser.add_argument("--intracellular-prior", required=True)
    parser.add_argument("--ligand-receptor-prior", required=True)
    parser.add_argument("--output-dir", default="Results")
    parser.add_argument("--min-cell-fraction", type=float, default=0.01)
    parser.add_argument("--min-gene-fraction", type=float, default=0.01)
    parser.add_argument("--target-sum", type=float, default=None)
    parser.add_argument(
        "--no-normalize", "--skip-normalize", dest="normalize", action="store_false"
    )
    parser.add_argument("--no-log1p", "--skip-log1p", dest="log1p", action="store_false")
    parser.set_defaults(normalize=True, log1p=True)
    parser.add_argument("--hvg-top-genes", type=int, default=5000)
    parser.add_argument(
        "--de-threshold", "--de-pvalue", dest="de_threshold", type=float,
        default=0.05,
        help="Raw one-vs-rest Wilcoxon p-value threshold (default: 0.05).",
    )
    parser.add_argument(
        "--alpha", type=float, default=0.05,
        help="Pearson/nonlinear edge BH threshold (default: 0.05).",
    )
    parser.add_argument("--min-abs-r", type=float, default=0.1)
    parser.add_argument("--min-delta-r2", type=float, default=0.01)
    parser.add_argument("--stability-resamples", type=int, default=20)
    parser.add_argument("--relation-weight", type=float, default=0.6)
    parser.add_argument("--stability-weight", type=float, default=0.2)
    parser.add_argument("--hub-penalty-weight", type=float, default=0.1)
    parser.add_argument("--hop-penalty", type=float, default=0.05)
    parser.add_argument(
        "--k-paths", type=int, default=5,
        help="Number of lowest-cost simple paths retained per R-TF pair (default: 5).",
    )
    parser.add_argument("--path-permutations", type=int, default=500)
    parser.add_argument("--path-alpha", type=float, default=0.05)
    parser.add_argument(
        "--data-edges-per-node", type=int, default=20, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--pair", action="append", default=[], metavar="SENDER:RECEIVER",
        help="Analyze one ordered pair; repeat as needed.",
    )
    parser.add_argument("--malignant-label", default="Malignant")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    if args.stage in {"all", "evidence"} and not args.expression:
        parser.error("--expression is required for --stage all or evidence.")
    for name in ("alpha", "de_threshold", "path_alpha"):
        if not 0 < getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 1.")
    for name in ("min_cell_fraction", "min_gene_fraction"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 1.")
    if args.hvg_top_genes <= 0:
        parser.error("--hvg-top-genes must be positive.")
    if args.k_paths <= 0:
        parser.error("--k-paths must be positive.")
    if min(args.stability_resamples, args.path_permutations) < 0:
        parser.error("Resample and permutation counts cannot be negative.")
    if not 0 <= args.min_abs_r <= 1 or not 0 <= args.min_delta_r2 <= 1:
        parser.error("--min-abs-r and --min-delta-r2 must be between 0 and 1.")
    weights = [
        args.relation_weight, args.stability_weight,
        args.hub_penalty_weight, args.hop_penalty,
    ]
    if min(weights) < 0 or args.relation_weight <= 0:
        parser.error("Cost weights must be nonnegative and --relation-weight must be positive.")
    return args


def _safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("._") or "cell_type"


def _direction_key(sender: str, receiver: str) -> str:
    return f"{sender}:{receiver}"


def _direction_stem(sender: str, receiver: str) -> str:
    return f"{_safe_name(sender)}_to_{_safe_name(receiver)}"


def _cell_type_stem(cell_type: str) -> str:
    return _safe_name(cell_type)


def _resolved_output_dir(output_dir: str | Path, expression: str | None) -> Path:
    base = Path(output_dir)
    if not expression:
        return base
    suffix = f"_{_safe_name(Path(expression).stem)}"
    if base.name.casefold().endswith(suffix.casefold()):
        return base
    return base.with_name(f"{base.name or 'Results'}{suffix}")


def _requested_pairs(
    specifications: list[str], cell_types: list[str], malignant_label: str,
) -> list[tuple[str, str]]:
    if not specifications:
        if malignant_label not in cell_types:
            raise ValueError(
                f"Malignant label {malignant_label!r} was not found. Available: {cell_types}"
            )
        pairs: list[tuple[str, str]] = []
        for cell_type in cell_types:
            if cell_type != malignant_label:
                pairs.extend([(cell_type, malignant_label), (malignant_label, cell_type)])
        return pairs
    pairs = []
    for item in specifications:
        if ":" not in item:
            raise ValueError(f"Invalid --pair {item!r}; expected SENDER:RECEIVER.")
        sender, receiver = item.split(":", 1)
        if sender not in cell_types or receiver not in cell_types:
            raise ValueError(f"Unknown cell type in --pair {item!r}. Available: {cell_types}")
        if sender == receiver:
            raise ValueError("Sender and receiver must be different cell types.")
        pairs.append((sender, receiver))
    return list(dict.fromkeys(pairs))


def _write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False, na_rep="NA")


def _software_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(), "anndata": package_version("anndata"),
        "networkx": package_version("networkx"), "numpy": np.__version__,
        "pandas": pd.__version__, "scanpy": package_version("scanpy"),
        "scipy": package_version("scipy"),
    }


def _parameter_subset(args: argparse.Namespace, names: list[str]) -> dict[str, Any]:
    return {name: getattr(args, name) for name in names}


def _base_manifest(args: argparse.Namespace) -> dict[str, Any]:
    manifest_path = Path(args.output_dir) / "other" / "run_manifest.json"
    if manifest_path.is_file():
        with manifest_path.open("r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    else:
        manifest = {}
    manifest.update({
        "cellsign_version": __version__,
        "interpretation": (
            "one-vs-rest-DE-filtered, data-driven, prior-annotated "
            "candidate signaling pathways"
        ),
        "software": _software_versions(),
    })
    return manifest


def _selected_cell_type_genes(
    de_table: pd.DataFrame,
    de_threshold: float,
) -> pd.DataFrame:
    selected = de_table.rename(columns={
        "Rank": "DE_Rank", "LogFC": "DE_LogFC",
    })
    selected = selected[
        np.isfinite(selected["DE_LogFC"])
        & (selected["DE_LogFC"] > 0)
        & np.isfinite(selected["P_Value"])
        & (selected["P_Value"] < de_threshold)
    ]
    selected = selected.sort_values(
        ["P_Value", "DE_Rank", "Gene"], kind="mergesort"
    ).reset_index(drop=True)
    columns = [
        "Gene", "CellType", "Reference", "DE_Rank", "DE_LogFC",
        "P_Value",
    ]
    return selected[columns]


def _run_evidence_stage(
    args: argparse.Namespace, gene_dir: Path, evidence_dir: Path, other_dir: Path,
) -> dict[str, Any]:
    print("[evidence] reading expression and preprocessing", flush=True)
    adata = read_expression(args.expression, args.cell_path, args.cell_type_column)
    if args.sample_column and args.sample_column not in adata.obs:
        raise ValueError(f"The h5ad file does not contain obs[{args.sample_column!r}].")
    input_cells, input_genes = int(adata.n_obs), int(adata.n_vars)
    sample_count = (
        int(adata.obs[args.sample_column].nunique()) if args.sample_column else None
    )
    adata = preprocess(
        adata, args.min_cell_fraction, args.min_gene_fraction, args.normalize,
        args.log1p, args.target_sum, args.hvg_top_genes,
    )
    filtered_cells = int(adata.uns.get("cellsign_filtered_cells", adata.n_obs))
    filtered_genes = int(adata.uns.get("cellsign_filtered_genes", adata.n_vars))
    matrix = dense_matrix(adata)
    hvg_genes = list(map(str, adata.var_names))
    with (other_dir / "processed_gene_list.txt").open("w", encoding="utf-8") as handle:
        handle.write("\n".join(hvg_genes) + "\n")

    prior = load_intracellular_prior(args.intracellular_prior, set(hvg_genes))
    ligand_receptor = load_ligand_receptor_prior(args.ligand_receptor_prior)
    cell_types = eligible_cell_types(adata, args.cell_type_column, minimum_cells=5)
    pairs = _requested_pairs(args.pair, cell_types, args.malignant_label)
    labels = adata.obs[args.cell_type_column].astype(str).to_numpy()
    hvg_positions = {gene: index for index, gene in enumerate(hvg_genes)}
    de_cell_types = cell_types
    de_tables = one_vs_rest_cell_type_de_tables(
        adata, args.cell_type_column, de_cell_types,
    )
    selected_by_type: dict[str, pd.DataFrame] = {}
    gene_paths: dict[str, Path] = {}
    deg_count_records: list[dict[str, object]] = []
    used_stems: dict[str, str] = {}
    for cell_type in de_cell_types:
        cell_stem = _cell_type_stem(cell_type)
        if cell_stem in used_stems and used_stems[cell_stem] != cell_type:
            raise ValueError(
                f"Cell-type names {used_stems[cell_stem]!r} and {cell_type!r} "
                "map to the same output filename."
            )
        used_stems[cell_stem] = cell_type
        selected = _selected_cell_type_genes(
            de_tables[cell_type], args.de_threshold
        )
        selected_by_type[cell_type] = selected
        gene_path = gene_dir / f"{cell_stem}_de_genes.txt"
        gene_paths[cell_type] = gene_path
        _write_table(selected, gene_path)
        full_de = de_tables[cell_type].copy()
        full_de["Selected"] = full_de["Gene"].isin(
            set(selected["Gene"].astype(str))
        )
        _write_table(
            full_de, evidence_dir / f"{cell_stem}_one_vs_rest_de.txt"
        )
        n_cells = int(np.sum(labels == cell_type))
        deg_count_records.append({
            "CellType": cell_type,
            "N_Cells": n_cells,
            "N_Other_Cells": int(len(labels) - n_cells),
            "N_HVG_Tests": len(hvg_genes),
            "N_DEGs": len(selected),
            "DE_Threshold": args.de_threshold,
        })
        print(
            f"[evidence] {cell_type} vs all other cells: {len(selected)} DEGs "
            f"(raw P < {args.de_threshold:g}, logFC > 0)",
            flush=True,
        )

    direction_tables: dict[str, dict[str, Any]] = {}
    count_records: list[dict[str, object]] = []
    receiver_networks: dict[str, dict[str, object]] = {}
    receivers = list(dict.fromkeys(receiver for _, receiver in pairs))

    for receiver_index, receiver in enumerate(receivers, start=1):
        print(
            f"[evidence] receiver network {receiver_index}/{len(receivers)}: {receiver}",
            flush=True,
        )
        receiver_genes = selected_by_type[receiver]
        receiver_gene_names = receiver_genes["Gene"].astype(str).tolist()
        receiver_stem = _cell_type_stem(receiver)
        gene_evidence_path = evidence_dir / f"{receiver_stem}_gene_evidence.txt"
        edge_path = evidence_dir / f"{receiver_stem}_edge_evidence.txt"
        candidate_count = 0
        retained_count = 0
        if len(receiver_gene_names) < 2:
            _write_table(pd.DataFrame(columns=GENE_EVIDENCE_COLUMNS), gene_evidence_path)
            _write_table(pd.DataFrame(columns=EDGE_EVIDENCE_COLUMNS), edge_path)
            print("[evidence] fewer than two receiver DE genes; empty network", flush=True)
        else:
            positions = np.asarray(
                [hvg_positions[gene] for gene in receiver_gene_names], dtype=int
            )
            try:
                result = build_cell_network(
                    receiver,
                    matrix[labels == receiver][:, positions],
                    receiver_gene_names,
                    receiver_genes,
                    prior,
                    alpha=args.alpha,
                    min_abs_r=args.min_abs_r,
                    min_delta_r2=args.min_delta_r2,
                    stability_resamples=args.stability_resamples,
                    relation_weight=args.relation_weight,
                    stability_weight=args.stability_weight,
                    hub_penalty_weight=args.hub_penalty_weight,
                    hop_penalty=args.hop_penalty,
                    seed=args.seed + receiver_index,
                )
                _write_table(result.gene_evidence, gene_evidence_path)
                _write_table(result.edge_evidence, edge_path)
                candidate_count = len(result.edge_evidence)
                retained_count = int(result.edge_evidence["Edge_Retained"].sum())
                prior_count = int(
                    result.edge_evidence.loc[
                        result.edge_evidence["Edge_Retained"], "Prior_Found"
                    ].sum()
                )
                print(
                    f"[evidence] {candidate_count} receiver DEG pairs tested; "
                    f"{retained_count} evidence-supported edges retained "
                    f"({prior_count} with prior direction)",
                    flush=True,
                )
            except ValueError as exc:
                _write_table(pd.DataFrame(columns=GENE_EVIDENCE_COLUMNS), gene_evidence_path)
                _write_table(pd.DataFrame(columns=EDGE_EVIDENCE_COLUMNS), edge_path)
                print(f"[evidence] {exc}; empty network", flush=True)
        receiver_networks[receiver] = {
            "gene_evidence": gene_evidence_path,
            "edges": edge_path,
            "candidate_count": candidate_count,
            "retained_count": retained_count,
        }

    for sender, receiver in pairs:
        sender_genes = selected_by_type[sender]
        receiver_genes = selected_by_type[receiver]
        network = receiver_networks[receiver]
        count_records.append({
            "Direction": _direction_key(sender, receiver), "Sender": sender,
            "Receiver": receiver, "N_Sender_Cells": int(np.sum(labels == sender)),
            "N_Receiver_Cells": int(np.sum(labels == receiver)),
            "N_Input_Genes": input_genes, "N_Filtered_Genes": filtered_genes,
            "N_HVG": len(hvg_genes),
            "N_DE_Tests_Per_Cell_Type": len(hvg_genes),
            "N_Sender_DEGs": len(sender_genes),
            "N_Receiver_DEGs": len(receiver_genes),
            "N_Candidate_Gene_Pairs": int(network["candidate_count"]),
            "N_Pearson": int(network["candidate_count"]),
            "N_Nonlinear_Edges": int(network["candidate_count"]),
            "N_Nonlinear_Model_Fits": 2 * int(network["candidate_count"]),
            "N_Retained_Edges": int(network["retained_count"]),
            "N_RTF_Pairs": 0, "N_Path_Permutations": 0,
            "DE_Threshold": args.de_threshold,
            "Pair_BH_Method": "Benjamini-Hochberg over all receiver DEG pairs",
            "Path_BH_Method": "Benjamini-Hochberg across R-TF pairs (diagnostic)",
            "Receptor_BH_Method": "Benjamini-Hochberg across receptors by direction",
        })
        direction_tables[_direction_key(sender, receiver)] = {
            "sender": sender,
            "receiver": receiver,
            "sender_genes": gene_paths[sender],
            "receiver_genes": gene_paths[receiver],
            "edges": network["edges"],
        }

    test_counts = pd.DataFrame(count_records, columns=TEST_COUNT_COLUMNS)
    _write_table(test_counts, other_dir / "test_counts.txt")
    deg_counts = pd.DataFrame(deg_count_records, columns=[
        "CellType", "N_Cells", "N_Other_Cells", "N_HVG_Tests", "N_DEGs",
        "DE_Threshold",
    ])
    _write_table(deg_counts, other_dir / "differential_gene_counts.txt")
    manifest = _base_manifest(args)
    manifest.update({
        "completed_stages": ["evidence"],
        "evidence_parameters": _parameter_subset(args, EVIDENCE_PARAMETER_NAMES),
        "inputs": {
            "expression": {"path": str(Path(args.expression).resolve()), "sha256": file_sha256(args.expression)},
            "cell_path": None if args.cell_path is None else {"path": str(Path(args.cell_path).resolve()), "sha256": file_sha256(args.cell_path)},
            "intracellular_prior": {"path": str(Path(args.intracellular_prior).resolve()), "sha256": file_sha256(args.intracellular_prior)},
            "ligand_receptor_prior": {"path": str(Path(args.ligand_receptor_prior).resolve()), "sha256": file_sha256(args.ligand_receptor_prior), "columns_used": "first two only", "score_filter": None},
        },
        "sample_metadata": {"column": args.sample_column, "n_samples": sample_count},
        "data_dimensions": {
            "input_cells": input_cells, "input_genes": input_genes,
            "filtered_cells": filtered_cells, "filtered_genes": filtered_genes,
            "selected_hvgs": adata.n_vars,
        },
        "directions": [list(pair) for pair in pairs],
        "de_design": "each eligible cell type versus all remaining cells, tested once",
        "de_cell_type_counts": {
            str(record["CellType"]): {
                "cells": int(record["N_Cells"]),
                "other_cells": int(record["N_Other_Cells"]),
                "tested_hvgs": int(record["N_HVG_Tests"]),
                "selected_degs": int(record["N_DEGs"]),
            }
            for record in deg_count_records
        },
        "candidate_edge_policy": (
            "all positive-logFC receiver genes from one-vs-rest DE passing "
            "the raw-p-value de_threshold; all unordered gene "
            "pairs are tested before intracellular prior directions are overlaid"
        ),
        "multiple_testing": {
            "one_vs_rest_de": "None; raw two-sided Wilcoxon p-values are used",
            "pearson": "BH over all DEG pairs within each receiver network",
            "nonlinear": "receiver-only directional p-values Bonferroni-combined per edge, then BH within receiver network",
            "path_pair_diagnostic": "BH across R-TF pairs within each direction",
            "receptor_omnibus": "BH across receptors within each direction",
        },
        "direction_counts": {}, "path_search": {},
    })
    _write_manifest(manifest, other_dir)
    return {
        "cell_types": cell_types, "pairs": pairs,
        "direction_tables": direction_tables, "prior": prior,
        "ligand_receptor": ligand_receptor, "test_counts": test_counts,
        "manifest": manifest,
    }


def _load_saved_context(
    args: argparse.Namespace, gene_dir: Path, evidence_dir: Path, other_dir: Path,
) -> dict[str, Any]:
    counts_path = other_dir / "test_counts.txt"
    genes_path = other_dir / "processed_gene_list.txt"
    if not counts_path.is_file() or not genes_path.is_file():
        raise FileNotFoundError("Saved evidence is incomplete. Run --stage evidence first.")
    test_counts = pd.read_csv(counts_path, sep="\t")
    required = {
        "Direction", "Sender", "Receiver", "DE_Threshold",
    }
    if not required.issubset(test_counts.columns):
        raise ValueError("Saved evidence is not from CellSigN 0.3.5; rerun --stage evidence.")
    available_pairs = list(zip(
        test_counts["Sender"].astype(str), test_counts["Receiver"].astype(str)
    ))
    cell_types = list(dict.fromkeys(
        [item for pair in available_pairs for item in pair]
    ))
    if args.pair:
        pairs = _requested_pairs(args.pair, cell_types, args.malignant_label)
        missing = [pair for pair in pairs if pair not in available_pairs]
        if missing:
            raise ValueError(f"No saved evidence for requested direction(s): {missing}")
    else:
        pairs = available_pairs
    direction_tables: dict[str, dict[str, Any]] = {}
    for sender, receiver in pairs:
        paths = {
            "sender": sender, "receiver": receiver,
            "sender_genes": gene_dir / f"{_cell_type_stem(sender)}_de_genes.txt",
            "receiver_genes": gene_dir / f"{_cell_type_stem(receiver)}_de_genes.txt",
            "edges": evidence_dir / f"{_cell_type_stem(receiver)}_edge_evidence.txt",
        }
        if not all(Path(value).is_file() for key, value in paths.items() if key.endswith("genes") or key == "edges"):
            raise FileNotFoundError(f"Missing saved evidence for {sender} -> {receiver}.")
        direction_tables[_direction_key(sender, receiver)] = paths
    hvg_genes = {
        line.strip() for line in genes_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    print(f"[paths] loaded saved evidence for {len(pairs)} directions", flush=True)
    return {
        "cell_types": cell_types, "pairs": pairs,
        "direction_tables": direction_tables,
        "prior": load_intracellular_prior(args.intracellular_prior, hvg_genes),
        "ligand_receptor": load_ligand_receptor_prior(args.ligand_receptor_prior),
        "test_counts": test_counts, "manifest": _base_manifest(args),
    }


def _run_path_stage(
    args: argparse.Namespace, context: dict[str, Any], output_dir: Path,
    evidence_dir: Path, other_dir: Path,
) -> None:
    pairs: list[tuple[str, str]] = context["pairs"]
    direction_tables: dict[str, dict[str, Any]] = context["direction_tables"]
    prior: pd.DataFrame = context["prior"]
    ligand_receptor: pd.DataFrame = context["ligand_receptor"]
    test_counts: pd.DataFrame = context["test_counts"].copy()
    direction_counts: dict[str, Any] = {}
    path_search: dict[str, Any] = {}

    for direction_index, (sender, receiver) in enumerate(pairs, start=1):
        key = _direction_key(sender, receiver)
        files = direction_tables[key]
        sender_table = pd.read_csv(files["sender_genes"], sep="\t")
        receiver_table = pd.read_csv(files["receiver_genes"], sep="\t")
        edge_table = pd.read_csv(files["edges"], sep="\t")
        sender_genes = set(sender_table.get("Gene", pd.Series(dtype=str)).astype(str))
        receiver_genes = set(receiver_table.get("Gene", pd.Series(dtype=str)).astype(str))
        lr = ligand_receptor[
            ligand_receptor["Ligand"].isin(sender_genes)
            & ligand_receptor["Receptor"].isin(receiver_genes)
        ]
        receptors = set(lr["Receptor"].astype(str))
        print(
            f"[paths] direction {direction_index}/{len(pairs)}: {sender} -> {receiver}; "
            f"{len(receptors)} eligible receptors",
            flush=True,
        )
        cache = prepare_receiver_path_cache(
            receiver, receiver_genes, edge_table, prior, receptors,
            path_permutations=args.path_permutations, k_paths=args.k_paths,
            data_edges_per_node=0, seed=args.seed + 10_000 + direction_index,
            progress=lambda message: print(f"[paths] {message}", flush=True),
        )
        result = attach_sender_paths(
            sender, receiver, sender_genes, receiver_genes, cache,
            ligand_receptor, path_alpha=args.path_alpha,
        )
        stem = _direction_stem(sender, receiver)
        _write_table(result.pathways, output_dir / f"{stem}_pathway.txt")
        _write_table(result.rmtf_paths, evidence_dir / f"{stem}_rmtf_paths.txt")
        _write_table(result.path_edges, evidence_dir / f"{stem}_path_edges.txt")
        _write_table(
            result.subnetwork_edges,
            evidence_dir / f"{stem}_rmtf_subnetwork_edges.txt",
        )
        row_mask = test_counts["Direction"].astype(str) == key
        test_counts.loc[row_mask, "N_RTF_Pairs"] = len(cache.pair_paths)
        test_counts.loc[row_mask, "N_Path_Permutations"] = (
            args.path_permutations if cache.pair_paths else 0
        )
        direction_counts[key] = {
            "rmtf_paths": len(result.rmtf_paths),
            "complete_pathways": len(result.pathways), **result.metadata,
        }
        path_search[key] = {
            "k_paths": args.k_paths,
            "retained_input_edges": cache.retained_input_edges,
            "path_graph_edges": cache.path_graph_edges,
            "cached_rtf_pairs": len(cache.pair_paths),
            "cached_receptors": len(cache.receptor_p),
        }
        print(
            f"[paths] {sender} -> {receiver}: {len(result.rmtf_paths)} R-M-TF paths",
            flush=True,
        )

    _write_table(test_counts[TEST_COUNT_COLUMNS], other_dir / "test_counts.txt")
    manifest = context["manifest"]
    completed = set(manifest.get("completed_stages", []))
    completed.add("paths")
    manifest.update({
        "completed_stages": sorted(completed),
        "path_parameters": _parameter_subset(args, PATH_PARAMETER_NAMES),
        "path_inputs": {
            "intracellular_prior": {"path": str(Path(args.intracellular_prior).resolve()), "sha256": file_sha256(args.intracellular_prior)},
            "ligand_receptor_prior": {"path": str(Path(args.ligand_receptor_prior).resolve()), "sha256": file_sha256(args.ligand_receptor_prior), "columns_used": "first two only", "score_filter": None},
        },
        "ligand_receptor_reading": "all rows; first two columns only; no score filter",
        "path_rule": "weighted k shortest simple directed paths; no maximum path length",
        "path_graph_policy": "intersection of evidence-supported edges and intracellular signaling prior",
        "path_null": "retained path-graph edge costs permuted over fixed topology",
        "direction_counts": direction_counts,
        "path_search": path_search,
    })
    _write_manifest(manifest, other_dir)


def _write_manifest(manifest: dict[str, Any], other_dir: Path) -> None:
    with (other_dir / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2, default=str)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    requested_output_dir = Path(args.output_dir)
    output_dir = _resolved_output_dir(args.output_dir, args.expression)
    args.output_dir = str(output_dir)
    if output_dir != requested_output_dir:
        print(f"[output] using directory: {output_dir}", flush=True)
    gene_dir = output_dir / "gene_lists"
    evidence_dir = output_dir / "evidence"
    other_dir = output_dir / "other"
    for directory in (output_dir, gene_dir, evidence_dir, other_dir):
        directory.mkdir(parents=True, exist_ok=True)
    if args.stage in {"all", "evidence"}:
        context = _run_evidence_stage(args, gene_dir, evidence_dir, other_dir)
    else:
        context = _load_saved_context(args, gene_dir, evidence_dir, other_dir)
    if args.stage in {"all", "paths"}:
        _run_path_stage(args, context, output_dir, evidence_dir, other_dir)
    print(f"CellSigN stage '{args.stage}' completed.", flush=True)


if __name__ == "__main__":
    main()
