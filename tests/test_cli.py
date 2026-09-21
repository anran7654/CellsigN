"""Tests for command-line direction and output-directory selection."""

from pathlib import Path
import subprocess
import sys

import pytest

import pandas as pd

from main.cli import (
    _requested_pairs,
    _resolved_output_dir,
    _selected_cell_type_genes,
    parse_args,
)


def test_python_module_entry_point():
    result = subprocess.run(
        [sys.executable, "-m", "main", "--help"],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--de-threshold" in result.stdout


def test_expression_stem_is_appended_to_output_directory():
    assert _resolved_output_dir("Results", "Data_Chung2017_Breast_all.h5ad") == Path(
        "Results_Data_Chung2017_Breast_all"
    )


def test_expression_suffix_is_not_duplicated():
    resolved = "Results_Data_Chung2017_Breast_all"
    assert _resolved_output_dir(resolved, "Data_Chung2017_Breast_all.h5ad") == Path(resolved)


def test_paths_without_expression_use_explicit_directory_literally():
    assert _resolved_output_dir("Results_Data_Chung2017_Breast_all", None) == Path(
        "Results_Data_Chung2017_Breast_all"
    )


def test_default_scope_is_malignant_other_bidirectional():
    cell_types = ["B_cell", "Malignant", "Myeloid", "T_cell"]
    assert _requested_pairs([], cell_types, "Malignant") == [
        ("B_cell", "Malignant"),
        ("Malignant", "B_cell"),
        ("Myeloid", "Malignant"),
        ("Malignant", "Myeloid"),
        ("T_cell", "Malignant"),
        ("Malignant", "T_cell"),
    ]


def test_explicit_pairs_override_default_scope():
    cell_types = ["B_cell", "Malignant", "Myeloid"]
    assert _requested_pairs(
        ["B_cell:Myeloid"], cell_types, "Malignant"
    ) == [("B_cell", "Myeloid")]


def test_default_scope_requires_configured_malignant_label():
    with pytest.raises(ValueError, match="Malignant label"):
        _requested_pairs([], ["B_cell", "Myeloid"], "Malignant")


def test_receptor_omnibus_path_defaults():
    args = parse_args([
        "--stage", "paths",
        "--intracellular-prior", "intracellular.tsv",
        "--ligand-receptor-prior", "ligand_receptor.tsv",
    ])
    assert args.path_permutations == 500
    assert args.path_alpha == 0.05
    assert args.k_paths == 5
    assert args.de_threshold == 0.05
    assert not hasattr(args, "de_correction")
    assert not hasattr(args, "de_pvalue")
    assert not hasattr(args, "de_alpha")
    assert not hasattr(args, "cell_top_genes")
    assert not hasattr(args, "node_weight")


def test_cell_type_gene_selection_can_use_raw_p_value():
    table = pd.DataFrame({
        "Gene": ["G1", "G2", "G3", "G4"],
        "CellType": ["Receiver"] * 4,
        "Reference": ["All_other_cells"] * 4,
        "Rank": [1, 2, 3, 4],
        "LogFC": [1.0, 0.5, 0.8, -1.0],
        "P_Value": [0.001, 0.049, 0.051, 0.001],
    })
    selected = _selected_cell_type_genes(table, 0.05)
    assert selected["Gene"].tolist() == ["G1", "G2"]
    assert selected.columns.tolist() == [
        "Gene", "CellType", "Reference", "DE_Rank", "DE_LogFC", "P_Value",
    ]
