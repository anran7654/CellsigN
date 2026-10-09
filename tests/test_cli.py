"""Tests for command-line direction and output-directory selection."""

from argparse import Namespace
import json
from pathlib import Path

import pytest

import pandas as pd

from main.cli import (
    _base_manifest,
    _requested_pairs,
    _resolved_output_dir,
    _selected_cell_type_genes,
    parse_args,
)


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


def test_path_permutation_defaults():
    args = parse_args([
        "--stage", "paths",
        "--intracellular-prior", "intracellular.tsv",
        "--ligand-receptor-prior", "ligand_receptor.tsv",
    ])
    assert args.path_permutations == 500
    assert args.k_paths == 1
    assert not hasattr(args, "path_alpha")
    assert args.de_threshold == 0.05
    assert args.de_gene_target == 500
    assert not hasattr(args, "de_correction")
    assert not hasattr(args, "de_pvalue")
    assert not hasattr(args, "de_alpha")
    assert not hasattr(args, "cell_top_genes")
    assert not hasattr(args, "node_weight")


def test_k_paths_can_still_be_overridden():
    args = parse_args([
        "--stage", "paths",
        "--intracellular-prior", "intracellular.tsv",
        "--ligand-receptor-prior", "ligand_receptor.tsv",
        "--k-paths", "3",
    ])
    assert args.k_paths == 3


def test_base_manifest_removes_legacy_path_tests_without_losing_evidence_settings(
    tmp_path, monkeypatch,
):
    manifest_path = tmp_path / "other" / "run_manifest.json"
    manifest_path.parent.mkdir()
    manifest_path.write_text(json.dumps({
        "completed_stages": ["evidence", "paths"],
        "path_parameters": {
            "k_paths": 3, "path_permutations": 500, "path_alpha": 0.05,
        },
        "evidence_parameters": {"stability_resamples": 100, "alpha": 0.05},
        "multiple_testing": {
            "pearson": "BH within receiver network",
            "nonlinear": "BH within receiver network",
            "path_pair_diagnostic": "old pair BH",
            "receptor_omnibus": "old receptor BH",
        },
    }), encoding="utf-8")
    monkeypatch.setattr("main.cli._software_versions", lambda: {"python": "test"})
    manifest = _base_manifest(Namespace(output_dir=str(tmp_path)))
    assert manifest["path_parameters"] == {
        "k_paths": 3, "path_permutations": 500,
    }
    assert manifest["evidence_parameters"] == {
        "stability_resamples": 100, "alpha": 0.05,
    }
    assert manifest["multiple_testing"] == {
        "pearson": "BH within receiver network",
        "nonlinear": "BH within receiver network",
    }
    assert manifest["completed_stages"] == ["evidence", "paths"]


def test_cell_type_gene_selection_can_use_raw_p_value():
    table = pd.DataFrame({
        "Gene": ["G1", "G2", "G3", "G4"],
        "CellType": ["Receiver"] * 4,
        "Reference": ["All_other_cells"] * 4,
        "Rank": [1, 2, 3, 4],
        "LogFC": [1.0, 0.5, 0.8, -1.0],
        "P_Value": [0.001, 0.049, 0.051, 0.001],
    })
    selected = _selected_cell_type_genes(table, 0.05, de_gene_target=0)
    assert selected["Gene"].tolist() == ["G1", "G2"]
    assert selected.columns.tolist() == [
        "Gene", "CellType", "Reference", "DE_Rank", "DE_LogFC", "P_Value",
        "Selection_Source",
    ]
    assert selected["Selection_Source"].eq("DE").all()


def test_de_gene_target_supplements_by_ascending_non_de_p_value():
    table = pd.DataFrame({
        "Gene": ["DE1", "NEG_SIG", "N3", "N1", "N2", "N4"],
        "CellType": ["Receiver"] * 6,
        "Reference": ["All_other_cells"] * 6,
        "Rank": [1, 2, 3, 4, 5, 6],
        "LogFC": [1.0, -1.0, 0.2, -0.1, 0.0, 0.3],
        "P_Value": [0.001, 0.002, 0.20, 0.051, 0.10, 0.30],
    })
    selected = _selected_cell_type_genes(table, 0.05, de_gene_target=4)
    assert selected["Gene"].tolist() == ["DE1", "N1", "N2", "N3"]
    assert selected["Selection_Source"].tolist() == [
        "DE", "p_value_supplement", "p_value_supplement", "p_value_supplement",
    ]


def test_de_gene_target_never_truncates_a_larger_de_set():
    table = pd.DataFrame({
        "Gene": ["G1", "G2", "G3", "G4"],
        "CellType": ["Receiver"] * 4,
        "Reference": ["All_other_cells"] * 4,
        "Rank": [1, 2, 3, 4],
        "LogFC": [1.0, 0.9, 0.8, 0.7],
        "P_Value": [0.001, 0.002, 0.003, 0.20],
    })
    selected = _selected_cell_type_genes(table, 0.05, de_gene_target=2)
    assert selected["Gene"].tolist() == ["G1", "G2", "G3"]


def test_de_gene_target_uses_all_genes_when_gene_universe_is_smaller():
    table = pd.DataFrame({
        "Gene": ["G1", "G2", "G3"],
        "CellType": ["Receiver"] * 3,
        "Reference": ["All_other_cells"] * 3,
        "Rank": [1, 2, 3],
        "LogFC": [1.0, -1.0, 0.2],
        "P_Value": [0.001, 0.002, 0.20],
    })
    selected = _selected_cell_type_genes(table, 0.05, de_gene_target=500)
    assert set(selected["Gene"]) == {"G1", "G2", "G3"}
    assert len(selected) == 3
    assert selected.loc[
        selected["Gene"] != "G1", "Selection_Source"
    ].eq("all_genes_below_target").all()
