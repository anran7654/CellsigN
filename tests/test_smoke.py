"""End-to-end test for the staged packaged example."""

import json
from pathlib import Path

import pandas as pd

from main.cli import main


def test_packaged_example_produces_a_path(tmp_path: Path):
    fixture_dir = Path(__file__).resolve().parent / "data"
    output_base = tmp_path / "Results"
    output_dir = tmp_path / "Results_expression"
    main([
        "--stage", "evidence",
        "--expression", str(fixture_dir / "expression.tsv"),
        "--cell-path", str(fixture_dir / "cell_types.tsv"),
        "--intracellular-prior", str(fixture_dir / "intracellular_prior.tsv"),
        "--ligand-receptor-prior", str(fixture_dir / "ligand_receptor.tsv"),
        "--output-dir", str(output_base),
        "--pair", "Sender:Receiver",
        "--alpha", "0.05",
        "--stability-resamples", "0",
    ])
    edge_path = output_dir / "evidence" / "Receiver_edge_evidence.txt"
    evidence_mtime = edge_path.stat().st_mtime_ns
    assert not (output_dir / "Sender_to_Receiver_pathway.txt").exists()

    # A saved evidence directory from the previous schema remains reusable.
    counts_path = output_dir / "other" / "test_counts.txt"
    legacy_counts = pd.read_csv(counts_path, sep="\t")
    legacy_counts["Path_BH_Method"] = "BH across R-TF pairs"
    legacy_counts["Receptor_BH_Method"] = "BH across receptors"
    legacy_counts.to_csv(counts_path, sep="\t", index=False)
    manifest_path = output_dir / "other" / "run_manifest.json"
    legacy_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    legacy_manifest["multiple_testing"]["path_pair_diagnostic"] = "old pair BH"
    legacy_manifest["multiple_testing"]["receptor_omnibus"] = "old receptor BH"
    legacy_manifest["path_parameters"] = {"path_alpha": 0.05}
    manifest_path.write_text(json.dumps(legacy_manifest), encoding="utf-8")

    main([
        "--stage", "paths",
        "--expression", str(fixture_dir / "expression.tsv"),
        "--intracellular-prior", str(fixture_dir / "intracellular_prior.tsv"),
        "--ligand-receptor-prior", str(fixture_dir / "ligand_receptor.tsv"),
        "--output-dir", str(output_base),
        "--pair", "Sender:Receiver",
        "--path-permutations", "5",
    ])
    assert edge_path.stat().st_mtime_ns == evidence_mtime
    output = pd.read_csv(output_dir / "Sender_to_Receiver_pathway.txt", sep="\t")
    assert not output.empty
    assert output.columns[:5].tolist() == ["Ligand", "Receptor", "Mediator", "TF", "Target"]
    assert {"LIG", "REC", "TF", "TG"}.issubset(set(output.iloc[0].astype(str)))
    assert "Path_P" in output.columns
    obsolete_columns = {
        "Path_Q", "Receptor_P", "Receptor_Q", "Receptor_Significant",
    }
    assert obsolete_columns.isdisjoint(output.columns)
    rmtf = pd.read_csv(
        output_dir / "evidence" / "Sender_to_Receiver_rmtf_paths.txt", sep="\t"
    )
    assert "Path_P" in rmtf.columns
    assert obsolete_columns.isdisjoint(rmtf.columns)
    assert rmtf["Path_Rank"].eq(1).all()
    assert rmtf["Alternative_Path_Count"].eq(1).all()
    assert (output_dir / "other" / "processed_gene_list.txt").is_file()
    assert (output_dir / "other" / "differential_gene_counts.txt").is_file()
    assert (output_dir / "other" / "test_counts.txt").is_file()
    assert (output_dir / "other" / "run_manifest.json").is_file()
    assert (output_dir / "gene_lists" / "Sender_de_genes.txt").is_file()
    assert (output_dir / "gene_lists" / "Receiver_de_genes.txt").is_file()
    assert (output_dir / "evidence" / "Sender_one_vs_rest_de.txt").is_file()
    assert (output_dir / "evidence" / "Receiver_one_vs_rest_de.txt").is_file()
    assert (output_dir / "evidence" / "Receiver_gene_evidence.txt").is_file()
    assert not list((output_dir / "gene_lists").glob("*_sender_de_genes.txt"))
    assert not list(output_dir.rglob("*.graphml"))
    assert not (output_dir / "all_pathways.txt").exists()
    manifest = json.loads((output_dir / "other" / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["ligand_receptor_reading"].endswith("no score filter")
    assert manifest["path_inputs"]["ligand_receptor_prior"]["columns_used"] == "first two only"
    assert manifest["path_inputs"]["ligand_receptor_prior"]["score_filter"] is None
    assert manifest["multiple_testing"]["one_vs_rest_de"].startswith("None")
    assert "all unordered gene pairs" in manifest["candidate_edge_policy"]
    assert set(manifest["completed_stages"]) == {"evidence", "paths"}
    assert manifest["path_parameters"]["k_paths"] == 1
    assert "path_alpha" not in manifest["path_parameters"]
    assert {
        "path_pair_diagnostic", "receptor_omnibus",
    }.isdisjoint(manifest["multiple_testing"])
    for direction_counts in manifest["direction_counts"].values():
        assert {
            "tested_receptors", "significant_receptors", "receptor_fdr_threshold",
        }.isdisjoint(direction_counts)
    counts = pd.read_csv(counts_path, sep="\t")
    assert {"Path_BH_Method", "Receptor_BH_Method"}.isdisjoint(counts.columns)
    assert "Pair_BH_Method" in counts.columns
    evidence = pd.read_csv(edge_path, sep="\t")
    assert {"P_Corr", "Q_Corr", "Q_Nonlinear", "Stability"}.issubset(evidence.columns)
    deg_counts = pd.read_csv(
        output_dir / "other" / "differential_gene_counts.txt", sep="\t"
    )
    assert set(deg_counts["CellType"]) == {"Sender", "Receiver"}
    assert deg_counts["CellType"].is_unique
