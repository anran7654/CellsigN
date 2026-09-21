"""End-to-end test using the small synthetic fixture in tests/data."""

import json
from pathlib import Path

import pandas as pd

from main.cli import main


def test_packaged_example_produces_a_path(tmp_path: Path):
    root = Path(__file__).resolve().parents[1]
    output_base = tmp_path / "Results"
    output_dir = tmp_path / "Results_expression"
    main([
        "--stage", "evidence",
        "--expression", str(root / "tests" / "data" / "expression.tsv"),
        "--cell-path", str(root / "tests" / "data" / "cell_types.tsv"),
        "--intracellular-prior", str(root / "tests" / "data" / "intracellular_prior.tsv"),
        "--ligand-receptor-prior", str(root / "tests" / "data" / "ligand_receptor.tsv"),
        "--output-dir", str(output_base),
        "--pair", "Sender:Receiver",
        "--alpha", "0.05",
        "--stability-resamples", "0",
    ])
    edge_path = output_dir / "evidence" / "Receiver_edge_evidence.txt"
    evidence_mtime = edge_path.stat().st_mtime_ns
    assert not (output_dir / "Sender_to_Receiver_pathway.txt").exists()

    main([
        "--stage", "paths",
        "--expression", str(root / "tests" / "data" / "expression.tsv"),
        "--intracellular-prior", str(root / "tests" / "data" / "intracellular_prior.tsv"),
        "--ligand-receptor-prior", str(root / "tests" / "data" / "ligand_receptor.tsv"),
        "--output-dir", str(output_base),
        "--pair", "Sender:Receiver",
        "--path-permutations", "5",
    ])
    assert edge_path.stat().st_mtime_ns == evidence_mtime
    output = pd.read_csv(output_dir / "Sender_to_Receiver_pathway.txt", sep="\t")
    assert not output.empty
    assert output.columns[:5].tolist() == ["Ligand", "Receptor", "Mediator", "TF", "Target"]
    assert {"LIG", "REC", "TF", "TG"}.issubset(set(output.iloc[0].astype(str)))
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
    deg_counts = pd.read_csv(
        output_dir / "other" / "differential_gene_counts.txt", sep="\t"
    )
    assert set(deg_counts["CellType"]) == {"Sender", "Receiver"}
    assert deg_counts["CellType"].is_unique
