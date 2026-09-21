"""Tests for prior-network parsing and direction handling."""

import gzip
from pathlib import Path

from main.priors import load_intracellular_prior, load_ligand_receptor_prior


def test_prior_removes_self_edges_and_expands_reversible_edges(tmp_path: Path):
    prior = tmp_path / "prior.tsv"
    prior.write_text(
        "Gene1\tGene2\tType\n"
        "A\tA\tinteracts-with\n"
        "A\tB\tinteracts-with\n"
        "B\tC\tcontrols-state-change-of\n",
        encoding="utf-8",
    )
    result = load_intracellular_prior(prior, {"A", "B", "C"})
    pairs = set(zip(result["source"], result["target"]))
    assert ("A", "A") not in pairs
    assert {("A", "B"), ("B", "A"), ("B", "C")}.issubset(pairs)
    assert ("C", "B") not in pairs


def test_intracellular_prior_reader_accepts_gzip(tmp_path: Path):
    prior = tmp_path / "prior.tsv.gz"
    with gzip.open(prior, "wt", encoding="utf-8") as handle:
        handle.write(
            "Gene1\tGene2\tType\n"
            "A\tB\tinteracts-with\n"
            "B\tC\tcontrols-state-change-of\n"
        )
    result = load_intracellular_prior(prior, {"A", "B", "C"})
    pairs = set(zip(result["source"], result["target"]))
    assert {("A", "B"), ("B", "A"), ("B", "C")}.issubset(pairs)


def test_ligand_receptor_reader_uses_all_rows_and_only_first_two_columns(tmp_path: Path):
    prior = tmp_path / "lr.tsv"
    prior.write_text(
        "Ligand\tReceptor\tScore\tUnused\n"
        "L1\tR1\t1\ta\n"
        "L2\tR2\t99\tb\n"
        "L2\tR2\t1\tc\n",
        encoding="utf-8",
    )
    result = load_ligand_receptor_prior(prior)
    assert result.columns.tolist() == ["Ligand", "Receptor"]
    assert len(result) == 3
    assert set(map(tuple, result.to_numpy())) == {("L1", "R1"), ("L2", "R2")}
