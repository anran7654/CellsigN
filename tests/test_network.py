"""Tests for receiver-only data-driven network construction and prior overlay."""

import numpy as np
import pandas as pd

from main.network import build_cell_network


def test_network_tests_all_receiver_gene_pairs_before_prior_annotation():
    rng = np.random.default_rng(3)
    gene_a = np.linspace(0, 2, 12)
    receiver = np.column_stack((
        gene_a,
        2 * gene_a + rng.normal(0, 0.01, 12),
        np.zeros(12),
        rng.normal(size=12),
    ))
    genes = ["A", "B", "C", "D"]
    ranking = pd.DataFrame({
        "Gene": genes,
        "DE_Rank": range(1, 5), "DE_LogFC": 1.0,
        "P_Value": 0.001,
    })
    prior = pd.DataFrame({
        "source": ["A", "C"], "target": ["B", "D"],
        "prior_type": ["controls-state-change-of", "controls-state-change-of"],
        "layer": ["signaling", "signaling"],
    })
    result = build_cell_network(
        "Receiver", receiver, genes, ranking, prior,
        min_abs_r=0.5, stability_resamples=0,
    )
    assert len(result.edge_evidence) == 6
    ab = result.edge_evidence[result.edge_evidence["Edge_ID"] == "A--B"].iloc[0]
    ac = result.edge_evidence[result.edge_evidence["Edge_ID"] == "A--C"].iloc[0]
    cd = result.edge_evidence[result.edge_evidence["Edge_ID"] == "C--D"].iloc[0]
    assert ab["Edge_Retained"] and ab["Prior_Found"]
    assert not ac["Prior_Found"]
    assert cd["Prior_Found"] and not cd["Edge_Retained"]
    assert "Node_Score" not in result.gene_evidence.columns
    assert "Edge_Source" not in result.edge_evidence.columns
