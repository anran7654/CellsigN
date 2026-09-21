"""Tests for weighted R-M-TF path selection."""

import pandas as pd

from main.pathways import (
    ReceiverPathCache,
    _build_graph,
    attach_sender_paths,
    infer_paths,
    select_path_graph_edges,
)


def _edge(gene1, gene2, cost, prior=True):
    return {
        "Edge_ID": f"{gene1}--{gene2}", "Receiver": "Receiver",
        "Gene1": gene1, "Gene2": gene2, "Edge_Retained": True,
        "Data_Accepted": True,
        "Relation_Type": "positive_linear",
        "Relation_Sign": "", "R_Receiver": 0.0, "P_Corr": 1.0, "Q_Corr": 1.0,
        "R_Background": 0.0, "P_Diff": 1.0, "Q_Diff": 1.0,
        "Specificity_Effect": 0.0, "Specificity_Label": "not_significant",
        "Nonlinear_F_UV": 0.0, "P_Nonlinear_UV": 1.0, "Delta_R2_UV": 0.0,
        "Nonlinear_F_VU": 0.0, "P_Nonlinear_VU": 1.0, "Delta_R2_VU": 0.0,
        "P_Nonlinear_Sym": 1.0, "Q_Nonlinear": 1.0, "Delta_R2": 0.0,
        "Stability": float("nan"), "Prior_Found": prior,
        "Prior_Direction": f"{gene1}->{gene2}" if prior else "",
        "Prior_Type": "known" if prior else "", "Data_Cost": cost,
        "Hub_Penalty": 0.0, "Hop_Penalty": 0.0, "Final_Cost": cost,
    }


def test_weighted_shortest_path_chooses_lowest_cost_route():
    edges = pd.DataFrame([
        _edge("R", "A", 0.1), _edge("A", "TF", 0.1),
        _edge("R", "B", 0.8), _edge("B", "TF", 0.8),
    ])
    prior = pd.DataFrame({
        "source": ["TF"], "target": ["TG"],
        "prior_type": ["controls-expression-of"], "layer": ["tf_target"],
    })
    lr = pd.DataFrame({"Ligand": ["L"], "Receptor": ["R"]})
    result = infer_paths(
        "Sender", "Receiver", {"L"}, {"R", "A", "B", "TF", "TG"},
        edges, prior, lr, path_permutations=0,
    )
    assert result.rmtf_paths.iloc[0]["Mediator"] == "A"
    assert result.pathways.columns[:5].tolist() == ["Ligand", "Receptor", "Mediator", "TF", "Target"]


def test_k_shortest_paths_are_retained_and_ranked_by_cost():
    edges = pd.DataFrame([
        _edge("R", "A", 0.1), _edge("A", "TF", 0.1),
        _edge("R", "B", 0.1), _edge("B", "TF", 0.1),
        _edge("R", "C", 0.1), _edge("C", "D", 0.1),
        _edge("D", "TF", 0.1),
    ])
    prior = pd.DataFrame({
        "source": ["TF"], "target": ["TG"],
        "prior_type": ["controls-expression-of"], "layer": ["tf_target"],
    })
    lr = pd.DataFrame({"Ligand": ["L"], "Receptor": ["R"]})
    result = infer_paths(
        "Sender", "Receiver", {"L"}, {"R", "A", "B", "C", "D", "TF", "TG"},
        edges, prior, lr, path_permutations=0,
    )
    assert set(result.rmtf_paths["Mediator"]) == {"A", "B", "C,D"}
    assert set(result.rmtf_paths["Path_Rank"]) == {1, 2, 3}
    assert set(result.rmtf_paths["Alternative_Path_Count"]) == {3}
    assert len(result.subnetwork_edges) == 7


def test_shortest_path_has_no_maximum_length_limit():
    mediators = [f"M{index}" for index in range(1, 12)]
    chain = ["R", *mediators, "TF"]
    edges = pd.DataFrame([
        _edge(source, target, 0.1)
        for source, target in zip(chain[:-1], chain[1:])
    ])
    prior = pd.DataFrame({
        "source": ["TF"], "target": ["TG"],
        "prior_type": ["controls-expression-of"], "layer": ["tf_target"],
    })
    lr = pd.DataFrame({"Ligand": ["L"], "Receptor": ["R"]})
    result = infer_paths(
        "Sender", "Receiver", {"L"}, set(chain) | {"TG"},
        edges, prior, lr, path_permutations=0,
    )
    assert result.rmtf_paths.iloc[0]["Mediator"] == ",".join(mediators)
    assert result.rmtf_paths.iloc[0]["Path_Length"] == 12


def test_receptor_that_is_also_a_tf_excludes_self_but_reaches_other_tf():
    edges = pd.DataFrame([
        _edge("R", "M", 0.1), _edge("M", "TF", 0.1),
    ])
    prior = pd.DataFrame({
        "source": ["R", "TF"], "target": ["TG1", "TG2"],
        "prior_type": ["controls-expression-of", "controls-expression-of"],
        "layer": ["tf_target", "tf_target"],
    })
    lr = pd.DataFrame({"Ligand": ["L"], "Receptor": ["R"]})
    result = infer_paths(
        "Sender", "Receiver", {"L"}, {"R", "M", "TF", "TG1", "TG2"},
        edges, prior, lr, path_permutations=0,
    )
    assert len(result.rmtf_paths) == 1
    assert result.rmtf_paths.iloc[0]["Receptor"] == "R"
    assert result.rmtf_paths.iloc[0]["TF"] == "TF"
    assert result.rmtf_paths.iloc[0]["Mediator"] == "M"


def test_direct_receptor_to_tf_path_is_allowed():
    edges = pd.DataFrame([_edge("R", "TF", 0.1)])
    prior = pd.DataFrame({
        "source": ["TF"], "target": ["TG"],
        "prior_type": ["controls-expression-of"], "layer": ["tf_target"],
    })
    lr = pd.DataFrame({"Ligand": ["L"], "Receptor": ["R"]})
    result = infer_paths(
        "Sender", "Receiver", {"L"}, {"R", "TF", "TG"},
        edges, prior, lr, path_permutations=0,
    )
    assert len(result.rmtf_paths) == 1
    assert pd.isna(result.rmtf_paths.iloc[0]["Mediator"]) or result.rmtf_paths.iloc[0]["Mediator"] == ""
    assert result.rmtf_paths.iloc[0]["Path_Length"] == 1


def test_number_of_paths_is_controlled_by_k_paths():
    mediators = [f"M{index:03d}" for index in range(101)]
    edges = pd.DataFrame([
        edge
        for mediator in mediators
        for edge in (_edge("R", mediator, 0.1), _edge(mediator, "TF", 0.1))
    ])
    prior = pd.DataFrame({
        "source": ["TF"], "target": ["TG"],
        "prior_type": ["controls-expression-of"], "layer": ["tf_target"],
    })
    lr = pd.DataFrame({"Ligand": ["L"], "Receptor": ["R"]})
    result = infer_paths(
        "Sender", "Receiver", {"L"}, {"R", *mediators, "TF", "TG"},
        edges, prior, lr, path_permutations=0, data_edges_per_node=0, k_paths=7,
    )
    assert len(result.rmtf_paths) == 7
    assert set(result.rmtf_paths["Alternative_Path_Count"]) == {7}
    assert set(result.rmtf_paths["Path_Rank"]) == set(range(1, 8))


def test_path_graph_keeps_all_retained_edges_with_or_without_prior_annotation():
    edges = pd.DataFrame([
        _edge("A", "B", 0.1, prior=False),
        _edge("A", "C", 0.2, prior=False),
        _edge("B", "C", 0.3, prior=False),
        _edge("D", "E", 0.9, prior=True),
    ])
    selected = select_path_graph_edges(edges, data_edges_per_node=1)
    assert set(selected["Edge_ID"]) == {"A--B", "A--C", "B--C", "D--E"}


def test_receptor_omnibus_statistics_are_reported_for_every_tf_path():
    edges = pd.DataFrame([
        _edge("R", "A", 0.1), _edge("A", "TF1", 0.1),
        _edge("R", "B", 0.4), _edge("B", "TF2", 0.4),
    ])
    prior = pd.DataFrame({
        "source": ["TF1", "TF2"], "target": ["TG1", "TG2"],
        "prior_type": ["controls-expression-of", "controls-expression-of"],
        "layer": ["tf_target", "tf_target"],
    })
    lr = pd.DataFrame({"Ligand": ["L"], "Receptor": ["R"]})
    result = infer_paths(
        "Sender", "Receiver", {"L"},
        {"R", "A", "B", "TF1", "TF2", "TG1", "TG2"},
        edges, prior, lr, path_permutations=9, seed=11,
    )
    assert len(result.rmtf_paths) == 2
    assert {
        "Receptor_P", "Receptor_Q", "Receptor_Significant"
    }.issubset(result.rmtf_paths.columns)
    assert result.rmtf_paths["Receptor_P"].nunique() == 1
    assert result.rmtf_paths["Receptor_Q"].nunique() == 1
    assert result.metadata["tested_receptors"] == 1


def test_receptor_bh_family_is_distinct_from_rtf_pair_family():
    edges = pd.DataFrame([
        _edge("R1", "TF1", 0.1),
        _edge("R1", "TF2", 0.2),
        _edge("R2", "TF1", 0.3),
    ])
    graph = _build_graph(edges)
    cache = ReceiverPathCache(
        receiver="Receiver",
        graph=graph,
        tf_target=pd.DataFrame({
            "source": ["TF1", "TF2"],
            "target": ["TG1", "TG2"],
        }),
        pair_paths={
            ("R1", "TF1"): [["R1", "TF1"]],
            ("R1", "TF2"): [["R1", "TF2"]],
            ("R2", "TF1"): [["R2", "TF1"]],
        },
        pair_p={("R1", "TF1"): 0.01, ("R1", "TF2"): 0.02, ("R2", "TF1"): 0.04},
        receptor_p={"R1": 0.01, "R2": 0.04},
        retained_input_edges=3,
        path_graph_edges=3,
        path_graph_unannotated_edges=0,
        path_graph_prior_annotated_edges=3,
    )
    result = attach_sender_paths(
        "Sender", "Receiver", {"L1", "L2"},
        {"R1", "R2", "TF1", "TF2", "TG1", "TG2"},
        cache,
        pd.DataFrame({"Ligand": ["L1", "L2"], "Receptor": ["R1", "R2"]}),
        path_alpha=0.03,
    )
    receptor_rows = result.rmtf_paths.groupby("Receptor", sort=True).first()
    assert receptor_rows.loc["R1", "Receptor_Q"] == 0.02
    assert receptor_rows.loc["R2", "Receptor_Q"] == 0.04
    assert bool(receptor_rows.loc["R1", "Receptor_Significant"])
    assert not bool(receptor_rows.loc["R2", "Receptor_Significant"])
    assert result.metadata["tested_receptors"] == 2
    assert result.metadata["significant_receptors"] == 1
