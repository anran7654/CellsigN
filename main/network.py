"""Build receiver-specific, data-driven association networks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .statistics import (
    benjamini_hochberg,
    cubic_increment_effect,
    cubic_nonlinear_receiver_pair_tests,
    pearson_pairs,
    q_to_score,
)


GENE_EVIDENCE_COLUMNS = [
    "Receiver", "Gene", "DE_Rank", "DE_LogFC", "P_Value",
]

EDGE_EVIDENCE_COLUMNS = [
    "Edge_ID", "Receiver", "Gene1", "Gene2", "Edge_Retained",
    "Relation_Type", "Relation_Sign", "R_Receiver",
    "P_Corr", "Q_Corr", "Nonlinear_F_UV", "P_Nonlinear_UV",
    "Delta_R2_UV", "Nonlinear_F_VU", "P_Nonlinear_VU", "Delta_R2_VU",
    "P_Nonlinear_Sym", "Q_Nonlinear", "Delta_R2", "Stability",
    "Prior_Found", "Prior_Direction", "Prior_Type", "Data_Cost",
    "Hub_Penalty", "Hop_Penalty", "Final_Cost",
]


@dataclass(frozen=True)
class CellNetworkResult:
    """Auditable gene and edge evidence for one receiver cell type."""

    gene_evidence: pd.DataFrame
    edge_evidence: pd.DataFrame


def _all_gene_pairs_with_prior(
    genes: list[str], prior_edges: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return every unordered gene pair with optional prior annotations."""
    source, target = np.triu_indices(len(genes), k=1)
    positions = {gene: index for index, gene in enumerate(genes)}
    annotations: dict[tuple[int, int], dict[str, set[str]]] = {}
    signaling = prior_edges[
        (prior_edges["layer"] == "signaling")
        & prior_edges["source"].isin(positions)
        & prior_edges["target"].isin(positions)
    ]
    for row in signaling.itertuples(index=False):
        a, b = positions[str(row.source)], positions[str(row.target)]
        if a == b:
            continue
        pair = (a, b) if a < b else (b, a)
        item = annotations.setdefault(pair, {"directions": set(), "types": set()})
        item["directions"].add(f"{row.source}->{row.target}")
        item["types"].update(
            token for token in str(row.prior_type).split(";") if token
        )
    found = np.asarray([
        (int(u), int(v)) in annotations for u, v in zip(source, target)
    ], dtype=bool)
    directions = np.asarray([
        ";".join(sorted(annotations.get((int(u), int(v)), {}).get("directions", set())))
        for u, v in zip(source, target)
    ], dtype=object)
    types = np.asarray([
        ";".join(sorted(annotations.get((int(u), int(v)), {}).get("types", set())))
        for u, v in zip(source, target)
    ], dtype=object)
    return source.astype(int), target.astype(int), found, directions, types


def _bootstrap_stability(
    receiver: np.ndarray,
    source_index: np.ndarray,
    target_index: np.ndarray,
    relation_type: np.ndarray,
    relation_sign: np.ndarray,
    min_abs_r: float,
    min_delta_r2: float,
    resamples: int,
    seed: int,
) -> np.ndarray:
    stability = np.full(source_index.size, np.nan, dtype=float)
    eligible = np.flatnonzero(np.isin(
        relation_type, ["positive_linear", "negative_linear", "nonlinear"]
    ))
    if resamples <= 0 or eligible.size == 0:
        return stability
    rng = np.random.default_rng(seed)
    successes = np.zeros(source_index.size, dtype=float)
    n_cells = receiver.shape[0]
    for _ in range(resamples):
        sampled = receiver[rng.integers(0, n_cells, size=n_cells)]
        linear = eligible[np.isin(
            relation_type[eligible], ["positive_linear", "negative_linear"]
        )]
        if linear.size:
            values, _ = pearson_pairs(
                sampled, source_index[linear], target_index[linear]
            )
            expected = np.where(relation_sign[linear] == "positive", 1.0, -1.0)
            successes[linear] += ((values * expected) > 0) & (
                np.abs(values) >= min_abs_r
            )
        nonlinear = eligible[relation_type[eligible] == "nonlinear"]
        for pair_index in nonlinear:
            u, v = source_index[pair_index], target_index[pair_index]
            delta = max(
                cubic_increment_effect(sampled[:, u], sampled[:, v]),
                cubic_increment_effect(sampled[:, v], sampled[:, u]),
            )
            successes[pair_index] += delta >= min_delta_r2
    stability[eligible] = successes[eligible] / resamples
    return stability


def build_cell_network(
    receiver: str,
    receiver_expression: np.ndarray,
    genes: list[str],
    receiver_gene_table: pd.DataFrame,
    prior_edges: pd.DataFrame,
    alpha: float = 0.05,
    min_abs_r: float = 0.1,
    min_delta_r2: float = 0.01,
    stability_resamples: int = 20,
    relation_weight: float = 0.6,
    stability_weight: float = 0.2,
    hub_penalty_weight: float = 0.1,
    hop_penalty: float = 0.05,
    seed: int = 0,
) -> CellNetworkResult:
    """Test all unordered DEG pairs within one receiver cell type.

    The statistical network is built before prior knowledge is overlaid. A
    Only retained edges found in the intracellular signaling prior enter path
    inference, where they follow their recorded prior directions. Unannotated
    retained edges remain in the statistical evidence table.
    """
    receiver_expression = np.asarray(receiver_expression, dtype=float)
    if receiver_expression.ndim != 2:
        raise ValueError("Receiver expression must be a cells-by-genes matrix.")
    if receiver_expression.shape[1] != len(genes):
        raise ValueError("Gene names and receiver expression columns do not align.")
    if receiver_expression.shape[0] < 5:
        raise ValueError("The receiver requires at least five cells.")
    if len(genes) < 2:
        raise ValueError("The receiver has fewer than two selected one-vs-rest DE genes.")
    if min(relation_weight, stability_weight, hub_penalty_weight, hop_penalty) < 0:
        raise ValueError("Cost weights and penalties must be nonnegative.")
    if relation_weight <= 0:
        raise ValueError("--relation-weight must be positive.")

    ranking = receiver_gene_table.set_index("Gene").reindex(genes)
    gene_evidence = pd.DataFrame({
        "Receiver": receiver,
        "Gene": genes,
        "DE_Rank": ranking["DE_Rank"].to_numpy(),
        "DE_LogFC": pd.to_numeric(ranking["DE_LogFC"], errors="coerce").to_numpy(),
        "P_Value": ranking["P_Value"].to_numpy(),
    }, columns=GENE_EVIDENCE_COLUMNS)

    u_index, v_index, prior_found, prior_direction, prior_type = (
        _all_gene_pairs_with_prior(genes, prior_edges)
    )
    gene_array = np.asarray(genes, dtype=object)
    gene1, gene2 = gene_array[u_index], gene_array[v_index]
    r_receiver, p_corr = pearson_pairs(receiver_expression, u_index, v_index)
    q_corr = benjamini_hochberg(p_corr)
    f_uv, p_uv, delta_uv, f_vu, p_vu, delta_vu = (
        cubic_nonlinear_receiver_pair_tests(
            receiver_expression, u_index, v_index
        )
    )
    p_nonlinear_sym = np.minimum(1.0, 2.0 * np.fmin(p_uv, p_vu))
    p_nonlinear_sym[~(np.isfinite(p_uv) | np.isfinite(p_vu))] = np.nan
    q_nonlinear = benjamini_hochberg(p_nonlinear_sym)
    delta_r2 = np.fmax(delta_uv, delta_vu)

    linear_pass = np.isfinite(q_corr) & (q_corr < alpha) & (
        np.abs(r_receiver) >= min_abs_r
    )
    nonlinear_pass = (
        np.isfinite(q_nonlinear) & (q_nonlinear < alpha)
        & np.isfinite(delta_r2) & (delta_r2 >= min_delta_r2)
    )
    data_accepted = linear_pass | nonlinear_pass
    relation_type = np.full(u_index.size, "not_supported", dtype=object)
    relation_sign = np.full(u_index.size, "", dtype=object)
    positive = linear_pass & (r_receiver > 0)
    negative = linear_pass & (r_receiver < 0)
    relation_type[positive] = "positive_linear"
    relation_type[negative] = "negative_linear"
    relation_sign[positive] = "positive"
    relation_sign[negative] = "negative"
    nonlinear_only = nonlinear_pass & ~linear_pass
    relation_type[nonlinear_only] = "nonlinear"
    relation_sign[nonlinear_only] = "nonlinear"

    stability = _bootstrap_stability(
        receiver_expression, u_index, v_index, relation_type, relation_sign,
        min_abs_r, min_delta_r2, stability_resamples, seed,
    )
    linear_score = (q_to_score(q_corr) + np.abs(r_receiver)) / 2.0
    nonlinear_score = (
        q_to_score(q_nonlinear) + np.nan_to_num(delta_r2, nan=0.0)
    ) / 2.0
    relation_score = np.where(
        relation_type == "nonlinear", nonlinear_score, linear_score
    )
    relation_score = np.clip(np.nan_to_num(relation_score), 0.0, 1.0)

    data_cost = np.empty(u_index.size, dtype=float)
    for pair_index in range(u_index.size):
        scores = [relation_score[pair_index]]
        weights = [relation_weight]
        if np.isfinite(stability[pair_index]) and stability_weight > 0:
            scores.append(stability[pair_index])
            weights.append(stability_weight)
        combined = float(np.dot(scores, weights) / sum(weights))
        data_cost[pair_index] = 1.0 - np.clip(combined, 0.0, 1.0)

    degrees = np.zeros(len(genes), dtype=float)
    np.add.at(degrees, u_index[data_accepted], 1)
    np.add.at(degrees, v_index[data_accepted], 1)
    max_degree = degrees.max(initial=0.0)
    normalized_degree = (
        np.log1p(degrees) / np.log1p(max_degree)
        if max_degree > 1 else degrees
    )
    hub_penalty = (
        normalized_degree[u_index] + normalized_degree[v_index]
    ) / 2.0
    final_cost = data_cost + hub_penalty_weight * hub_penalty + hop_penalty
    edge_id = np.asarray([
        f"{min(str(u), str(v))}--{max(str(u), str(v))}"
        for u, v in zip(gene1, gene2)
    ], dtype=object)

    edge_evidence = pd.DataFrame({
        "Edge_ID": edge_id,
        "Receiver": receiver,
        "Gene1": gene1,
        "Gene2": gene2,
        "Edge_Retained": data_accepted,
        "Relation_Type": relation_type,
        "Relation_Sign": relation_sign,
        "R_Receiver": r_receiver,
        "P_Corr": p_corr,
        "Q_Corr": q_corr,
        "Nonlinear_F_UV": f_uv,
        "P_Nonlinear_UV": p_uv,
        "Delta_R2_UV": delta_uv,
        "Nonlinear_F_VU": f_vu,
        "P_Nonlinear_VU": p_vu,
        "Delta_R2_VU": delta_vu,
        "P_Nonlinear_Sym": p_nonlinear_sym,
        "Q_Nonlinear": q_nonlinear,
        "Delta_R2": delta_r2,
        "Stability": stability,
        "Prior_Found": prior_found,
        "Prior_Direction": prior_direction,
        "Prior_Type": prior_type,
        "Data_Cost": data_cost,
        "Hub_Penalty": hub_penalty,
        "Hop_Penalty": hop_penalty,
        "Final_Cost": final_cost,
    }, columns=EDGE_EVIDENCE_COLUMNS)
    return CellNetworkResult(gene_evidence, edge_evidence)
