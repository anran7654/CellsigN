"""Sparse, cached weighted R-M-TF path inference and endpoint attachment."""

from __future__ import annotations

from dataclasses import dataclass, field
import heapq
from itertools import count as counter
from typing import Callable

import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

from .statistics import benjamini_hochberg


MAIN_COLUMNS = [
    "Ligand", "Receptor", "Mediator", "TF", "Target", "Path_ID", "Path_Rank",
    "Path_Length", "Path_Cost", "Path_P", "Path_Q", "Receptor_P",
    "Receptor_Q", "Receptor_Significant", "Edge_Types", "Direction_Status",
]

RMTF_COLUMNS = [
    "Path_ID", "Sender", "Receiver", "Receptor", "Mediator", "TF", "Path_Rank",
    "Path_Length", "Path_Cost", "Path_P", "Path_Q", "Receptor_P",
    "Receptor_Q", "Receptor_Significant", "Alternative_Path_Count",
    "Min_Stability", "Direction_Status",
]

PATH_EDGE_COLUMNS = [
    "Path_ID", "Path_Rank", "Edge_Order", "Gene1", "Gene2", "Relation_Type",
    "Relation_Sign", "Prior_Found", "Prior_Direction",
    "Stability", "Data_Cost", "Hub_Penalty", "Hop_Penalty", "Final_Cost",
]

SUBNETWORK_EDGE_COLUMNS = [
    "Sender", "Receiver", "Gene1", "Gene2", "Best_Path_Rank",
    "Path_Count", "Relation_Type", "Relation_Sign", "Prior_Direction",
    "Min_Final_Cost",
]

@dataclass(frozen=True)
class PathResult:
    """Path, edge, merged-subnetwork, and metadata outputs for one direction."""

    pathways: pd.DataFrame
    rmtf_paths: pd.DataFrame
    path_edges: pd.DataFrame
    subnetwork_edges: pd.DataFrame
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass
class ReceiverPathCache:
    """Reusable R-TF paths calculated once for one receiver network."""

    receiver: str
    graph: nx.DiGraph
    tf_target: pd.DataFrame
    pair_paths: dict[tuple[str, str], list[list[str]]]
    pair_p: dict[tuple[str, str], float]
    receptor_p: dict[str, float]
    retained_input_edges: int
    path_graph_edges: int
    path_graph_unannotated_edges: int
    path_graph_prior_annotated_edges: int


def _empty_result(metadata: dict[str, object] | None = None) -> PathResult:
    return PathResult(
        pd.DataFrame(columns=MAIN_COLUMNS),
        pd.DataFrame(columns=RMTF_COLUMNS),
        pd.DataFrame(columns=PATH_EDGE_COLUMNS),
        pd.DataFrame(columns=SUBNETWORK_EDGE_COLUMNS),
        metadata or {},
    )


def select_path_graph_edges(
    edge_evidence: pd.DataFrame,
    data_edges_per_node: int,
) -> pd.DataFrame:
    """Return all retained receiver-association edges.

    ``data_edges_per_node`` is accepted only for command-line compatibility
    with 0.3.3 and has no effect.
    """
    retained = edge_evidence[edge_evidence["Edge_Retained"]].copy()
    return retained.drop_duplicates("Edge_ID").reset_index(drop=True)


def _build_graph(edge_table: pd.DataFrame) -> nx.DiGraph:
    graph = nx.DiGraph()
    for row in edge_table.itertuples(index=False):
        attributes = row._asdict()
        attributes["weight"] = float(row.Final_Cost)
        arcs: set[tuple[str, str]] = set()
        if row.Prior_Found:
            for token in str(row.Prior_Direction).split(";"):
                if "->" in token:
                    source, target = token.split("->", 1)
                    arcs.add((source, target))
        if not row.Prior_Found:
            arcs.add((str(row.Gene1), str(row.Gene2)))
            arcs.add((str(row.Gene2), str(row.Gene1)))
        if row.Prior_Found and not arcs:
            arcs.add((str(row.Gene1), str(row.Gene2)))
            arcs.add((str(row.Gene2), str(row.Gene1)))
        for source, target in arcs:
            if graph.has_edge(source, target):
                if attributes["weight"] < graph[source][target]["weight"]:
                    graph[source][target].update(attributes)
            else:
                graph.add_edge(source, target, **attributes)
    return graph


def _receiver_roles(
    receiver_genes: set[str],
    prior_edges: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str], set[str]]:
    tf_target = prior_edges[
        (prior_edges["layer"] == "tf_target")
        & prior_edges["source"].isin(receiver_genes)
        & prior_edges["target"].isin(receiver_genes)
    ].copy()
    tf_candidates = set(tf_target["source"].astype(str))
    signaling = prior_edges[prior_edges["layer"] == "signaling"]
    signaling_nodes = set(signaling["source"].astype(str)) | set(
        signaling["target"].astype(str)
    )
    terminal_only = (
        set(tf_target["target"].astype(str)) - signaling_nodes - tf_candidates
    )
    return tf_target, tf_candidates, terminal_only


def _shortest_path_predecessors(
    graph: nx.DiGraph,
    receptor: str,
    targets: set[str],
    all_receptors: set[str],
    all_tfs: set[str],
    terminal_only: set[str],
) -> tuple[dict[str, list[str]], dict[str, float], int]:
    """Run one unrestricted Dijkstra search and retain every tied predecessor."""
    if receptor not in graph or not targets:
        return {}, {}, 0
    order = counter()
    queue: list[tuple[float, int, str]] = [(0.0, next(order), receptor)]
    distances: dict[str, float] = {receptor: 0.0}
    predecessors: dict[str, list[str]] = {receptor: []}
    expansions = 0
    while queue:
        cost, _, node = heapq.heappop(queue)
        if cost > distances.get(node, np.inf):
            continue
        if node in all_tfs and node != receptor:
            continue
        expansions += 1
        for neighbor in sorted(graph.successors(node)):
            if neighbor in terminal_only:
                continue
            if (
                neighbor in all_receptors and neighbor != receptor
                and neighbor not in targets
            ):
                continue
            weight = float(graph[node][neighbor]["weight"])
            if weight < 0:
                raise ValueError("Dijkstra path costs must be nonnegative.")
            next_cost = cost + float(graph[node][neighbor]["weight"])
            previous = distances.get(neighbor, np.inf)
            if next_cost < previous:
                distances[neighbor] = next_cost
                predecessors[neighbor] = [node]
                heapq.heappush(queue, (next_cost, next(order), neighbor))
            elif next_cost == previous and node not in predecessors.get(neighbor, []):
                predecessors.setdefault(neighbor, []).append(node)
    return predecessors, distances, expansions


def _k_shortest_paths(
    graph: nx.DiGraph,
    receptor: str,
    target: str,
    all_receptors: set[str],
    all_tfs: set[str],
    terminal_only: set[str],
    k_paths: int,
) -> list[list[str]]:
    """Return up to k lowest-cost simple paths with valid R-M-TF roles."""
    if receptor == target or receptor not in graph or target not in graph:
        return []
    blocked = set(terminal_only)
    blocked.update(all_receptors - {receptor, target})
    blocked.update(all_tfs - {receptor, target})
    view = nx.subgraph_view(graph, filter_node=lambda node: node not in blocked)
    try:
        generator = nx.shortest_simple_paths(
            view, receptor, target, weight="weight"
        )
        paths: list[list[str]] = []
        for path in generator:
            paths.append(list(path))
            if len(paths) >= k_paths:
                break
        return paths
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return []


def _grouped_shortest_costs(
    graph: nx.DiGraph,
    receptor: str,
    targets: set[str],
    all_receptors: set[str],
    all_tfs: set[str],
    terminal_only: set[str],
) -> dict[str, float]:
    """Find unrestricted minimum costs to multiple TFs with one Dijkstra run."""
    valid_targets = set(targets) - {receptor}
    _, distances, _ = _shortest_path_predecessors(
        graph, receptor, valid_targets, all_receptors, all_tfs, terminal_only
    )
    return {
        target: distances[target]
        for target in valid_targets
        if target in distances
    }


def _path_cost(graph: nx.DiGraph, path: list[str]) -> float:
    return float(sum(graph[u][v]["weight"] for u, v in zip(path[:-1], path[1:])))


def _direction_status(attributes: list[dict[str, object]]) -> str:
    prior_flags = [bool(item["Prior_Found"]) for item in attributes]
    if prior_flags and not any(prior_flags):
        return "unresolved"
    if prior_flags and not all(prior_flags):
        return "mixed"
    return "known"


def _permutation_p_values(
    graph: nx.DiGraph,
    pair_paths: dict[tuple[str, str], list[list[str]]],
    receptors: set[str],
    all_tfs: set[str],
    terminal_only: set[str],
    permutations: int,
    seed: int,
    progress: Callable[[str], None] | None,
) -> tuple[dict[tuple[str, str], float], dict[str, float]]:
    """Return pair-level and receptor-level fixed-topology permutation p-values.

    The receptor statistic is the minimum path cost from that receptor to any
    reachable candidate TF. The same minimum-over-TFs selection is repeated in
    every permutation, so the receptor p-value accounts for choosing its best
    TF endpoint.
    """
    pairs = sorted(pair_paths)
    receptor_ids = sorted({receptor for receptor, _ in pairs})
    if permutations <= 0 or not pairs:
        return (
            {pair: np.nan for pair in pairs},
            {receptor: np.nan for receptor in receptor_ids},
        )
    observed = {
        pair: _path_cost(graph, pair_paths[pair][0]) for pair in pairs
    }
    edge_ids = sorted({str(data["Edge_ID"]) for _, _, data in graph.edges(data=True)})
    original_cost = {
        str(data["Edge_ID"]): float(data["weight"])
        for _, _, data in graph.edges(data=True)
    }
    costs = np.asarray([original_cost[edge_id] for edge_id in edge_ids], dtype=float)
    # A receptor can leave its own node, but no other receptor or TF may be
    # used as an internal node. Give each source receptor a private entrance
    # whose outgoing arcs share the original edge's permuted cost. In the
    # common graph, receptor and TF nodes are endpoints (no outgoing arcs).
    nodes = sorted(graph.nodes)
    node_index = {node: index for index, node in enumerate(nodes)}
    edge_index = {edge_id: index for index, edge_id in enumerate(edge_ids)}
    source_index = {
        receptor: len(nodes) + index
        for index, receptor in enumerate(receptor_ids)
    }
    receptor_index = {
        receptor: index for index, receptor in enumerate(receptor_ids)
    }
    arcs: dict[tuple[int, int], int] = {}
    for source, target, data in graph.edges(data=True):
        if target in terminal_only:
            continue
        cost_index = edge_index[str(data["Edge_ID"])]
        if source not in receptors and source not in all_tfs:
            arcs[(node_index[source], node_index[target])] = cost_index
        if source in source_index:
            arcs[(source_index[source], node_index[target])] = cost_index
    arc_keys = sorted(arcs)
    rows = np.fromiter((row for row, _ in arc_keys), dtype=np.intp)
    columns = np.fromiter((column for _, column in arc_keys), dtype=np.intp)
    arc_cost_index = np.fromiter(
        (arcs[key] for key in arc_keys), dtype=np.intp
    )
    size = len(nodes) + len(receptor_ids)
    sparse_graph = csr_matrix(
        (costs[arc_cost_index], (rows, columns)), shape=(size, size)
    )
    # Arcs were sorted by row and column; the CSR data order is identical.
    pair_sources = np.asarray(
        [receptor_index[receptor] for receptor, _ in pairs], dtype=np.intp
    )
    pair_targets = np.asarray(
        [node_index[tf] for _, tf in pairs], dtype=np.intp
    )
    observed_costs = np.asarray([observed[pair] for pair in pairs])
    receptor_starts = np.unique(pair_sources, return_index=True)[1]
    observed_receptor = np.minimum.reduceat(observed_costs, receptor_starts)
    extreme = np.zeros(len(pairs), dtype=np.int64)
    receptor_extreme = np.zeros(len(receptor_ids), dtype=np.int64)
    rng = np.random.default_rng(seed)
    report_every = max(1, permutations // 10)
    for permutation in range(1, permutations + 1):
        shuffled_costs = rng.permutation(costs)
        sparse_graph.data[:] = shuffled_costs[arc_cost_index]
        distances = dijkstra(
            sparse_graph, directed=True,
            indices=list(source_index.values()), return_predecessors=False,
        )
        null_costs = distances[pair_sources, pair_targets]
        extreme += null_costs <= observed_costs
        receptor_extreme += (
            np.minimum.reduceat(null_costs, receptor_starts) <= observed_receptor
        )
        if progress and (permutation % report_every == 0 or permutation == permutations):
            progress(f"path permutations {permutation}/{permutations}")
    denominator = permutations + 1.0
    return (
        {pair: (extreme[index] + 1.0) / denominator for index, pair in enumerate(pairs)},
        {
            receptor: (receptor_extreme[index] + 1.0) / denominator
            for index, receptor in enumerate(receptor_ids)
        },
    )


def prepare_receiver_path_cache(
    receiver: str,
    receiver_genes: set[str],
    receiver_edges: pd.DataFrame,
    prior_edges: pd.DataFrame,
    receptors: set[str],
    path_permutations: int = 500,
    k_paths: int = 5,
    data_edges_per_node: int = 20,
    seed: int = 0,
    progress: Callable[[str], None] | None = None,
) -> ReceiverPathCache:
    """Calculate the weighted k shortest simple paths for one direction."""
    if k_paths <= 0:
        raise ValueError("k_paths must be a positive integer.")
    path_edges = select_path_graph_edges(receiver_edges, data_edges_per_node)
    graph = _build_graph(path_edges)
    tf_target, tf_candidates, terminal_only = _receiver_roles(
        receiver_genes, prior_edges
    )
    receptors = set(receptors) & receiver_genes & set(graph.nodes)
    tf_candidates &= set(graph.nodes)
    retained_count = int(receiver_edges["Edge_Retained"].sum()) if not receiver_edges.empty else 0
    unannotated_count = int((~path_edges["Prior_Found"].astype(bool)).sum()) if not path_edges.empty else 0
    prior_count = int(path_edges["Prior_Found"].sum()) if not path_edges.empty else 0
    if progress:
        progress(
            f"Receiver {receiver}: {retained_count} retained evidence edges -> "
            f"{len(path_edges)} path-graph edges "
            f"({prior_count} prior-annotated, {unannotated_count} unannotated); "
            f"{len(receptors)} receptors, {len(tf_candidates)} TFs"
        )

    pair_paths: dict[tuple[str, str], list[list[str]]] = {}
    receptor_list = sorted(receptors)
    for index, receptor in enumerate(receptor_list, start=1):
        found: dict[str, list[list[str]]] = {}
        for tf in sorted(tf_candidates - {receptor}):
            paths = _k_shortest_paths(
                graph, receptor, tf, receptors, tf_candidates, terminal_only,
                k_paths,
            )
            found[tf] = paths
            if paths:
                pair_paths[(receptor, tf)] = paths
        if progress:
            progress(
                f"Receiver {receiver}: receptor {index}/{len(receptor_list)} "
                f"({receptor}), found {sum(bool(v) for v in found.values())}/"
                f"{len(tf_candidates)} TF endpoints; up to {k_paths} paths each"
            )

    pair_p, receptor_p = _permutation_p_values(
        graph, pair_paths, receptors, tf_candidates, terminal_only,
        path_permutations, seed, progress,
    )
    return ReceiverPathCache(
        receiver=receiver,
        graph=graph,
        tf_target=tf_target,
        pair_paths=pair_paths,
        pair_p=pair_p,
        receptor_p=receptor_p,
        retained_input_edges=retained_count,
        path_graph_edges=len(path_edges),
        path_graph_unannotated_edges=unannotated_count,
        path_graph_prior_annotated_edges=prior_count,
    )


def attach_sender_paths(
    sender: str,
    receiver: str,
    sender_genes: set[str],
    receiver_genes: set[str],
    cache: ReceiverPathCache,
    ligand_receptor: pd.DataFrame,
    path_alpha: float = 0.05,
) -> PathResult:
    """Attach sender ligands and receiver targets to a cached R-TF network."""
    lr = ligand_receptor[
        ligand_receptor["Ligand"].isin(sender_genes)
        & ligand_receptor["Receptor"].isin(receiver_genes)
    ].copy()
    available_receptors = set(lr["Receptor"].astype(str))
    pair_paths = {
        pair: paths for pair, paths in cache.pair_paths.items()
        if pair[0] in available_receptors
    }
    receptor_ids = sorted({receptor for receptor, _ in pair_paths})
    metadata = {
        "eligible_lr_pairs": int(len(lr)),
        "eligible_receptors": int(len(available_receptors)),
        "tested_rtf_pairs": int(len(pair_paths)),
        "tested_receptors": len(receptor_ids),
        "significant_receptors": 0,
        "receptor_fdr_threshold": path_alpha,
        "retained_input_edges": cache.retained_input_edges,
        "path_graph_edges": cache.path_graph_edges,
        "path_graph_unannotated_edges": cache.path_graph_unannotated_edges,
        "path_graph_prior_annotated_edges": cache.path_graph_prior_annotated_edges,
    }
    if not pair_paths:
        return _empty_result(metadata)

    pairs = list(pair_paths)
    p_values = np.asarray([cache.pair_p.get(pair, np.nan) for pair in pairs])
    q_values = benjamini_hochberg(p_values)
    pair_statistics = {
        pair: (float(p_values[index]), float(q_values[index]))
        for index, pair in enumerate(pairs)
    }
    receptor_p_values = np.asarray([
        cache.receptor_p.get(receptor, np.nan) for receptor in receptor_ids
    ])
    receptor_q_values = benjamini_hochberg(receptor_p_values)
    receptor_statistics = {
        receptor: (
            float(receptor_p_values[index]),
            float(receptor_q_values[index]),
        )
        for index, receptor in enumerate(receptor_ids)
    }
    metadata.update({
        "tested_receptors": len(receptor_ids),
        "significant_receptors": int(np.sum(receptor_q_values <= path_alpha)),
        "receptor_fdr_threshold": path_alpha,
    })
    ligands_by_receptor = lr.groupby("Receptor")["Ligand"].apply(
        lambda values: ",".join(sorted(set(map(str, values))))
    ).to_dict()
    targets_by_tf = cache.tf_target.groupby("source")["target"].apply(
        lambda values: ",".join(sorted(set(map(str, values))))
    ).to_dict()

    rmtf_records: list[dict[str, object]] = []
    edge_records: list[dict[str, object]] = []
    main_records: list[dict[str, object]] = []
    path_counter = 0
    for (receptor, tf), paths in pair_paths.items():
        path_p, path_q = pair_statistics[(receptor, tf)]
        receptor_p, receptor_q = receptor_statistics[receptor]
        receptor_significant: bool | object = (
            bool(receptor_q <= path_alpha) if np.isfinite(receptor_q) else pd.NA
        )
        for rank, path in enumerate(paths, start=1):
            path_counter += 1
            path_id = f"{sender}_to_{receiver}_P{path_counter:06d}"
            attributes = [
                cache.graph[u][v] for u, v in zip(path[:-1], path[1:])
            ]
            cost = _path_cost(cache.graph, path)
            stabilities = [
                float(item["Stability"]) for item in attributes
                if pd.notna(item["Stability"])
            ]
            status = _direction_status(attributes)
            mediator = ",".join(path[1:-1])
            edge_types = ";".join(str(item["Relation_Type"]) for item in attributes)
            rmtf_records.append({
                "Path_ID": path_id, "Sender": sender, "Receiver": receiver,
                "Receptor": receptor, "Mediator": mediator, "TF": tf,
                "Path_Rank": rank, "Path_Length": len(path) - 1,
                "Path_Cost": cost, "Path_P": path_p, "Path_Q": path_q,
                "Receptor_P": receptor_p, "Receptor_Q": receptor_q,
                "Receptor_Significant": receptor_significant,
                "Alternative_Path_Count": len(paths),
                "Min_Stability": min(stabilities) if stabilities else np.nan,
                "Direction_Status": status,
            })
            for order, ((u, v), item) in enumerate(
                zip(zip(path[:-1], path[1:]), attributes), start=1
            ):
                edge_records.append({
                    "Path_ID": path_id, "Path_Rank": rank, "Edge_Order": order,
                    "Gene1": u, "Gene2": v,
                    "Relation_Type": item["Relation_Type"],
                    "Relation_Sign": item["Relation_Sign"],
                    "Prior_Found": item["Prior_Found"],
                    "Prior_Direction": item["Prior_Direction"],
                    "Stability": item["Stability"], "Data_Cost": item["Data_Cost"],
                    "Hub_Penalty": item["Hub_Penalty"],
                    "Hop_Penalty": item["Hop_Penalty"],
                    "Final_Cost": item["Final_Cost"],
                })
            ligand = ligands_by_receptor.get(receptor, "")
            target = targets_by_tf.get(tf, "")
            if ligand and target:
                main_records.append({
                    "Ligand": ligand, "Receptor": receptor, "Mediator": mediator,
                    "TF": tf, "Target": target, "Path_ID": path_id,
                    "Path_Rank": rank, "Path_Length": len(path) - 1,
                    "Path_Cost": cost, "Path_P": path_p, "Path_Q": path_q,
                    "Receptor_P": receptor_p, "Receptor_Q": receptor_q,
                    "Receptor_Significant": receptor_significant,
                    "Edge_Types": edge_types,
                    "Direction_Status": status,
                })
    path_edge_table = pd.DataFrame(edge_records, columns=PATH_EDGE_COLUMNS)
    subnetwork_records: list[dict[str, object]] = []
    if not path_edge_table.empty:
        for (gene1, gene2), group in path_edge_table.groupby(
            ["Gene1", "Gene2"], sort=True
        ):
            first = group.sort_values(
                ["Path_Rank", "Final_Cost", "Path_ID"], kind="mergesort"
            ).iloc[0]
            subnetwork_records.append({
                "Sender": sender,
                "Receiver": receiver,
                "Gene1": gene1,
                "Gene2": gene2,
                "Best_Path_Rank": int(group["Path_Rank"].min()),
                "Path_Count": int(group["Path_ID"].nunique()),
                "Relation_Type": first["Relation_Type"],
                "Relation_Sign": first["Relation_Sign"],
                "Prior_Direction": first["Prior_Direction"],
                "Min_Final_Cost": float(group["Final_Cost"].min()),
            })
    return PathResult(
        pd.DataFrame(main_records, columns=MAIN_COLUMNS),
        pd.DataFrame(rmtf_records, columns=RMTF_COLUMNS),
        path_edge_table,
        pd.DataFrame(subnetwork_records, columns=SUBNETWORK_EDGE_COLUMNS),
        metadata,
    )


def infer_paths(
    sender: str,
    receiver: str,
    sender_genes: set[str],
    receiver_genes: set[str],
    receiver_edges: pd.DataFrame,
    prior_edges: pd.DataFrame,
    ligand_receptor: pd.DataFrame,
    path_permutations: int = 500,
    k_paths: int = 5,
    data_edges_per_node: int = 20,
    path_alpha: float = 0.05,
    seed: int = 0,
) -> PathResult:
    """Compatibility wrapper for one sender-receiver analysis."""
    lr = ligand_receptor[
        ligand_receptor["Ligand"].isin(sender_genes)
        & ligand_receptor["Receptor"].isin(receiver_genes)
    ]
    receptors = set(lr["Receptor"].astype(str))
    cache = prepare_receiver_path_cache(
        receiver, receiver_genes, receiver_edges, prior_edges, receptors,
        path_permutations=path_permutations,
        k_paths=k_paths,
        data_edges_per_node=data_edges_per_node,
        seed=seed,
    )
    return attach_sender_paths(
        sender, receiver, sender_genes, receiver_genes, cache, ligand_receptor,
        path_alpha=path_alpha,
    )
