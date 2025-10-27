from scipy.sparse import coo_matrix
import networkx as nx
import numpy as np


def find_shortest_path(G, source, target, max_length=10):
    """Find a path between source and target nodes using NetworkX."""
    try:
        shortest_paths = list(nx.all_shortest_paths(G, source=source, target=target))
        if len(shortest_paths) > 0 and len(shortest_paths[0]) <= max_length:
            return shortest_paths
        else:
            return None
    except nx.NetworkXNoPath:
        return None


def find_all_path(G, source, target, cutoff=5):
    """Find a path between source and target nodes using NetworkX."""
    try:
        paths = list(nx.all_simple_paths(G, source=source, target=target, cutoff=cutoff))
        return paths
    except nx.NetworkXNoPath:
        return None


def find_paths_fixed_length(G, source, target, length):
    """
    Find all paths of a given length in a directed graph.

    Parameters:
    G (networkx.DiGraph): Directed graph.
    source (int): Source node.
    target (int): Target node.
    length (int): Desired path length.

    Returns:
    List[List[int]]: All paths of the specified length from source to target.
    """
    paths = []  # To store all valid paths

    def dfs(current_node, current_path):
        # If path length is met, check if the last node is the target
        if len(current_path) == length + 1:
            if current_node == target:
                paths.append(current_path[:])  # Add a copy of the current path
            return

        # Explore neighbors
        for neighbor in G.successors(current_node):
            if neighbor not in current_path:  # Avoid cycles
                current_path.append(neighbor)
                dfs(neighbor, current_path)
                current_path.pop()  # Backtrack

    # Start DFS from the source node
    dfs(source, [source])
    return paths


def infer_pathway(gene_set, lg_rp_dict, A3_sparse, A4_sparse, max_length):
    G = nx.from_numpy_array(A3_sparse.toarray(), create_using=nx.DiGraph)
    lg_rp_pairs = [(k, v) for k, s in lg_rp_dict.items() for v in s]

    pathways_lg_rp = [
        [ligand, receptor]
        for ligand in gene_set
        for _, receptor in lg_rp_pairs if ligand == _
    ]

    source_set = set(p[-1] for p in pathways_lg_rp)  # set of receptors
    target_set = set(A4_sparse.row)  # Genes with outgoing edges in A4

    shortest_pathways = []
    for source in source_set:
        for target in target_set:
            if source != target:
                shortest_path = find_shortest_path(G, source=source, target=target, max_length=max_length)
                if shortest_path:
                    for path in shortest_path:
                        for p in pathways_lg_rp:
                            if p[-1] == source:  # Match source gene in [a, b, c]
                                shortest_pathways.append(p + path[1:])
    return shortest_pathways

