import numpy as np
import networkx as nx
from collections import deque


def find_path_with_networkx(A, source, target):
    """Find a path between source and target nodes using NetworkX."""
    G = nx.from_numpy_array(A, create_using=nx.DiGraph)  # Use from_numpy_array
    try:
        path = nx.shortest_path(G, source=source, target=target)
        return path
    except nx.NetworkXNoPath:
        return None  # No path found


def bfs_shortest_path(A, source, target, max_length=10):
    """
    Find the shortest path from source to target using BFS with a max_length limit.
    If the path length exceeds max_length, terminate the search and return None.
    """
    n = len(A)
    queue = deque([(source, [source])])  # Queue stores (node, path) tuples
    visited = set()

    while queue:
        node, path = queue.popleft()

        # If path length exceeds max_length, stop search
        if len(path) > max_length:
            return None

        # If target node is found, return the path
        if node == target:
            return path

        # Explore neighbors if the node is not visited
        if node not in visited:
            visited.add(node)
            for neighbor, connected in enumerate(A[node]):
                if connected and neighbor not in visited:
                    queue.append((neighbor, path + [neighbor]))
    return None  # No path found within max_length

# Example usage


A = np.array([
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0]
])

source = 0
target = 3

path = bfs_shortest_path(A, source, target, 3)
print(path)  # Output: [0, 1, 2, 3]
