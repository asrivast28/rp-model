import igraph as ig
from itertools import izip
import numpy as np


def datatype(value_max):
    """
    @brief  Returns the unsigned NumPy datatype suitable for storing value_max.
    """
    import numpy as np
    from bisect import bisect
    all_dtypes = (np.uint8, np.uint16, np.uint32, np.uint64)
    dtype_max = [np.iinfo(dtype).max for dtype in all_dtypes]
    dtype_max.insert(0, 0)
    return all_dtypes[bisect(dtype_max, value_max) - 1]

def count_simple_paths(G, sources, predecessor, successor, paths):
    """
    @brief  Counts the number of simple paths from the source vertices in the network.

    @param G            ig.Graph representation of the network.
    @param sources      Source vertices in the network.
    @param predecessor  Mode for finding predecessor vertices in the network.
    @param successor    Mode for finding successor vertices in the network.
    @param paths        Array of the number of simple paths from the sources to every vertex.
    """
    weights = G['weights']
    while sources:
        next_level = set()
        for u in sources:
            paths[u] = np.dot(paths[G.neighbors(u, mode=predecessor)], weights[G.incident(u, mode=predecessor)])
            next_level.update(G.neighbors(u, mode=successor))
        sources = next_level
