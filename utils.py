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

def count_simple_paths(G, G_T, sources, paths):
    """
    @brief  Counts the number of simple paths from the source vertices in the network.

    @param G         np.array representation of the network's weighted adjacency matrix.
    @param G_T       Transpose of G.
    @param sources   Source vertices in the network.
    @param paths     np.array containing the number of simple paths from the sources to every vertex.
    """
    while sources:
        next_level = set()
        for u in sources:
            paths[u] = np.dot(paths, G_T[u])
            next_level.update(np.nonzero(G[u])[0])
        sources = next_level
