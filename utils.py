import igraph as ig

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

def forward_iter(G, n, weight=True):
    """
    @brief  Forward neighbor iterator for a node in ig.Graph
    """
    if weight:
        return [(G.es[e].target, G.es[e]['weight']) for e in G.incident(n, mode=ig.OUT)]
    else:
        return [G.es[e].target for e in G.incident(n, mode=ig.OUT)]

def reverse_iter(G, n, weight=True):
    """
    @brief  Reverse neighbor iterator for a node in ig.Graph
    """
    if weight:
        return [(G.es[e].source, G.es[e]['weight']) for e in G.incident(n, mode=ig.IN)]
    else:
        return [G.es[e].source for e in G.incident(n, mode=ig.IN)]

def count_simple_paths(G, predecessors_iter, successors_iter, sources, paths):
    """
    @brief  Counts the number of simple paths from the source vertices in the network.
    
    @param G                  ig.Graph representation of the network.
    @param predecessors_iter  Iterator provider over predecessors of a vertex in the network.
    @param successors_iter    Iterator provider over successors of a vertex in the network.
    @param sources            Source vertices in the network.
    @param paths              Array of the number of simple paths from the sources to every vertex.
    """
    while sources:
        next_level = set()
        for u in sources:
            paths[u] = sum((paths[p] * w) for p, w in predecessors_iter(G, u, weight=True))
            next_level.update(set(s for s in successors_iter(G, u, weight=False)))
        sources = next_level
