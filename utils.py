def forward_iter(G, n):
    """
    @brief  Forward neighbor iterator for a node in nx.DiGraph or nx.MultiDiGraph
    """
    for u, v in G.out_edges_iter(n):
        yield v

def reverse_iter(G, n):
    """
    @brief  Reverse neighbor iterator for a node in nx.DiGraph or nx.MultiDiGraph
    """
    for u, v in G.in_edges_iter(n):
        yield u

def count_simple_paths(G, predecessors_iter, successors_iter, sources, paths):
    """
    @brief  Counts the number of simple paths from the source vertices in the network
    
    @param G                  nx.DiGraph or nx.MultiDiGraph representation of the network
    @param predecessors_iter  Iterator provider over predecessors of a vertex in the network
    @param successors_iter    Iterator provider over successors of a vertex in the network
    @param sources            Source vertices in the network
    @param paths              Array of the number of simple paths from the sources to every vertex
    """
    next_level = set()
    for u in sources:
        paths[u] = max(paths[u], sum(paths[p] for p in predecessors_iter(G, u)))
        next_level.update(set(s for s in successors_iter(G, u)))
    if next_level:
        count_simple_paths(G, predecessors_iter, successors_iter, next_level, paths)
