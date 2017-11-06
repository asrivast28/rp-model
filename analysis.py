from itertools import combinations
import numpy as np


def _compute_num_paths(G, this_level, predecessors_iter, successors_iter, P_v):
    # Compute number of paths for all the vertices
    next_level = set()
    for u in this_level:
        P_v[u] = max(P_v[u], sum(P_v[p] for p in predecessors_iter(G, u)))
        next_level.update(set(s for s in successors_iter(G, u)))
    if next_level:
        _compute_num_paths(G, np.array(list(next_level)), predecessors_iter, successors_iter, P_v)

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

def path_centrality(G, source, target, pathtype=np.uint64):
    """
    @brief  Computes path centrality values for all the vertices in a network.

    @param G           nx.DiGraph or nx.MultiDiGraph representation of the network.
    @param source      NumPy array of type bool with 1 for every source vertex.
    @param target      NumPy array of type bool with 1 for every target vertex.
    @param pathtype    NumPy datatype provided as a hint for storing centrality values.

    @return  A NumPy array containg centrality value for every vertex in the given network.
    """

    # Compute the number of paths from sources to every vertex: complexity
    # Initialize complexity for the source vertices to 1
    P_s = source.astype(pathtype)
    _compute_num_paths(G, np.where(source)[0], reverse_iter, forward_iter, P_s)

    # Compute the number of paths from every vertex to targets: generality
    # Initialize generality for the target vertices to 1
    P_t = target.astype(pathtype)
    _compute_num_paths(G, np.where(target)[0], forward_iter, reverse_iter, P_t)

    # Multiply complexity and generality to get the path centrality
    centrality = P_s * P_t
    return centrality

def remove_vertex(G, vertex, source, target, in_degree, out_degree):
    """
    @brief  Removes all the edges to and from the given vertex in the given network.

    @param G           nx.DiGraph representation of the network.
    @param vertex      Vertex to be removed from the network.
    @param in_degree   NumPy array containing in-degree for every vertex.
    @param out_degree  NumPy array containing out-degree for every vertex.
    """
    for p, v in G.in_edges_iter(vertex):
        out_degree[p] -= 1
    for v, s in G.out_edges_iter(vertex):
        in_degree[s] -= 1
    in_degree[vertex] = 0
    out_degree[vertex] = 0
    source[source & (out_degree == 0)] = False
    target[target & (in_degree == 0)] = False
    G.remove_node(vertex)
    G.add_node(vertex)

def core_vertices(G, source, target, tau, vertextype=np.uint32):
    """
    @brief  Greedily finds core vertices in the given network.

    @param G           nx.DiGraph or nx.MultiDiGraph representation of the network.
    @param source      NumPy array of type bool with 1 for every source vertex.
    @param target      NumPy array of type bool with 1 for every target vertex.
    @param tau         Fraction of S-T paths to be covered by the core vertices.
    @param vertextype  NumPy datatype provided as a hint for storing degrees.

    @return  List containing all the core vertices.
    """
    # Copy objects that are going to be modified
    G = G.copy()
    source = np.copy(source)
    target = np.copy(target)

    # Compute initial in and out degree for the network
    in_degree = np.array(list(G.in_degree(n) for n in xrange(G.number_of_nodes())), dtype=vertextype)
    out_degree = np.array(list(G.out_degree(n) for n in xrange(G.number_of_nodes())), dtype=vertextype)

    # Compute initial path centralities
    centrality = path_centrality(G, source, target)

    # Compute the initial number of S-T paths
    # Same as the number of paths exiting from the sources
    P = np.sum(source * centrality)

    C = []
    P_R = 0
    while float(P_R) < (tau * P):
        max_centrality = np.amax(centrality)
        candidate_vertices = np.where(centrality == max_centrality)[0]
        # Remove the first vertex with maximum centrality
        candidate_vertex = set()
        if candidate_vertices.size == 1:
            # Remove the vertex from the network
            candidate_vertex.add(candidate_vertices[0])
        else:
            PES = []
            # Check if two or more vertices are on the same path
            for vertex in candidate_vertices:
                _G = G.copy()
                _source = np.copy(source)
                _target = np.copy(target)
                _in_degree = np.copy(in_degree)
                _out_degree = np.copy(out_degree)
                remove_vertex(_G, vertex, _source, _target, _in_degree, _out_degree)
                _centrality = path_centrality(_G, _source, _target)
                neighbors = set([vertex])
                for other in candidate_vertices:
                    if other != vertex and _centrality[other] == 0:
                        neighbors.add(other)
                PES.append(neighbors)
            updates = True
            while updates:
                for i, j in combinations(range(len(PES)), 2):
                    if PES[i].intersection(PES[j]):
                        PES[i].update(PES[j])
                        PES.pop(j)
                        break
                else:
                    updates = False
            candidate_vertex = PES[0]
        for vertex in candidate_vertex:
            remove_vertex(G, vertex, source, target, in_degree, out_degree)
        C.append(candidate_vertex if len(candidate_vertex) > 1 else candidate_vertex.pop())

        # Recompute path centralities
        centrality = path_centrality(G, source, target)

        P_R = P - np.sum(source * centrality)
    return C
