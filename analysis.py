from itertools import combinations
import numpy as np

import utils


def path_centrality(G, P_s, P_t, centrality):
    """
    @brief  Computes path centrality values for all the vertices in a network.

    @param G           nx.DiGraph or nx.MultiDiGraph representation of the network.
    @param P_s         NumPy array with 1 for every source vertex and 0 for other vertices.
    @param P_s         NumPy array with 1 for every target vertex and 0 for other vertices.
    @param centrality  NumPy array containg centrality value for every vertex in the given network.
    """
    # Compute the number of paths from sources to every vertex: complexity
    first_level = set()
    for s in np.where(P_s)[0]:
        first_level.update(n for n in utils.forward_iter(G, s))
    utils.count_simple_paths(G, utils.reverse_iter, utils.forward_iter, first_level, P_s)

    # Compute the number of paths from every vertex to targets: generality
    first_level = set()
    for t in np.where(P_t)[0]:
        first_level.update(n for n in utils.reverse_iter(G, t))
    utils.count_simple_paths(G, utils.forward_iter, utils.reverse_iter, first_level, P_t)

    # Multiply complexity and generality to get the path centrality
    np.multiply(P_s, P_t, out=centrality)

def remove_vertex(G, vertex, source, target, in_degree, out_degree):
    """
    @brief  Removes all the edges to and from the given vertex in the given network.

    @param G           nx.DiGraph representation of the network.
    @param vertex      Vertex to be removed from the network.
    @param in_degree   NumPy array containing in-degree for every vertex.
    @param out_degree  NumPy array containing out-degree for every vertex.

    @return  List containing all the deleted edges.
    """
    for p, v in G.in_edges_iter(vertex):
        out_degree[p] -= 1
    for v, s in G.out_edges_iter(vertex):
        in_degree[s] -= 1
    in_degree[vertex] = 0
    out_degree[vertex] = 0
    source[source & (out_degree == 0)] = False
    target[target & (in_degree == 0)] = False
    vertex_edges = G.out_edges(vertex)
    vertex_edges.extend(G.in_edges(vertex))
    G.remove_node(vertex)
    G.add_node(vertex)
    return vertex_edges

def core_vertices(G, source, target, tau, datatype=np.uint64):
    """
    @brief  Greedily finds core vertices in the given network.

    @param G         nx.DiGraph or nx.MultiDiGraph representation of the network.
    @param source    NumPy array of type bool with 1 for every source vertex.
    @param target    NumPy array of type bool with 1 for every target vertex.
    @param tau       Fraction of S-T paths to be covered by the core vertices.
    @param datatype  NumPy datatype provided as a hint for storage.

    @return  List containing all the core vertices.
    """
    # Copy objects that are going to be modified
    G = G.copy()
    source = np.copy(source)
    target = np.copy(target)

    # Compute initial in and out degree for the network
    in_degree = np.array(list(G.in_degree(n) for n in xrange(G.number_of_nodes())), dtype=datatype)
    out_degree = np.array(list(G.out_degree(n) for n in xrange(G.number_of_nodes())), dtype=datatype)

    # Compute initial path centralities
    # Initialize complexity for the source vertices to 1
    P_s = source.astype(datatype)
    # Initialize generality for the target vertices to 1
    P_t = target.astype(datatype)
    centrality = np.empty(source.size, datatype)
    path_centrality(G, P_s, P_t, centrality)

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
                _source = np.copy(source)
                _target = np.copy(target)
                vertex_edges = remove_vertex(G, vertex, _source, _target, np.copy(in_degree), np.copy(out_degree))
                P_s[source] = 1
                P_s[~source] = 0
                P_t[target] = 1
                P_t[~target] = 0
                path_centrality(G, P_s, P_t, centrality)
                neighbors = set([vertex])
                for other in candidate_vertices:
                    if other != vertex and centrality[other] == 0:
                        neighbors.add(other)
                PES.append(neighbors)
                G.add_edges_from(vertex_edges)
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
        P_s[source] = 1
        P_s[~source] = 0
        P_t[target] = 1
        P_t[~target] = 0
        path_centrality(G, P_s, P_t, centrality)

        P_R = P - np.sum(source * centrality)
    return C
