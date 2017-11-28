from itertools import combinations
import numpy as np

import utils


def update_path_centrality(G, P_s, P_t, centrality, vertex=None):
    """
    @brief  Updates path centrality values for the vertices in a network.

    @param G           np.array representation of the network's weighted adjacency matrix.
    @param P_s         np.array with 1 for every source vertex and 0 for other vertices.
    @param P_s         np.array with 1 for every target vertex and 0 for other vertices.
    @param centrality  np.array containg centrality value for every vertex in the given network.
    @param vertex      Only update successors and predecessors of the vertex, if it is not None.
    """
    # Compute the number of paths from sources to every vertex: complexity
    first_level = set()
    for s in np.where(P_s)[0] if vertex is None else vertex:
        first_level.update(np.nonzero(G[s])[0])
    utils.count_simple_paths(G, G.T, first_level, P_s)

    # Compute the number of paths from every vertex to targets: generality
    first_level = set()
    for t in np.where(P_t)[0] if vertex is None else vertex:
        first_level.update(np.nonzero(G.T[t])[0])
    utils.count_simple_paths(G.T, G, first_level, P_t)

    # Multiply complexity and generality to get the path centrality
    np.multiply(P_s, P_t, out=centrality)

def remove_vertex(G, vertex, source, target, in_degree, out_degree):
    """
    @brief  Sets weights for all the edges to and from the given vertices to zero.

    @param G           np.array representation of the network's weighted adjacency matrix.
    @param vertex      Vertex to be removed from the network.
    @param in_degree   np.array containing in-degree for every vertex.
    @param out_degree  np.array containing out-degree for every vertex.
    """
    for v in vertex:
        in_degree[np.nonzero(G[v])[0]] -= 1
        G[v] = 0
        out_degree[np.nonzero(G.T[v])[0]] -= 1
        G.T[v] = 0
    in_degree[vertex] = 0
    out_degree[vertex] = 0
    source[source & (out_degree == 0)] = False
    target[target & (in_degree == 0)] = False

def core_vertices(G, source, target, tau, datatype=np.float64):
    """
    @brief  Greedily finds core vertices in the given network.

    @param G         np.array representation of the network's weighted adjacency matrix.
    @param source    np.array of type bool with 1 for every source vertex.
    @param target    np.array of type bool with 1 for every target vertex.
    @param tau       Fraction of S-T paths to be covered by the core vertices.
    @param datatype  np.dtype provided as a hint for storage.

    @return  List containing all the core vertices.
    """
    # Copy objects that are going to be modified
    source = np.copy(source)
    target = np.copy(target)

    # Compute initial in and out degree for the network
    vertextype = utils.datatype(G.shape[0])
    out_degree = np.count_nonzero(G, axis=1)
    in_degree = np.count_nonzero(G, axis=0)

    # Initialize complexity for the source vertices to 1
    P_s = source.astype(datatype)
    # Initialize generality for the target vertices to 1
    P_t = target.astype(datatype)
    # Compute initial path centralities
    centrality = np.empty(source.size, datatype)
    update_path_centrality(G, P_s, P_t, centrality)

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
                _P_s = np.copy(P_s)
                _P_t = np.copy(P_t)
                _P_s[vertex] = 0
                _P_t[vertex] = 0
                update_path_centrality(G, _P_s, _P_t, centrality, [vertex])
                neighbors = set([vertex])
                for other in candidate_vertices:
                    if other != vertex and centrality[other] == 0:
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
        candidate_vertex = sorted(list(candidate_vertex))
        P_s[candidate_vertex] = 0
        P_t[candidate_vertex] = 0
        update_path_centrality(G, P_s, P_t, centrality, candidate_vertex)
        remove_vertex(G, candidate_vertex, source, target, in_degree, out_degree)
        C.append(candidate_vertex if len(candidate_vertex) > 1 else candidate_vertex[0])

        P_R = P - np.sum(source * centrality)

    return P, C
