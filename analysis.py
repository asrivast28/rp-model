from itertools import combinations
import numpy as np

def _compute_num_paths(vertex, neighbor_iter, degree, pathtype):
    # Initialize complexity for the initial vertices to 1
    P_v = vertex.astype(pathtype)
    # Compute number of paths for all the vertices 
    this_level = np.where(P_v)[0]
    while this_level.size > 0:
        next_level = set()
        for u in this_level:
            for n in neighbor_iter(u):
                degree[n] -= 1
                if degree[n] == 0:
                    next_level.add(n)
                P_v[n] += P_v[u]
        this_level = np.array(list(next_level))
    return P_v

def path_centrality(G, source, target, in_degree, out_degree, pathtype=np.uint64):
    """
    @brief  Computes path centrality values for all the vertices in a network.

    @param G           nx.DiGraph or nx.MultiDiGraph representation of the network.
    @param source      NumPy array of type bool with 1 for every source vertex.
    @param target      NumPy array of type bool with 1 for every target vertex.
    @param in_degree   NumPy array containing in-degree for every vertex.
    @param out_degree  NumPy array containing out-degree for every vertex.
    @param pathtype    NumPy datatype provided as a hint for storing centrality values.

    @return  A NumPy array containg centrality value for every vertex in the given network.
    """
    source = source & (out_degree > 0)
    target = target & (in_degree > 0)

    # Compute the number of paths from sources to every vertex: complexity
    P_s = _compute_num_paths(source, G.successors_iter, np.copy(in_degree), pathtype)

    # Compute the number of paths from every vertex to targets: generality
    P_t = _compute_num_paths(target, G.predecessors_iter, np.copy(out_degree), pathtype)

    # Multiply complexity and generality to get the path centrality
    centrality = P_s * P_t
    return centrality

def remove_vertex(G, vertex, in_degree, out_degree):
    """
    @brief  Removes all the edges to and from the given vertex in the given network.

    @param G           nx.DiGraph representation of the network.
    @param vertex      Vertex to be removed from the network.
    @param in_degree   NumPy array containing in-degree for every vertex.
    @param out_degree  NumPy array containing out-degree for every vertex.
    """
    for u, v in G.in_edges_iter(vertex):
        out_degree[u] -= 1
    for u, v in G.out_edges_iter(vertex):
        in_degree[v] -= 1
    in_degree[vertex] = 0
    out_degree[vertex] = 0
    G.remove_node(vertex)
    G.add_node(vertex)

def core_vertices(G, tau, vertextype=np.uint32):
    """
    @brief  Greedily finds core vertices in the given network.

    @param G           nx.DiGraph or nx.MultiDiGraph representation of the network.
    @param tau         Fraction of S-T paths to be covered by the core vertices.
    @param vertextype  NumPy datatype provided as a hint for storing degrees.

    @return  List containing all the core vertices.
    """
    G = G.copy()

    # Compute initial in and out degree for the network
    in_degree = np.array(list(G.in_degree(n) for n in xrange(G.number_of_nodes())), dtype=vertextype)
    out_degree = np.array(list(G.out_degree(n) for n in xrange(G.number_of_nodes())), dtype=vertextype)

    # Compute boolean array specifying if every vertex is source or target
    source = (in_degree == 0)
    target = (out_degree == 0)

    # Compute initial path centralities
    centrality = path_centrality(G, source, target, in_degree, out_degree)

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
                _in_degrees = np.copy(in_degree)
                _out_degrees = np.copy(out_degree)
                remove_vertex(_G, vertex, _in_degrees, _out_degrees)
                _centrality = path_centrality(_G, source, target, _in_degrees, _out_degrees)
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
            remove_vertex(G, vertex, in_degree, out_degree)
        C.append(candidate_vertex if len(candidate_vertex) > 1 else candidate_vertex.pop())

        # Recompute path centralities
        centrality = path_centrality(G, source, target, in_degree, out_degree)

        P_R = P - np.sum(source * centrality)
    return C
