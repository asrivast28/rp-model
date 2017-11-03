import networkx as nx
import numpy as np


def read(filename):
    """
    @brief  Reads a network from the file in edge-list format.

    @param filename  Name of the file to be read.

    @return  A nx.DiGraph instance of the read network.
    """
    # initialize a graph
    G = nx.DiGraph()
    # with open('centrality_test,txt', 'rb') as f:
    with open(filename, 'rb') as f:
        for edge in f.xreadlines():
            if not edge.startswith('%'):
                G.add_edge(*tuple(int(u) for u in edge.split()))
    return G

def rp_model(S, M, T, alpha, d_in):
    """
    @brief  Creates a network using the RP-model.

    @param S      Number of source vertices.
    @param M      Number of intermediate vertices.
    @param T      Number of target vertices.
    @param alpha  Preference for reuse.
    @param d_in   Function that generates in-degrees.

    @return  A nx.DiGraph instance of the generated network.
    """
    # Create an empty directed network
    G = nx.DiGraph()
    # Add all the vertices to the network
    G.add_nodes_from(xrange(S + M + T))

    # Source ranks array
    source_ranks = np.arange(S, dtype=np.uint32)
    # Intermediate ranks array
    inter_ranks = np.zeros(M, dtype=np.uint32)
    inter_indices = np.arange(M)
    # Create connections for the M intermediates
    for m in xrange(M):
        # Increase source ranks by one
        # for calculating the probabilities
        source_ranks = source_ranks + 1
        # Randomly assign ranks m through
        # m + (S-1) to the sources
        np.random.shuffle(source_ranks)
        # Increase intermediate ranks by one
        # for calculating the probabilities
        inter_ranks = inter_ranks + (inter_indices < m)
        # Ranks of all the vertices for creating this vertex
        all_ranks = np.concatenate((source_ranks, inter_ranks)).astype(np.float)
        with np.errstate(divide='ignore'):
            numerators = np.power(all_ranks, -1.0 * alpha)
            numerators[~np.isfinite(numerators)] = 0.0
        # Probability of connecting to every older vertex
        probabilities = numerators / np.sum(numerators)
        # Pick unique source vertices and add incoming edges from them
        for u in np.random.choice(S + M, size=d_in(), replace=False, p=probabilities):
            G.add_edge(u, S + m)

    # Increase ranks by one for calculating the probabilities
    source_ranks = source_ranks + 1
    inter_ranks = inter_ranks + 1
    # Randomly assign ranks M through
    # M + (S-1) to the sources
    np.random.shuffle(source_ranks)
    all_ranks = np.concatenate((source_ranks, inter_ranks)).astype(np.float)
    numerators = np.power(all_ranks, -1.0 * alpha)
    # Probability of connecting to every older vertex
    probabilities = numerators / np.sum(numerators)

    # Create connections for the T targets in a batch
    for t in xrange(T):
        # Pick unique source vertices and add incoming edges from them
        for u in np.random.choice(S + M, size=d_in(), replace=False, p=probabilities):
            G.add_edge(u, S + M + t)

    return G

def count_simple_paths(G, s, t, npath=None, pathtype=np.uint32):
    """
    @brief  Counts the number of simple s-t paths in the network.
            Modified the code from https://cs.stackexchange.com/q/3087

    @param G         nx.DiGraph representation of the network.
    @param s         Source vertex in the network.
    @param t         Target vertex in the network.
    @param npath     Number of paths from every vertex to t.
    @param pathtype  NumPy datatype provided as a hint for storing the path counts.

    @return  Number of paths from s to t in the given network.
    """
    if npath is None:
        npath = np.zeros(G.number_of_nodes(), dtype=pathtype)
    if s == t:
        return 1
    else:
        # assume sum returns 0 if s has no children
        if npath[s] == 0:
            npath[s] = sum(count_simple_paths(G, u, t, npath) for u in G.successors(s))
        return npath[s]

def flatten(G, vertextype=np.uint32):
    """
    @brief  Flattens the given dependency network.

    @param G           nx.DiGraph representation of the network.
    @param vertextype  NumPy datatype given as hint for storing simple paths.

    @return  nx.MultiDiGraph representation of the flattened dependency network.
    """
    # Create a directed flat network which allows parallel edges
    G_f = nx.MultiDiGraph()
    # Add the same nodes as the original network
    G_f.add_nodes_from(xrange(G.number_of_nodes()))

    # Find sources and targets in the original network
    source = np.array(list(G.in_degree(n) == 0 for n in xrange(G.number_of_nodes())), dtype=np.bool)
    target = np.array(list(G.out_degree(n) == 0 for n in xrange(G.number_of_nodes())), dtype=np.bool)

    # Create the flat dependency network
    for s in np.where(source)[0]:
        for t in np.where(target)[0]:
            # Compute the number of simple s-t paths
            st_paths = count_simple_paths(G, s, t)
            # Create multiple parallel edges between s and t
            # one for every s-t path
            G_f.add_edges_from((s, t) for p in xrange(st_paths))
    return G_f
