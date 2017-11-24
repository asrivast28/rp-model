import igraph as ig
import numpy as np

import utils


def read(filename):
    """
    @brief  Reads a network from the file in edge-list format.

    @param filename  Name of the file to be read.

    @return  A ig.Graph instance of the read network.
    """
    # initialize a graph
    # with open('centrality_test,txt', 'rb') as f:
    with open(filename, 'rb') as elf:
        G = ig.Graph.Read_Edgelist(elf)
    vertex = np.arange(G.vcount(), dtype=utils.datatype(G.vcount()))
    source = np.vectorize(lambda v : G.degree(v, mode=ig.IN) == 0, otypes=[np.bool])(vertex)
    target = np.vectorize(lambda v : G.degree(v, mode=ig.OUT) == 0, otypes=[np.bool])(vertex)
    return G, source, target

def rp_model(S, M, T, alpha, d_in, out):
    """
    @brief  Creates a network using the RP-model.

    @param S      Number of source vertices.
    @param M      Number of intermediate vertices.
    @param T      Number of target vertices.
    @param alpha  Preference for reuse.
    @param d_in   Function that generates in-degrees.
    @param out    Prefix of the files to which the network is written.

    @return  A ig.Graph instance of the generated network.
    """
    # Add all the vertices to the network
    V = S + M + T
    vertextype = utils.datatype(V)
    # Create a directed network
    G = ig.Graph(V, directed=True)

    # Source ranks array
    source_ranks = np.arange(S, dtype=vertextype)
    # Intermediate ranks array
    inter_ranks = np.zeros(M, dtype=vertextype)
    inter_indices = np.arange(M)

    k = float(d_in())
    G.add_edges((s, S) for s in xrange(S))
    G.es['weight'] = k/S

    # Create connections for rest of the (M-1) intermediates
    for m in xrange(1, M):
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
        if not np.allclose(alpha, 0.0):
            with np.errstate(divide='ignore'):
                numerators = np.power(all_ranks, -1.0 * alpha)
                numerators[~np.isfinite(numerators)] = 0.0
        else:
            numerators = (all_ranks > 0).astype(np.float)
        # Probability of connecting to every older vertex
        probabilities = numerators / np.sum(numerators)
        # Pick unique source vertices and add incoming edges from them
        for u in np.random.choice(S + M, size=d_in(), replace=False, p=probabilities):
            G.add_edge(u, S + m, weight=1.0)

    # Increase ranks by one for calculating the probabilities
    source_ranks = source_ranks + 1
    inter_ranks = inter_ranks + 1
    # Randomly assign ranks M through
    # M + (S-1) to the sources
    np.random.shuffle(source_ranks)
    all_ranks = np.concatenate((source_ranks, inter_ranks)).astype(np.float)
    if not np.allclose(alpha, 0.0):
        numerators = np.power(all_ranks, -1.0 * alpha)
    else:
        numerators = (all_ranks > 0).astype(np.float)
    # Probability of connecting to every older vertex
    probabilities = numerators / np.sum(numerators)

    # Create connections for the T targets in a batch
    for t in xrange(T):
        # Pick unique source vertices and add incoming edges from them
        for u in np.random.choice(S + M, size=d_in(), replace=False, p=probabilities):
            G.add_edge(u, S + M + t, weight=1.0)

    vertex = np.arange(V, dtype=vertextype)
    source = (vertex < S)
    source[list(s for s in xrange(S) if G.degree(s, mode=ig.OUT) == 0)] = False
    target = (vertex >= (S + M))
    if out is not None:
        with open('%s_links.txt'%out, 'wb') as elf:
            G.write_edgelist(elf)
        with open('%s_sources.txt'%out, 'wb') as sf:
            sf.write('\n'.join(str(s) for s in np.where(source)[0]))
        with open('%s_targets.txt'%out, 'wb') as tf:
            tf.write('\n'.join(str(t) for t in np.where(target)[0]))
    return G, source, target

def flatten(G, source, target, datatype=np.uint64):
    """
    @brief  Flattens the given dependency network.

    @param G         ig.Graph representation of the network.
    @param source    NumPy array of type bool with 1 for every source vertex.
    @param target    NumPy array of type bool with 1 for every target vertex.
    @param datatype  NumPy datatype provided as a hint for storage.

    @return  Weighted ig.Graph representation of the flattened dependency network.
    """
    # Create a weighted directed network
    G_f = ig.Graph(directed=True, edge_attrs={'weight' : []})
    # Add the same nodes as the original network
    G_f.add_vertices(G.vcount())

    # Create the flat dependency network
    P_s = np.empty(G.vcount(), dtype=datatype)
    targets = np.where(target)[0]
    weight = np.array([])
    for s in np.where(source)[0]:
        P_s.fill(0)
        P_s[s] = 1
        utils.count_simple_paths(G, utils.reverse_iter, utils.forward_iter, set(n for n in utils.forward_iter(G, s, weight=False)), P_s)
        G_f.add_edges((s, t) for t in targets)
        weight = np.append(weight, [P_s[t] for t in targets])
    G_f.es['weight'] = weight
    return G_f
