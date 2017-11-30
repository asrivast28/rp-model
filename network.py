import numpy as np

import utils


def read(filename):
    """
    @brief  Reads a network from the file in edge-list format.

    @param filename  Name of the file to be read.

    @return  np.array representation of the network's weighted adjacency matrix.
    """
    # initialize a graph
    # with open('centrality_test,txt', 'rb') as f:
    with open(filename, 'rb') as elf:
        G = ig.Graph.Read_Edgelist(elf)
    vertex = np.arange(G.vcount(), dtype=utils.datatype(G.vcount()))
    source = np.vectorize(lambda v : G.degree(v, mode=ig.IN) == 0, otypes=[np.bool])(vertex)
    target = np.vectorize(lambda v : G.degree(v, mode=ig.OUT) == 0, otypes=[np.bool])(vertex)
    return G, source, target

def rp_model(S, M, T, alpha, d_in, out, datatype=np.float64):
    """
    @brief  Creates a network using the RP-model.

    @param S         Number of source vertices.
    @param M         Number of intermediate vertices.
    @param T         Number of target vertices.
    @param alpha     Preference for reuse.
    @param d_in      Function that generates in-degrees.
    @param out       Prefix of the files to which the network is written.
    @param datatype  np.dtype provided as a hint for storage.

    @return  np.array representation of the network's weighted adjacency matrix.
    """
    # Add all the vertices to the network
    V = S + M + T
    vertextype = utils.datatype(V)
    # Create a directed network
    G = np.zeros((V, V), dtype=datatype)

    # Source ranks array
    source_ranks = np.arange(S, dtype=vertextype)
    # Intermediate ranks array
    inter_ranks = np.zeros(M, dtype=vertextype)
    inter_indices = np.arange(M)

    vertex = np.arange(V, dtype=vertextype)
    source = (vertex < S)
    target = (vertex >= (S + M))

    k = float(d_in())
    G.T[S][source] = k/S

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
        k = d_in()
        if k > S + m:
            k = S + m
        # Pick unique source vertices and add incoming edges from them
        for u in np.random.choice(S + M, size=k, replace=False, p=probabilities):
            G[u][S + m] = 1.0

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
        k = d_in()
        if k > S + M:
            k = S + M
        # Pick unique source vertices and add incoming edges from them
        for u in np.random.choice(S + M, size=k, replace=False, p=probabilities):
            G[u][S + M + t] = 1.0

    if out is not None:
        with open('%s_links.txt'%out, 'wb') as elf:
            for u in vertex:
                for v in np.nonzero(G[u])[0]:
                    elf.write('%d\t%d\t%f\n'%(u, v, G[u][v]))
        with open('%s_sources.txt'%out, 'wb') as sf:
            sf.write('\n'.join(str(s) for s in np.where(source)[0]))
        with open('%s_targets.txt'%out, 'wb') as tf:
            tf.write('\n'.join(str(t) for t in np.where(target)[0]))
    return G, source, target

def flatten(G, G_T, source, target, datatype=np.float64):
    """
    @brief  Flattens the given dependency network.

    @param G         np.array representation of the network's weighted adjacency matrix.
    @param G_T       Transpose of G.
    @param source    np.array of type bool with 1 for every source vertex.
    @param target    np.array of type bool with 1 for every target vertex.
    @param datatype  np.dtype provided as a hint for storage.

    @return  np.array representation of the flattened network's weighted adjacency matrix.
    """
    # Create a weighted directed network
    G_f = np.zeros(G.shape, dtype=datatype)

    # Create the flat dependency network
    targets = np.where(target)[0]
    for s in np.where(source)[0]:
        G_f[s][s] = 1
        utils.count_simple_paths(G, G_T, np.where(G[s])[0], G_f[s])
        G_f[s][~target] = 0
    return G_f
