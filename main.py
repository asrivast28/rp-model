import argparse
import numpy as np

import network
import analysis


def parse_args():
    """
    @brief  Parses command line arguments. 
    
    @return  Structure containing parsed arguments. 
    """
    parser = argparse.ArgumentParser(description = 'Create a network as per RP-model and analyze it.')
    parser.add_argument('--file', '-f', metavar = 'FILE', type = str, help = 'name of the file from which the network is read')
    parser.add_argument('--vertices', '-v', metavar = 'V', type = int, default = 1000, help = 'total number of vertices')
    parser.add_argument('--sources', '-s', metavar = 'S', type = int, help = 'number of sources')
    parser.add_argument('--targets', '-t', metavar = 'T', type = int, help = 'number of targets')
    parser.add_argument('--degreedist', '-d', metavar = 'D', type = int, default = 1, help = 'degree distribution: 1- 1 + Poisson(1), 2- Constant(2), 3- 1 + Poisson(3)')
    parser.add_argument('--alpha', '-a', metavar = 'A', type = float, default = 0.0, help = 'preference for reuse')
    parser.add_argument('--tau', '-u', metavar = 'U', type = float, default = 0.9, help = 'path coverage threshold')
    parser.add_argument('--seed', '-n', metavar = 'N', type = int, default = 0, help = 'seed for the PRNG')

    args = parser.parse_args()
    # In-degree distribution for the network
    if args.degreedist == 1:
        args.d_in = lambda : 1 + np.random.poisson(1.0)
    elif args.degreedist == 2:
        args.d_in = lambda : 2
    elif args.degreedist == 3:
        args.d_in = lambda : 1 + np.random.poisson(3.0)
    else:
        raise argparse.ArgumentTypeError, 'degree distribution should be one of {1, 2, 3}'

    return args


def main():
    """
    @brief  Main function. 
    """
    # Parse the command line arguments.
    args = parse_args()
    # Seed the PRNG
    np.random.seed(args.seed)

    # Read the network if a file name is provided
    # Otherwise, create a network using the RP-model 
    if args.file:
        G = network.read(args.file)
    else:
        intermediates = args.vertices - (args.sources + args.targets)
        G = network.rp_model(args.sources, intermediates, args.targets, args.alpha, args.d_in)

    # Get the core vertices for the network
    C = analysis.core_vertices(G, args.tau)

    # Get the flattened network corresponding to the original network
    G_f = network.flatten(G)
    # Get the core vertices for the flattened network
    C_f = analysis.core_vertices(G_f, args.tau)

    # H-score of the original network
    H = 1 - (float(len(C)) / len(C_f))
    print '%d\t%f'%(len(C), H)

if __name__ == '__main__':
    main()