import networkx as nx
import sys, os, re

def main():
    """
    python /home/quim/PHD/Projects/DIANA/diana/scripts/create_network_without_loops.py
    """

    input_network = '/home/quim/PHD/Projects/DIANA/diana/data/network_cheng.txt'
    output_network = '/home/quim/PHD/Projects/DIANA/diana/data/network_cheng_without_loops.sif'

    # Read the input network
    network = nx.Graph()
    network_without_loops = nx.Graph()
    num_loops = 0
    with open(input_network, 'r') as net_fd:
        for line in net_fd:
            fields = line.strip().split('\t')
            node1 = fields[0]
            node2 = fields[1]
            network.add_edge(node1, node2)
            if node1 != node2:
                network_without_loops.add_edge(node1, node2)
            else:
                num_loops+=1
    print('\nTHE NETWORK HAS {} EDGES AND {} NODES (GENEIDS)'.format(network.number_of_edges(), network.number_of_nodes()))
    print('\nTHE NETWORK FILTERING LOOPS HAS {} EDGES AND {} NODES (GENEIDS)'.format(network_without_loops.number_of_edges(), network_without_loops.number_of_nodes()))
    print('\nNUMBER OF LOOPS: {}'.format(num_loops))

    # Write the new network in sif format and without self-interactions
    with open(output_network, 'w') as out_fd:
        for u,v in network_without_loops.edges():
            out_fd.write('{}\t1\t{}\n'.format(u,v))

    return

def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


if  __name__ == "__main__":
    main()