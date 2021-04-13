import networkx as nx
import sys, os, re

def main():
    """
    python /home/quim/PHD/Projects/DIANA/diana/scripts/create_connected_network.py
    """

    input_network = '/home/quim/PHD/Projects/DIANA/diana/data/network_cheng.txt'
    output_network = '/home/quim/PHD/Projects/DIANA/diana/data/network_cheng_lcc.sif'
    output_network_without_loops = '/home/quim/PHD/Projects/DIANA/diana/data/network_cheng_lcc_without_loops.sif'

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

    # Network LCC
    components = nx.connected_components(network)
    lcc_nodes = set()
    for component in components:
        if len(component) > len(lcc_nodes):
            lcc_nodes = component

    lcc = network.subgraph(lcc_nodes)

    # Network without loops LCC
    components = nx.connected_components(network_without_loops)
    lcc_without_loops_nodes = set()
    for component in components:
        if len(component) > len(lcc_without_loops_nodes):
            lcc_without_loops_nodes = component

    lcc_without_loops = network_without_loops.subgraph(lcc_without_loops_nodes)

    print('\nTHE NETWORK HAS {} EDGES AND {} NODES (GENEIDS)'.format(network.number_of_edges(), network.number_of_nodes()))
    print('\nNUMBER OF LOOPS: {}'.format(num_loops))
    print('\nTHE LONGEST CONNECTED COMPONENT HAS {} EDGES AND {} NODES (GENEIDS)'.format(lcc.number_of_edges(), lcc.number_of_nodes()))
    print('\nTHE LONGEST CONNECTED COMPONENT WITHOUT LOOPS HAS {} EDGES AND {} NODES (GENEIDS)'.format(lcc_without_loops.number_of_edges(), lcc_without_loops.number_of_nodes()))

    # Write the new network in sif format
    with open(output_network, 'w') as out_fd:
        for u,v in lcc.edges():
            out_fd.write('{}\t1\t{}\n'.format(u,v))

    # Write the new network in sif format and without loops
    with open(output_network_without_loops, 'w') as out_fd:
        for u,v in lcc_without_loops.edges():
            out_fd.write('{}\t1\t{}\n'.format(u,v))

    return


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