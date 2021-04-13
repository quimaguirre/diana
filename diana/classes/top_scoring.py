import sys, os, re
import math
import networkx as nx

import diana.classes.network_translation as NT
import diana.classes.functional_analysis as FA


def node_top_scoring(node_to_score, threshold, threshold_type, output_file, verbose=False):
    """
    Creates a profile with the most relevant nodes using a percentage threshold of the score provided by the user.

    Parameters:
        @node_to_score:         Dictionary of all the nodes to their GUILD score
        @threshold:             Top percentage / number of nodes of the scores file in which we will cut to obtain the most relevant nodes
        @threshold_type:        Type of threshold used to define the top scoring nodes (it can either be 'percentage' or 'number_of_nodes')
        @output_file:           Resulting file which will contain the most relevant nodes.
    """

    # Get number of top scoring nodes
    if threshold_type == 'percentage':
        # By calculating the percentage over the total
        ntop=round(float(threshold)*len([x for x in node_to_score.items()])/100.0)
    elif threshold_type == 'number_of_nodes':
        # By exact number of nodes
        ntop=int(threshold)
    else:
        raise TOP.IncorrectThresholdType(wrong_top_type=threshold_type, top_types=['percentage', 'number_of_nodes'])

    # Now, write the profile with the top scoring nodes
    top_nodes = set()
    last_score = ''
    ii=0
    with open(output_file, 'w') as f:
        f.write('#Id\tScore\n')
        for node, score in sorted(node_to_score.items(),key=lambda x: x[1],reverse=True):
            if ii < ntop:
                top_nodes.add(node)
                last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last node, it is included as well                    
                f.write('{}\t{}\n'.format(node, score))
                ii=ii+1
            else:
                if score == last_score: # Here, if a score above the threshold has the same score as the last, it is also recorded
                    top_nodes.add(node)
                    f.write('{}\t{}\n'.format(node, score))
                else:
                    break

    if verbose:
        print('  DIANA INFO:\t{} file created.\n'.format(output_file))
    return


def edge_top_scoring(network_file, node_to_score, threshold, threshold_type, output_file, verbose=False):
    """
    Creates profiles with the most relevant edges using the thresholds provided by the user.
    The profile is created selecting the top most relevant nodes, and creating the subsequent scored subnetwork.

    Parameters:
        @network_file:          File containing the scored network
        @node_to_score:         Dictionary of all the nodes to their GUILD score
        @threshold:             Top percentage of the network in which we will cut to obtain the most relevant edges
        @threshold_type:        Type of threshold used to define the top scoring nodes (it can either be 'percentage' or 'number_of_nodes')
        @output_file:           Resulting file which will contain the most relevant edges
    """

    # Create a network from the network file
    g = create_network_from_sif_file(network_file, use_edge_data=True)

    # Get number of top scoring nodes
    if threshold_type == 'percentage':
        # By calculating the percentage over the total
        ntop=round(float(threshold)*len([x for x in node_to_score.items()])/100.0)
    elif threshold_type == 'number_of_nodes':
        # By exact number of nodes
        ntop=int(threshold)
    else:
        raise TOP.IncorrectThresholdType(wrong_top_type=threshold_type, top_types=['percentage', 'number_of_nodes'])

    top_nodes = set()
    ii=0
    last_score = ''

    # Filter the nodes by top scoring nodes
    for node, score in sorted(node_to_score.items(),key=lambda x: x[1],reverse=True):
        if ii < ntop:
            top_nodes.add(node)
            last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last node, it is included as well
            ii=ii+1
        else:
            if score == last_score: # Here, if a score above the threshold has the same score as the last, it is also recorded
                top_nodes.add(node)
            else:
                break

    # Get subnetwork induced by top scoring nodes
    g_sub = get_subgraph(g, top_nodes)

    with open(output_file, 'w') as f:

        # The subnetwork induced by the top scoring nodes is scored and recorded in a file
        for u, v in g_sub.edges():
            score_u = node_to_score[u]
            score_v = node_to_score[v]
            score = (score_u + score_v) / 2
            f.write('{}\t{:f}\t{}\n'.format(u, score, v))

    if verbose:
        print('  DIANA INFO:\t{} file created.\n'.format(output_file))
    return


def functional_top_scoring(top_geneids, type_correction, associations_file, output_file, verbose=False):
    """
    Creates profiles with the most relevant functions using the thresholds provided by the user.
    The profile is created selecting the top most relevant nodes, and computing a functional enrichment analysis of them.

    Parameters:
        @top_geneids:           A list containing the selected Entrez gene IDs to do the enrichment analysis
        @type_correction:       Type of p-value multiple test correction (fda_bh or bonferroni)
        @output_file:           Resulting file which will contain the functions enriched.
        @associations_file:     File containing the function-gene associations
    """
    FA.calculate_functional_enrichment_profile(top_geneids=top_geneids, type_correction=type_correction, associations_file=associations_file, output_file=output_file)
    if verbose:
        print('  DIANA INFO:\t{} file created.\n'.format(output_file))
    return


def score_network(network, node_to_score, output_file):
    """
    Scores a network
    """

    # Get all nodes
    all_nodes = node_to_score.keys()

    # Here, the network is scored and stored in a file
    with open(output_file, 'w') as f:
        for u, v in network.edges():
            score_u = node_to_score[u]
            score_v = node_to_score[v]
            score = (float(score_u) + float(score_v)) / 2
            f.write('{}\t{:f}\t{}\n'.format(u, score, v))
 
    return


def create_network_from_sif_file(network_file_in_sif, use_edge_data = False, delim = None, include_unconnected=True):
    """
    Creates a NetworkX graph object from a sif file
    """
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type = use_edge_data, delim = delim, data_to_float = False )
    g=nx.Graph()
    if include_unconnected:
        g.add_nodes_from(setNode)
    if use_edge_data:
        for e,w in dictEdge.items():
            u,v = e
            g.add_edge(u, v, w=w)
    else:
        g.add_edges_from(setEdge)
    return g


def get_nodes_and_edges_from_sif_file(file_name, store_edge_type = False, delim=None, data_to_float=True):
    """
    Parse sif file into node and edge sets and dictionaries
    returns setNode, setEdge, dictNode, dictEdge
    store_edge_type: if True, dictEdge[(u,v)] = edge_value
    delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    f=open(file_name)
    for line in f:
        if delim is None:
            words = line[:-1].split()
        else:
            words = line[:-1].split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            if data_to_float:
                score = float(words[1])
            else:
                score = words[1]
            dictNode[id1] = score
        elif len(words) == 3: 
            id2 = words[2]
            setNode.add(id2)
            setEdge.add((id1, id2))
            if store_edge_type:
                if data_to_float:
                    dictEdge[(id1, id2)] = float(words[1])
                else:
                    dictEdge[(id1, id2)] = words[1]
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    return setNode, setEdge, dictNode, dictEdge


def get_subgraph(G, nodes):
    """
    NetworkX subgraph method wrapper
    """
    return G.subgraph(nodes)


class IncorrectThresholdType(NameError):
    """
    Subclass exception of the NameError which raises when a network format is
    not provided correctly
    """
    def __init__(self, wrong_top_type, top_types):
        self.wrong_top_type = wrong_top_type
        self.top_types = top_types

    def __str__(self):
        return 'The threshold type {} is not valid. It must be one of the following threshold types: {}'.format(self.wrong_top_type, self.top_types)
