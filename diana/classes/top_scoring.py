import sys, os, re
import math
import networkx as nx

import network_translation as NT
from goatools.go_enrichment import GOEnrichmentStudy


def node_top_scoring(node_to_vals, threshold, output_file):
    """
    Creates a profile with the most relevant nodes using a percentage threshold of the score provided by the user.

    Parameters:
        @node_to_vals:          Dictionary of all the nodes to a tuple containing (score, pvalue)
        @threshold:             Top percentage of the pvalue_file in which we will cut to obtain the most relevant nodes
        @output_file:           Resulting file which will contain the most relevant nodes
    """

    # Get top scoring, i.e. nodes that are in a given percentage of top score
    ntop=float(threshold)*len([x for x in node_to_vals.iteritems()])/100.0

    top_nodes = set()
    last_score = ''
    ii=0

    f = open(output_file, 'w')
    f.write('#Id Score P-value\n')

    # Now, write the profile with the top scoring nodes
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
        score, pval = vals
        if ii < ntop:
            top_nodes.add(node)
            last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last node, it is included as well                    
            f.write('{} {} {}\n'.format(node, score, pval))
            ii=ii+1
        else:
            if score == last_score: # Here, if a score above the threshold has the same score as the last, it is also recorded
                top_nodes.add(node)
                f.write('{} {} {}\n'.format(node, score, pval))
            else:
                break
    f.close()

    print('  DIANA INFO:\t{} file created.\n'.format(output_file))

    return


def edge_top_scoring(network_file, node_to_vals, threshold, output_file):
    """
    Creates profiles with the most relevant edges using the thresholds provided by the user.
    The profile is created selecting the top most relevant nodes, and creating the subsequent scored subnetwork.

    Parameters:
        @network_file:          File containing the scored network
        @node_to_vals:          Dictionary of all the nodes to a tuple containing (score, pvalue)
        @threshold:             Top percentage of the network in which we will cut to obtain the most relevant edges
        @output_file:           Resulting file which will contain the most relevant edges
    """

    # Create a network from the network file
    g = create_network_from_sif_file(network_file, use_edge_data=True)

    # Get top scoring, i.e. nodes that are in a given percentage of top score
    ntop=float(threshold)*len([x for x in node_to_vals.iteritems()])/100.0

    top_nodes = set()
    ii=0
    last_score = ''

    # Filter the nodes by top scoring nodes
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
        score, pval = vals
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

    f = open(output_file, 'w')

    # The subnetwork induced by the top scoring nodes is scored and recorded in a file
    for u, v in g_sub.edges():
        score_u = node_to_vals[u][0]
        score_v = node_to_vals[v][0]
        score = (score_u + score_v) / 2
        f.write('{}\t{:f}\t{}\n'.format(u, score, v))

    f.close()

    print('  DIANA INFO:\t{} file created.\n'.format(output_file))

    return


def functional_top_scoring(obodag, geneid2gos, all_nodes, top_nodes, output_file, temp_file):
    """
    Creates profiles with the most relevant functions using the thresholds provided by the user.
    The profile is created selecting the top most relevant nodes, and computing a functional enrichment analysis of them.

    Parameters:
        @obodag:                Dictionary containing the Gene Ontology
        @geneid2gos:            Dictionary containing the equivalences from geneID to GOs
        @all_nodes:             A list containing all the nodes of the network
        @top_nodes:             A list containing the selected nodes to do the enrichment analysis
        @output_file:           Resulting file which will contain the most relevant edges
        @temp_file:             A file where the functional enrichment will be temporary calculated
    """

    # Get background genes
    all_nodes = [ int(x) for x in all_nodes ]

    # Get top genes
    top_nodes = [ int(x) for x in top_nodes ]

    calculate_functional_enrichment_profile(obodag, geneid2gos, top_nodes, all_nodes, temp_file, output_file)
    #GOA.calculate_functional_enrichment_profile(obodag, geneid2gos, top_nodes, all_nodes, temp_file, output_file)

    print("  DIANA INFO:\t%s file created.\n" %(output_file))

    return


def score_network(network_file, node_to_vals, output_file):
    """
    Scores a network
    """

    # Create a network from the network file
    g = create_network_from_sif_file(network_file, use_edge_data=True)

    # Get all nodes
    all_nodes = node_to_vals.keys()

    f  = open(output_file, 'w')
    # Here, the network is scored and stored in a file
    for u, v in g.edges():
        score_u = node_to_vals[u][0]
        score_v = node_to_vals[v][0]
        score = (float(score_u) + float(score_v)) / 2
        f.write('{}\t{:f}\t{}\n'.format(u, score, v))
    f.close()
 

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
        for e,w in dictEdge.iteritems():
            u,v = e
            g.add_edge(u,v,{'w':w})
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


def calculate_functional_enrichment_profile(obodag, geneid2go, top_nodes, all_nodes, temp_file, output_file):

    print("{N:,} annotated human genes".format(N=len(geneid2go)))

    # Create a dictionary of GOs to their geneids (it will be used for the Log of Odds calculus)
    go2geneids = {}
    for geneid in geneid2go:
        for function in geneid2go[geneid]:
            go2geneids.setdefault(function, set())
            go2geneids[function].add(geneid)

    # Define the GOEnrichmentStudy object
    goeaobj = GOEnrichmentStudy(
            all_nodes, # List of background genes
            geneid2go, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.9999, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method

    # Writing the temorary file (without the Log of odds ratio)
    goea_results_all = goeaobj.run_study(top_nodes)
    goeaobj.wr_txt(temp_file, goea_results_all)

    results = []
    ft = open(temp_file, 'r')
    k = len(top_nodes) # num genes in top
    t = len(all_nodes) # num genes in network

    # Calculating the Log of odds ratio
    for line in ft:
        words = line.strip().split()
        go = words[0]
        type_func = words[1]
        pvalue = words[2]
        name = ' '.join(words[4:])
        q = int(words[3]) # num genes assoc with function in top
        m = len(go2geneids[go]) # num genes assoc with func in network
        log_of_odds = calculate_log_of_odds_ratio(q, k, m, t)
        results.append([q, m, log_of_odds, '-', pvalue, go, name, type_func])

    ft.close()

    # Writing the output file
    fo = open(output_file, 'w')
    fo.write('# num of genes\tnum of total genes\tLog of odds ratio\tP-value\tAdjusted p-value\tGO term ID\tGO term name\tGO type\n')

    for line in sorted(results, key=lambda x: x[2], reverse=True):
        if line[0] <= 0: # Skip the functions that do not have at leaast 1 top gene associated (the log of odds ratio is -infinite!!)
            continue
        if line[2] < 0: # Skip the log of odds ratio that are negative (we will only use functions with log of odds ratio from 0 to inf)
            continue
        if line[7].upper() == 'BP' or line[7].upper() == 'MF':
            new_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7])
            fo.write(new_line)

    fo.close()

    # Removing the temporary file
    command = "rm {}".format(temp_file)
    os.system(command)

    return


def calculate_log_of_odds_ratio(q, k, m, t):
    """Calculates the log of the odds ratio"""
    # print("q: {}".format(q))
    # print("k: {}".format(k))
    # print("m: {}".format(m))
    # print("t: {}".format(t))
    odds_ratio = ( (float(q)/float(k)) / (float(m)/float(t)) )
    #odds_ratio = ( (float(q)/(float(k)-float(q))) / (float(m)/((float(t)-float(m)))) )
    if odds_ratio == 0:
        return -float('inf')
    else:
        return float(math.log(odds_ratio, 2))
