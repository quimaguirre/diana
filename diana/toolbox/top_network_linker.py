#!/usr/bin/env python
import sys
import os
import imp
guild_utilities = imp.load_source('guild_utilities','.')
network_utilities = imp.load_source('network_utilities', '.')
stat_utilities = imp.load_source('stat_utilities', '.')
from toolbox import guild_utilities
from toolbox import network_utilities
from toolbox import stat_utilities

def main():
    """
	Get nodes that are top scoring w.r.t. GUILD scores.
	Assumes that GUILD scores have been calculated already (i.e. python hello_world.py).
    """

    # Set the seed and network files

    if len(sys.argv)>4:
     pval_threshold=float(sys.argv[4])
     scoring_folder = sys.argv[3]+"/"
     pvalue_file=sys.argv[2]
     data_dir=sys.argv[1]+"/"
    elif len(sys.argv)==4:
     pval_threshold=0.05
     scoring_folder = sys.argv[3]+"/"
     data_dir=sys.argv[1]+"/"
     pvalue_file=sys.argv[2]
    elif len(sys.argv)==3:
     pval_threshold=0.05
     data_dir=sys.argv[1]+"/"
     pvalue_file=sys.argv[2]
     scoring_folder =  data_dir
    elif len(sys.argv)==2:
     pval_threshold=0.05
     data_dir=sys.argv[1]+"/"
     scoring_folder =  data_dir
     pvalue_file=scoring_folder + "output_scores.sif.netcombo.pval"
    else:
     pval_threshold=0.05
     data_dir = "./"
     scoring_folder =  data_dir
     pvalue_file=scoring_folder + "output_scores.sif.netcombo.pval"

    seed_file = data_dir + "seeds.txt"
    network_file = data_dir + "interactions.sif"
    enrichment_file = data_dir + "enrichment.txt"
    subnetwork_file = scoring_folder + "subnetwork.sif"
    subnetwork_linker_file = scoring_folder + "subnetwork_linkers.sif"
    linker_file = scoring_folder + "linkers.txt"


    print("P-value threshold= %e\n"%(pval_threshold))

    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Get top scoring, i.e. nodes that have p-value <= 0.05
    top_nodes = set()
    for node, vals in node_to_vals.iteritems():
	score, pval = vals
	if pval <= pval_threshold:
	    top_nodes.add(node)

    # Load interaction network
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)

    # Get subnetwork induced by top scoring nodes
    g_sub = network_utilities.get_subgraph(g, top_nodes)

    # Get subnetwork linkers induced by only seeds
    seeds=set()
    fs = open(seed_file,"r")
    for line in fs:
      seeds.add(line.strip().split()[0]) # Introduce the seeds into a set
    fs.close()
    ld_nodes=seeds.intersection(top_nodes) # Performs the intersection between the seeds and the top nodes, obtaining the seeds inside the top nodes
    node_ld = network_utilities.get_node_linker_degrees(g_sub, seeds) # Finds if the nodes in the network are linkers or not and give them a value (+1 for the each seed that is neighbor of the node)
    # node_ld is a dictionary containing the nodes of the network, and a value corresponding to the degree of linker
    sorted_nodes_by_ld = sorted(node_ld.items(), key=lambda x: x[1],reverse=True) # This fantastic line of code produces a list of tuples (node, linker-degree) sorted by degree
    for node,ld in node_ld.iteritems(): # This loop adds the linkers found previously (degree > 0) into the ld_nodes set where the seeds were introduced before
      if ld>0:
       ld_nodes.add(node)
    g_link = network_utilities.get_subgraph(g_sub, ld_nodes) # Now, a subnetwork of the subnetwork is created using the linkers, obtaining a NETWORK of the LINKERS!!

    # Output subnetwork along with the inverted p-value scores (z-scores) calculated for edges
    f = open(subnetwork_file, 'w')
    fl = open(subnetwork_linker_file, 'w')
    # Here, the INITIAL SUBNETWORK and the LINKERS NETWORK are displayed and SCORED
    for u, v in g_sub.edges():
	zscore_u = stat_utilities.convert_p_values_to_z_scores([node_to_vals[u][1]])[0]
	zscore_v = stat_utilities.convert_p_values_to_z_scores([node_to_vals[v][1]])[0]
	score = (zscore_u + zscore_v) / 2
	f.write("%s\t%f\t%s\n" % (u, score, v))
        if (u,v) in g_link.edges():
          fl.write("%s\t%f\t%s\n" % (u, score, v))
    f.close()
    fl.close()
    # Here we write a LINKER FILE containing: 'node', 'linker-degree', 'score', 'p-value'
    fn = open(linker_file, 'w')
    for node,ld in sorted_nodes_by_ld:
      fn.write("%s\t%d\t%f\t%e\n" % (node,ld,node_to_vals[node][0],node_to_vals[node][1]))
    fn.close()

    return

if __name__ == "__main__":
    main()

