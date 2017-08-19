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

    if len(sys.argv)>6:
     pval_threshold=float(sys.argv[6])
     top_threshold=float(sys.argv[5])
     scoring_folder = sys.argv[4]+"/"
     network_file = sys.argv[3]
     pvalue_file=sys.argv[2]
     data_dir=sys.argv[1]+"/"
    if len(sys.argv)==6:
     pval_threshold=float(sys.argv[5])
     top_threshold=float(sys.argv[4])
     scoring_folder = sys.argv[3]+"/"
     pvalue_file=sys.argv[2]
     data_dir=sys.argv[1]+"/"
    if len(sys.argv)==5:
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
    # else:
    #  pval_threshold=0.05
    #  data_dir = "./"
    #  scoring_folder =  data_dir
    #  pvalue_file=scoring_folder + "output_scores.sif.netcombo.pval"

    complete_network_scored = data_dir + "complete_network_scored.sif"
    edges_profile = scoring_folder + "edges_profile.sif"
    print("P-value threshold= %e\n"%(pval_threshold))

    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Get top scoring, i.e. nodes that are in a given percentage of top score
    top_nodes = set()
    top_nodes_filtered = set()
    ntop=float(top_threshold)*len([x for x in node_to_vals.iteritems()])/100.0
    ii=0

    # Now, filter top_nodes by p-value <= (0.05). Include them in the profile if entry_type == "uniprotentry"
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
        score, pval = vals
        if ii < ntop:
            top_nodes.add(node)
            if pval <= pval_threshold:
                top_nodes_filtered.add(node)
            ii=ii+1
        else:
            break

    # Load interaction network
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)

    # # Output the complete network along with the inverted p-value scores (z-scores) calculated for edges
    # f = open(complete_network_scored, 'w')
    # for u, v in g.edges():
    #     zscore_u = stat_utilities.convert_p_values_to_z_scores([node_to_vals[u][1]])[0]
    #     zscore_v = stat_utilities.convert_p_values_to_z_scores([node_to_vals[v][1]])[0]
    #     score = (zscore_u + zscore_v) / 2
    #     f.write("%s\t%f\t%s\n" % (u, score, v))
    # f.close()

    # Get subnetwork induced by top scoring nodes
    g_sub = network_utilities.get_subgraph(g, top_nodes_filtered)

    # Output subnetwork along with the inverted p-value scores (z-scores) calculated for edges
    f = open(edges_profile, 'w')
    for u, v in g_sub.edges():
        zscore_u = stat_utilities.convert_p_values_to_z_scores([node_to_vals[u][1]])[0]
        zscore_v = stat_utilities.convert_p_values_to_z_scores([node_to_vals[v][1]])[0]
        score = (zscore_u + zscore_v) / 2
        f.write("%s\t%f\t%s\n" % (u, score, v))
    f.close()

    return

if __name__ == "__main__":
    main()

