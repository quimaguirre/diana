#!/usr/bin/env python
import sys
import os
import imp
sys.path.append("/home/boliva/PROJECTS/GUILD/guild_tutorial")
guild_utilities = imp.load_source('guild_utilities.py','/home/boliva/PROJECTS/GUILD/guild_tutorial/toolbox/guild_utilities.py')
network_utilities = imp.load_source('network_utilities', '/home/boliva/PROJECTS/GUILD/guild_tutorial/toolbox/network_utilities.py')
stat_utilities = imp.load_source('stat_utilities', '/home/boliva/PROJECTS/GUILD/guild_tutorial/toolbox/stat_utilities.py')

def main():
    """
	Get nodes that are top scoring w.r.t. GUILD scores but use only seeds to generate the subgraph
	Assumes that GUILD scores have been calculated already (i.e. python hello_world.py).
    """

    # Set the seed and network files

    if len(sys.argv)>4:
     top_threshold=float(sys.argv[4])
     scoring_folder = sys.argv[3]+"/"
     pvalue_file=sys.argv[2]
     data_dir=sys.argv[1]+"/"
    elif len(sys.argv)==4:
     top_threshold=1.0
     scoring_folder = sys.argv[3]+"/"
     data_dir=sys.argv[1]+"/"
     pvalue_file=sys.argv[2]
    elif len(sys.argv)==3:
     top_threshold=1.0
     data_dir=sys.argv[1]+"/"
     pvalue_file=sys.argv[2]
     scoring_folder =  data_dir
    elif len(sys.argv)==2:
     top_threshold=1.0
     data_dir=sys.argv[1]+"/"
     scoring_folder =  data_dir
     pvalue_file=scoring_folder + "output_scores.sif.netcombo.pval"
    else:
     top_threshold=1.0
     data_dir = "./"
     scoring_folder =  data_dir
     pvalue_file=scoring_folder + "output_scores.sif.netcombo.pval"

    seed_file = data_dir + "seeds.txt"
    network_file = data_dir + "interactions.sif"
    enrichment_file = data_dir + "enrichment.txt"
    subnetwork_file = scoring_folder + "subnetwork.sif"
    subnetwork_linker_file = scoring_folder + "subnetwork_linkers.sif"
    linker_file = scoring_folder + "linkers.txt"

    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Get top scoring, i.e. nodes that have p-value <= 0.05
    top_nodes = set()
    #for node, vals in node_to_vals.iteritems():
    ntop=float(top_threshold)*len([x for x in node_to_vals.iteritems()])/100.0
    print("Top percentage threshold= %e are first %e nodes\n"%(top_threshold,ntop))  
    ii=0
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
	score, pval = vals
	#if pval <= top_threshold:
        if ii < ntop:
	    top_nodes.add(node)
            print("Add Node %s Score %e P-value: %e\n"%(node,score,pval))
            ii=ii+1

    # Load interaction network
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)

    # Get subnetwork induced by top scoring nodes
    g_sub = network_utilities.get_subgraph(g, top_nodes)

    # Get subnetwork linkers induced by only seeds
    seeds=set()
    fs = open(seed_file,"r")
    for line in fs:
      seeds.add(line.strip().split()[0])
    fs.close()
    ld_nodes=seeds.intersection(top_nodes)
    node_ld = network_utilities.get_node_linker_degrees(g_sub, seeds)
    sorted_nodes_by_ld = sorted(node_ld.items(), key=lambda x: x[1],reverse=True)
    
    for node,ld in node_ld.iteritems():
      if ld>0:
       ld_nodes.add(node)
    g_link = network_utilities.get_subgraph(g_sub, ld_nodes)

    # Output subnetwork along with the inverted p-value scores (z-scores) calculated for edges
    f  = open(subnetwork_file, 'w')
    fl = open(subnetwork_linker_file, 'w')
    for u, v in g_sub.edges():
	zscore_u = stat_utilities.convert_p_values_to_z_scores([node_to_vals[u][1]])[0]
	zscore_v = stat_utilities.convert_p_values_to_z_scores([node_to_vals[v][1]])[0]
	score = (zscore_u + zscore_v) / 2
	f.write("%s\t%f\t%s\n" % (u, score, v))
        if (u,v) in g_link.edges():
          fl.write("%s\t%f\t%s\n" % (u, score, v))
    f.close()
    fl.close()
    fn = open(linker_file, 'w')
    for node,ld in sorted_nodes_by_ld:
      fn.write("%s\t%d\t%f\t%e\n" % (node,ld,node_to_vals[node][0],node_to_vals[node][1]))
    fn.close()

    return
    

if __name__ == "__main__":
    main()

