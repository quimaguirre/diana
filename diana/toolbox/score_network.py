import guild_utilities
import network_utilities
import sys
from os.path import isdir, join

def main():
    """
    Scores a network
    """

    if len(sys.argv)==7:
        network_file=sys.argv[1]
        pvalue_file=sys.argv[2]
        seed_file=sys.argv[3]
        top_threshold=float(sys.argv[4])
        pval_threshold=float(sys.argv[5])
        scoring_folder=sys.argv[6] + '/'
    else:
        print("Please, introduce as first argument the network file, and as second argument the P-value file, and as third argument the folder where you want to introduce the output files")
        sys.exit(10)


    # Output files
    network_scored = join(scoring_folder + "network_scored.sif")
    linker_network_scored = join(scoring_folder + "network_linkers.sif")
    subnetwork_file = join(scoring_folder + "edge_profile.sif")
    linker_file = join(scoring_folder + "linkers.txt")
    linker_node_profile = join(scoring_folder + "linker_node_profile.txt")


    # Get GUILD scores and p-values
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Create a network from the network file
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)

    # Get the seeds from the seed file
    seeds=set()
    fs = open(seed_file,"r")
    for line in fs:
        seeds.add(line.strip().split()[0]) # Introduce the seeds into a set
    fs.close()

    # Get all nodes
    all_nodes = node_to_vals.keys()

    # Get network linkers induced by seeds
    ld_nodes=seeds.intersection(all_nodes) # Performs the intersection between the seeds and the nodes, obtaining the seeds inside the nodes. Introduces them into the linkers list
    node_ld, seeds_to_d = network_utilities.get_node_linker_degrees(g, seeds) # Finds if the nodes in the network are linkers or not and give them a value (+1 for the each seed that is neighbor of the node)
    # node_ld is a dictionary containing the nodes of the network, and a value corresponding to the degree of linker
    sorted_nodes_by_ld = sorted(node_ld.items(), key=lambda x: x[1],reverse=True) # This fantastic line of code produces a list of tuples (node, linker-degree) sorted by degree

    for node,ld in node_ld.iteritems(): # This loop adds the linkers found previously (degree > 0) into the ld_nodes set where the seeds were introduced before
        if ld>0:
            ld_nodes.add(node)
    g_link = network_utilities.get_subgraph(g, ld_nodes) # Now, a subnetwork is created using the linkers, obtaining a NETWORK of the LINKERS!!


    # Here we write a LINKER FILE containing: 'node', 'linker-degree', 'score', 'p-value'
    # And we write a LINKER NODE PROFILE containing: 'node' 'score' 'p-value' (only if degree > 0)
    linker2ld = {}
    fn = open(linker_file, 'w')
    fp = open(linker_node_profile, 'w')
    linker_node_profile_list = []
    for node in seeds_to_d:
        ld = seeds_to_d[node]
        linker2ld[node] = ld
        fp.write("%s %d %e\n" % (node, int(ld), 0)) # Add the seeds inside the linker profile with their corresponding degree (number of connections + 1)
    for node,ld in sorted_nodes_by_ld:
        fn.write("%s\t%d\t%f\t%e\n" % (node,ld,node_to_vals[node][0],node_to_vals[node][1]))
        if ld>0:
            linker_node_profile_list.append(node)
            linker2ld[str(node)] = ld
            fp.write("%s %d %e\n" % (node, int(ld), node_to_vals[node][1])) # Add the linker inside the linker profile with its degree (number of times connected with a seed)
    fn.close()
    fp.close()


    f  = open(network_scored, 'w')
    fl = open(linker_network_scored, 'w')

    # Here, the network and the network of linkers are scored
    # They are recorded into files
    for u, v in g.edges():
        score_u = node_to_vals[u][0]
        score_v = node_to_vals[v][0]
        score = (score_u + score_v) / 2
        f.write("%s\t%f\t%s\n" % (u, score, v))
        if (u,v) in g_link.edges():
            # Here we multiply the score of the node by the linker degree (if it is a linker)
            if str(u) in linker2ld:
                score_u = float(linker2ld[u])
            if str(v) in linker2ld:
                score_v = float(linker2ld[v])
            score = (score_u + score_v) / 2
            fl.write("%s\t%f\t%s\n" % (u, score, v))
    f.close()
    fl.close()


    # Get top scoring, i.e. nodes that are in a given percentage of top score
    top_nodes = set()
    top_nodes_filtered = set()
    ntop=float(top_threshold)*len([x for x in node_to_vals.iteritems()])/100.0
    ii=0
    last_score = ''

    # Now, filter top_nodes by p-value <= (0.05). Include them in the profile if entry_type == "uniprotentry"
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
        score, pval = vals
        if ii < ntop:
            top_nodes.add(node)
            last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last node, it is included as well
            if pval <= pval_threshold:
                top_nodes_filtered.add(node)
            ii=ii+1
        else:
            if score == last_score: # Here, if a score above the threshold has the same score as the last, it is also recorded
                top_nodes.add(node)
                if pval <= pval_threshold:
                    top_nodes_filtered.add(node)
            else:
                break

    # Get subnetwork induced by top scoring nodes
    g_sub = network_utilities.get_subgraph(g, top_nodes_filtered)

    fs  = open(subnetwork_file, 'w')
    # The subnetwork induced by the top scoring nodes is scored and recorded in a file
    for u, v in g_sub.edges():
        score_u = node_to_vals[u][0]
        score_v = node_to_vals[v][0]
        score = (score_u + score_v) / 2
        fs.write("%s\t%f\t%s\n" % (u, score, v))
    fs.close()


    return

if __name__ == "__main__":
    main()