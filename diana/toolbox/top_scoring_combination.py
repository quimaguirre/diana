#!/usr/bin/env python
import sys
import os
import imp
from toolbox import guild_utilities
from toolbox import functional_enrichment
#guild_utilities = imp.load_source('guild_utilities', '/home/boliva/PROJECTS/GUILD/guild_tutorial/toolbox/guild_utilities.py'')
#functional_enrichment = imp.load_source('functional_enrichment', '/home/boliva/PROJECTS/GUILD/guild_tutorial/toolbox/functional_enrichment.py')



def main():
    """
	Get nodes that are top scoring w.r.t. GUILD scores.
	Assumes that GUILD scores have been calculated already (i.e. python hello_world.py).
    """
    if len(sys.argv)>6:
     analysis_type=sys.argv[6]
     pval_threshold=float(sys.argv[5])
     top_threshold=float(sys.argv[4])
     scoring_folder = sys.argv[3]+"/"
     pvalue_file=sys.argv[2]
     data_dir=sys.argv[1]+"/"
    if len(sys.argv)==6:
     analysis_type=sys.argv[5]
     pval_threshold=float(sys.argv[4])
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

    print("\t\tINFO:\tTop threshold= %d%%\n"%(top_threshold))
    print("\t\tINFO:\tP-value threshold= %f\n"%(pval_threshold))

    # Set enrichment file
    enrichment_file = scoring_folder + "enrichment.txt"

    # Set uniprotentry profile file
    node_profile_file = scoring_folder + "node_profile.txt"

    #### IF FINALLY GENESPACE APPROACH WITH R IMPLEMENTED ####
    # # Set gene space file
    # gene_space_weighted = scoring_folder + "gene_space_weighted.tsv"
    # fgs = open(gene_space_weighted, 'w')
    # fgs.write("#Gene with GO attribute\tGene weight\n")
    
    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Get all nodes
    all_nodes = node_to_vals.keys()

    # Get top scoring, i.e. nodes that are in a given percentage of top score
    top_nodes = set()
    top_nodes_filtered = set()
    ntop=float(top_threshold)*len([x for x in node_to_vals.iteritems()])/100.0
    #print("\t\tINFO:\tTop percentage threshold= %d%% top scoring nodes are first %d nodes\n"%(top_threshold,ntop))  
    ii=0

    # If analysis_type == "node" or "functional", start writing a node profile
    if (analysis_type == "node"):
        f = open(node_profile_file, 'w')
        delim = " "
        f.write("#Id%sScore%sP-value\n" % (delim, delim))
    
    # Now, filter top_nodes by p-value <= (0.05). Include them in the profile if analysis_type == "node"
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
        score, pval = vals
        #### IF FINALLY GENESPACE APPROACH WITH R IMPLEMENTED ####
        # # Writing the genespace weighted
        # fgs.write("%s\t%s\n" %(node, int(score*1000)))
        if ii < ntop:
            top_nodes.add(node)
            last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last node, it is included as well
            if pval <= pval_threshold:
                top_nodes_filtered.add(node)
                if (analysis_type == "node"):
                    f.write("%s%s%f%s%s\n" % (node, delim, score, delim, pval))
            ii=ii+1
        else:
            if score == last_score: # Here, if a score above the threshold has the same score as the last, it is also recorded
                top_nodes.add(node)
                if pval <= pval_threshold:
                    top_nodes_filtered.add(node)
                    if (analysis_type == "node"):
                        f.write("%s%s%f%s%s\n" % (node, delim, score, delim, pval))
            #break

    #### IF FINALLY GENESPACE APPROACH WITH R IMPLEMENTED ####
    # # Closing the genespace weighted
    # fgs.close()

    if (analysis_type == "node"):
        f.close()

    print ("\t\tINFO:\t%d proteins in network, %d in top %d%% scoring proteins, and %d under the p-value threshold of %f\n" %
    (len(all_nodes), len(top_nodes), top_threshold, len(top_nodes_filtered), pval_threshold))


    #### IF FINALLY GENESPACE APPROACH WITH R IMPLEMENTED ####
#     query_file = scoring_folder + "query_file.txt"
#     results_file = scoring_folder + "R_functional_enrichment.txt"
#     fq = open(query_file, 'w')
#     i = 0
#     for node in top_nodes_filtered:
#         i+=1
#         if i < len(top_nodes_filtered):
#             fq.write("%s "%(node))
#         else:
#             fq.write("%s"%(node))
#     fq.close()
#     import os
#     command = '/soft/R/R-3.3.0/bin/Rscript toolbox/Funcassociate_client_example_geneSpace_weighted.R'
# #    command = '/soft/R/R-3.3.0/bin/Rscript toolbox/Funcassociate_client_example_geneSpace_weighted.R %s %s %s' % (query_file, gene_space_weighted, results_file)
#     os.system(command)


    # Get functional enrichment for these nodes
    if (analysis_type == "functional"):
        f = open(enrichment_file, 'w')
        functional_enrichment.check_functional_enrichment(list(top_nodes_filtered), list(all_nodes), "geneid",
                                                          f.write, species = "Homo sapiens")
        f.close()


    return

if __name__ == "__main__":
    main()


