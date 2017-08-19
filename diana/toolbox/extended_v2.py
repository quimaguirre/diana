from toolbox import functional_enrichment
from toolbox import guild_utilities
from toolbox import network_utilities
import toolbox.comparison as comp
import toolbox.run_goatools as GOA


def extended_analysis_comparison(profiles_list_node, profiles_list_edge, profiles_list_functional, pvalue_files, scored_networks, top_thresholds_list, output_file):
    """
    Performs the comparison between all the different profiles obtained from the extended analysis
    """
    f = open(output_file, 'w')
    f.write("##### EXTENDED ANALYSIS #####\n")
    extended = {}

    # We will perform two extended analyses, one for the node profile and the other for the edge profile
    for type_of_profile in ('node', 'edge', 'functional'):

        f.write("\n# %s PROFILE\n"%(type_of_profile.upper()))

        if type_of_profile == 'node':
            profiles_drug1 = profiles_list_node[0]
            profiles_drug2 = profiles_list_node[1]
            pvalue_file_drug1 = pvalue_files[0]
            pvalue_file_drug2 = pvalue_files[1]
            # Create a list of lists containing the profile_drug1_top_x, the profile_drug2_top_x and the top x
            profiles_comparison = list( zip(profiles_list_node[0], profiles_list_node[1], top_thresholds_list) )
        elif type_of_profile == 'edge':
            profiles_drug1 = profiles_list_edge[0]
            profiles_drug2 = profiles_list_edge[1]
            network_scored_drug1 = scored_networks[0]
            network_scored_drug2 = scored_networks[1]
            # Create a list of lists containing the profile_drug1_top_x, the profile_drug2_top_x and the top x
            profiles_comparison = list( zip(profiles_list_edge[0], profiles_list_edge[1], top_thresholds_list) )
        elif type_of_profile == 'functional':
            profiles_drug1 = profiles_list_functional[0]
            profiles_drug2 = profiles_list_functional[1]
            # Create a list of lists containing the profile_drug1_top_x, the profile_drug2_top_x and the top x
            profiles_comparison = list( zip(profiles_list_functional[0], profiles_list_functional[1], top_thresholds_list) )
        data_frame = []

        for (file1, file2, top) in profiles_comparison:

            # f.write("\n## TOP THRESHOLD: %s%% ##\n" %(top))

            if type_of_profile == 'node':
                top_node_to_vals_drug1 = guild_utilities.get_values_from_pvalue_file(file1)
                top_node_to_vals_drug2 = guild_utilities.get_values_from_pvalue_file(file2)
                node_to_vals_drug1 = guild_utilities.get_values_from_pvalue_file(pvalue_file_drug1)
                node_to_vals_drug2 = guild_utilities.get_values_from_pvalue_file(pvalue_file_drug2)
            elif type_of_profile == 'edge':
                top_node_to_vals_drug1 = comp.parse_edge_profile(file1)
                top_node_to_vals_drug2 = comp.parse_edge_profile(file2)
                node_to_vals_drug1 = comp.parse_edge_profile(network_scored_drug1)
                node_to_vals_drug2 = comp.parse_edge_profile(network_scored_drug2)
            elif type_of_profile == 'functional':
                top_node_to_vals_drug1 = comp.parse_enrichment_file(file1)
                top_node_to_vals_drug2 = comp.parse_enrichment_file(file2)

            common_items = len(set(top_node_to_vals_drug1.keys()).intersection(set(top_node_to_vals_drug2.keys())))

            # Calculating the differences in keys between network 1 and network 2
            diff1 = set(top_node_to_vals_drug1.keys()) - set(top_node_to_vals_drug2.keys())
            # Calculating the differences in keys between network 2 and network 1
            diff2 = set(top_node_to_vals_drug2.keys()) - set(top_node_to_vals_drug1.keys())

            # For diff2 set, add as values a set (0, 1) --> score = 0, pval = 1
            diff2_dict = dict.fromkeys(diff2, (0, 1))
            # If we sum the first network with the difference between 2 and 1, we obtain all the nodes from both networks
            top_node_to_vals_drug1.update(diff2_dict)

            # For diff1 set, add as values a set (0, 1) --> score = 0, pval = 1
            diff1_dict = dict.fromkeys(diff1, (0, 1))
            # If we sum the second network with the difference between 1 and 2, we obtain all the nodes from both networks
            top_node_to_vals_drug2.update(diff1_dict)

            ratio_common_total = float(common_items) / float(len(top_node_to_vals_drug1))

            if (len(top_node_to_vals_drug1) == len(top_node_to_vals_drug2)):
                pass
            else:
                print("  DIANA INFO:\tThe vectors from both networks are not equal, so comparison cannot continue.\n\t\tSorry for the inconvenience.\n")
                sys.exit(10)

            if type_of_profile == 'functional':
                (scores_drug1, scores_drug2) = comp.create_vectors_from_dicts(top_node_to_vals_drug1, top_node_to_vals_drug2)
            else:
                # Get the top scoring items for both dictionaries
                filt_node_to_vals_drug1 = comp.filter_top_scoring_nodes(top_node_to_vals_drug1, node_to_vals_drug1)
                filt_node_to_vals_drug2 = comp.filter_top_scoring_nodes(top_node_to_vals_drug2, node_to_vals_drug2)

                (scores_drug1, scores_drug2) = comp.create_vectors_from_dicts(filt_node_to_vals_drug1, filt_node_to_vals_drug2)

            # Compute the dot product for the two vectors, multiplying the scores at each position
            dot_product = comp.calculate_norm_dot_product(scores_drug1, scores_drug2)

            (spearman, spearman_tie) = comp.calculate_spearman(scores_drug1, scores_drug2)

            data_frame.append([top, spearman[0], spearman[1], dot_product, common_items, len(top_node_to_vals_drug1), ratio_common_total])

        extended[type_of_profile] = data_frame
        # Output the data frame of the results
        f.write("# Top %\tSpearman\tP-value\tDot_Product\tcommon_items\ttotal_items\tratio_common_total\n")
        for row in data_frame:
            (top, spearman, pvalue, dot_product, common, total, ratio) = row
            f.write("%s\t%0.3f\t%0.3f\t%0.3f\t%d\t%d\t%0.3f\n" %(top, spearman, pvalue, dot_product, common, total, ratio))
        f.write("# END DATAFRAME\n")

    f.close()
    return extended


def run_extended_analysis_node(scoring_folder, pvalue_file, top_thresholds_list):
    """
    Creates a node profile for each of the top thesholds in the list.
    Returns a list with the names of the node profiles created
    """

    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Set files list to store the names of the files created
    profiles_list = []

    for top in top_thresholds_list:

        node_profile_file = scoring_folder + "/node_profile_"+str(top)+".txt"

        # Add the name of the file created inside the list
        profiles_list.append(node_profile_file)

        # Get top scoring, i.e. nodes that are in a given percentage of top score
        top_nodes = set()
        top_nodes_filtered = set()
        ntop=float(top)*len([x for x in node_to_vals.iteritems()])/100.0
        ii=0

        f = open(node_profile_file, 'w')
        delim = " "
        f.write("#Id%sScore%sP-value\n" % (delim, delim))

        # Now, write the profile with the top scoring nodes
        for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
            score, pval = vals
            if ii < ntop:
                top_nodes.add(node)
                f.write("%s%s%f%s%s\n" % (node, delim, score, delim, pval))
                ii=ii+1
            else:
                break

        f.close()

        print("  DIANA INFO:\t%s file created.\n" %(node_profile_file))


    return profiles_list


def run_extended_analysis_edge(scoring_folder, network_file, pvalue_file, top_thresholds_list):
    """
    Create a series of scored networks on basis of the thresholds defined in the extended analysis
    """

    # Get GUILD scores and p-values
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Create a network from the network file
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)

    # Set files list to store the names of the files created
    profiles_list = []

    for top_threshold in top_thresholds_list:

        edge_profile_file = scoring_folder + "/edge_profile_top_"+str(top_threshold)+".txt"

        # Add the name of the file created inside the list
        profiles_list.append(edge_profile_file)

        # Get top scoring, i.e. nodes that are in a given percentage of top score
        top_nodes = set()
        ntop=float(top_threshold)*len([x for x in node_to_vals.iteritems()])/100.0
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
        g_sub = network_utilities.get_subgraph(g, top_nodes)

        f = open(edge_profile_file, 'w')

        # The subnetwork induced by the top scoring nodes is scored and recorded in a file
        for u, v in g_sub.edges():
            score_u = node_to_vals[u][0]
            score_v = node_to_vals[v][0]
            score = (score_u + score_v) / 2
            f.write("%s\t%f\t%s\n" % (u, score, v))

        f.close()

        print("  DIANA INFO:\t%s file created.\n" %(edge_profile_file))

    return profiles_list


def run_extended_analysis_functional(scoring_folder, top_thresholds_list):
    """
    Using the functional profile previously created using a threshold, the edge profile is cut in
    different subprofiles using the top thesholds in the list.
    Returns a list with the names of the profiles created
    """

    # Set files list to store the names of the files created
    profiles_list = []

    for top in top_thresholds_list:

        node_profile_file = scoring_folder + "/node_profile_"+str(top)+".txt"
        functional_profile_file = scoring_folder + "/functional_profile_top_"+str(top)+".txt"

        # Add the name of the file created inside the list
        profiles_list.append(functional_profile_file)

        # Get values from node profile
        node_to_vals = guild_utilities.get_values_from_pvalue_file(node_profile_file)
        # Get only the nodes
        all_nodes = node_to_vals.keys()
        # Calculate the enrichment analysis using the nodes
        f = open(functional_profile_file, 'w')
        functional_enrichment.check_functional_enrichment(all_nodes, all_nodes, "geneid",
                                                              f.write, species = "Homo sapiens")
        f.close()

        print("  DIANA INFO:\t%s file created.\n" %(functional_profile_file))

    return profiles_list


def run_extended_analysis_functional_GOA(scoring_folder, top_thresholds_list, obodag, geneid2gos):
    """
    Using the functional profile previously created using a threshold, the edge profile is cut in
    different subprofiles using the top thesholds in the list.
    Returns a list with the names of the profiles created
    """

    # Set files list to store the names of the files created
    profiles_list = []

    for top in top_thresholds_list:

        node_profile_file = scoring_folder + "/node_profile_"+str(top)+".txt"
        functional_profile_file = scoring_folder + "/functional_profile_top_"+str(top)+".txt"
        pvalue_file = scoring_folder+"/output_scores.sif.netcombo.pval"

        # Add the name of the file created inside the list
        profiles_list.append(functional_profile_file)

        # Get background genes
        node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)
        all_nodes = node_to_vals.keys()
        all_nodes = [ int(x) for x in all_nodes ]

        # Get top genes
        top_node_to_vals = guild_utilities.get_values_from_pvalue_file(node_profile_file)
        top_nodes = top_node_to_vals.keys()
        top_nodes = [ int(x) for x in top_nodes ]

        temp_file = scoring_folder+"/temp_enrichment_goatools.txt"

        GOA.calculate_functional_enrichment_profile(obodag, geneid2gos, top_nodes, all_nodes, temp_file, functional_profile_file)

        print("  DIANA INFO:\t%s file created.\n" %(functional_profile_file))

    return profiles_list
