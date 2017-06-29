from toolbox import guild_utilities
import toolbox.comparison as comp


def extended_analysis_comparison(profiles_list_node, profiles_list_edge, profiles_list_functional, top_thresholds_list, output_file):
    """
    Performs the comparison between all the different profiles obtained from the extended analysis
    """
    f = open(output_file, 'w')
    f.write("##### EXTENDED ANALYSIS #####\n")

    # We will perform two extended analyses, one for the node profile and the other for the edge profile
    for type_of_profile in ('node', 'edge', 'functional'):

        # Skip if speed == 2 (so that edge profile has not being calculated)
        if profiles_list_edge == []:
            continue

        f.write("\n# %s PROFILE\n"%(type_of_profile.upper()))

        if type_of_profile == 'node':
            profiles_drug1 = profiles_list_node[0]
            profiles_drug2 = profiles_list_node[1]
            # Create a list of lists containing the profile_drug1_top_x, the profile_drug2_top_x and the top x
            profiles_comparison = list( zip(profiles_list_node[0], profiles_list_node[1], top_thresholds_list) )
        elif type_of_profile == 'edge':
            profiles_drug1 = profiles_list_edge[0]
            profiles_drug2 = profiles_list_edge[1]
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
                node_to_vals_drug1 = guild_utilities.get_values_from_pvalue_file(file1)
                node_to_vals_drug2 = guild_utilities.get_values_from_pvalue_file(file2)
            elif type_of_profile == 'edge':
                node_to_vals_drug1 = comp.parse_edge_profile(file1)
                node_to_vals_drug2 = comp.parse_edge_profile(file2)
            elif type_of_profile == 'functional':
                node_to_vals_drug1 = comp.parse_enrichment_top_file(file1)
                node_to_vals_drug2 = comp.parse_enrichment_top_file(file2)

            common_items = len(set(node_to_vals_drug1.keys()).intersection(set(node_to_vals_drug2.keys())))

            # f.write("Initially, the length of network 1 is: %d\n" %(len(node_to_vals_drug1)))
            # f.write("Initially, the length of network 2 is: %d\n" %(len(node_to_vals_drug2)))
            # f.write("The common items between 1 and 2 are: %d\n" %(common_items))

            # Calculating the differences in keys between network 1 and network 2
            diff1 = set(node_to_vals_drug1.keys()) - set(node_to_vals_drug2.keys())
            # Calculating the differences in keys between network 2 and network 1
            diff2 = set(node_to_vals_drug2.keys()) - set(node_to_vals_drug1.keys())
            # For diff2 set, add as values a set (0, 1) --> score = 0, pval = 1
            diff2_dict = dict.fromkeys(diff2, (0, 1))
            # If we sum the first network with the difference between 2 and 1, we obtain all the nodes from both networks
            node_to_vals_drug1.update(diff2_dict)
            # For diff1 set, add as values a set (0, 1) --> score = 0, pval = 1
            diff1_dict = dict.fromkeys(diff1, (0, 1))
            # If we sum the second network with the difference between 1 and 2, we obtain all the nodes from both networks
            node_to_vals_drug2.update(diff1_dict)

            # f.write("Updating network 1 with diff between 2 and 1, the new length of network 1 is: %d\n" %(len(node_to_vals_drug1)))

            ratio_common_total = float(common_items) / float(len(node_to_vals_drug1))

            if (len(node_to_vals_drug1) == len(node_to_vals_drug2)):
                pass
            else:
                print("  DIANA INFO:\tThe vectors from both networks are not equal, so comparison cannot continue.\n\t\tSorry for the inconvenience.\n")
                sys.exit(10)

            (scores_drug1, scores_drug2) = comp.create_vectors_from_dicts(node_to_vals_drug1, node_to_vals_drug2)

            # Compute the dot product for the two vectors, multiplying the scores at each position
            dot_product = comp.calculate_norm_dot_product(scores_drug1, scores_drug2)

            (spearman, spearman_tie) = comp.calculate_spearman(scores_drug1, scores_drug2)

            # Output the results in the extended analysis file
            # f.write("# Spearman: %0.3f (P-value: %0.3f)\n" %(spearman[0],spearman[1]))
            # f.write("# Dot product: %0.3f\n" %(dot_product))

            data_frame.append([top, spearman[0], spearman[1], dot_product, common_items, len(node_to_vals_drug1), ratio_common_total])

        # Output the data frame of the results
        f.write("# Top %\tSpearman\tP-value\tDot_Product\tcommon_items\ttotal_items\tratio_common_total\n")
        for row in data_frame:
            (top, spearman, pvalue, dot_product, common, total, ratio) = row
            f.write("%s\t%0.3f\t%0.3f\t%0.3f\t%d\t%d\t%0.3f\n" %(top, spearman, pvalue, dot_product, common, total, ratio))
        f.write("# END DATAFRAME\n")

    f.close()


def run_extended_analysis_node(scoring_folder, pvalue_file, top_thresholds_list, nodes_or_linkers):
    """
    Creates a node profile for each of the top thesholds in the list.
    Returns a list with the names of the node profiles created
    """
    # Set p-value threshold
    pval_threshold = 1

    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Set files list to store the names of the files created
    profiles_list = []

    for top in top_thresholds_list:

        # Set uniprotentry profile file
        if nodes_or_linkers == "nodes":
            uniprotentry_profile_file = scoring_folder + "/node_profile_"+str(pval_threshold)+"_"+str(top)+".txt"
        elif nodes_or_linkers == "linkers":
            uniprotentry_profile_file = scoring_folder + "/linker_node_profile_"+str(pval_threshold)+"_"+str(top)+".txt"

        # Add the name of the file created inside the list
        profiles_list.append(uniprotentry_profile_file)

        # Get top scoring, i.e. nodes that are in a given percentage of top score
        top_nodes = set()
        top_nodes_filtered = set()
        ntop=float(top)*len([x for x in node_to_vals.iteritems()])/100.0
        ii=0

        f = open(uniprotentry_profile_file, 'w')
        delim = " "
        f.write("#Id%sScore%sP-value\n" % (delim, delim))

        # Now, filter top_nodes by p-value <= (0.05). Include them in the profile if entry_type == "uniprotentry"
        for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
            score, pval = vals
            if ii < ntop:
                top_nodes.add(node)
                if pval <= pval_threshold:
                    top_nodes_filtered.add(node)
                    f.write("%s%s%f%s%s\n" % (node, delim, score, delim, pval))
                ii=ii+1
            else:
                break

        f.close()

        print("  DIANA INFO:\t%s file created.\n" %(uniprotentry_profile_file))


    return profiles_list


def run_extended_analysis_edge(scoring_folder, scored_profile, top_thresholds_list, edges_or_linkers):
    """
    Using the edge profile previously created using a threshold, the edge profile is cut in different 
    subprofiles using the top thesholds in the list.
    Returns a list with the names of the profiles created.
    It can be used with the edge profile or the linkers subnetwork.
    """

    # Get values from the scored edge profile
    edge_profile_list = guild_utilities.get_values_from_edges_file_scored(scored_profile)

    # Set files list to store the names of the files created
    profiles_list = []

    for top in top_thresholds_list:

        # Set uniprotentry profile file
        if edges_or_linkers == "edges":
            edge_profile_file = scoring_folder + "/edge_profile_top_"+str(top)+".txt"
        elif edges_or_linkers == "linkers":
            edge_profile_file = scoring_folder + "/linker_edge_profile_top_"+str(top)+".txt"

        # Add the name of the file created inside the list
        profiles_list.append(edge_profile_file)

        # Calculate the number of top edges
        ntop=float(top)*len(edge_profile_list)/100.0
        ii=0

        f = open(edge_profile_file, 'w')

        for edge1, score, edge2 in sorted(edge_profile_list, key=lambda x: x[1], reverse=True):
            if ii < ntop:
                f.write("%s\t%f\t%s\n" % (edge1, float(score), edge2))
                ii=ii+1
            else:
                break

        f.close()

        print("  DIANA INFO:\t%s file created.\n" %(edge_profile_file))

    return profiles_list


def run_extended_analysis_functional(scoring_folder, functional_profile, top_thresholds_list, functional_or_linkers):
    """
    Using the functional profile previously created using a threshold, the edge profile is cut in
    different subprofiles using the top thesholds in the list.
    Returns a list with the names of the profiles created
    """

    # Get values from the scored edge profile
    functional_profile_parsed = comp.parse_enrichment_file(functional_profile)

    # Set files list to store the names of the files created
    profiles_list = []

    for top in top_thresholds_list:

        # Set uniprotentry profile file
        if functional_or_linkers == "functional":
            functional_profile_file = scoring_folder + "/functional_profile_top_"+str(top)+".txt"
        elif functional_or_linkers == "linkers":
            functional_profile_file = scoring_folder + "/functional_linker_profile_top_"+str(top)+".txt"

        # Add the name of the file created inside the list
        profiles_list.append(functional_profile_file)

        # Calculate the number of top edges
        ntop=float(top)*len([x for x in functional_profile_parsed.iteritems()])/100.0
        ii=0

        f = open(functional_profile_file, 'w')

        for go, vals in sorted(functional_profile_parsed.items(),key=lambda x: x[1][0],reverse=True):
            odds, pval = vals
            if ii < ntop:
                f.write("%s\t%f\t%s\n" % (go, float(odds), pval))
                ii=ii+1
            else:
                break

        f.close()

        print("  DIANA INFO:\t%s file created.\n" %(functional_profile_file))

    return profiles_list
