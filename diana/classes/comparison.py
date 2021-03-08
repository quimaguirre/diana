import sys, os, re
import Bio.Cluster
import numpy as np
import scipy.stats
import diana.classes.indigo_toolkit.indigo as indigo


class ComparisonResult(object):
    """ 
    Class defining a the results of the comparison between two profiles
    """

    def __init__(self):
        """ 
        """

        self.target_results = {'target':[], 'pfam':[], 'gobp-fdr_bh':[], 'gobp-bonferroni':[], 'gomf-fdr_bh':[], 'gomf-bonferroni':[], 'reactome-fdr_bh':[], 'reactome-bonferroni':[]}
        self.guild_results = {'node':{}, 'edge':{}, 'gobp-fdr_bh':{}, 'gobp-bonferroni':{}, 'gomf-fdr_bh':{}, 'gomf-bonferroni':{}, 'reactome-fdr_bh':{}, 'reactome-bonferroni':{}}
        self.structure_result = None
        self.atc_result = []
        self.se_result = []
        self.example_threshold = None
        self.node_examples = []
        self.edge_examples = []
        self.function_examples = []

    ###########
    # METHODS #
    ###########

    def add_target_result(self, data_type, result):
        """
        Adds the result of a comparison from Target of a certain type of data
        (target, pfam, function)
        """
        if data_type in self.target_results:
            self.target_results[data_type] = result
        else:
            print('Incorrect type of data for Target! It must be target, pfam or function.\n')
            sys.exit(10)
        return

    def add_guild_result(self, data_type, threshold, result):
        """
        Adds the result of a comparison from guild of a certain type of data
        (node, edge, function) and a certain top threshold
        """
        if data_type in self.guild_results:
            self.guild_results[data_type][threshold] = result
        else:
            print('Incorrect type of data for guild! It must be node, edge or function.\n')
            sys.exit(10)
        return

    def get_median_threshold(self, threshold_list):
        """
        Obtains the threshold in the middle of the threshold list.
        It is used as example threshold, which is the top threshold used
        in the profiles used as examples
        """
        if len(threshold_list) % 2 == 0:
            self.example_threshold = threshold_list[int(len(threshold_list)/2)-1]
        else:
            self.example_threshold = threshold_list[int(len(threshold_list)/2)]
        return

    def output_results_table(self, results_table, threshold_list):
        """
        Outputs a results table containing the results of all the comparisons
        """
        with open(results_table, 'w') as results_table_fd:
            results_table_fd.write('method\tdata_type\tthreshold\tdot_product\tspearman\tjaccard\n')

            results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'target', 'target', '-', self.target_results['target'][0], self.target_results['target'][1], self.target_results['target'][2] ))
            results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'target', 'pfam', '-', self.target_results['pfam'][0],  self.target_results['pfam'][1],  self.target_results['pfam'][2] ))
            results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'target', 'function', '-',  self.target_results['function'][0], self.target_results['function'][1], self.target_results['function'][2] ))

            for top_threshold in threshold_list:
                results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'guild', 'node', top_threshold, self.guild_results['node'][top_threshold][0], self.guild_results['node'][top_threshold][1], self.guild_results['node'][top_threshold][2] ))
                results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'guild', 'edge', top_threshold, self.guild_results['edge'][top_threshold][0], self.guild_results['edge'][top_threshold][1], self.guild_results['edge'][top_threshold][2] ))
                results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'guild', 'function', top_threshold, self.guild_results['function'][top_threshold][0], self.guild_results['function'][top_threshold][1], self.guild_results['function'][top_threshold][2] ))

            results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'structure', 'structure', '-', self.structure_result, self.structure_result, self.structure_result ))
            results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'atc', 'atc', '-', self.atc_result[0], self.atc_result[1], self.atc_result[2] ))
            results_table_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( 'se', 'se', '-', self.se_result[0], self.se_result[1], self.se_result[2] ))

    def output_venn_diagram_guild(self, drug1_name, drug2_name, output_venn_nodes, output_venn_edges, output_venn_functions):
        """
        Outputs three Venn diagrams with the overlap between the nodes, edges and functions in one of the thresholds used
        """
        nodes_drug1 = set(self.node_examples[0].node_to_values.keys())
        nodes_drug2 = set(self.node_examples[1].node_to_values.keys())
        shared_nodes = nodes_drug1 & nodes_drug2
        print(shared_nodes)
        plot_venn_2(drug1_name, drug2_name, nodes_drug1, nodes_drug2, shared_nodes, output_venn_nodes)

        edges_drug1 = set(self.edge_examples[0].edge_to_values.keys())
        edges_drug2 = set(self.edge_examples[1].edge_to_values.keys())
        shared_edges = edges_drug1 & edges_drug2
        plot_venn_2(drug1_name, drug2_name, edges_drug1, edges_drug2, shared_edges, output_venn_edges)

        functions_drug1 = set(self.function_examples[0].go_id_to_values.keys())
        functions_drug2 = set(self.function_examples[1].go_id_to_values.keys())
        shared_functions = functions_drug1 & functions_drug2
        plot_venn_2(drug1_name, drug2_name, functions_drug1, functions_drug2, shared_functions, output_venn_functions)

        return



#############
#############
# FUNCTIONS #
#############
#############


def calculate_comparison(node_to_vals_drug1, node_to_vals_drug2, n_total):
    """
    Calculates the comparison between two drugs from a dict containing as keys the nodes/edges/GOs and as values
    a set (score, pval). Returns a list with the results.
    n_total = total number of nodes/edges in the network (or total number of GOs in the case of functions)
    """

    ### Compute the Jaccard Index ###
    jaccard = calculate_jaccard_index(set(node_to_vals_drug1.keys()), set(node_to_vals_drug2.keys()))

    print('  DIANA INFO:\tThe Jaccard Index is {:.3f}\n'.format(jaccard))


    ### Compute the Fisher's test ###
    n_common, fisher_oddsratio, fisher_pvalue = calculate_fisher_test(set(node_to_vals_drug1.keys()), set(node_to_vals_drug2.keys()), n_total)

    print('  DIANA INFO:\tFisher\'s test results. Number of common: {}. Odds ratio: {:.3f}. P-value: {:.3f}.\n'.format(n_common, fisher_oddsratio, fisher_pvalue))


    ### Compute the Dot Product ###

    #### First, the keys of both vectors need to be equal ####

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

    # Check that the two dictionaries have now the same nodes!
    if sorted(node_to_vals_drug1.keys()) != sorted(node_to_vals_drug2.keys()):
        print('  DIANA INFO:\tThe vectors from both networks are not equal, so comparison cannot continue.\n\t\tSorry for the inconvenience.\n')
        sys.exit(10)

    # Create two lists of vectors containing the scores of their corresponding keys in each position
    (scores_drug1, scores_drug2) = create_vectors_from_dicts(node_to_vals_drug1, node_to_vals_drug2)

    # Compute the dot product for the two vectors, multiplying the scores at each position
    dot_product = calculate_norm_dot_product(scores_drug1, scores_drug2)

    print('  DIANA INFO:\tThe dot product between the scores of the profiles is: {:.3f}\n'.format(dot_product))


    ### Compute the Spearman correlation coefficient ###

    #(spearman_result, spearman_pval) = calculate_spearman_without_tie(scores_drug1, scores_drug2)
    #print('  DIANA INFO:\tThe Spearman rank correlation coefficient is {:.3f} with P-value {:.3f}'.format(spearman_result, spearman_pval))

    spearman_tie = calculate_spearman_with_tie(scores_drug1, scores_drug2)
    print('  DIANA INFO:\tThe Spearman rank correlation coefficient with tie correction is {:.3f}\n'.format(spearman_tie))

    return dot_product, spearman_tie, jaccard, n_common, fisher_oddsratio, fisher_pvalue


def calculate_comparison_top_scoring(top_node_to_vals_drug1, node_to_vals_drug1, top_node_to_vals_drug2, node_to_vals_drug2):
    """
    Calculates the comparison between two drugs from a dict containing as keys the nodes/edges and as values
    a set (score, pval). 
    If the value of a top scoring profile is not in the profile of the other drug, the score for the latter is obtained from the complete profile
    Returns a list with the results.
    """

    # Compute the Jaccard Index
    jaccard = calculate_jaccard_index(set(top_node_to_vals_drug1.keys()), set(top_node_to_vals_drug2.keys()))

    print('  DIANA INFO:\tThe Jaccard Index is {:.3f}\n'.format(jaccard))


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


    # Check that the two dictionaries have now the same nodes!
    if set(top_node_to_vals_drug1.keys()) != set(top_node_to_vals_drug2.keys()):
        print('  DIANA INFO:\tThe vectors from both networks are not equal, so comparison cannot continue.\n\t\tSorry for the inconvenience.\n')
        sys.exit(10)

    # Get the top scoring items for both dictionaries
    filt_node_to_vals_drug1 = filter_top_scoring_values(top_node_to_vals_drug1, node_to_vals_drug1)
    filt_node_to_vals_drug2 = filter_top_scoring_values(top_node_to_vals_drug2, node_to_vals_drug2)

    # Check that the two dictionaries have now the same items!
    if set(filt_node_to_vals_drug1.keys()) != set(filt_node_to_vals_drug2.keys()):
        print('  DIANA INFO:\tThe vectors from the filtered networks are not equal, so comparison cannot continue.\n\t\tSorry for the inconvenience.\n')
        sys.exit(10)

    # Create two lists of vectors containing the scores of their corresponding keys in each position
    (scores_drug1, scores_drug2) = create_vectors_from_dicts(filt_node_to_vals_drug1, filt_node_to_vals_drug2)

    # Compute the dot product for the two vectors, multiplying the scores at each position
    dot_product = calculate_norm_dot_product(scores_drug1, scores_drug2)

    print('  DIANA INFO:\tThe dot product between the scores of the profiles is: {:.3f}\n'.format(dot_product))


    # Compute the Spearman correlation coefficient
    #(spearman_result, spearman_pval) = calculate_spearman_without_tie(scores_drug1, scores_drug2)
    #print('  DIANA INFO:\tThe Spearman rank correlation coefficient is {:.3f} with P-value {:.3f}'.format(spearman_result, spearman_pval))
    spearman_tie = calculate_spearman_with_tie(scores_drug1, scores_drug2)
    print('  DIANA INFO:\tThe Spearman rank correlation coefficient with tie correction is {:.3f}\n'.format(spearman_tie))

    return dot_product, spearman_tie, jaccard


def create_vectors_from_dicts(dict1, dict2):
    """
       From two dicts that contain as values a tuple = (score, pval), two lists are created in which
       each position contains the score of a concrete key which has to be the same for both lists
       Example: dict1 = {"C":(1,0), "B":(5,0), "A":(9,0)}
                dict2 = {"C":(5,0), "A":(7,0), "B":(0.01,0)}
                The lists returned will be:
                list1 = [9, 1, 5]
                list2 = [7, 5, 0.01]
    """

    list1 = []
    list2 = []
    for key in dict1:
        list1.append(dict1[key][0])
        list2.append(dict2.get(key)[0])
    return(list1, list2)


def calculate_norm_dot_product(list1, list2):
    """
    Calculates the normalized dot product of two lists
    """
    # Obtains the length of the list, which is the square root of the list multiplied by itself
    # Or what is the same, the square root of all the components of the list squared
    length_list1 = np.sqrt(np.dot(list1, list1))
    length_list2 = np.sqrt(np.dot(list2, list2))
    # Divides all the components of the list by its length to obtain the unitary vector
    list1 = [x / length_list1 for x in list1]
    list2 = [x / length_list2 for x in list2]
    # Performs the dot product of the unitary vectors
    return np.dot(list1, list2)


def calculate_spearman_without_tie(list1, list2):
    """
    Calculates the Spearman rank correlation coefficient (without tie correction)
    More info: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rank_correlations/#Spearman
    """
    (spearman_result, spearman_pval) = scipy.stats.stats.spearmanr(list1, list2)
    return spearman_result, spearman_pval


def calculate_spearman_with_tie(list1, list2):
    """
    Calculates the Spearman rank correlation coefficient (with tie correction)
    More info: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rank_correlations/#Spearman
    """
    spearman_tie = 1 - Bio.Cluster.distancematrix((list1,list2), dist="s")[1][0]
    return spearman_tie


def calculate_jaccard_index(set1, set2):
    """
    Calculates the Jaccard index of two sets
    """
    size_intersection = float(len(set1.intersection(set2)))
    size_union = float(len(set1.union(set2)))
    return size_intersection / size_union


def calculate_simpson_index(set1, set2):
    """
    Calculates the Simpson index of two sets
    """
    size_intersection = float(len(set1.intersection(set2)))
    size_smaller_set = min(float(len(set1)), float(len(set2)))
    if size_smaller_set <= 0:
        return 0
    else:
        return size_intersection / size_smaller_set


def calculate_fisher_test(set1, set2, n_total):
    """
    Calculates the Fisher's test between two sets
    set1 = proteins of drug 1
    set2 = proteins of drug 2
    n_total = total number of nodes/edges in the network (or total number of GOs in case of functions)
    """
    n1 = len(set1)
    n2 = len(set2)
    # Get the common proteins between the two sessions
    common = set1 & set2
    n_common = len(common)

    # Create the contingency table and calculate the fisher test

    #             |---------------|------------------------------|
    #             | Top 2         | Non-Top 2                    |
    # |-----------|---------------|------------------------------|
    # | Top 1     | n_common      | n1 - n_common                |
    # |-----------|---------------|------------------------------|
    # | Non-Top 1 | n2 - n_common | n_total - n1 - n2 + n_common |
    # |-----------|---------------|------------------------------|

    contigency = [[n_common, n1 - n_common], [n2 - n_common, n_total - n1 - n2 + n_common]]
    oddsratio, pvalue = scipy.stats.fisher_exact(contigency, alternative="greater")
    return n_common, oddsratio, pvalue


def generate_targets_dict_for_comparison(targets_list):
    """
       Generates a dict which contains the target name as key, and the score '1.00' and a p-value 'foo'
       as a set in the value of the dict. This dict will be used for the targets comparison
    """

    targets_dict = {}

    for target in targets_list:
        targets_dict[target] = (1.0, "foo")

    return targets_dict


def filter_top_scoring_values(top_node_to_vals, node_to_vals):
    """
    From the dict that contains all the nodes in common in both drugs scored, we get the top ones and their values.
    """
    filtered_node_to_vals = {}

    for item in top_node_to_vals:
        if item in node_to_vals:
            filtered_node_to_vals[item] = node_to_vals[item]
        else:
            filtered_node_to_vals[item] = top_node_to_vals[item]

    return filtered_node_to_vals


def get_smiles_similarity(smiles1, smiles2, fp_type = "sub", metric = "tanimoto"):
    """
    Calculates the structural similarity between the SMILES of two compounds by means of the Indigo package.
    fp_type: sim | sub
      - "sim" = similarity fingerprints: shorter
      - "sub" = substructure fingerprints: more descriptive
    metric: tanimoto | tversky
    """
    if len(smiles1) == 0 or len(smiles2) == 0:
        return None
    ind = indigo.Indigo()

    # Sometimes the SMILES are not correct, and they fail when loading. This is why I have introduced this exception
    try:
        m = ind.loadMolecule(smiles1)
        m2 = ind.loadMolecule(smiles2)
    except indigo.IndigoException: 
        return None
    m.aromatize()
    fp = m.fingerprint(fp_type)
    m2.aromatize() # Aromatize molecules in case they are not in aromatic form
    fp2 = m2.fingerprint(fp_type) # Calculate similarity between "similarity" fingerprints
    d = ind.similarity(fp, fp2, metric)
    return d


def plot_venn_2(drug1_name, drug2_name, items_drug1, items_drug2, shared_items, output_plot):
    """
    Plot a venn diagram of 2 variables
    """
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2, venn2_circles

    # Create a figure of size 8x6 inches, 80 dots per inch
    fig = plt.figure(figsize=(8, 8), dpi=100)

    # Create a new subplot from a grid of 2x1 --> bar plot for the lost PDIs
    ax = fig.add_subplot(111)
    # s will contain the subset sizes #
    s = (len(items_drug1)- len(shared_items), len(items_drug2)- len(shared_items), len(shared_items))

    v = venn2(subsets=s, set_labels=(drug1_name, drug2_name))

    # Subset labels
    v.get_label_by_id('10').set_text(len(items_drug1) - len(shared_items))
    v.get_label_by_id('01').set_text(len(items_drug2) - len(shared_items))
    if len(shared_items) > 0:
        v.get_label_by_id('11').set_text(len(shared_items))

    # Subset fontsize
    for text in v.set_labels:
        if text:
            text.set_fontsize(22)
    for text in v.subset_labels:
        if text:
            text.set_fontsize(22)

    # Subset colors
    v.get_patch_by_id('10').set_color('c')
    v.get_patch_by_id('10').set_alpha(1.0)
    v.get_patch_by_id('01').set_color('#993333')
    v.get_patch_by_id('01').set_alpha(0.6)

    # Border styles
    c = venn2_circles(subsets=s, linestyle='solid')

    #plt.title("Consistence between different methods", fontsize=22)

    fig.savefig(output_plot, bbox_inches='tight')

    return


