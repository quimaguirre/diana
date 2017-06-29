import numpy as np
import Bio.Cluster
import scipy.stats


def calculate_comparison(node_to_vals_drug1, node_to_vals_drug2):
    """
    Calculates the comparison between two drugs from a dict containing as keys the nodes/edges/GOs and as values
    a set (score, pval). Returns a list with the results.
    """

    summary = []

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

    if (len(node_to_vals_drug1) == len(node_to_vals_drug2)):
        pass
    else:
        print("  DIANA INFO:\tThe vectors from both networks are not equal, so comparison cannot continue.\n\t\tSorry for the inconvenience.\n")
        sys.exit(10)

    # Create two lists of vectors containing the scores of their corresponding keys in each position
    (scores_drug1, scores_drug2) = create_vectors_from_dicts(node_to_vals_drug1, node_to_vals_drug2)

    # Compute the dot product for the two vectors, multiplying the scores at each position
    dot_product = calculate_norm_dot_product(scores_drug1, scores_drug2)
    summary.append(dot_product)

    print("  DIANA INFO:\tThe dot product between the scores of the profiles is: %0.3f\n" %(dot_product))

    # Compute the Spearman correlation coefficient
    (spearman, spearman_tie) = calculate_spearman(scores_drug1, scores_drug2)
    summary.append(spearman[0])
    summary.append(spearman[1])

    print("  DIANA INFO:\tThe Spearman rank correlation coefficient is %0.3f with P-value %0.3f" %(spearman[0], spearman[1]))
    print("\t\tThe Spearman rank correlation coefficient with tie correction is %0.3f\n" %(spearman_tie))

    return summary

def calculate_comparison_top_scoring(top_node_to_vals_drug1, node_to_vals_drug1, top_node_to_vals_drug2, node_to_vals_drug2, type_of_profile):
    """
    Calculates the comparison between two drugs from a dict containing as keys the nodes/edges/GOs and as values
    a set (score, pval). 
    If the value of a top scoring profile is not in the profile of the other drug, the score for the latter is obtained from the complete profile
    Returns a list with the results.
    """

    summary = []

    # Calculating the differences in keys between network 1 and network 2
    diff1 = set(top_node_to_vals_drug1.keys()) - set(top_node_to_vals_drug2.keys())
    # Calculating the differences in keys between network 2 and network 1
    diff2 = set(top_node_to_vals_drug2.keys()) - set(top_node_to_vals_drug1.keys())

    # For diff2 set, add as values a set (0, 1) --> score = 0, pval = 1
    diff2_dict = dict.fromkeys(diff2, (0, 1))
    #print("  DIANA INFO:\tInitially, the length of network 1 is: %d" %(len(top_node_to_vals_drug1)))
    #print("\t\tAnd the difference between 2 and 1 is: %d" %(len(diff2_dict)))
    # If we sum the first network with the difference between 2 and 1, we obtain all the nodes from both networks
    top_node_to_vals_drug1.update(diff2_dict)
    #print("\t\tUpdating network 1 with diff between 2 and 1, the new length of network 1 is: %d\n" %(len(top_node_to_vals_drug1)))

    # For diff1 set, add as values a set (0, 1) --> score = 0, pval = 1
    diff1_dict = dict.fromkeys(diff1, (0, 1))
    #print("  DIANA INFO:\tInitially, the length of network 2 is: %d" %(len(top_node_to_vals_drug2)))
    #print("\t\tAnd the difference between 1 and 2 is: %d" %(len(diff1_dict)))
    # If we sum the second network with the difference between 1 and 2, we obtain all the nodes from both networks
    top_node_to_vals_drug2.update(diff1_dict)
    #print("\t\tUpdating network 2 with diff between 1 and 2, the new length of network 2 is: %d\n" %(len(top_node_to_vals_drug2)))

    if (len(top_node_to_vals_drug1) == len(top_node_to_vals_drug2)):
        #print("  DIANA INFO:\tAs the vectors for both networks are equal, the comparison will continue.\n")
        pass
    else:
        print("  DIANA INFO:\tThe vectors from both networks are not equal, so comparison cannot continue.\n\t\tSorry for the inconvenience.\n")
        sys.exit(10)

    # Get the top scoring items for both dictionaries
    if type_of_profile == 'node':
        filt_node_to_vals_drug1 = filter_top_scoring_nodes(top_node_to_vals_drug1, node_to_vals_drug1)
        filt_node_to_vals_drug2 = filter_top_scoring_nodes(top_node_to_vals_drug2, node_to_vals_drug2)
    elif type_of_profile == 'edge':
        filt_node_to_vals_drug1 = filter_top_scoring_edges(top_node_to_vals_drug1, node_to_vals_drug1)
        filt_node_to_vals_drug2 = filter_top_scoring_edges(top_node_to_vals_drug2, node_to_vals_drug2)
        filt1 = set(filt_node_to_vals_drug1.keys())
        filt2 = set(filt_node_to_vals_drug2.keys())
        if filt1 != filt2:
            print("  DIANA INFO:\tThere is an error in the comparison of edges. The edges of both vectors are not equal")
            print('filt1: {}'.format(filt1))
            print('filt2: {}'.format(filt2))
            sys.exit(10)

    # Create two lists of vectors containing the scores of their corresponding keys in each position
    (scores_drug1, scores_drug2) = create_vectors_from_dicts(filt_node_to_vals_drug1, filt_node_to_vals_drug2)

    # Compute the dot product for the two vectors, multiplying the scores at each position
    dot_product = calculate_norm_dot_product(scores_drug1, scores_drug2)
    summary.append(dot_product)

    print("  DIANA INFO:\tThe dot product between the scores of the profiles is: %0.3f\n" %(dot_product))

    # Compute the Spearman correlation coefficient
    (spearman, spearman_tie) = calculate_spearman(scores_drug1, scores_drug2)
    summary.append(spearman[0])
    summary.append(spearman[1])

    print("  DIANA INFO:\tThe Spearman rank correlation coefficient is %0.3f with P-value %0.3f" %(spearman[0], spearman[1]))
    print("\t\tThe Spearman rank correlation coefficient with tie correction is %0.3f\n" %(spearman_tie))

    return summary

def calculate_dot_product(list1, list2):
    """Calculates the dot product of two lists"""
    return np.dot(list1, list2)

def calculate_jaccard_index(set1, set2):
    """Calculates the Jaccard index of two sets"""
    size_intersection = float(len(set1.intersection(set2)))
    size_union = float(len(set1.union(set2)))
    return size_intersection / size_union

def calculate_norm_dot_product(list1, list2):
    """Calculates the normalized dot product of two lists"""
    # Obtains the length of the list, which is the square root of the list multiplied by itself
    # Or what is the same, the square root of all the components of the list squared
    length_list1 = np.sqrt(np.dot(list1, list1))
    length_list2 = np.sqrt(np.dot(list2, list2))
    # Divides all the components of the list by its length to obtain the unitary vector
    list1 = [x / length_list1 for x in list1]
    list2 = [x / length_list2 for x in list2]
    # Performs the dot product of the unitary vectors
    return np.dot(list1, list2)

def calculate_spearman(list1, list2):
    """
    Calculates the Spearman rank correlation coefficient (without tie correction)
    More info: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rank_correlations/#Spearman
    """
    spearman = scipy.stats.stats.spearmanr(list1, list2)
    spearman_tie = 1 - Bio.Cluster.distancematrix((list1,list2), dist="s")[1][0]
    return (spearman, spearman_tie)

def compare_profiles(type_profile, profile1, profile2, pval=0.05):
    """Obtains a list with the common elements in both profiles"""

    elements1 = obtain_elements_from_profile(type_profile, profile1, pval)
    elements2 = obtain_elements_from_profile(type_profile, profile2, pval)
    intersection = list(elements1 & elements2)

    return intersection

def create_vectors_from_dicts(dict1, dict2):
    """
       From two dicts that contain as values a set = (score, pval), two lists are created in which
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

def filter_top_scoring_nodes(top_node_to_vals, node_to_vals):
    """
    From the dict that contains all the nodes in common in both drugs scored, we only take the top ones.
    """
    filt_node_to_vals = {}

    for node in top_node_to_vals:
        if node in node_to_vals:
            filt_node_to_vals[node] = node_to_vals[node]
        else:
            filt_node_to_vals[node] = top_node_to_vals[node]

    return filt_node_to_vals

def filter_top_scoring_edges(top_node_to_vals, node_to_vals):
    """
    From the dict that contains all the nodes in common in both drugs scored, we only take the top ones.
    """
    filt_node_to_vals = {}

    for nodes1 in top_node_to_vals:
        (node1, node2) = nodes1.split(' ')
        nodes2 = '{} {}'.format(node2, node1)
        if nodes1 in node_to_vals:
            filt_node_to_vals[nodes1] = node_to_vals[nodes1]
        elif nodes2 in node_to_vals:
            filt_node_to_vals[nodes1] = node_to_vals[nodes2]
        elif nodes1 not in node_to_vals and nodes2 not in node_to_vals:
            filt_node_to_vals[nodes1] = top_node_to_vals[nodes1]

    return filt_node_to_vals

def get_top_functions(enrichment_file, top):
    """
       Obtains from the top x functions from the enrichment profile.
       If the last one has the same score as the following, the following is included.
    """
    f = open(enrichment_file,"r")
    enrichment_dict = {}

    num_line = 0
    last_score = ''

    for line in f:
        if (line[0] == "#"):
            continue
        num_line += 1
        words = line.strip().split("\t")
        score = float(words[2])
        if num_line > int(top):
            if score == last_score:
                enrichment_dict[words[5]] = (score, words[4])
            else:
                break
        else:
            enrichment_dict[words[5]] = (score, words[4])
            last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last function, it is included as well
    f.close()
    return enrichment_dict

def obtain_elements_from_profile(type_profile, profile, pval = 0.05):
    """Obtains a set with the items of the profile"""

    f = open(profile,"r")
    elements = set()
    for line in f:
        if line[0] == "#":
            continue
        if type_profile == 'node':
            words = line.split()
            elements.add(words[0])
        elif type_profile == 'edge':
            words = line.strip().split("\t")
            nodes = "%s-%s" %(words[0], words[2])
            elements.add(nodes)
        elif type_profile == 'functional':
            words = line.strip().split("\t")
            if words[4][0] == '<':
                words[4] = 0
            else:
                words[4] = float(words[4])
            if words[4] < float(pval):
                elements.add(words[5])
    f.close()

    return elements

def parse_edge_profile(edge_profile):
    """
       Obtains from the edges profile a dict containing:
       - Key = "node1 node2"
       - Value = Set --> (Score, "useless")
    """
    f = open(edge_profile,"r")
    edges_dict = {}

    for line in f:
        words = line.strip().split("\t")
        nodes = "%s %s" %(words[0], words[2])
        edges_dict[nodes] = (float(words[1]), "useless")
    f.close()
    return edges_dict

def parse_edge_profiles_2_drugs(edge_profile_drug1, edge_profile_drug2):
    """
       Obtains from the edges profile a dict containing:
       - Key = "node1 node2"
       - Value = Set --> (Score, "useless")
       It cares about the position of the edges!!
       So, if drug 1 has the nodes in position (node1 score node2)
       and drug 2 has the nodes in position (node2 score node 1),
       the function puts them in the same order
    """
    edges_dict_drug1 = {}
    edges_dict_drug2 = {}
    f1 = open(edge_profile_drug1,"r")
    f2 = open(edge_profile_drug2,"r")

    for line in f1:
        words = line.strip().split("\t")
        nodes = "%s %s" %(words[0], words[2])
        edges_dict_drug1[nodes] = (float(words[1]), "useless")
    f1.close()

    for line in f2:
        words = line.strip().split("\t")
        nodes1 = "%s %s" %(words[0], words[2])
        nodes2 = "%s %s" %(words[2], words[0])
        if nodes2 in edges_dict_drug1:
            edges_dict_drug2[nodes2] = (float(words[1]), "useless")
        else:
            edges_dict_drug2[nodes1] = (float(words[1]), "useless")
    f2.close()

    return edges_dict_drug1, edges_dict_drug2

def parse_enrichment_file(enrichment_file):
    """
       Obtains from the functional enrichment profile a dict containing:
       - Key = GO term id
       - Value = Set --> (Log of odds ratio, Adjusted p-value)
    """
    f = open(enrichment_file,"r")
    enrichment_dict = {}

    for line in f:
        if (line[0] == "#"):
            continue
        words = line.strip().split("\t")
        enrichment_dict[words[5]] = (float(words[2]), words[4])
    f.close()
    return enrichment_dict

def parse_enrichment_top_file(enrichment_top_file):
    """
    Parses the enrichement file created by the Extended Analysis
    """
    f = open(enrichment_top_file,"r")
    enrichment_dict = {}

    for line in f:
        words = line.split()
        enrichment_dict[words[0]] = (float(words[1]), words[2])
    f.close()
    return enrichment_dict