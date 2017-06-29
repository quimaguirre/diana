import cPickle
import sys, os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import mysql.connector


def main():

    analyze_example()

    return


def analyze_example():

    ################
    #### INPUTS ####
    ################

    dcdb_drug1 = 'DCC0073'
    dcdb_drug2 = 'DCC0290'

    threshold = 10
    type_threshold = 'percentage' # percentage or item

    # Folder containing everything about the analysis
    current_dir = '/home/quim/project/diana_results'
    results_dir = current_dir + '/3_targets_analysis'
    drug_dir1 = results_dir + '/data/' + dcdb_drug1.upper() + '/guild_results_using_sif' 
    drug_dir2 = results_dir + '/data/' + dcdb_drug2.upper() + '/guild_results_using_sif'

    # Output plots
    comb_name = '{}-{}'.format(dcdb_drug1, dcdb_drug2)
    targets_venn = results_dir + '/venn_targets_{}.png'.format(comb_name)
    nodes_venn = results_dir + '/venn_nodes_{}.png'.format(comb_name)
    edges_venn = results_dir + '/venn_edges_{}.png'.format(comb_name)
    functions_venn = results_dir + '/venn_functions_{}.png'.format(comb_name)

    # Analyze targets
    dump_file = current_dir + "/dcdb2targets.pcl"
    dcdb2targets = cPickle.load(open(dump_file))

    targets_drug1 = dcdb2targets[dcdb_drug1]
    targets_drug2 = dcdb2targets[dcdb_drug2]
    shared_targets = targets_drug1 & targets_drug2
    print(len(targets_drug1 & targets_drug2))
    print(len(targets_drug1 | targets_drug2))

    plot_venn_2(dcdb_drug1, dcdb_drug2, targets_drug1, targets_drug2, shared_targets, targets_venn)


    # Analyze nodes
    node_prof_drug1 = os.path.join(drug_dir1, 'node_profile_{}_{}.txt'.format(str(threshold), type_threshold))
    node_prof_drug2 = os.path.join(drug_dir2, 'node_profile_{}_{}.txt'.format(str(threshold), type_threshold))
    nodes_drug1 = set(parse_node_profile(node_prof_drug1).keys())
    nodes_drug2 = set(parse_node_profile(node_prof_drug2).keys())
    shared_nodes = nodes_drug1 & nodes_drug2
    print(len(nodes_drug1 & nodes_drug2))
    print(len(nodes_drug1 | nodes_drug2))

    plot_venn_2(dcdb_drug1, dcdb_drug2, nodes_drug1, nodes_drug2, shared_nodes, nodes_venn)


    # Analyze edges
    edge_prof_drug1 = os.path.join(drug_dir1, 'edge_profile_top_{}_{}.txt'.format(str(threshold), type_threshold))
    edge_prof_drug2 = os.path.join(drug_dir2, 'edge_profile_top_{}_{}.txt'.format(str(threshold), type_threshold))
    edges_drug1 = set(parse_edge_profile(edge_prof_drug1).keys())
    edges_drug2 = set(parse_edge_profile(edge_prof_drug2).keys())
    shared_edges = edges_drug1 & edges_drug2
    print(len(edges_drug1 & edges_drug2))
    print(len(edges_drug1 | edges_drug2))

    plot_venn_2(dcdb_drug1, dcdb_drug2, edges_drug1, edges_drug2, shared_edges, edges_venn)


    # Analyze functions
    func_prof_drug1 = os.path.join(drug_dir1, 'functional_profile_top_{}_{}.txt'.format(str(threshold), type_threshold))
    func_prof_drug2 = os.path.join(drug_dir2, 'functional_profile_top_{}_{}.txt'.format(str(threshold), type_threshold))
    functions_drug1 = set(parse_functional_profile(func_prof_drug1).keys())
    functions_drug2 = set(parse_functional_profile(func_prof_drug2).keys())
    shared_functions = functions_drug1 & functions_drug2
    print(len(functions_drug1 & functions_drug2))
    print(len(functions_drug1 | functions_drug2))

    plot_venn_2(dcdb_drug1, dcdb_drug2, functions_drug1, functions_drug2, shared_functions, functions_venn)


    return


def obtain_dcdb2name(cnx):
    """
    Obtain dictionary DCDB_drugID : drug_name
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, N.value FROM externalEntityDCDB_drugID D, externalEntityName N WHERE D.externalEntityID = N.externalEntityID
             ''')

    cursor.execute(query)

    dcdb2name = {}
    for dcdb, name in cursor:
        dcdb2name[dcdb] = name

    cursor.close()

    return dcdb2name


def parse_node_profile(profile):
    """
    Parse the node profile
    """

    node_profile = {}

    profile_fd = open(profile, 'r')

    for line in profile_fd:

        if line[0] == '#':
            continue

        node, score, pval = line.strip().split()
        node_profile[node] = [float(score), float(pval)]

    profile_fd.close()

    return node_profile

def parse_edge_profile(profile):
    """
    Parse the edge profile
    """

    edge_profile = {}

    profile_fd = open(profile, 'r')

    for line in profile_fd:

        if line[0] == '#':
            continue

        node1, score, node2 = line.strip().split()
        edge = frozenset((node1,node2))
        edge_profile[edge] = float(score)

    profile_fd.close()

    return edge_profile

def parse_functional_profile(profile):
    """
    Parse the functional profile
    """

    functional_profile = {}

    profile_fd = open(profile, 'r')

    for line in profile_fd:

        if line[0] == '#':
            continue

        fields = line.strip().split()
        odds = fields[2]
        pval = fields[4]
        go = fields[5]
        functional_profile[go] = [float(odds), float(pval)]

    profile_fd.close()

    return functional_profile

def plot_venn_2(dcdb_drug1, dcdb_drug2, items_drug1, items_drug2, shared_items, output_plot):
    """
    Plot a venn diagram of 2 variables
    """

    cnx = mysql.connector.connect(user='quim', password='',
                                  host='localhost',
                                  database='BIANA_JAN_2017')
    dcdb2name = obtain_dcdb2name(cnx)


    # Create a figure of size 8x6 inches, 80 dots per inch
    fig = plt.figure(figsize=(8, 8), dpi=100)

    # Create a new subplot from a grid of 2x1 --> bar plot for the lost PDIs
    ax = fig.add_subplot(111)
    # s will contain the subset sizes #
    s = (len(items_drug1)- len(shared_items), len(items_drug2)- len(shared_items), len(shared_items))

    v = venn2(subsets=s, set_labels=(dcdb2name[dcdb_drug1], dcdb2name[dcdb_drug2]))

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

    plt.show()
    fig.savefig(output_plot, bbox_inches='tight')

    return

if  __name__ == "__main__":
    main()