import argparse
import copy
import mysql.connector
import numpy as np
import time
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug
import diana.classes.comparison as comparison
import diana.classes.network_analysis as network_analysis




def main():

    options = parse_user_arguments()
    compare_profiles(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Compare the profiles of the input drugs",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-d1','--drug_name1',dest='drug_name1',action = 'store',
                        help = """ Name of the drug number 1. If you do not provide targets for this drug or the number of targets is not large enough,
                        the program will use this name to search for targets in BIANA database. If targets are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between 
                        double quotes. """)
    parser.add_argument('-d2','--drug_name2',dest='drug_name2',action = 'store',
                        help = """ Name of the drug number 2. If you do not provide targets for this drug or the number of targets is not large enough,
                        the program will use this name to search for targets in BIANA database. If targets are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between 
                        double quotes. """)
    parser.add_argument('-t1','--targets1',dest='targets1',action = 'store',
                        help = 'Input file with the targets of the drug 1. Each target must be separated by a newline character.')
    parser.add_argument('-t2','--targets2',dest='targets2',action = 'store',
                        help = 'Input file with the targets of the drug 2. Each target must be separated by a newline character.')
    parser.add_argument('-pt','--proteins_type_id',dest='proteins_type_id',action = 'store', default='geneid',
                        help = 'Input the type of ID of the targets introduced / proteins of the network. It must be the same! (default is geneid).')
    parser.add_argument('-th','--threshold_list',dest='threshold_list',action = 'store',
                        help = """List of percentages that will be used as cut-offs to define the profiles of the drugs. It has to be a file containing:
                        - Different numbers that will be the threshold values separated by newline characters. 
                        For example, a file called "top_threshold.list" containing:
                        0.1
                        0.5
                        1
                        5
                        10
                        """)
    parser.add_argument('-ws','--worspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory will be created""")
    parser.add_argument('-db','--database',dest='database',action = 'store',default='BIANA_JUN_2017',
                        help = """Define the database to use for the search of targets: 
                        (default is BIANA_JUN_2017)""")
    parser.add_argument('-dbu','--db_user',dest='db_user',action = 'store',default='quim',
                        help = """Define the MySQL user to access to the database: 
                        (default is quim)""")
    parser.add_argument('-dbp','--db_pass',dest='db_pass',action = 'store',default='',
                        help = """Define the MySQL password to access to the database: 
                        (default is '')""")
    parser.add_argument('-dbh','--db_host',dest='db_host',action = 'store',default='localhost',
                        help = """Define the MySQL host to access to the database: 
                        (default is localhost)""")
    parser.add_argument('-up','--unification',dest='unification_protocol',action = 'store',default='geneid_seqtax_v1',
                        help = """Define the unification protocol used in BIANA database (default is BIANA_JUN_2017)""")

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def compare_profiles(options):
    """
    Compares the profiles of two input drugs
    """

    # Start marker for time measure
    start = time.time()

    print("\n\t\t------------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Second part: Comparison of drug profiles\n")
    print("\t\t------------------------------------------------------------------------------------------------------------------------\n")

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')

    # Check that the profiles directory exists
    data_dir = os.path.join(options.workspace, "profiles")
    check_directory(data_dir)

    # Create a directory for the results of the comparison
    results_dir = os.path.join(options.workspace, "comparisons")
    create_directory(results_dir)

    # Create drug instances
    drug1_instance = diana_drug.Drug(options.drug_name1)
    drug2_instance = diana_drug.Drug(options.drug_name2)

    # Create a ComparisonResult instance
    comparison_instance = comparison.ComparisonResult()



    #------------------------#
    #   TARGETS CONTROLLER   #
    #------------------------#

    # TARGETS CONTROLLER: Checks the targets provided by the user. If necessary, performs a search 
    # in BIANA database to obtain more targets
    # Gets the ID of the drug using the drug_name and the targets provided. 
    # The ID is used to find the folder containing the results

    for drug_instance, targets, num_drug in [(drug1_instance, options.targets1, 1), (drug2_instance, options.targets2, 2)]:

        # Check if the targets file is provided
        if targets and os.path.isfile(targets):
            drug_instance.obtain_targets_from_file(targets, options.proteins_type_id)
        else:
            # Create a connection to BIANA database
            biana_cnx = mysql.connector.connect(user=options.db_user, password=options.db_pass,
                                                host=options.db_host,
                                                database=options.database)
            # Obtain the targets from BIANA
            drug_instance.obtain_targets_from_BIANA(biana_cnx, options.proteins_type_id, options.unification_protocol)
            biana_cnx.close()

        print( "  DIANA INFO:\tThe targets provided for the drug {} are:\n\t\t{}.\n".format( drug_instance.drug_name, ', '.join([ str(x) for x in drug_instance.targets]) ) )

        # Create a directory for the drug
        drug_id = diana_drug.generate_drug_id(drug_instance.drug_name, drug_instance.targets)
        print('  DIANA INFO:\tThe ID for drug {} is: {}\n'.format(drug_instance.drug_name, drug_id))
        if num_drug == 1:
            drug1_dir = drug_dir = os.path.join(data_dir, drug_id)
            check_directory(drug1_dir)
        elif num_drug == 2:
            drug2_dir = drug_dir = os.path.join(data_dir, drug_id)
            check_directory(drug2_dir)



    #--------------------------------------#
    #   COMPARISON OF DCTARGETS PROFILES   #
    #--------------------------------------#

    print('  DIANA INFO:\tSTARTING COMPARISON OF dcTARGETS PROFILES\n')

    # Directories for the dcTargets results
    dctargets_dir1 = os.path.join(drug1_dir, 'dctargets_profiles')
    dctargets_dir2 = os.path.join(drug2_dir, 'dctargets_profiles')

    pfam_file1 = os.path.join(dctargets_dir1, 'pfam_profile.txt')
    pfam_file2 = os.path.join(dctargets_dir2, 'pfam_profile.txt')

    targets_functional_file1 = os.path.join(dctargets_dir1, 'targets_functional_profile.txt')
    targets_functional_file2 = os.path.join(dctargets_dir2, 'targets_functional_profile.txt')

    # Get PFAMs from PFAM profiles
    check_file(pfam_file1)
    drug1_instance.obtain_pfams_from_file(pfam_file1)

    check_file(pfam_file2)
    drug2_instance.obtain_pfams_from_file(pfam_file2)

    # Obtain functional profile instances
    check_file(targets_functional_file1)
    targets_functions_instance1 = network_analysis.FunctionalProfile(targets_functional_file1, 'targets', 'targets')

    check_file(targets_functional_file2)
    targets_functions_instance2 = network_analysis.FunctionalProfile(targets_functional_file2, 'targets', 'targets')

    # Compare targets
    targets1_dict = comparison.generate_targets_dict_for_comparison(drug1_instance.targets)
    targets2_dict = comparison.generate_targets_dict_for_comparison(drug2_instance.targets)
    summary_targets = comparison.calculate_comparison(targets1_dict, targets2_dict)
    comparison_instance.add_dctargets_result('target', summary_targets)
    print(summary_targets)

    # Compare PFAMs
    pfams1_dict = comparison.generate_targets_dict_for_comparison(drug1_instance.pfams)
    pfams2_dict = comparison.generate_targets_dict_for_comparison(drug2_instance.pfams)
    summary_pfams = comparison.calculate_comparison(pfams1_dict, pfams2_dict)
    comparison_instance.add_dctargets_result('pfam', summary_pfams)
    print(summary_pfams)

    # Compare functional profiles
    summary_target_functions = comparison.calculate_comparison(targets_functions_instance1.go_id_to_values, targets_functions_instance2.go_id_to_values)
    comparison_instance.add_dctargets_result('function', summary_target_functions)
    print(summary_target_functions)



    #------------------------------------#
    #   COMPARISON OF DCGUILD PROFILES   #
    #------------------------------------#

    print('  DIANA INFO:\tSTARTING COMPARISON OF dcGUILD PROFILES\n')

    # Directories for the dcGUILD results
    dcguild_dir1 = os.path.join(drug1_dir, 'dcguild_profiles')
    dcguild_dir2 = os.path.join(drug2_dir, 'dcguild_profiles')

    # Get the complete node profile (the pvalue_file)
    pvalue_file1 = os.path.join(dcguild_dir1, 'output_scores.sif.netcombo.pval')
    check_file(pvalue_file1)
    guild_profile_instance1 = network_analysis.GUILDProfile(pvalue_file1, options.proteins_type_id, 100)

    pvalue_file2 = os.path.join(dcguild_dir2, 'output_scores.sif.netcombo.pval')
    check_file(pvalue_file2)
    guild_profile_instance2 = network_analysis.GUILDProfile(pvalue_file2, options.proteins_type_id, 100)

    # Get the complete edge profile (the network_scored)
    scored_network_file1 = os.path.join(dcguild_dir1, 'network_scored.txt')
    check_file(scored_network_file1)
    scored_network_instance1 = network_analysis.EdgeProfile(scored_network_file1, None, options.proteins_type_id, 'sif', 100)

    scored_network_file2 = os.path.join(dcguild_dir2, 'network_scored.txt')
    check_file(scored_network_file2)
    scored_network_instance2 = network_analysis.EdgeProfile(scored_network_file2, None, options.proteins_type_id, 'sif', 100)



    # Get the list of thresholds to create the profiles
    if options.threshold_list and fileExist(options.threshold_list):
        threshold_list = get_values_from_threshold_file(options.threshold_list)
    else:
        threshold_list = [1, 5, 10, 20, 50]
    print('  DIANA INFO:\tList of percentages used to define the drug profiles: {}\n'.format(', '.join([str(x) for x in threshold_list])))

    # We will use the median threshold to output an example in the results
    comparison_instance.get_median_threshold(threshold_list)

    for top_threshold in threshold_list:

        # Obtain node profiles
        node_file1 = os.path.join(dcguild_dir1, 'node_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id.lower()))
        node_file2 = os.path.join(dcguild_dir2, 'node_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id.lower()))

        check_file(node_file1)
        node_profile_instance1 = network_analysis.GUILDProfile(node_file1, options.proteins_type_id, top_threshold)

        check_file(node_file2)
        node_profile_instance2 = network_analysis.GUILDProfile(node_file2, options.proteins_type_id, top_threshold)

        # Obtain edge profiles
        edge_file1 = os.path.join(dcguild_dir1, 'edge_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id.lower()))
        edge_file2 = os.path.join(dcguild_dir2, 'edge_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id.lower()))

        check_file(edge_file1)
        edge_profile_instance1 = network_analysis.EdgeProfile(edge_file1, None, options.proteins_type_id, 'sif', top_threshold)

        check_file(edge_file2)
        edge_profile_instance2 = network_analysis.EdgeProfile(edge_file2, None, options.proteins_type_id, 'sif', top_threshold)

        # Obtain functional profiles
        functional_file1 = os.path.join(dcguild_dir1, 'functional_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id.lower()))
        functional_file2 = os.path.join(dcguild_dir2, 'functional_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id.lower()))

        check_file(functional_file1)
        functional_profile_instance1 = network_analysis.FunctionalProfile(functional_file1, top_threshold, node_file1)

        check_file(functional_file2)
        functional_profile_instance2 = network_analysis.FunctionalProfile(functional_file2, top_threshold, node_file2)


        # Compare node profiles
        print('NODE PROFILES, THRESHOLD: {}'.format(top_threshold))
        node_profile1_values = copy.copy(node_profile_instance1.node_to_values)
        node_profile2_values = copy.copy(node_profile_instance2.node_to_values)
        guild_profile1_values = copy.copy(guild_profile_instance1.node_to_values)
        guild_profile2_values = copy.copy(guild_profile_instance2.node_to_values)
        summary_nodes = comparison.calculate_comparison_top_scoring(node_profile1_values, guild_profile1_values , node_profile2_values, guild_profile2_values)
        comparison_instance.add_dcguild_result('node', top_threshold, summary_nodes)
        print(summary_nodes)

        # Compare edge profiles
        print('EDGE PROFILES, THRESHOLD: {}'.format(top_threshold))
        edge_profile1_values = copy.copy(edge_profile_instance1.edge_to_values)
        edge_profile2_values = copy.copy(edge_profile_instance2.edge_to_values)
        scored_network1_values = copy.copy(scored_network_instance1.edge_to_values)
        scored_network2_values = copy.copy(scored_network_instance2.edge_to_values)
        summary_edges = comparison.calculate_comparison_top_scoring(edge_profile1_values, scored_network1_values , edge_profile2_values, scored_network2_values)
        comparison_instance.add_dcguild_result('edge', top_threshold, summary_edges)
        print(summary_edges)

        # Compare functional profiles
        print('FUNCTIONAL PROFILES, THRESHOLD: {}'.format(top_threshold))
        functional_profile1_values = copy.copy(functional_profile_instance1.go_id_to_values)
        functional_profile2_values = copy.copy(functional_profile_instance2.go_id_to_values)
        summary_functions = comparison.calculate_comparison(functional_profile1_values, functional_profile2_values)
        comparison_instance.add_dcguild_result('function', top_threshold, summary_functions)
        print(summary_functions)

        # Get the profiles to represent an example
        if top_threshold == comparison_instance.example_threshold:
            comparison_instance.node_examples = [node_profile_instance1, node_profile_instance2]
            comparison_instance.edge_examples = [edge_profile_instance1, edge_profile_instance2]
            comparison_instance.function_examples = [functional_profile_instance1, functional_profile_instance2]


    #----------------------------------------#
    #   COMPARISON OF DCSTRUCTURE PROFILES   #
    #----------------------------------------#

    print('  DIANA INFO:\tSTARTING COMPARISON OF dcSTRUCTURE PROFILES\n')

    # Directories for the dcStructure results
    dcstructure_dir1 = os.path.join(drug1_dir, 'dcstructure_profiles')
    dcstructure_dir2 = os.path.join(drug2_dir, 'dcstructure_profiles')

    structure_file1 = os.path.join(dcstructure_dir1, 'structure_profile.txt')
    structure_file2 = os.path.join(dcstructure_dir2, 'structure_profile.txt')

    # Get SMILES from structure profiles
    check_file(structure_file1)
    drug1_instance.obtain_SMILES_from_file(structure_file1)

    check_file(structure_file2)
    drug2_instance.obtain_SMILES_from_file(structure_file2)

    similarity_results = []
    if len(drug1_instance.smiles) > 0 and len(drug2_instance.smiles) > 0:
        for smiles1 in drug1_instance.smiles:
            for smiles2 in drug2_instance.smiles:
                similarity_result = comparison.get_smiles_similarity(smiles1, smiles2, fp_type = "sub", metric = "tanimoto")
                if similarity_result:
                    similarity_results.append(similarity_result)
        if len(similarity_results) > 0:
            similarity_results = np.mean(similarity_results)
            comparison_instance.dcstructure_result = similarity_results
            print('  DIANA INFO:\tStructural similarity between the two drugs: {:.3f}\n'.format(similarity_results))
        else:
            similarity_results = None
            print('  DIANA INFO:\tStructural similarity unavailable\n')
    else:
        similarity_results = None
        print('  DIANA INFO:\tThe SMILES of the drugs are missing! It is not possible to compute the structural similarity\n')



    #-------------------#
    #   WRITE RESULTS   #
    #-------------------#

    # Create a directory for the comparison of the two drugs
    comparison_id =  '{}_{}'.format(diana_drug.generate_drug_id(drug1_instance.drug_name, drug1_instance.targets), diana_drug.generate_drug_id(drug2_instance.drug_name, drug2_instance.targets))
    comparison_dir = os.path.join(results_dir, comparison_id)
    create_directory(comparison_dir)

    # Write the results table
    results_table = os.path.join(comparison_dir, 'results_table.tsv')
    comparison_instance.output_results_table(results_table, threshold_list)

    # Output the venn diagram of common targets
    venn_plot_targets = os.path.join(comparison_dir, 'venn_plot_targets.png')
    shared_targets = set(drug1_instance.targets) & set(drug2_instance.targets)
    comparison.plot_venn_2(drug1_instance.drug_name, drug2_instance.drug_name, set(drug1_instance.targets), set(drug2_instance.targets), shared_targets, venn_plot_targets)

    # Output the venn diagram of common nodes/edges/functions using median_threshold
    venn_plot_nodes = os.path.join(comparison_dir, 'venn_plot_nodes.png')
    venn_plot_edges = os.path.join(comparison_dir, 'venn_plot_edges.png')
    venn_plot_functions = os.path.join(comparison_dir, 'venn_plot_functions.png')
    comparison_instance.output_venn_diagram_dcguild(drug1_instance.drug_name, drug2_instance.drug_name, venn_plot_nodes, venn_plot_edges, venn_plot_functions)



    # End marker for time
    end = time.time()
    print('\n  DIANA INFO:\tTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))

    return



#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


def check_file(file):
    """
    Checks if a file exists and if not, raises FileNotFound exception
    """
    if not fileExist(file):
        raise FileNotFound(file)


def check_directory(directory):
    """
    Checks if a directory exists and if not, raises DirNotFound exception
    """
    try:
        os.stat(directory)
    except:
        raise DirNotFound(directory)


class FileNotFound(Exception):
    """
    Exception raised when a file is not found.
    """
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return 'The file {} has not been found.\nTherefore, the comparison cannot be performed. Please, check that all the profiles have been correctly generated.\n'.format(self.file)


class DirNotFound(Exception):
    """
    Exception raised when a directory is not found.
    """
    def __init__(self, directory):
        self.directory = directory

    def __str__(self):
        return 'The directory {} has not been found.\nTherefore, the comparison cannot be performed. Please, check that all the parameters have been correctly introduced and the profiles have been correctly generated.\n'.format(self.directory)


if  __name__ == "__main__":
    main()
