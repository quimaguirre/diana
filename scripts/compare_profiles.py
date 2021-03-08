import argparse
import configparser
import copy
import ntpath
import numpy as np
import pandas as pd
import time
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug
import diana.classes.comparison as comparison
import diana.classes.network_analysis as network_analysis
import diana.classes.functional_analysis as functional_analysis




def main():
    """
    Generate profiles for drugs using GUILD.
    Optimized for Python 3.
    python /home/quim/PHD/Projects/DIANA/diana/scripts/compare_profiles.py -d1 DB11699 -d2 DB00177 -sif /home/quim/PHD/Projects/DIANA/diana/data/network_cheng.sif
    """

    options = parse_user_arguments()
    compare_profiles(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Compare the profiles of the input drugs",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-j1','--job_id_drug1',dest='job_id_drug1',action = 'store',
                        help = """ Identifier of the drug number 1. 
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between
                        double quotes. """)
    parser.add_argument('-j2','--job_id_drug2',dest='job_id_drug2',action = 'store',
                        help = """ Identifier of the drug number 2.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between
                        double quotes. """)
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = 'Input file with a protein-protein interaction network in SIF format.')
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
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory will be created""")

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

    # Get the script path and define directories used
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    scripts_dir = os.path.join(main_path, 'scripts')
    mappings_dir = os.path.join(main_path, 'mappings')
    data_dir = os.path.join(main_path, 'data')
    workspace_dir = options.workspace
    create_directory(workspace_dir)
    profiles_dir = os.path.join(workspace_dir, 'profiles')

    # Create a directory for the results of the comparison
    results_dir = os.path.join(workspace_dir, "comparisons")
    create_directory(results_dir)

    # Create directories for additional data
    other_data_dir = os.path.join(workspace_dir, 'additional_data')
    create_directory(other_data_dir)
    random_networks_dir = os.path.join(other_data_dir, 'random_networks')
    create_directory(random_networks_dir)
    associations_dir = os.path.join(other_data_dir, 'gene_function_associations')
    create_directory(associations_dir)
    numbers_dir = os.path.join(other_data_dir, 'numbers')
    create_directory(numbers_dir)


    # Create a ComparisonResult instance
    comparison_instance = comparison.ComparisonResult()


    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_file)


    #--------------------#
    #   SIF CONTROLLER   #
    #--------------------#

    # SIF CONTROLLER: Checks the network in SIF format provided by the user.

    # Check if the network file is provided
    if options.sif and fileExist(options.sif):
        # Get the network name
        network_filename = ntpath.basename(options.sif)
    else:
        # If not, we output an error
        print('  DIANA INFO:\tThe network SIF file is missing. Please, introduce the parameter -sif.\n\t\tIf you do not have a network, use one of the networks in the sif folder.\n')
        sys.exit(10)


    #------------------------------#
    #   CREATE/READ NUMBERS FILE   #
    #------------------------------#

    # Define parameters for the functional enrichment
    type_functions = ['gobp', 'gomf', 'reactome']
    type_corrections = ['fdr_bh', 'bonferroni']

    # Numbers file associated to all targets
    target_numbers_file = os.path.join(numbers_dir, 'target_numbers.txt')
    if not fileExist(target_numbers_file):
        with open(target_numbers_file, 'w') as num_fd:
            num_fd.write('#feature\tnumber\n')
            # Get targets
            drugbank_geneid_mapping_file = os.path.join(mappings_dir, 'drugbank_geneid_drug_target_interactions.txt')
            targets = diana_drug.get_all_targets_from_mappings(drugbank_geneid_mapping_file)
            num_fd.write('target\t{}\n'.format(len(targets)))
            # Get PFAMs
            geneid_target_mapping_file = os.path.join(mappings_dir, 'geneid_target_mappings.txt')
            pfams = diana_drug.get_all_pfams_from_mappings(geneid_target_mapping_file)
            num_fd.write('pfam\t{}\n'.format(len(pfams)))
            # Get functions
            target_associations_dir = os.path.join(associations_dir, 'targets')
            for type_function in type_functions:
                associations_file = os.path.join(target_associations_dir, '{}_to_gene.txt'.format(type_function))
                functions = functional_analysis.get_functions_from_associations_file(associations_file)
                num_fd.write('{}\t{}\n'.format(type_function, len(functions)))
            # Get ATCs
            drugbank_atc_file = os.path.join(mappings_dir, 'drugbank_drug_atc.txt')
            atcs = diana_drug.get_all_atcs_from_mappings(drugbank_atc_file)
            num_fd.write('atc\t{}\n'.format(len(atcs)))
            # Get SEs
            drugbank_side_effects_file = os.path.join(mappings_dir, 'drugbank_drug_side_effects.txt')
            ses = diana_drug.get_all_ses_from_mappings(drugbank_side_effects_file)
            num_fd.write('se\t{}\n'.format(len(ses)))
    else:
        target_numbers_df = pd.read_csv(target_numbers_file, sep='\t', index_col=None)

    # Numbers file associated to network
    network_numbers_file = os.path.join(numbers_dir, '{}_numbers.txt'.format(network_filename))
    if not fileExist(network_numbers_file):
        with open(network_numbers_file, 'w') as num_fd:
            num_fd.write('#feature\tnumber\n')
            # We create a Network instance
            network_instance = network_analysis.Network(network_file=options.sif, type_id='geneid', network_format='sif')
            # We keep the number of nodes
            num_fd.write('node\t{}\n'.format(len(network_instance.get_nodes())))
            num_fd.write('edge\t{}\n'.format(len(network_instance.get_edges())))
            # Get functions
            network_associations_dir = os.path.join(associations_dir, network_filename)
            for type_function in type_functions:
                associations_file = os.path.join(network_associations_dir, '{}_to_gene.txt'.format(type_function))
                functions = functional_analysis.get_functions_from_associations_file(associations_file)
                num_fd.write('{}\t{}\n'.format(type_function, len(functions)))
    else:
        network_numbers_df = pd.read_csv(network_numbers_file, sep='\t', index_col=None)


    #-------------------#
    #   READ PROFILES   #
    #-------------------#

    print('  DIANA INFO:\tREADING PROFILES\n')

    # Get the list of thresholds to create the profiles
    if options.threshold_list and fileExist(options.threshold_list):
        threshold_list = get_values_from_threshold_file(options.threshold_list)
    else:
        threshold_list = [1, 2, 5, 'functions']
    print('  DIANA INFO:\tList of percentages used to define the drug profiles: {}\n'.format(', '.join([str(x) for x in threshold_list])))

    # Check if the directories of the drugs exist
    if options.job_id_drug1:
        drug_dir1 = os.path.join(profiles_dir, options.job_id_drug1)
        check_directory(drug_dir1)
    else:
        print('  DIANA INFO:\tjob_id_drug1 parameter is missing. Please, introduce the parameter -j1 with the job identifier of the drug.\n')
        sys.exit(10)

    if options.job_id_drug2:
        drug_dir2 = os.path.join(profiles_dir, options.job_id_drug2)
        check_directory(drug_dir2)
    else:
        print('  DIANA INFO:\tjob_id_drug2 parameter is missing. Please, introduce the parameter -j2 with the job identifier of the drug.\n')
        sys.exit(10)

    # Read profiles for drug 1
    drug_instance1, guild_profile_instance1, scored_network_instance1, target_function_results1, guild_results1 = read_drug_profiles(drug_dir=drug_dir1, threshold_list=threshold_list)
    # Read profiles for drug 2
    drug_instance2, guild_profile_instance2, scored_network_instance2, target_function_results2, guild_results2 = read_drug_profiles(drug_dir=drug_dir2, threshold_list=threshold_list)


    #----------------------#
    #   COMPARE PROFILES   #
    #----------------------#

    # Compare targets
    targets_dict1 = comparison.generate_targets_dict_for_comparison(drug_instance1.targets)
    targets_dict2 = comparison.generate_targets_dict_for_comparison(drug_instance2.targets)
    num_targets = int(target_numbers_df[target_numbers_df['#feature'] == 'target']['number'])
    summary_targets = comparison.calculate_comparison(targets_dict1, targets_dict2, num_targets)
    comparison_instance.add_target_result('target', summary_targets)
    print(summary_targets)

    # Compare PFAMs
    pfams_dict1 = comparison.generate_targets_dict_for_comparison(drug_instance1.pfams)
    pfams_dict2 = comparison.generate_targets_dict_for_comparison(drug_instance2.pfams)
    num_pfams = int(target_numbers_df[target_numbers_df['#feature'] == 'pfam']['number'])
    summary_pfams = comparison.calculate_comparison(pfams_dict1, pfams_dict2, num_pfams)
    comparison_instance.add_target_result('pfam', summary_pfams)
    print(summary_pfams)

    # Compare functional profiles
    for type_function in type_functions:
        num_functions = int(target_numbers_df[target_numbers_df['#feature'] == type_function]['number'])
        for type_correction in type_corrections:
            targets_functions_instance1 = target_function_results1['{}-{}'.format(type_function, type_correction)]
            targets_functions_instance2 = target_function_results2['{}-{}'.format(type_function, type_correction)]
            summary_target_functions = comparison.calculate_comparison(targets_functions_instance1.term_id_to_pvalue, targets_functions_instance2.term_id_to_pvalue, num_functions)
            comparison_instance.add_target_result('{}-{}'.format(type_function, type_correction), summary_target_functions)
            print('Target: {} - {}'.format(type_function, type_correction))
            print(summary_target_functions)

    # Compare GUILD profiles
    for top_threshold in threshold_list:

        if top_threshold == 'functions':
            for type_function in type_functions:
                num_functions = int(target_numbers_df[target_numbers_df['#feature'] == type_function]['number'])
                for type_correction in type_corrections:
                    # Compare node profiles
                    print('NODE PROFILES, THRESHOLD: {} - {} - {}'.format(top_threshold, type_function, type_correction))
                    node_profile1_values = copy.copy(guild_results1['node-{}-{}-{}'.format(top_threshold, type_function, type_correction)].node_to_score)
                    node_profile2_values = copy.copy(guild_results2['node-{}-{}-{}'.format(top_threshold, type_function, type_correction)].node_to_score)
                    guild_profile1_values = copy.copy(guild_profile_instance1.node_to_score)
                    guild_profile2_values = copy.copy(guild_profile_instance2.node_to_score)
                    summary_nodes = comparison.calculate_comparison_top_scoring(node_profile1_values, guild_profile1_values, node_profile2_values, guild_profile2_values)
                    comparison_instance.add_dcguild_result('node', top_threshold, summary_nodes)
                    print(summary_nodes)
        else:
            # Compare node profiles
            print('NODE PROFILES, THRESHOLD: {}'.format(top_threshold))
            node_profile1_values = copy.copy(guild_results1['node-{}'.format(top_threshold)].node_to_score)
            node_profile2_values = copy.copy(guild_results2['node-{}'.format(top_threshold)].node_to_score)
            guild_profile1_values = copy.copy(guild_profile_instance1.node_to_score)
            guild_profile2_values = copy.copy(guild_profile_instance2.node_to_score)
            summary_nodes = comparison.calculate_comparison_top_scoring(node_profile1_values, guild_profile1_values, node_profile2_values, guild_profile2_values)
            comparison_instance.add_dcguild_result('node', top_threshold, summary_nodes)
            print(summary_nodes)

    sys.exit(0)

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

    # Compare structures
    similarity_results = []
    if len(drug_instance1.smiles) > 0 and len(drug_instance2.smiles) > 0:
        for smiles1 in drug_instance1.smiles:
            for smiles2 in drug_instance2.smiles:
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

    # Compare ATC profiles
    drug_instance1.ATCs = set([x[0] for x in drug_instance1.ATCs])
    drug_instance2.ATCs = set([x[0] for x in drug_instance2.ATCs])
    ATCs1_dict = comparison.generate_targets_dict_for_comparison(drug_instance1.ATCs)
    ATCs2_dict = comparison.generate_targets_dict_for_comparison(drug_instance2.ATCs)
    print('  DIANA INFO:\tComparing ATCs {} of drug1 and ATcs {} of drug2\n'.format(', '.join(drug_instance1.ATCs), ', '.join(drug_instance2.ATCs)))
    summary_ATCs = comparison.calculate_comparison(ATCs1_dict, ATCs2_dict)
    comparison_instance.dcatc_result = summary_ATCs
    print(summary_ATCs)

    # Compare SE profiles
    SEs1_dict = comparison.generate_targets_dict_for_comparison(drug_instance1.SEs)
    SEs2_dict = comparison.generate_targets_dict_for_comparison(drug_instance2.SEs)
    summary_SEs = comparison.calculate_comparison(SEs1_dict, SEs2_dict)
    comparison_instance.dcse_result = summary_SEs
    print(summary_SEs)


    #-------------------#
    #   WRITE RESULTS   #
    #-------------------#

    # Create a directory for the comparison of the two drugs
    comparison_id =  '{}---{}'.format(diana_drug.generate_drug_id(drug_instance1.drug_name, drug_instance1.targets, network_filename), diana_drug.generate_drug_id(drug_instance2.drug_name, drug_instance2.targets, network_filename))
    comparison_dir = os.path.join(results_dir, comparison_id)
    create_directory(comparison_dir)

    # Write the results table
    results_table = os.path.join(comparison_dir, 'results_table.tsv')
    comparison_instance.output_results_table(results_table, threshold_list)

    if use_server == 'false':

        # Output the venn diagram of common targets
        venn_plot_targets = os.path.join(comparison_dir, 'venn_plot_targets.png')
        shared_targets = set(drug_instance1.targets) & set(drug_instance2.targets)
        comparison.plot_venn_2(drug_instance1.drug_name, drug_instance2.drug_name, set(drug_instance1.targets), set(drug_instance2.targets), shared_targets, venn_plot_targets)
        # Output the venn diagram of common target-functions
        venn_plot_func_tar = os.path.join(comparison_dir, 'venn_plot_functions-targets.png')
        targets_functions_instance1 = network_analysis.FunctionalProfile(targets_functional_file1, 'targets', 'targets')
        targets_functions_instance2 = network_analysis.FunctionalProfile(targets_functional_file2, 'targets', 'targets')
        functions_targets_drug1 = targets_functions_instance1.go_id_to_values
        functions_targets_drug2 = targets_functions_instance2.go_id_to_values
        shared_func_tar = set(functions_targets_drug1) & set(functions_targets_drug2)
        comparison.plot_venn_2(drug_instance1.drug_name, drug_instance2.drug_name, set(functions_targets_drug1), set(functions_targets_drug2), shared_func_tar, venn_plot_func_tar)

        # Output the venn diagram of common nodes/edges/functions using median_threshold
        venn_plot_nodes = os.path.join(comparison_dir, 'venn_plot_nodes.png')
        venn_plot_edges = os.path.join(comparison_dir, 'venn_plot_edges.png')
        venn_plot_functions = os.path.join(comparison_dir, 'venn_plot_functions.png')
        comparison_instance.output_venn_diagram_dcguild(drug_instance1.drug_name, drug_instance2.drug_name, venn_plot_nodes, venn_plot_edges, venn_plot_functions)



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


def get_targets_in_sif_file(sif_file, targets):
    """
    Get the targets that are inside the network given by the user
    """
    targets_in_network = set()
    str_tar = [str(x) for x in targets]
    with open(sif_file, 'r') as sif_fd:
        for line in sif_fd:
            node1, score, node2 = line.strip().split('\t')
            if node1 in str_tar:
                targets_in_network.add(node1)
            if node2 in str_tar:
                targets_in_network.add(node2)
    return list(targets_in_network)


def read_parameters_file(parameters_file):
    """
    Reads the parameters file of a drug profile
    """
    with open(parameters_file, 'r') as parameters_fd:
        header = parameters_fd.readline()
        content = parameters_fd.readline()
        fields = content.strip().split('\t')
    return fields


def read_drug_profiles(drug_dir, threshold_list=[1, 2, 5, 'functions']):
    """
    Read the profiles of a drug.
    """
    # Check/Read parameters file
    output_parameters_file = os.path.join(drug_dir, 'parameters.txt')
    check_file(output_parameters_file)
    parameters = read_parameters_file(output_parameters_file)
    drugname = parameters[1]

    # Create drug instance
    drug_instance = diana_drug.Drug(drugname)

    # Read target profile
    target_dir = os.path.join(drug_dir, 'target_profiles')
    target_file = os.path.join(target_dir, '{}_targets.txt'.format(drugname))
    check_file(target_file)
    drug_instance.obtain_targets_from_file(target_file, target_type_id='geneid')
    print('  DIANA INFO:\tTARGETS OF {}: {}\n'.format(drugname, ', '.join(drug_instance.targets)))

    # Read PFAM profile
    pfam_file = os.path.join(target_dir, 'pfam_profile.txt')
    check_file(pfam_file)
    drug_instance.obtain_pfams_from_file(pfam_file)

    # Read functional profiles
    type_functions = ['gobp', 'gomf', 'reactome']
    type_corrections = ['fdr_bh', 'bonferroni']
    target_function_results = {}
    for type_function in type_functions:
        for type_correction in type_corrections:
            targets_functional_file = os.path.join(target_dir, 'targets_functional_profile_{}_{}.txt'.format(type_function, type_correction))
            check_file(targets_functional_file)
            targets_functions_instance = network_analysis.FunctionalProfile(targets_functional_file, 'targets', 'targets')
            target_function_results['{}-{}'.format(type_function, type_correction)] = targets_functions_instance

    # Read GUILD node scores
    guild_dir = os.path.join(drug_dir, 'guild_profiles')
    scores_file = os.path.join(guild_dir, 'output_scores.sif.netcombo')
    check_file(scores_file)
    guild_profile_instance = network_analysis.GUILDProfile(scores_file, type_id='geneid', top=100, top_type='percentage')
    drug_instance.targets_in_network = set([target for target in drug_instance.targets if target in guild_profile_instance.node_to_score.keys()])

    # Read GUILD edge scores
    scored_network_file = os.path.join(guild_dir, 'network_scored.txt')
    check_file(scored_network_file)
    scored_network_instance = network_analysis.EdgeProfile(network_file=scored_network_file, type_id='geneid', network_format='sif', top=100)

    guild_results = {}
    for top_threshold in threshold_list:
        if top_threshold == 'functions':
            # Read profiles associated to a functions threshold
            for type_function in type_functions:
                for type_correction in type_corrections:
                    # Obtain cut-off
                    output_sliding_window_file = os.path.join(guild_dir, 'sliding_window_{}_{}.txt'.format(type_function, type_correction))
                    cutoff_central_position, cutoff_right_interval = functional_analysis.read_sliding_window_file(output_sliding_window_file=output_sliding_window_file, num_seeds=len(drug_instance.targets_in_network))
                    print('{} - {} - {} - {}: Cut-off central position: {}. Cut-off right interval position: {}'.format(drug_instance.drug_name, top_threshold, type_function, type_correction, cutoff_central_position, cutoff_right_interval))
                    # Obtain node profile
                    node_file = os.path.join(guild_dir, 'node_profile_top_{}_{}_{}_{}.txt'.format('functions', type_function, type_correction, guild_profile_instance.type_id))
                    check_file(node_file)
                    node_profile_instance = network_analysis.GUILDProfile(scores_file=node_file, type_id=guild_profile_instance.type_id, top=cutoff_right_interval, top_type='number_of_nodes')
                    guild_results['node-{}-{}-{}'.format(top_threshold, type_function, type_correction)] = node_profile_instance
                    # Obtain edge profile
                    edge_file = os.path.join(guild_dir, 'edge_profile_top_{}_{}_{}_{}.txt'.format('functions', type_function, type_correction, guild_profile_instance.type_id))
                    check_file(edge_file)
                    edge_profile_instance = network_analysis.EdgeProfile(network_file=edge_file, type_id=guild_profile_instance.type_id, network_format='sif', top=cutoff_right_interval, top_type='number_of_nodes')
                    guild_results['edge-{}-{}-{}'.format(top_threshold, type_function, type_correction)] = edge_profile_instance
                    # Obtain functional profile
                    function_file = os.path.join(guild_dir, 'functional_profile_top_{}_{}_{}.txt'.format('functions', type_function, type_correction))
                    functional_profile_instance = network_analysis.FunctionalProfile(functional_file=function_file, top=cutoff_right_interval, node_file=node_file)
                    guild_results['functional-{}-{}-{}'.format(top_threshold, type_function, type_correction)] = functional_profile_instance
        else:
            # Obtain node profile
            node_file = os.path.join(guild_dir, 'node_profile_top_{}_{}.txt'.format(str(top_threshold), guild_profile_instance.type_id))
            check_file(node_file)
            node_profile_instance = network_analysis.GUILDProfile(node_file, type_id=guild_profile_instance.type_id, top=top_threshold, top_type='percentage')
            guild_results['node-{}'.format(top_threshold)] = node_profile_instance

            # Obtain edge profiles
            edge_file = os.path.join(guild_dir, 'edge_profile_top_{}_{}.txt'.format(str(top_threshold), guild_profile_instance.type_id))
            check_file(edge_file)
            edge_profile_instance = network_analysis.EdgeProfile(network_file=edge_file, type_id=guild_profile_instance.type_id, network_format='sif', top=top_threshold, top_type='percentage')
            guild_results['edge-{}'.format(top_threshold)] = edge_profile_instance

            # Obtain functional profiles
            for type_function in type_functions:
                for type_correction in type_corrections:
                    functional_file = os.path.join(guild_dir, 'functional_profile_top_{}_{}_{}.txt'.format(str(top_threshold), type_function, type_correction))
                    check_file(functional_file)
                    functional_profile_instance = network_analysis.FunctionalProfile(functional_file=functional_file, top=top_threshold, node_file=node_file)
                    guild_results['functional-{}-{}-{}'.format(top_threshold, type_function, type_correction)] = functional_profile_instance

    # Read structures
    structure_file = os.path.join(drug_dir, 'structure_profiles/structure_profile.txt')
    check_file(structure_file)
    drug_instance.obtain_SMILES_from_file(structure_file)

    # Read ATCs
    atc_file = os.path.join(drug_dir, 'atc_profiles/ATC_profile.txt')
    check_file(atc_file)
    drug_instance.obtain_ATCs_from_file(atc_file)

    # Read side effects
    se_file = os.path.join(drug_dir, 'se_profiles/SE_profile.txt')
    check_file(se_file)
    drug_instance.obtain_SE_from_file(se_file)

    return drug_instance, guild_profile_instance, scored_network_instance, target_function_results, guild_results


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
