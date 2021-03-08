import argparse
import configparser
import ntpath
import shutil
import time
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug
import diana.classes.network_analysis as network_analysis
import diana.classes.top_scoring as top_scoring
import diana.classes.functional_analysis as functional_analysis


def main():
    """
    Generate profiles for drugs using GUILD.
    Optimized for Python 3.
    python /home/quim/PHD/Projects/DIANA/diana/scripts/generate_profiles.py -j entresto -d entresto -sif /home/quim/PHD/Projects/DIANA/diana/data/network_cheng.sif
    python /home/quim/PHD/Projects/DIANA/diana/scripts/generate_profiles.py -d DB11699 -sif /home/quim/PHD/Projects/DIANA/diana/data/network_cheng.sif
    """

    options = parse_user_arguments()
    generate_profiles(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of the input drug",
        epilog      = "@oliva's lab 2020")
    parser.add_argument('-j','--job_id',dest='job_id',action = 'store',
                        help = 'Identifier of the job. It will be used to create a directory with the name of the identifier to store the results')
    parser.add_argument('-d','--drug_name',dest='drug_name',action = 'store',
                        help = """ Name of the drug. If you do not provide targets for this drug or the number of targets is not large enough,
                        the program will use this name to search for targets in BIANA database. If targets are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between
                        double quotes. """)
    parser.add_argument('-t','--targets',dest='targets',action = 'store',
                        help = 'Input file with the targets of the drug. Each target must be separated by a newline character.')
    parser.add_argument('-pt','--proteins_type_id',dest='proteins_type_id',action = 'store', default='geneid',
                        help = 'Input the type of ID of the targets introduced / proteins of the network. It must be the same! (default is geneid).')
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
                        functions
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

def generate_profiles(options):
    """
    Generates the profiles of the input drug
    """

    # Start marker for time measure
    start = time.time()

    print("\n\t\t-----------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. First part: Generation of drug profiles\n")
    print("\t\t-----------------------------------------------------------------------------------------------------------------------\n")

    # Get the script path and define directories used
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    scripts_dir = os.path.join(main_path, 'scripts')
    mappings_dir = os.path.join(main_path, 'mappings')
    data_dir = os.path.join(main_path, 'data')
    workspace_dir = options.workspace
    create_directory(workspace_dir)

    # Create a directory for the data
    profiles_dir = os.path.join(workspace_dir, 'profiles')
    create_directory(profiles_dir)

    # Create directories for additional data
    other_data_dir = os.path.join(workspace_dir, 'additional_data')
    create_directory(other_data_dir)
    random_networks_dir = os.path.join(other_data_dir, 'random_networks')
    create_directory(random_networks_dir)
    associations_dir = os.path.join(other_data_dir, 'gene_function_associations')
    create_directory(associations_dir)

    # Create a drug instance
    if options.drug_name:
        drug_instance = diana_drug.Drug(options.drug_name)
    else:
        print('  DIANA INFO:\tdrug_name parameter is missing. Please, introduce the parameter -d with the name of the drug.\n')
        sys.exit(10)


    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_file)


    #------------------------#
    #   TARGETS CONTROLLER   #
    #------------------------#

    # TARGETS CONTROLLER: Checks the targets provided by the user. If necessary, performs a search
    # in BIANA database to obtain more targets

    # Check if the targets file is provided
    if options.targets and fileExist(options.targets):
        drug_instance.obtain_targets_from_file(options.targets, options.proteins_type_id)
    else:
        # Obtain the targets from a table
        drug_to_targets_file = os.path.join(mappings_dir, 'drugbank_geneid_drug_to_targets.txt')
        drug_mapping_file = os.path.join(mappings_dir, 'drugbank_drug_mappings.txt')
        if fileExist(drug_to_targets_file) and fileExist(drug_mapping_file):
            # First, translate the drug input name to drugbankid (if necessary)
            drugbankids = drug_instance.obtain_drugbankids_from_table(drug_mapping_file)
            # Then, get the targets (in geneid) from the table
            drug_instance.obtain_targets_from_table(drugbankids, drug_to_targets_file, target_type_id='geneid')
        else:
            if not fileExist(drug_to_targets_file):
                print("  DIANA INFO:\tMissing drug to targets file: {}.\n".format(drug2targets_file))
            if not fileExist(drug_mapping_file):
                print("  DIANA INFO:\tMissing drug mappings file: {}.\n".format(drug_mapping_file))
            sys.exit(10)

    print( "  DIANA INFO:\tThe targets provided for the drug {} are:\n\t\t{}.\n".format( options.drug_name, ', '.join([ str(x) for x in drug_instance.targets]) ) )


    #--------------------#
    #   SIF CONTROLLER   #
    #--------------------#

    # SIF CONTROLLER: Checks the network in SIF format provided by the user.

    # Check if the network file is provided
    if options.sif and fileExist(options.sif):
        # If the network file is provided, we create a Network instance
        network_instance = network_analysis.Network(network_file=options.sif, type_id=options.proteins_type_id, network_format='sif')
        # We search for the targets in the network
        drug_instance.targets_in_network = network_instance.get_targets_in_network(drug_instance.targets)
        # We create a directory in the random networks directory for this network
        network_filename = ntpath.basename(options.sif)
        random_networks_dir = os.path.join(random_networks_dir, network_filename)
        create_directory(random_networks_dir)
        # We create a directory in gene function associations directory for this network
        network_associations_dir = os.path.join(associations_dir, network_filename)
        create_directory(network_associations_dir)
        # We create a directory in gene function associations directory for targets
        target_associations_dir = os.path.join(associations_dir, 'targets')
        create_directory(target_associations_dir)
    else:
        # If not, we output an error
        print('  DIANA INFO:\tThe network SIF file is missing. Please, introduce the parameter -sif.\n\t\tIf you do not have a network, use one of the networks in the sif folder.\n')
        sys.exit(10)


    # Check if the number of targets provided is sufficient for the analysis
    if len(drug_instance.targets_in_network) < 1:
        raise diana_drug.InsufficientTargets(drug_instance.targets_in_network, 1)
    else:
        print( "  DIANA INFO:\tThe targets found in the network are:\n\t\t{}.\n".format( ', '.join([ str(x) for x in drug_instance.targets_in_network]) ) )


    #------------------------------------------#
    #   CREATE DIRECTORIES AND GENERAL FILES   #
    #------------------------------------------#

    # Create a directory for the drug
    if options.job_id and options.job_id != '':
        drug_id = options.job_id
    else:
        drug_id = diana_drug.generate_diana_id(drug_instance.drug_name, drug_instance.targets, network_filename)
    drug_dir = os.path.join(profiles_dir, drug_id)
    create_directory(drug_dir)
    print('  DIANA INFO:\tThe ID given to the drug, which will be used to create a directory and store the results, is: {}\n'.format(drug_id))


    # Create a directory for the Target results
    targets_dir = os.path.join(drug_dir, 'target_profiles')
    create_directory(targets_dir)

    # Create a directory for the GUILD results
    guild_dir = os.path.join(drug_dir, 'guild_profiles')
    create_directory(guild_dir)

    # Create a directory for the Structure results
    structure_dir = os.path.join(drug_dir, 'structure_profiles')
    create_directory(structure_dir)

    # Create a directory for the ATCs results
    atc_dir = os.path.join(drug_dir, 'atc_profiles')
    create_directory(atc_dir)

    # Create a directory for the dcse results
    se_dir = os.path.join(drug_dir, 'se_profiles')
    create_directory(se_dir)

    # Create a targets file
    targets_file = os.path.join(targets_dir, '{}_targets.txt'.format(drug_instance.drug_name))
    diana_drug.create_targets_file(drug_instance.targets, targets_file)


    #------------------------------#
    #   CREATE ASSOCIATION FILES   #
    #------------------------------#

    # Define parameters for the functional enrichment
    type_functions = ['gobp', 'gomf', 'reactome']
    type_corrections = ['fdr_bh', 'bonferroni']

    # Check if the gene-function association files are created
    functions_data_dir = os.path.join(data_dir, 'functions_data')
    for type_function in type_functions:
        # Association files for the network
        associations_file = os.path.join(network_associations_dir, '{}_to_gene.txt'.format(type_function))
        if not fileExist(associations_file):
            # Create associations file for GUILD (using the geneids of the network)
            functional_analysis.create_association_file(all_geneids=network_instance.get_nodes(), type_function=type_function, taxID=9606, output_file=associations_file, functions_data_dir=functions_data_dir)
        # Association files for all targets
        associations_file = os.path.join(target_associations_dir, '{}_to_gene.txt'.format(type_function))
        if not fileExist(associations_file):
            # Get all geneids associated to targets
            drugbank_geneid_mapping_file = os.path.join(mappings_dir, 'drugbank_geneid_drug_target_interactions.txt')
            targets = diana_drug.get_all_targets_from_mappings(drugbank_geneid_mapping_file)
            # Create associations file for targets (using all geneids)
            functional_analysis.create_association_file(all_geneids=targets, type_function=type_function, taxID=9606, output_file=associations_file, functions_data_dir=functions_data_dir)


    #--------------------------------#
    #  SCORING OF NETWORKS (GUILD)   #
    #--------------------------------#

    # Run GUILD
    print("  DIANA INFO:\tRunning GUILD (network scoring program).\n")

    # Create a directory for GUILD results
    guild_output_dir = os.path.join(drug_dir, 'guild_output')
    create_directory(guild_output_dir)

    # Create targets file for the targets in the network (it will be used by GUILD)
    network_targets_file = os.path.join(guild_output_dir, '{}_targets_in_network.txt'.format(drug_instance.drug_name))
    diana_drug.create_targets_file(drug_instance.targets_in_network, network_targets_file)

    # Run GUILD
    scores_file = os.path.join(guild_output_dir, 'output_scores.sif.netcombo')
    if not fileExist(scores_file):

        guild_command = 'python {} {} {} {} {} {} {}'.format( os.path.join(scripts_dir, 'run_guild.py'), drug_dir, network_targets_file, options.sif, guild_output_dir, random_networks_dir, config.get('Paths', 'guild_path') )
        os.system(guild_command)
        print('  DIANA INFO:\tGUILD has finished.\n')

        # Remove files not needed
        files_to_remove = ['edge_scores_netshort.sif', 'edge_scores.sif', 'node_scores_background.sif', 'node_scores.sif', 'output_scores.sif.netscore.log', 'output_scores.sif.netshort.log', 'output_scores.sif.netzcore.log', 'seed_scores_background.sif', 'seed_scores.sif']
        for file_to_remove in files_to_remove:
            if fileExist(os.path.join(guild_output_dir, file_to_remove)):
                command = 'rm {}'.format(os.path.join(guild_output_dir, file_to_remove))
                os.system(command)

    else:
        print('  DIANA INFO:\tThe scoring of the network with GUILD for {} was already done and it has been skipped.\n'.format(options.drug_name))

    # Creating an instance of the file generated by GUILD
    guild_profile_instance = network_analysis.GUILDProfile(scores_file, type_id=network_instance.type_id, top=100, top_type='percentage')


    #-----------------------------#
    #   GENERATE GUILD PROFILES   #
    #-----------------------------#

    print('  DIANA INFO:\tSTARTING GENERATION OF GUILD PROFILES\n')

    # Copy the scores file at the guild directory
    new_scores_file = os.path.join(guild_dir, 'output_scores.sif.netcombo')
    shutil.copyfile(scores_file, new_scores_file)

    # Score the network
    scored_network_file = os.path.join(guild_dir, 'network_scored.txt')
    if not fileExist(scored_network_file):
        scored_network_instance = network_instance.score_network(guild_profile_instance.node_to_score, scored_network_file)
    else:
        scored_network_instance = network_analysis.EdgeProfile(network_file=scored_network_file, type_id=network_instance.type_id, network_format=network_instance.network_format, top=100)

    # Get the list of thresholds to create the profiles
    if options.threshold_list and fileExist(options.threshold_list):
        threshold_list = get_values_from_threshold_file(options.threshold_list)
    else:
        threshold_list = [1, 2, 5, 'functions']
    print('  DIANA INFO:\tList of percentages used to define the drug profiles: {}\n'.format(', '.join([str(x) for x in threshold_list])))

    # Calculate functional threshold and derived profiles
    functional_thresholds = {}
    if 'functions' in threshold_list:
        # Calculate the seeds functional enrichment (using targets of the network, because the cutoff is for GUILD profiles)
        for type_function in type_functions:
            associations_file = os.path.join(network_associations_dir, '{}_to_gene.txt'.format(type_function))
            for type_correction in type_corrections:
                print(type_function, type_correction)
                output_seeds_enrichment_file = os.path.join(guild_dir, 'functional_profile_targets_network_{}_{}.txt'.format(type_function, type_correction))
                output_sliding_window_file = os.path.join(guild_dir, 'sliding_window_{}_{}.txt'.format(type_function, type_correction))
                if not fileExist(output_sliding_window_file):
                    cutoff_central_position, cutoff_right_interval = functional_analysis.calculate_functions_threshold(seed_geneids=drug_instance.targets_in_network, geneid_to_score=guild_profile_instance.node_to_score, type_correction=type_correction, associations_file=associations_file, output_seeds_enrichment_file=output_seeds_enrichment_file, output_sliding_window_file=output_sliding_window_file, seed_functional_enrichment=False)
                else:
                    cutoff_central_position, cutoff_right_interval = functional_analysis.read_sliding_window_file(output_sliding_window_file=output_sliding_window_file, num_seeds=len(drug_instance.targets_in_network))
                print('Cut-off central position: {}. Cut-off right interval position: {}'.format(cutoff_central_position, cutoff_right_interval))
                functional_thresholds.setdefault(type_function, {})
                functional_thresholds[type_function][type_correction] = cutoff_right_interval

                # Plot sliding window
                #output_plot_file = os.path.join(guild_dir, 'sliding_window_{}_{}.png'.format(type_function, type_correction))
                #functional_analysis.plot_sliding_window(output_sliding_window_file, output_plot_file, maximum_rank=cutoff_right_interval)

                # Generate the NODE profile from the top % scoring nodes
                output_file = os.path.join(guild_dir, 'node_profile_top_{}_{}_{}_{}.txt'.format('functions', type_function, type_correction, guild_profile_instance.type_id))
                if not fileExist(output_file):
                    node_profile_instance = guild_profile_instance.create_node_profile(threshold=cutoff_right_interval, threshold_type='number_of_nodes', output_file=output_file)
                else:
                    node_profile_instance = network_analysis.GUILDProfile(scores_file=output_file, type_id=guild_profile_instance.type_id, top=cutoff_right_interval, top_type='number_of_nodes')

                # Generate the FUNCTIONAL profile from the top % scoring nodes
                output_file = os.path.join(guild_dir, 'functional_profile_top_{}_{}_{}.txt'.format('functions', type_function, type_correction))
                if not fileExist(output_file):
                    functional_profile_instance = node_profile_instance.create_functional_profile(type_correction, output_file, associations_file)

                # Generate the EDGE profile from the top % scoring nodes
                output_file = os.path.join(guild_dir, 'edge_profile_top_{}_{}_{}_{}.txt'.format('functions', type_function, type_correction, guild_profile_instance.type_id))
                if not fileExist(output_file):
                    edge_profile_instance = scored_network_instance.create_edge_profile(node_to_score=guild_profile_instance.node_to_score, threshold=cutoff_right_interval, threshold_type='number_of_nodes', output_file=output_file)


    # Create GUILD profiles for each non-functional threshold defined 
    for top_threshold in threshold_list:

        if top_threshold != 'functions':

            # Generate the NODE profile from the top % scoring nodes
            output_file = os.path.join(guild_dir, 'node_profile_top_{}_{}.txt'.format(str(top_threshold), guild_profile_instance.type_id))
            if not fileExist(output_file):
                node_profile_instance = guild_profile_instance.create_node_profile(threshold=top_threshold, threshold_type='percentage', output_file=output_file)
            else:
                node_profile_instance = network_analysis.GUILDProfile(scores_file=output_file, type_id=guild_profile_instance.type_id, top=top_threshold, top_type='percentage')

            # Generate the FUNCTIONAL profile from the top % scoring nodes
            for type_function in type_functions:
                associations_file = os.path.join(network_associations_dir, '{}_to_gene.txt'.format(type_function))
                for type_correction in type_corrections:
                    output_file = os.path.join(guild_dir, 'functional_profile_top_{}_{}_{}.txt'.format(str(top_threshold), type_function, type_correction))
                    if not fileExist(output_file):
                        functional_profile_instance = node_profile_instance.create_functional_profile(type_correction, output_file, associations_file)

            # Generate the EDGE profile from the top % scoring nodes
            output_file = os.path.join(guild_dir, 'edge_profile_top_{}_{}.txt'.format(str(top_threshold), guild_profile_instance.type_id))
            if not fileExist(output_file):
                edge_profile_instance = scored_network_instance.create_edge_profile(node_to_score=guild_profile_instance.node_to_score, threshold=top_threshold, threshold_type='percentage', output_file=output_file)


    #------------------------------#
    #   GENERATE TARGET PROFILES   #
    #------------------------------#

    print('  DIANA INFO:\tSTARTING GENERATION OF TARGET PROFILES\n')

    # Create PFAM profile from targets
    pfam_file = os.path.join(targets_dir, 'pfam_profile.txt')

    if not fileExist(pfam_file):
        target_mapping_file = os.path.join(mappings_dir, 'geneid_target_mappings.txt')
        if fileExist(target_mapping_file):
            # Obtain the PFAMs from a table
            drug_instance.obtain_pfams_from_geneid_target_table(drug_instance.targets, target_mapping_file)
            # Create a PFAM file
            diana_drug.create_targets_file(drug_instance.pfams, pfam_file)
        else:
            print("  DIANA INFO:\tMissing target mappings file: {}.\n".format(target_mapping_file))
            sys.exit(10)
    else:
        drug_instance.obtain_pfams_from_file(pfam_file)

    print( "  DIANA INFO:\tThe PFAMs obtained from the targets are\n\t\t{}.\n".format( ', '.join([ str(x) for x in drug_instance.pfams]) ) )

    # Create the FUNCTIONAL profile from targets
    for type_function in type_functions:
        associations_file = os.path.join(target_associations_dir, '{}_to_gene.txt'.format(type_function))
        for type_correction in type_corrections:
            output_file = os.path.join(targets_dir, 'targets_functional_profile_{}_{}.txt'.format(type_function, type_correction))
            if not fileExist(output_file):
                top_scoring.functional_top_scoring(top_geneids=drug_instance.targets, type_correction=type_correction, output_file=output_file, associations_file=associations_file)


    #---------------------------#
    #   GENERATE ATC PROFILES   #
    #---------------------------#

    print('  DIANA INFO:\tSTARTING GENERATION OF ATCs PROFILES\n')

    #inside TARGETS CONTROLLER insert function that creates directory. By now: ATC_dir
    ATC_file = os.path.join(atc_dir,'ATC_profile.txt')

    if not fileExist(ATC_file):
        # Obtain ATCs from a table
        drugbank_atc_file = os.path.join(mappings_dir, 'drugbank_drug_atc.txt')
        if fileExist(drugbank_atc_file):
            # First, translate the drug input name to drugbankid (if necessary)
            drugbankids = drug_instance.obtain_drugbankids_from_table(drug_mapping_file)
            # Search ATCs in the table
            drug_instance.obtain_ATCs_from_table(drugbankids, drugbank_atc_file)
            # Create an ATC file
            diana_drug.create_targets_file(drug_instance.ATCs, ATC_file)
        else:
            print("  DIANA INFO:\tMissing DrugBank - ATCs file: {}.\n".format(drugbank_atc_file))
            sys.exit(10)
    else:
        drug_instance.obtain_ATCs_from_file(ATC_file)

    print( "  DIANA INFO:\tThe ATCs obtained from the targets are\n\t\t{}.\n".format( ', '.join([ str(x) for x in drug_instance.ATCs])))


    #---------------------------------#
    #   GENERATE STRUCTURE PROFILES   #
    #---------------------------------#

    print('  DIANA INFO:\tSTARTING GENERATION OF STRUCTURE PROFILES\n')

    # Obtain the SMILES of the compound from the database
    structure_file = os.path.join(structure_dir, 'structure_profile.txt')

    if not fileExist(structure_file):
        # Obtain SMILES from a table
        drugbank_smiles_file = os.path.join(mappings_dir, 'drugbank_drug_smiles.txt')
        if fileExist(drugbank_smiles_file):
            # First, translate the drug input name to drugbankid (if necessary)
            drugbankids = drug_instance.obtain_drugbankids_from_table(drug_mapping_file)
            # Search SMILES in the table
            drug_instance.obtain_SMILES_from_table(drugbankids, drugbank_smiles_file)
            # Create the SMILES file
            diana_drug.create_targets_file(drug_instance.smiles, structure_file)
        else:
            print("  DIANA INFO:\tMissing DrugBank - SMILES file: {}.\n".format(drugbank_smiles_file))
            sys.exit(10)
    else:
        drug_instance.obtain_SMILES_from_file(structure_file)

    print( "  DIANA INFO:\tThe SMILES obtained are\n\t\t{}.\n".format( ', '.join([ str(x) for x in drug_instance.smiles]) ) )


    #------------------------------------#
    #   GENERATE SIDE EFFECTS PROFILES   #
    #------------------------------------#

    print('  DIANA INFO:\tSTARTING GENERATION OF SIDE EFFECTS PROFILES\n')

    # Obtain the side effects of the compound from the database
    SE_file = os.path.join(se_dir, 'SE_profile.txt')

    if not fileExist(SE_file):
        # Obtain side effects from a table
        drugbank_side_effects_file = os.path.join(mappings_dir, 'drugbank_drug_side_effects.txt')
        if fileExist(drugbank_side_effects_file):
            # First, translate the drug input name to drugbankid (if necessary)
            drugbankids = drug_instance.obtain_drugbankids_from_table(drug_mapping_file)
            # Search side effects in the table
            drug_instance.obtain_SE_from_table(drugbankids, drugbank_side_effects_file)
            # Create the side effects file
            diana_drug.create_targets_file(drug_instance.SEs, SE_file)
        else:
            print("  DIANA INFO:\tMissing DrugBank - side effects file: {}.\n".format(drugbank_side_effects_file))
            sys.exit(10)
    else:
        drug_instance.obtain_SE_from_file(SE_file)

    print( "  DIANA INFO:\tThe side effects obtained are\n\t\t{}.\n".format( ', '.join([ str(x) for x in drug_instance.SEs]) ) )


    #--------------------------------#
    #   GENERATE A PARAMETERS FILE   #
    #--------------------------------#

    # Generate parameters file
    parameters_file = os.path.join(drug_dir, 'parameters.txt')
    if not fileExist(parameters_file):
        create_parameters_file(drug_id, drug_instance.drug_name, drug_instance.targets, network_filename, parameters_file)


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

def create_parameters_file(diana_id, drug_name, targets, network_name, output_file):
    """
    Creates a file containing the parameters used in the analysis.
    """
    targets = [str(x) for x in targets] # Transform the targets to strings
    with open(output_file, 'w') as output_fd:
        output_fd.write('#diana_id\tdrug\ttargets\tnetwork\n')
        output_fd.write('{}\t{}\t{}\t{}\n'.format( diana_id, drug_name, ';'.join(targets), network_name ))
    return

def get_targets_in_network_of_expansion(node_file):
    """
    In the network of expansion, the targets are scored 1.00.
    This function gets the targets in the network by identifying the nodes with score 1.
    """
    targets_in_network = set()
    with open(node_file, 'r') as node_fd:
        for line in node_fd:
            node, score = line.strip().split('\t')
            score = float(score)
            if score == 1:
                targets_in_network.add(node)
    return list(targets_in_network)

# def get_targets_in_sif_file(sif_file, targets):
#     """
#     Get the targets that are inside the network given by the user
#     """
#     targets_in_network = set()
#     str_tar = [str(x) for x in targets]
#     with open(sif_file, 'r') as sif_fd:
#         for line in sif_fd:
#             node1, score, node2 = line.strip().split('\t')
#             if node1 in str_tar:
#                 targets_in_network.add(node1)
#             if node2 in str_tar:
#                 targets_in_network.add(node2)
#     return list(targets_in_network)

def get_values_from_threshold_file(threshold_file):
    """
    Get the targets that are inside the network given by the user
    """
    threshold_list = []
    with open(threshold_file, 'r') as threshold_fd:
        for line in threshold_fd:
            threshold = line.strip()
            try:
                threshold_list.append(int(threshold))
            except:
                try:
                    threshold_list.append(float(threshold))
                except:
                    threshold_list.append(str(threshold))
    return threshold_list

def translate_targets_to_type_id(targets_translation_file, targets):
    """
    Obtain the translated targets using the targets_translation_file
    """
    # Parse the targets_translation_file
    # Fill new with the translations
    new={}
    with open(targets_translation_file, 'r') as targets_trans_fd:
        for line in targets_trans_fd:
            word=line.strip().split('\t')
            node=word[0]
            for transet in word[1].split("'"):
                for trans in transet.split(","):
                    if trans != "" and trans != ",":
                        new.setdefault(trans,set())
                        new[trans].add(node)

    # Get the translation from the targets
    translated_targets = set()
    for target in targets:
        target = str(target)
        if target in new:
            for translation in new[target]:
                translated_targets.add(translation)

    return list(translated_targets)

if  __name__ == "__main__":
    main()
