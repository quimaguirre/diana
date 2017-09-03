import argparse
import mysql.connector
import ntpath
import time
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug
import diana.classes.network_analysis as network_analysis
import diana.classes.network_generation as network_generation
import diana.classes.top_scoring as top_scoring
from diana.classes.goatools.obo_parser import GODag
from diana.classes.goatools.associations import read_ncbi_gene2go




def main():

    options = parse_user_arguments()
    generate_profiles(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of the input drug",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-d','--drug_name',dest='drug_name',action = 'store',
                        help = """ Name of the drug number 1. If you do not provide targets for this drug or the number of targets is not large enough,
                        the program will use this name to search for targets in BIANA database. If targets are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between 
                        double quotes. """)
    parser.add_argument('-t','--targets',dest='targets',action = 'store',
                        help = 'Input file with the targets of the drug. Each target must be separated by a newline character.')
    parser.add_argument('-pt','--proteins_type_id',dest='proteins_type_id',action = 'store', default='geneid',
                        help = 'Input the type of ID of the targets introduced / proteins of the network. It must be the same! (default is geneid).')
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = """" Input file with a protein-protein interaction network in SIF format.
                        If not introduced, the program will create a network of expansion using the targets as center and expanding as many neighbors
                        as specified in the parameter radius. """)
    parser.add_argument('-rad','--radius',dest='radius',action = 'store',default='3',
                        help = """ Define the radius of expansion for the creation of the network from targets (default is 3).
                        Only used if no network file provided. """)
    parser.add_argument('-tax','--taxid',dest='taxid',action = 'store',default='9606',
                        help = """Define the restriction of species for the creation of the network from targets using a Taxonomy ID (default is '9606' (human))""")
    parser.add_argument('-res','--restriction',dest='restriction',action = 'store',
                        help = """Define an experiment restriction for the creation of the network from targets.\n
                        Options:\n
                        - AFF: Use interactions at least described by affinity methods (i.e. Tandem Affinity Purification)\n
                        - Y2H: Use interactions at least described by yeast two hybrid methods (Y2H)\n
                        - eAFF: Use all interactions except those described by affinity methods (i.e. Tandem Affinity Purification)\n
                        - eY2H: Use all interactions except those described by yeast two hybrid methods (Y2H)\n
                        - None: Not use experiment restrictions
                        """)
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
                        help = """Define the database to use for the generation of the network of expansion / search of targets: 
                        (default is BIANA_JUN_2017)""")
    parser.add_argument('-up','--unification',dest='unification_protocol',action = 'store',default='geneid_seqtax_v1',
                        help = """Define the unification protocol used in BIANA database (default is BIANA_JUN_2017)""")

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

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')

    # Create a directory for the data
    data_dir = os.path.join(options.workspace, "profiles")
    create_directory(data_dir)

    # Create a directory for the random networks
    other_data_dir = os.path.join(options.workspace, "other_data")
    create_directory(other_data_dir)

    random_networks_dir = os.path.join(other_data_dir, "random_networks")
    create_directory(random_networks_dir)

    # Create a drug instance
    drug_instance = diana_drug.Drug(options.drug_name)



    #------------------------#
    #   TARGETS CONTROLLER   #
    #------------------------#

    # TARGETS CONTROLLER: Checks the targets provided by the user. If necessary, performs a search 
    # in BIANA database to obtain more targets

    # Check if the targets file is provided
    if options.targets and os.path.isfile(options.targets):
        drug_instance.obtain_targets_from_file(options.targets, options.proteins_type_id)
    else:
        # Create a connection to BIANA database
        biana_cnx = mysql.connector.connect(user='quim', password="",
                                            host='localhost',
                                            database=options.database)
        # Obtain the targets from BIANA
        drug_instance.obtain_targets_from_BIANA(biana_cnx, options.proteins_type_id, options.unification_protocol)
        biana_cnx.close()

    print( "  DIANA INFO:\tThe targets provided for the drug {} are:\n\t\t{}.\n".format( options.drug_name, ', '.join([ str(x) for x in drug_instance.targets]) ) )

    # Create a directory for the drug
    drug_dir = os.path.join(data_dir, diana_drug.generate_drug_id(drug_instance.drug_name, drug_instance.targets))
    create_directory(drug_dir)

    # Create a directory for the dcTargets results
    dctargets_dir = os.path.join(drug_dir, 'dctargets_profiles')
    create_directory(dctargets_dir)

    # Create a directory for the dcGUILD results
    dcguild_dir = os.path.join(drug_dir, 'dcguild_profiles')
    create_directory(dcguild_dir)

    # Create a targets file
    targets_file = os.path.join(dctargets_dir, '{}_targets.txt'.format(drug_instance.drug_name))
    diana_drug.create_targets_file(drug_instance.targets, targets_file)



    #--------------------#
    #   SIF CONTROLLER   #
    #--------------------#

    # SIF CONTROLLER: Checks the network in SIF format provided by the user.
    # If no file provided, generates a network expanding it from the targets provided.
    # The network is expanded as many neighbors as indicated in the parameter 'radius'.

    # Check if the network file is provided
    if options.sif and os.path.isfile(options.sif):
        # If the network file is provided, we create a Network instance
        network_instance = network_analysis.Network(options.sif, None, options.proteins_type_id, 'sif')
        targets_in_network = get_targets_in_sif_file(options.sif, drug_instance.targets)
        translation_file = None # Just defining the variable to know that there is no translation file when sif is provided

        # We create a directory in the random networks directory for this network
        network_filename = ntpath.basename(options.sif).split('.')[0]
        random_networks_dir = os.path.join(random_networks_dir, network_filename)
        create_directory(random_networks_dir)

    else:
        # If not, we create a network from the targets
        print('  DIANA INFO:\tGenerating network for {}. This can take a few minutes...\n'.format(options.drug_name))

        network_dir = os.path.join(drug_dir, 'network_of_expansion')
        create_directory(network_dir)
        node_file = os.path.join(network_dir, '{}_network.nodes'.format(drug_instance.drug_name))
        edge_file = os.path.join(network_dir, '{}_network.edges'.format(drug_instance.drug_name))
        translation_file = os.path.join(network_dir, '{}_network_trans_to_{}.trans'.format(drug_instance.drug_name, options.proteins_type_id))
        targets_translation_file = os.path.join(network_dir, 'targets_to_BIANA.trans')

        if not fileExist(edge_file):

            # By direct command (it works)
            restriction = process_restriction(options.restriction)
            command = 'python {} -iseed {} -radius {} -taxid {} -stype {} -ttype {} -trans {} -strans {} -node {} -edge {} -db {} -up {} {}'.format(
                      os.path.join(toolbox_dir, 'generate_netscore_files_vapr2017.py'),
                      targets_file, options.radius, options.taxid, options.proteins_type_id, options.proteins_type_id, translation_file, targets_translation_file, node_file, edge_file,
                      options.database, options.unification_protocol, restriction
                      )
            os.system(command)

            # By importing the module (it does not work!!!)
            # restricted_to_TAP, restricted_to_Y2H, restricted_to_user, except_TAP, except_Y2H, except_user = network_generation.check_restriction(options.restriction)
            # network_generation.generate_network(drug_instance.targets, drug_instance.type_id, options.radius, options.taxid, translation_file, options.proteins_type_id, node_file, edge_file,
            #              restricted_to_TAP = restricted_to_TAP, restricted_to_Y2H = restricted_to_Y2H, restricted_to_user = restricted_to_user,
            #              except_TAP = except_TAP, except_Y2H = except_Y2H, except_user = except_user,
            #              database = options.database, unification_protocol = options.unification_protocol,
            #              output_format = 'sif', verbose = False)

        else:
            print('  DIANA INFO:\tThe network of expansion for {} was already done and it has been skipped.\n'.format(options.drug_name))

        network_instance = network_analysis.Network(edge_file, None, 'biana', 'sif')
        options.sif = edge_file
        targets_in_network = get_targets_in_network_of_expansion(node_file)

        # Create a directory of random networks corresponding to this network
        random_networks_dir = os.path.join(random_networks_dir, '{}_network'.format(drug_instance.drug_name))
        create_directory(random_networks_dir)


    # Check if the number of targets provided is sufficient for the analysis
    if len(targets_in_network) < 3:
        raise diana_drug.InsufficientTargets(targets_in_network)



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
    diana_drug.create_targets_file(targets_in_network, network_targets_file)

    # Run GUILD
    pvalue_file = os.path.join(guild_output_dir, 'output_scores.sif.netcombo.pval')
    if not fileExist(pvalue_file):

        guild_command = 'python {} {} {} {} {} {}'.format( os.path.join(toolbox_dir, 'run_guild.py'), drug_dir, network_targets_file, options.sif, guild_output_dir, random_networks_dir )
        os.system(guild_command)
        print('  DIANA INFO:\tGUILD has finished.\n')

    else:
        print('  DIANA INFO:\tThe scoring of the network with GUILD for {} was already done and it has been skipped.\n'.format(options.drug_name))

    # Creating an instance of the file generated by GUILD
    guild_profile_instance = network_analysis.GUILDProfile(pvalue_file, network_instance.type_id, 100)

    # Translate the NODE profile to protein_type_id if the type of id is 'biana'
    if guild_profile_instance.type_id == 'biana' and translation_file:
        output_file = os.path.join(guild_output_dir, 'node_profile_top_100_{}.txt'.format(options.proteins_type_id))
        guild_profile_geneid = guild_profile_instance.translate_pvalue_file(translation_file, options.proteins_type_id, output_file, verbose=False)



    #-------------------------------#
    #   GENERATE DCGUILD PROFILES   #
    #-------------------------------#

    print('  DIANA INFO:\tSTARTING GENERATION OF dcGUILD PROFILES\n')

    # Score the network
    scored_network_file = os.path.join(dcguild_dir, 'network_scored.txt')
    scored_network_instance = network_instance.score_network(guild_profile_instance.node_to_values, scored_network_file)

    # Get the list of thresholds to create the profiles
    if options.threshold_list and fileExist(options.threshold_list):
        threshold_list = get_values_from_threshold_file(options.threshold_list)
    else:
        threshold_list = [1, 5, 10, 20, 50]
    print('  DIANA INFO:\tList of percentages used to define the drug profiles: {}\n'.format(', '.join([str(x) for x in threshold_list])))

    # Load the files go and gene2go files for the functional profiles
    print("  DIANA INFO:\tLoading GOATOOLS for the functional analysis. Remember to update the files go-basic.obo and gene2go.")
    print("\t\tpython /path/to/diana/scripts/download_go_files.py\n")
    obodag = GODag(os.path.join(main_path, "diana/toolbox/go-basic.obo"))
    geneid2gos_human = read_ncbi_gene2go(os.path.join(main_path, "diana/toolbox/gene2go"), taxids=[9606])

    for top_threshold in threshold_list:

        output_file = os.path.join(dcguild_dir, 'node_profile_top_{}_{}.txt'.format(str(top_threshold), guild_profile_instance.type_id))
        if not fileExist(output_file):

            # Generate the NODE profile from the top % scoring nodes
            node_profile_instance = guild_profile_instance.create_node_profile(top_threshold, output_file)

            # Translate the NODE profile to protein_type_id if the type of id is 'biana'
            if node_profile_instance.type_id == 'biana' and translation_file:
                output_file = os.path.join(dcguild_dir, 'node_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id))
                node_profile_geneid = node_profile_instance.translate_pvalue_file(translation_file, options.proteins_type_id, output_file, verbose=False)


            output_file = os.path.join(dcguild_dir, 'functional_profile_top_{}_{}.txt'.format(str(top_threshold), guild_profile_instance.type_id))
            if not fileExist(output_file):

                # Generate the FUNCTIONAL profile from the top % scoring nodes of the pvalue file
                if node_profile_instance.type_id == 'biana' and translation_file:
                    functional_profile_instance = node_profile_geneid.create_functional_profile(obodag, geneid2gos_human, guild_profile_geneid.node_to_values, output_file)
                else:
                    functional_profile_instance = node_profile_instance.create_functional_profile(obodag, geneid2gos_human, guild_profile_instance.node_to_values, output_file)


        output_file = os.path.join(dcguild_dir, 'edge_profile_top_{}_{}.txt'.format(str(top_threshold), guild_profile_instance.type_id))
        if not fileExist(output_file):

            # Generate the EDGE profile from the top % scoring nodes of the pvalue file
            edge_profile_instance = scored_network_instance.create_edge_profile(guild_profile_instance.node_to_values, top_threshold, output_file)

            # Translate the EDGE profile to protein_type_id if the type of id is 'biana'
            if edge_profile_instance.type_id == 'biana' and translation_file:
                output_file = os.path.join(dcguild_dir, 'edge_profile_top_{}_{}.txt'.format(str(top_threshold), options.proteins_type_id))
                edge_profile_geneid = edge_profile_instance.translate_network(translation_file, options.proteins_type_id, output_file, verbose=False)



    #---------------------------------#
    #   GENERATE DCTARGETS PROFILES   #
    #---------------------------------#

    print('  DIANA INFO:\tSTARTING GENERATION OF dcTARGETS PROFILES\n')

    # Create PFAM profile from targets
    pfam_file = os.path.join(dctargets_dir, 'pfam_profile.txt')

    if not fileExist(pfam_file):
        # Create a connection to BIANA database
        biana_cnx = mysql.connector.connect(user='quim', password="",
                                            host='localhost',
                                            database=options.database)
        # Obtain the PFAMs from BIANA
        drug_instance.obtain_pfams_from_targets(biana_cnx, pfam_file, options.unification_protocol)
        biana_cnx.close()
    else:
        drug_instance.obtain_pfams_from_file(pfam_file)

    print( "  DIANA INFO:\tThe PFAMs obtained from the targets are\n\t\t{}.\n".format( ', '.join([ str(x) for x in drug_instance.pfams]) ) )


    # Create the FUNCTIONAL profile from targets
    output_file = os.path.join(dctargets_dir, 'targets_functional_profile.txt')
    if not fileExist(output_file):

        if guild_profile_instance.type_id == 'biana' and translation_file:
            #targets_in_network_geneid = translate_targets_to_type_id(targets_translation_file, targets_in_network)
            #targets_functional_profile_instance = top_scoring.functional_top_scoring(obodag, geneid2gos_human, guild_profile_geneid.node_to_values.keys(), targets_in_network_geneid, output_file)
            # First we add the targets that are not in the network among all the nodes in the network, to use them as background
            # Because if one of them is not among the background genes, the program raises an error!
            all_nodes_geneid = set(guild_profile_geneid.node_to_values.keys())
            for target in drug_instance.targets:
                all_nodes_geneid.add(target)
            targets_functional_profile_instance = top_scoring.functional_top_scoring(obodag, geneid2gos_human, list(all_nodes_geneid), drug_instance.targets, output_file)
        else:
            # Here we also add the targets that are not in the network among all the nodes in the network, to use them as background
            all_nodes_geneid = set(guild_profile_instance.node_to_values.keys())
            for target in drug_instance.targets:
                all_nodes_geneid.add(target)
            targets_functional_profile_instance = top_scoring.functional_top_scoring(obodag, geneid2gos_human, list(all_nodes_geneid), drug_instance.targets, output_file)


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

def process_restriction(restriction):
    """
    Checks if the restriction has been correctly introduced.
    Returns the value that will be introduced in the command to generate the network of expansion.
    """
    if not restriction:
        return ''
    else:
        res = restriction.lower()
        if res == 'eaff':
            return '-eAFF'
        elif res == 'ey2h':
            return '-eY2H'
        elif res == 'y2h':
            return '-rY2H'
        elif res == 'aff':
            return '-rAFF'
        else:
            raise network_generation.IncorrectRestrictionType(res)

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

def get_targets_in_sif_file(sif_file, targets):
    """
    Get the targets that are inside the network given by the user
    """
    targets_in_network = set()
    with open(sif_file, 'r') as sif_fd:
        for line in sif_fd:
            node1, score, node2 = line.strip().split('\t')
            if node1 in targets:
                targets_in_network.add(node1)
            if node2 in targets:
                targets_in_network.add(node2)
    return list(targets_in_network)

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
                threshold_list.append(float(threshold))
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

