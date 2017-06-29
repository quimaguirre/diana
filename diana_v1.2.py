import sys
import os
import argparse
import re
import toolbox
import toolbox.comparison as comp
import toolbox.drugbank_searcher as dbs
import toolbox.guild_utilities as GU
from toolbox import functional_enrichment
import time

def main():

    options = parse_user_arguments()
    run_diana(options)

def parse_user_arguments(*args, **kwds):
    """Parses the arguments of the program"""

    parser = argparse.ArgumentParser(
        description = "Create and compare the profiles of two drugs",
        epilog      = "@oliva's lab 2016")
    parser.add_argument('-d1','--drug1_name',dest='drug1',action = 'store',
                        help = """Name of the drug number 1. If you do not provide seeds for this drug or the number of seeds is not large enough,
                        the program will look for seeds of this drug at the DrugBank database. If seeds are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has two words, introduce the two words joined by an underscore (i.e. acetylsalicylic_acid)""")
    parser.add_argument('-s1','--seeds_drug1',dest='seeds1',action = 'store',default='drug1.seeds',
                        help = 'Input file with seeds of drug 1. Use UniprotEntries separated by newline (default is drug1.seeds)')
    parser.add_argument('-d2','--drug2_name',dest='drug2',action = 'store',
                        help = """Name of the drug number 2. If you do not provide seeds for this drug or the number of seeds is not large enough,
                        the program will look for seeds of this drug at the DrugBank database. If seeds are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has two words, introduce the two words joined by an underscore (i.e. acetylsalicylic_acid)""")
    parser.add_argument('-s2','--seeds_drug2',dest='seeds2',action = 'store',default='drug2.seeds',
                        help = 'Input file with seeds of drug 2. Use UniprotEntries separated by newline (default is drug2.seeds)')
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = 'Input file with a protein-protein interaction network in SIF format (if not introduced, networks of expansion are used)')
    parser.add_argument('-t','--top',dest='top_threshold',action = 'store',default='10',
                        help = 'Define the percentage of top scoring nodes to create the node profile (default is 10)')
    parser.add_argument('-p','--pvalue',dest='pval_threshold',action = 'store',default='1',
                        help = 'Define the pvalue threshold to filter the top nodes and create the node profile (default is 0.05)')
    parser.add_argument('-te','--top_ed',dest='top_threshold_edges',action = 'store',default='10',
                        help = 'Define the percentage of top scoring nodes to create the edge profile (default is 10)')
    parser.add_argument('-pe','--pvalue_ed',dest='pval_threshold_edges',action = 'store',default='1',
                        help = 'Define the pvalue threshold to filter the top nodes and create the edge profile (default is 1)')
    parser.add_argument('-tf','--top_fun',dest='top_threshold_functional',action = 'store',default='10',
                        help = 'Define the percentage of top scoring nodes to create the functional profile (default is 10)')
    parser.add_argument('-pf','--pvalue_fun',dest='pval_threshold_functional',action = 'store',default='1',
                        help = 'Define the pvalue threshold to filter the top nodes and create the functional profile (default is 1)')
    parser.add_argument('-rad','--radius',dest='radius',action = 'store',default='3',
                        help = 'Define the radius of expansion for the creation of the network from seeds (default is 3)')
    parser.add_argument('-tax','--taxid',dest='taxid',action = 'store',default='9606',
                        help = 'Define the restriction of taxid for the creation of the network from seeds (default is 9606 (human))')
    parser.add_argument('-res','--restriction',dest='restriction',action = 'store',default='Y2H',
                        help = """Define an experiment restriction for the creation of the network from seeds (default is Y2H). 
                        Options:
                        - AFF: Use interactions at least described by affinity methods (i.e. Tandem Affinity Purification)
                        - Y2H: Use interactions at least described by yeast two hybrid methods (Y2H)
                        - eAFF: Use all interactions except those described by affinity methods (i.e. Tandem Affinity Purification)
                        - eY2H: Use all interactions except those described by yeast two hybrid methods (Y2H)
                        - None: Not use experiment restrictions
                        """)
    parser.add_argument('-sp','--speed',dest='speed',action = 'store', type=int,default=1,
                        help = """Speed of the analysis: 
                        - 1: Normal speed, it performs the complete analysis
                        - 2: Quick speed, it skips the comparison of edges
                        - Default: 1
                        """)
    parser.add_argument('-ex','--extended',dest='extended',action = 'store', default=None,
                        help = """Extended analysis. It performs a comparison of the Spearman coefficients using several Top/P-val cut-offs. The input consists in:
                        - A file with different top threshold values separated by spaces. For example, a file called "top_thresholds.list" containing:
                        0.5 1 2 3 4 5 10 15 20 25 30 40 50
                        It can only be used if the speed of the analysis is 1!!
                        Default: None
                        """)


    options=parser.parse_args()

    return options

#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def check_two_words_drug(drug):
    """
    If a drug is composed of two words joined by underscore, the function detects it and returns the words separated by space.
    If it is not composed of two words, returns the initial drug name
    """
    processed_name = drug.split("_")
    if len(processed_name) > 1:
        return " ".join(processed_name)
    else:
        return drug

def check_type_of_seed(seeds_list):
    """
    Checks which kind of seed has been provided. Returns 'uniprotentry', 'geneid', 'other' (if it 
    is neither uniprotentry or geneid), 'different' (if it has different types of seeds)
    """

    # UniprotEntry regular expression pattern
    uniprot_pattern = re.compile('[a-zA-Z0-9]{,10}_[a-zA-Z0-9]{,5}')
    #uniprot_pattern = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')

    # GeneId pattern
    geneid_pattern = re.compile('[0-9]+')

    # In these boolean lists we store the results of the type checking
    uniprot_boolean_list = []
    geneid_boolean_list = []
    other_boolean_list = []

    for seed in seeds_list:
        if uniprot_pattern.match(seed) != None:
            uniprot_boolean_list.append(True)
            geneid_boolean_list.append(False)
            other_boolean_list.append(False)
        elif geneid_pattern.match(seed) != None:
            uniprot_boolean_list.append(False)
            geneid_boolean_list.append(True)
            other_boolean_list.append(False)
        else:
            uniprot_boolean_list.append(False)
            geneid_boolean_list.append(False)
            other_boolean_list.append(True)

    if all(item == True for item in uniprot_boolean_list) == True:
        return 'uniprotentry'
    elif all(item == True for item in geneid_boolean_list) == True:
        return 'geneid'
    elif all(item == True for item in other_boolean_list) == True:
        return 'other'
    else:
        return 'different'

def fileExist(file):
    """Checks if a file exists AND is a file"""
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else: 
        return False

def find_seeds_in_nodes_file(drug, nodes_file, new_seeds_file):
    """
       Finds seeds (scored 1.00000) in a nodes file and creates a new seeds file. Returns the seeds
       in a list
    """

    fn = open(nodes_file,"r")
    fs = open(new_seeds_file,"w")
    new_seeds = []

    for line in fn:
        words = line.split()
        if (words[1] == "1.00000"):
            fs.write("%s %s\n" %(words[0], words[1]))
            new_seeds.append(words[0])

    fn.close()
    fs.close()

    # Produce a copy of the seeds file in the results directory
    command = "cp "+new_seeds_file+" "+drug+"/guild_results/seeds.txt"
    os.system(command)

    return new_seeds

def find_seeds_in_sif_file(drug, sif_file, initial_seeds, new_seeds_file):
    """Finds seeds in a sif file and creates a new seeds file. Returns the seeds in a list"""

    fn = open(sif_file,"r")
    fs = open(new_seeds_file,"w")
    new_seeds = []

    for line in fn:
        words = line.split()
        if (words[0] in initial_seeds or words[2] in initial_seeds):
            if words[0] in initial_seeds:
                word = words[0]
            else:
                word = words[2]
            if word not in new_seeds:
                fs.write("%s 1.00000\n" %(word))
                new_seeds.append(word)

    fn.close()
    fs.close()

    # Produce a copy of the seeds file in the results directory
    command = "cp "+new_seeds_file+" "+drug+"/guild_results_using_sif/seeds.txt"
    os.system(command)

    return new_seeds

def generate_seeds_dict_for_comparison(seeds_list):
    """
       Generates a dict which contains the seed name as key, and the score '1.00' and a p-value 'foo'
       as a set in the value of the dict. This dict will be used for the seeds comparison
    """

    seeds_dict = {}

    for seed in seeds_list:
    	seeds_dict[seed] = (1.0, "foo")

    return seeds_dict


def process_seeds(seeds_file):
    """Takes the seeds from the input seeds file and stores them into a list"""

    fd = open(seeds_file,"r")
    seeds = []

    for line in fd:
        seeds.append(line.strip())

    fd.close()

    return seeds

def run_generate_netscore_files(drug, seeds_file, radius, taxid, restriction):
    """Runs the program generate_netscore_files.py that generates a network from a seeds file"""

    nodes = drug+"/human_nodes_BIANA_"+drug+".txt"
    loc = drug+"/noLoc_bPPI_"+drug+".txt"
    edges = drug+"/human_edges_BIANA_"+drug+".txt"
    trans = drug+"/translation_BIANA_"+drug+".txt"
    radius = int(radius)
    taxid = int(taxid)
    if (restriction == 'eAFF'):
        restriction = '-eAFF'
    elif (restriction == 'eY2H'):
        restriction = '-eY2H'
    elif (restriction == 'Y2H'):
        restriction = '-r'+restriction
    elif (restriction == 'AFF'):
        restriction = '-r'+restriction
    elif (restriction == 'None'):
        restriction = ''

    command = "/soft/devel/python-2.7/bin/python2.7 toolbox/generate_netscore_files_geneID.py -iseed %s -radius %d -stype geneid -ttype geneid -score 0.0 -taxid %d -node %s -loc %s -edge %s -trans %s -format guild %s " % (seeds_file, radius, taxid, nodes, loc, edges, trans, restriction)
    os.system(command)

    return

def run_top_network(data_dir, pvalue_file, network_file, seed_file, top_threshold, pval_threshold):
    """
       Runs the program top_network.py which creates:
       - A sif file which contains the edges of the filtered nodes scored
    """

    top_threshold = float(top_threshold)
    pval_threshold = float(pval_threshold)

    command = "/soft/devel/python-2.7/bin/python toolbox/top_network.py %s %s %s %s %s %f %f " % (data_dir, pvalue_file, network_file, seed_file, data_dir, top_threshold, pval_threshold)
    #command = "/soft/devel/python-2.7/bin/python toolbox/top_network.py %s %s %s %f %f >& %s/top_network.log" % (data_dir, pvalue_file, data_dir, top_threshold, pval_threshold, data_dir)
    os.system(command)

    return

def run_top_scoring(data_dir, pvalue_file, top_threshold, pval_threshold, entry_type):
    """
       Runs the program top_scoring_combination.py that has the following characteristics:
       - Takes the entries which are at a given percentage of the top scoring ones
       - Then, filters the top scoring entries by a given threshold of p-value
       - If entry_type = "uniprotentry": Creates a uniprotentry profile containing Id/Score/P-value
       - If entry_type = "refseq": Creates a functional enrichment profile
    """

    top_threshold = float(top_threshold)
    pval_threshold = float(pval_threshold)

    command = "/soft/devel/python-2.7/bin/python  toolbox/top_scoring_combination.py %s %s %s %f %f %s " % (data_dir, pvalue_file, data_dir, top_threshold, pval_threshold, entry_type)
    #command = "/soft/devel/python-2.7/bin/python  toolbox/top_scoring_combination.py %s %s %s %f %f %s >& %s/top_scoring_%s.log" % (data_dir, pvalue_file, data_dir, top_threshold, pval_threshold, entry_type, data_dir, entry_type)
    os.system(command)

    return

def run_transform_network(network_input_file, nodes_input_file, translation_file, output_file):
    """Runs the program TransformNetwork.py that translates a file in BIANA codes to a desired notation, using the previously obtained translation file from generate_netscore_files"""

    output_edges_file = output_file + ".edges"
    output_nodes_file = output_file + ".nodes"

    command = "/soft/devel/python-2.7/bin/python toolbox/TransformNetwork.py -i %s -n %s -trans %s -oe %s -on %s" % (network_input_file, nodes_input_file, translation_file, output_edges_file, output_nodes_file)
    os.system(command)

    return

def run_translate_network(network_input_file, nodes_input_file, taxid, translation_type, output_file, translation_of_nodes_file = None):
    """Runs the program TranslateNetwork.py that translates a file in BIANA codes to a desired notation"""

    taxid = int(taxid)

    if (translation_of_nodes_file == None):
        command = "/soft/devel/python-2.7/bin/python toolbox/TranslateNetwork.py -i %s -n %s -x -taxid %d -ttype %s -o %s" % (network_input_file, nodes_input_file, taxid, translation_type, output_file)
    else:
        command = "/soft/devel/python-2.7/bin/python toolbox/TranslateNetwork.py -i %s -n %s -trans %s -x -taxid %d -ttype %s -o %s" % (network_input_file, nodes_input_file, translation_of_nodes_file, taxid, translation_type, output_file)
    os.system(command)

    return

def translate_seeds(seeds_list, taxid):
    """
    Uses TranslateNetwork.py to translate a list of seeds from UniprotEntry to GeneID.
    Returns the translated seeds and the final seed file created
    """

    # Checking if a directory for seeds exists. If not, we will create it so that we have the seeds ordered in the same directory
    try:
        os.stat("seeds")
    except:
        print("  DIANA INFO:\tCreating directory for seeds.\n")
        os.mkdir("seeds")

    # Write a file similar to a p-value file from GUILD with the UniprotEntry seeds
    seeds_file = "seeds/"+drug+".seeds.foo"
    fs = open(seeds_file, 'w')
    for seed in seeds_list:
        fs.write("%s 0 0 0\n" %(seed))
    fs.close()

    # Translation of the uniprotentry file created to 
    input_format = "uniprotentry"
    translation_type = "geneid"
    output_file = "seeds/"+drug+".seeds.geneid"
    taxid = int(taxid)

    command = "/soft/devel/python-2.7/bin/python toolbox/TranslateNetwork.py -i %s -n %s -iformat %s -x -taxid %d -ttype %s -o %s" % (seeds_file, seeds_file, input_format, taxid, translation_type, output_file)
    os.system(command)
    # Remove the edge file created, which is not useful
    command = "rm %s.edges" %(output_file)
    os.system(command)

    # Get the uniprotaccession seeds from the node file created
    new_seeds = set()
    created_file = output_file+".nodes"
    f = open(created_file, 'r')
    for line in f:
        words = line.split()
        geneid = words[0]
        new_seeds.add(geneid)
    f.close()

    # Write a final seeds file with the seeds correctly written
    final_seeds_file = "seeds/"+drug+".seeds.final"
    fs = open(final_seeds_file, 'w')
    for item in new_seeds:
        fs.write("%s\n" %(item))
    fs.close()

    # Remove the node file created, which is not useful anymore
    command = "rm %s" %(created_file)
    os.system(command)

    return (new_seeds, final_seeds_file)



#################
#################
# MAIN FUNCTION #
#################
#################

def run_diana(options):
    """Runs the program DIANA"""

    print("\n\t\t------------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlyzer (DIANA), a program created by @OLIVA'S LAB. Thank you for your trust in our software\n")
    print("\t\t------------------------------------------------------------------------------------------------------------------------\n")

    # Start marker for time measure
    start = time.time()

    # Inform about the SPEED of the analysis
    print("  DIANA INFO:\tThe ANALYSIS OF SPEED %s option has been activated\n" %(options.speed))

    # Inform if EXTENDED analysis selected
    if options.extended != None:
        print("  DIANA INFO:\tEXTENDED ANALYSIS ACTIVATED\n")    

    # If no names for drug1 or drug2 are provided, the analysis will be performed without looking for drug-drug interactions. The analysis will just consider two groups of seeds
    # The skip_ddi variable will be used to skip the drug-drug interactions checker if this situation occurs
    skip_ddi = False
    if options.drug1 == None and options.drug2 == None:
        print("  DIANA INFO:\tAs there are no drug names, the analysis will be performed without looking for a drug-drug interaction.\n\t\tDIANA will only compare and look for a relation between the two groups of seeds provided.\n")
        options.drug1 = "seeds1"
        options.drug2 = "seeds2"
        skip_ddi = True

    # If just one name for a drug is provided, the other will be called drug1 or drug2
    if options.drug1 == None and options.drug2 != None:
        print(" DIANA INFO:\tNo name has been provided for drug 1, so it will be called as 'drug1'\n")
        options.drug1 = "drug1"

    if options.drug1 != None and options.drug2 == None:
        print(" DIANA INFO:\tNo name has been provided for drug 2, so it will be called as 'drug2'\n")
        options.drug2 = "drug2"


    # Taking drugs and seeds in a list
    list_seeds_drugs = [ (options.drug1, options.seeds1), (options.drug2, options.seeds2) ]

    # Looping through each drug
    for (drug, seeds) in list_seeds_drugs:

        # Create a directory for each drug
        try:
            os.stat(drug)
            print("  DIANA INFO:\tIt already exists a directory for %s.\n" %(drug))
        except:
            print("  DIANA INFO:\tCreating directory for %s.\n" %(drug))
            os.mkdir(drug)

	    #----------------------#
	    #   SEEDS CONTROLLER   #
	    #----------------------#

	    # SEEDS CONTROLLER: Checks the seeds provided by the user. If necessary, performs a search 
	    # in DrugBank to obtain more seeds

        # Check if the seeds files are provided
        if not fileExist(seeds):

            # If the seeds file is not provided but the name of the drug is known, seeds will be searched in DrugBank
            if drug != 'drug1' and drug != 'drug2' and drug != None:
                print("  DIANA INFO:\tFile with input seeds for %s is missing. Seeds will be searched in DrugBank database.\n" % (drug))
                name = check_two_words_drug(drug)
                seeds = dbs.search_seeds_in_drugbank(name, provided_seeds = None)
                print("  DIANA INFO:\tThe seeds found in DrugBank for %s are:\n\t\t%s\n" % (name, seeds))

                # If the number of seeds is still small, the program is stopped
                if len(seeds) < 3:
                    print("  DIANA INFO:\tThe number of seeds in the network is below 3.\n\t\t\tThe analysis will not be reliable with this small number of seeds.\n\t\t\tTry to find a drug with more reported targets.\n\t\t\tSorry for the inconvenience.\n")
                    sys.exit(10)
                # If the number of seeds is 3 or 4, the analysis will not be considered highly reliable but will carry on
                elif len(seeds) > 2 and len(seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 3 and 5.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
                # If the number of seeds is above 4, the analysis will be considered highly reliable
                else:
                    print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")

            # If the seeds file is not provided and the name of the drug is unknown, the program finishes with error
            else:
                print("  DIANA INFO:\tFile with input seeds for %s is missing and the name of the drug is unknown.\n\t\tPlease, provide a file with seeds for the drug AND/OR the name of the drug in order to search for seeds in DrugBank.\n" % (drug))
                sys.exit(10)

        # If the seeds file is provided, the seeds are processed
        else:
            if drug == list_seeds_drugs[0][0]:
                seeds1 = process_seeds(seeds)
                seeds = seeds1
            else:
                seeds2 = process_seeds(seeds)
                seeds = seeds2

            print("  DIANA INFO:\tThe seeds provided for %s are:\n\t\t%s\n" % (drug, seeds))

            #------------------------#
            #   SEEDS TYPE CHECKER   #
            #------------------------#

            # SEEDS TYPE CHECKER: Checks the type of seeds provided by the user. Returns 'uniprotentry', 'geneid', 'other' or 'different'
            seeds_type = check_type_of_seed(seeds)
            if seeds_type == 'uniprotentry':
                # If the seeds are in UniprotEntry, they are translated to GeneID using the function translate_seeds()
                (seeds, seeds_file) = translate_seeds(seeds, options.taxid)
                # The new seed files and seed list are reassigned
                if drug == list_seeds_drugs[0][0]:
                    seeds1 = process_seeds(seeds)
                    seeds = seeds1
                    options.seeds1 = seeds_file
                else:
                    seeds2 = process_seeds(seeds)
                    seeds = seeds2
                    options.seeds2 = seeds_file
            # If the seeds are other types or more than one type (different), the program is stopped
            elif seeds_type == 'other':
                print("  DIANA INFO:\tFile with input seeds for %s contains an inappropriate type of seed.\n\t\tPlease, provide a file containing GeneID or UniprotEntry seeds.\n" % (drug))
                sys.exit(10)
            elif seeds_type == 'different':
                print("  DIANA INFO:\tFile with input seeds for %s contains different types of seed in the same file.\n\t\tPlease, provide a file containing only GeneID or UniprotEntry seeds.\n" % (drug))
                sys.exit(10)


            # If the number of seeds is below 3, new seeds will be searched in DrugBank
            if len(seeds) < 3:
                print("  DIANA INFO:\tThe number of seeds in the network is below 3.\n\t\tThe analysis will not be reliable with this small number of seeds.\n\t\tWe will try to find more seeds in DrugBank.\n")
                name = check_two_words_drug(drug)
                seeds = dbs.search_seeds_in_drugbank(name, provided_seeds = seeds)
                print("  DIANA INFO:\tThe new list of seeds for %s are:\n\t\t%s\n" % (name, seeds))

                # If the number of seeds is still small, the program is stopped
                if len(seeds) < 3:
                    print("  DIANA INFO:\tThe number of seeds in the network is below 3.\n\t\t\tThe analysis will not be reliable with this small number of seeds.\n\t\t\tTry to find a drug with more reported targets.\n\t\t\tSorry for the inconvenience.\n")
                    sys.exit(10)
                # If the number of seeds is 3 or 4, the analysis will not be considered highly reliable but will carry on
                elif len(seeds) > 2 and len(seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 3 and 5.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
                # If the number of seeds is above 4, the analysis will be considered highly reliable
                else:
                    print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")

            # If the number of seeds is 3 or 4, the analysis will not be highly reliable, so more seeds will be searched in DrugBank
            elif len(seeds) > 2 and len(seeds) < 5:
                print("  DIANA INFO:\tThe number of seeds in the network is between 3 and 5.\n\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n\t\tMore seeds will be searched in DrugBank\n")
                name = check_two_words_drug(drug)
                seeds = dbs.search_seeds_in_drugbank(name, provided_seeds = seeds)
                print("  DIANA INFO:\tThe new list of seeds for %s are:\n\t\t%s\n" % (name, seeds))

                # If the number of seeds is 3 or 4, the analysis will not be considered highly reliable but will carry on
                if len(seeds) > 2 and len(seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 3 and 5.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
                # If the number of seeds is above 4, the analysis will be considered highly reliable
                else:
                    print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")

            # If the number of seeds is above 4, the analysis will be considered highly reliable
            else:
                print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")
                continue

        # Assign the new list of seeds to its corresponding variable
        if drug == list_seeds_drugs[0][0]:
            seeds1 = seeds
        else:
            seeds2 = seeds

        # Checking if a directory for seeds exists
        # If not, we will create it so that we have the seeds ordered in the same directory
        try:
            os.stat("seeds")
        except:
            print("  DIANA INFO:\tCreating directory for seeds.\n")
            os.mkdir("seeds")

        # Write a file with the new seeds (if there have been new seeds added) 
        # If not, the loop has been skipped and this will not be written
        seeds_file = 'seeds/'+drug+'.newseeds'
        f = open(seeds_file,"w")
        for seed in seeds:
            f.write("%s\n" %(seed))
        f.close()


	#----------------------#
	#   SEEDS COMPARISON   #
	#----------------------#

    print("  DIANA INFO:\t#### STARTING SEEDS COMPARISON ####\n")
    print("  DIANA INFO:\tStarting comparison between seed node profiles.\n")
    seeds1_dict = generate_seeds_dict_for_comparison(seeds1)
    seeds2_dict = generate_seeds_dict_for_comparison(seeds2)
    summary_seed_node = comp.calculate_comparison(seeds1_dict, seeds2_dict)
    
    print("  DIANA INFO:\tStarting comparison between seed functional profiles.\n")

    for (drug, seeds) in list_seeds_drugs:
        if drug == list_seeds_drugs[0][0]:
            seeds = seeds1
        else:
            seeds = seeds2

        # Perform a functional enrichment analysis using the seeds
        seed_functional_profile = drug+"/seed_functional_profile.txt"
        if not fileExist(seed_functional_profile):
            sf = open(seed_functional_profile, 'w')
            functional_enrichment.check_functional_enrichment(seeds, seeds, "geneid",
                                                                  sf.write, species = "Homo sapiens")
            sf.close()
        else:
            print("  DIANA INFO:\tFunctional profile was already created for %s, so the enrichment analysis has been skipped\n" %(drug))

    enrichment_drug1 = options.drug1+"/seed_functional_profile.txt"
    enrichment_drug2 = options.drug2+"/seed_functional_profile.txt"
    # From the enrichment file, obtain a dict with GO terms as keys and a set of Log of odds and Adjusted pval as values
    enrichment_dict_drug1 = comp.parse_enrichment_file(enrichment_drug1)
    enrichment_dict_drug2 = comp.parse_enrichment_file(enrichment_drug2)
    summary_seed_functional = comp.calculate_comparison(enrichment_dict_drug1, enrichment_dict_drug2)


    # Check if the SIF file is provided
    ## If the SIF file is not provided, a network of expansion will be generated for each of the drugs using their seeds
    if options.sif is None:
        print("  DIANA INFO:\tThere is no SIF input file, so networks will be created using seeds provided.\n")

	#--------------------------------------------#
	#   GENERATE NETWORKS OF EXPANSION (BIANA)   #
	#--------------------------------------------#

        for drug in (options.drug1, options.drug2):

            print("  DIANA INFO:\t#### STARTING NETWORK GENERATION AND NETWORK SCORING PROCESS FOR %s ####\n" %(drug.upper()))

            # Assign to each drug its corresponding seeds file
            if (drug == options.drug1):
                seeds_list = seeds1
                new_seeds_file = 'seeds/'+drug+'.newseeds'
                if fileExist(new_seeds_file):
                    seeds = new_seeds_file
                else:
                    seeds = options.seeds1
            else:
                seeds_list = seeds2
                new_seeds_file = 'seeds/'+drug+'.newseeds'
                if fileExist(new_seeds_file):
                    seeds = new_seeds_file
                else:
                    seeds = options.seeds2

            nodes_file = drug+"/human_nodes_BIANA_"+drug+".txt"
            edges_file = drug+"/human_edges_BIANA_"+drug+".txt"

            # If there are no nodes or edges file, the network will be generated. Else, this part will be skipped
            if not fileExist(nodes_file) and not fileExist(edges_file):

                # Generate network from seeds provided running the script "generate_netscore_files.py"
                print("  DIANA INFO:\tGenerating network for %s. This can take a few minutes...\n" %(drug))
                run_generate_netscore_files(drug, seeds, options.radius, options.taxid, options.restriction)
                print("  DIANA INFO:\tNetwork finished for %s.\n" %(drug))

                # Change the name and folder of the translation of seeds
                translation = drug+'/translation_'+drug+'_seeds2BIANAcodes.txt'
                command = 'mv translation_seeds_to_BIANA_codes.txt '+translation
                os.system(command)

                # Create the folder where the results will be stored
                try:
                    os.stat(drug+"/guild_results")
                except:
                    print("  DIANA INFO:\tCreating directory for GUILD results of %s.\n" %(drug))
                    os.mkdir(drug+"/guild_results")

                # Obtain the seeds which have been included inside the network
                # There could be seeds with no interactions which do not appear in the network
                # If the new number of seeds is too small (< 3) the analysis is stopped
                print("  DIANA INFO:\tFinding the seeds which appear in the network of %s.\n" %(drug))
                nodes_file = drug+"/human_nodes_BIANA_"+drug+".txt"
                new_seeds_file = drug+"/"+drug+"_seeds_guild.txt"
                new_seeds = find_seeds_in_nodes_file(drug, nodes_file, new_seeds_file)
                print("  DIANA INFO:\tThere are %d seeds in the network from the %d initial ones.\n" %(len(new_seeds), len(seeds_list)) )
                if len(new_seeds) < 3:
                    print("  DIANA INFO:\tThe number of seeds in the network is below 3.\n\t\t\tThe analysis will not be reliable with this small number of seeds.\n\t\t\tTry to find a drug with more reported targets.\n\t\t\tSorry for the inconvenience.\n")
                    sys.exit(10)
                elif len(new_seeds) > 2 and len(new_seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 3 and 5.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
                else:
                    print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")


                # Transform the edges file to SIF format (network format)
                print("  DIANA INFO:\tTransforming the edges file to SIF format.\n")
                edges_file = drug+"/human_edges_BIANA_"+drug+".txt"
                sif_file = drug+"/network_"+drug+"_guild.sif"
                command = 'cat '+edges_file+'|cut -f 1,2|sed -e "s/\t//" > '+sif_file
                os.system(command)

            #--------------------------------#
            #  SCORING OF NETWORKS (GUILD)   #
            #--------------------------------#

                # Run GUILD
                print("  DIANA INFO:\tRunning GUILD (network scoring program). This can take a few minutes...\n")
                data_dir = drug
                results_dir = drug+'/guild_results'

                guild_command = "/soft/devel/python-2.7/bin/python toolbox/run_guild.py %s %s %s %s" % (data_dir, new_seeds_file, sif_file, results_dir)
                os.system(guild_command)
                print("\n  DIANA INFO:\tGUILD has finished.\n")

            else:
                print("  DIANA INFO:\tThe network of expansion for %s was already done and it has been skipped.\n" %(drug))


            #---------------------------------#
            #   TRANSLATION PROCESS (BIANA)   #
            #---------------------------------#

            print("  DIANA INFO:\t#### STARTING TRANSLATION PROCESS FOR %s ####\n" %(drug.upper()))

            # Translate from BIANA code to GeneID

            ## Translation of the NODES in the network to GeneID codes
            translation_type = "geneid"
            network_input_file = results_dir+"/edge_scores.sif"
            nodes_input_file = results_dir+"/output_scores.sif.netcombo.pval"
            translation_file = drug+"/translation_BIANA_%s.txt" %(drug)
            output_file = results_dir+"/interactions_"+translation_type
            # If the output file does not exist, the translation will be performed. Else, it will be skipped
            if not fileExist(output_file+".nodes") and not fileExist(output_file+".edges"):
                print("  DIANA INFO:\tTranslation of network\n")
                run_transform_network(network_input_file, nodes_input_file, translation_file, output_file)
                #run_translate_network(network_input_file, nodes_input_file, translation_file, options.taxid, translation_type, output_file)
            else:
                print("  DIANA INFO:\tThe translation of the network to %s was already done and it has been skipped\n" %(translation_type))

            ## Translation of the SEEDS in the network to GeneID codes
            print("  DIANA INFO:\tTranslation of seeds\n")
            seeds_input_file = results_dir+"/seeds.txt"
            output_file = results_dir+"/seeds_"+translation_type

            run_transform_network(network_input_file, seeds_input_file, translation_file, output_file)

            # Remove the node file created, which is not useful anymore
            created_file = output_file+".edges"
            command = "rm %s" %(created_file)
            os.system(command)

            #### NOTE: The seeds translation has been skipped for purpose of speeding the analysis, as it is not an important part for the analysis
            # ### Translation of the seeds file
            # seeds_file = results_dir+"/seeds.txt"
            # output_file = results_dir+"/seeds_"+translation_type
            # # If the output file does not exist, the translation will be performed. Else, it will be skipped
            # if not fileExist(output_file+".txt"):
            #     run_translate_network(network_input_file, seeds_file, options.taxid, translation_type, output_file, translation_of_nodes_file = translation_file)
            #     command = "mv %s.nodes %s.txt" %(output_file, output_file)
            #     os.system(command)
            #     command = "rm %s.edges" %(output_file)
            #     os.system(command)
            # else:
            #     print("  DIANA INFO:\tThe translation of the seeds to %s was already done and it has been skipped" %(translation_type))

            ## Translation into RefSeq
            # translation_type = "refseq"
            # output_file = results_dir+"/interactions_"+translation_type
            # # If the output file does not exist, the translation will be performed. Else, it will be skipped
            # if not fileExist(output_file+".nodes") and not fileExist(output_file+".edges"):
            #     run_translate_network(network_input_file, nodes_input_file, options.taxid, translation_type, output_file)
            # else:
            #     print("  DIANA INFO:\tThe translation of the network to %s was already done and it has been skipped" %(translation_type))


            print("  DIANA INFO:\tTranslation process finished for %s.\n" %(drug.upper()))



    ## If the SIF file is provided, the analysis will directly use GUILD to score the network
    else:

	#-------------------------------#
	#   USING A SIF FILE PROVIDED   #
	#-------------------------------#

        print("  DIANA INFO:\tThere is a SIF file provided, so %s will be used for profiles creation.\n" % (options.sif))

        for drug in (options.drug1, options.drug2):

            print("  DIANA INFO:\t#### STARTING NETWORK SCORING PROCESS FOR %s ####\n" %(drug.upper()))

            # Assign to each drug its corresponding seeds file
            if (drug == options.drug1):
                seeds = options.seeds1
                seeds_list = seeds1
            else:
                seeds = options.seeds2
                seeds_list = seeds2

            # Create a directory for the drug
            try:
                os.stat(drug)
            except:
                print("  DIANA INFO:\tCreating directory for %s.\n" %(drug))
                os.mkdir(drug)

            # Create the folder where the results will be stored
            results_dir = drug+"/guild_results_using_sif"
            try:
                os.stat(results_dir)
            except:
                print("  DIANA INFO:\tCreating directory for GUILD results of %s.\n" %(drug))
                os.mkdir(results_dir)

            # Obtain the seeds which have been included inside the network
            # There could be seeds with no interactions which do not appear in the network
            print("  DIANA INFO:\tFinding the seeds which appear in the network %s.\n" %(options.sif))
            new_seeds_file = drug+"/"+drug+"_seeds_guild_using_sif.txt"
            new_seeds = find_seeds_in_sif_file(drug, options.sif, seeds_list, new_seeds_file)
            print("  DIANA INFO:\tThere are %d seeds in the network from the %d initial ones.\n" %(len(new_seeds), len(seeds_list)) )

            #--------------------------------#
            #  SCORING OF NETWORKS (GUILD)   #
            #--------------------------------#

            # Run GUILD
            print("  DIANA INFO:\tRunning GUILD (network scoring program).\n")
            data_dir = drug

            guild_command = "/soft/devel/python-2.7/bin/python toolbox/run_guild.py %s %s %s %s" % (data_dir, new_seeds_file, options.sif, results_dir)
            os.system(guild_command)
            print("  DIANA INFO:\tGUILD has finished.\n")


    profiles_list_node = []
    profiles_list_edge = []
    profiles_list_functional = []
    profiles_list_linker_node = []
    profiles_list_linker_edge = []
    profiles_list_linker_functional = []

    for drug in (options.drug1, options.drug2):

        if options.sif is None:
            results_dir = drug+"/guild_results"
        else:
            results_dir = drug+"/guild_results_using_sif"

		#-------------------------#
		#   TOP SCORING PROCESS   #
		#-------------------------#

        print("  DIANA INFO:\tStarting Top Scoring process for %s.\n" %(drug.upper()))

        # TOP SCORE procedure: obtain a file of the uniprotentries under a given threshold and a 
        # second file of the functional enrichment using refseq

        entry_type = "geneid"
        pvalue_file = results_dir+"/interactions_"+entry_type+".nodes"
        top_threshold = options.top_threshold
        pval_threshold = options.pval_threshold
        run_top_scoring(results_dir, pvalue_file, top_threshold, pval_threshold, entry_type)

        # entry_type = "refseq"
        # pvalue_file = results_dir+"/interactions_"+entry_type+".nodes"
        # top_threshold = options.top_threshold_functional
        # pval_threshold = options.pval_threshold_functional
        # run_top_scoring(results_dir, pvalue_file, top_threshold, pval_threshold, entry_type)

        print("  DIANA INFO:\tTop scoring process finished for %s.\n" %(drug.upper()))

		#-------------------------#
		#   TOP NETWORK PROCESS   #
		#-------------------------#

        if options.speed == 1:

            print("  DIANA INFO:\tStarting Top network process for %s.\n" %(drug.upper()))

            # TOP NETWORK procedure (edge profile creation): obtain a file of scored edges and a
            # subnetwork created using nodes ander a given threshold
            output_file = results_dir+"/edge_profile.sif"
            if not fileExist(output_file):
                entry_type = "geneid"
                pvalue_file = results_dir+"/interactions_"+entry_type+".nodes"
                network_file = results_dir+"/interactions_"+entry_type+".edges"
                seed_file = results_dir+"/seeds_%s.nodes"%(entry_type)
                top_threshold = options.top_threshold_edges
                pval_threshold = options.pval_threshold_edges
                run_top_network(results_dir, pvalue_file, network_file, seed_file, top_threshold, pval_threshold)

                # Translation of LINKER NODE PROFILE to REFSEQ
                network_input_file = results_dir+"/subnetwork_linkers.sif"
                nodes_input_file = results_dir+"/linker_node_profile.txt"
                # input_format = "geneid"
                # taxid = int(options.taxid)
                # translation_type = "refseq"
                # output_file = results_dir+"/linker_node_profile_"+translation_type

                # command = "/soft/devel/python-2.7/bin/python toolbox/TranslateNetwork.py -i %s -n %s -iformat %s -x -taxid %d -ttype %s -o %s" % (network_input_file, nodes_input_file, input_format, taxid, translation_type, output_file)
                # os.system(command)

                # Parse the nodes of the translated linker node profile to create the linker functional profile
                linker_node_profile = nodes_input_file
                linker_node_profile_list = []
                fl = open(linker_node_profile, 'r')
                for line in fl:
                    linker_node = line.strip().split()[0]
                    linker_node_profile_list.append(linker_node)
                fl.close()

                # Creation of the LINKER FUNCTIONAL PROFILE
                linker_functional_profile = results_dir+"/linker_functional_profile.txt"
                ff = open(linker_functional_profile, 'w')
                functional_enrichment.check_functional_enrichment(linker_node_profile_list, linker_node_profile_list, "geneid",
                                                                      ff.write, species = "Homo sapiens")
                ff.close()


                print("  DIANA INFO:\tTop network process finished for %s.\n" %(drug.upper()))
            else:
                print("  DIANA INFO:\tThe top network procedure was already done and it will be skipped\n")


		#---------------------------------#
		#   EXTENDED ANALYSIS PROCEDURE   #
		#---------------------------------#

        # EXTENDED ANALYSIS PROCEDURE: Build the profiles defined by different tops

        if options.extended != None:

            # Read the list of top thresholds
            f = open(options.extended,"r")
            line = f.readline()
            top_thresholds_list = line.split()
            f.close()

            print("  DIANA INFO:\tStarting EXTENDED ANALYSIS process for %s.\n" %(drug.upper()))

            import toolbox.extended as ext

            # Run the Extended Analysis for the NODE PROFILE

            pvalue_file = results_dir+"/interactions_geneid.nodes"
            profiles_list_node.append(ext.run_extended_analysis_node(results_dir, pvalue_file, top_thresholds_list, "nodes"))

            # Run the Extended Analysis for the EDGE PROFILE

            edge_profile_scored = results_dir+"/edge_profile.sif"
            profiles_list_edge.append(ext.run_extended_analysis_edge(results_dir, edge_profile_scored, top_thresholds_list, "edges"))

            # Run the Extended Analysis for the FUNCTIONAL PROFILE

            functional_profile = results_dir+"/enrichment.txt"
            profiles_list_functional.append(ext.run_extended_analysis_functional(results_dir, functional_profile, top_thresholds_list, "functional"))
            
            print("  DIANA INFO:\tExtended analysis finished for %s.\n" %(drug.upper()))

            # Run the Extended Analysis for the LINKER NODE PROFILE

            linker_node_profile = results_dir+"/linker_node_profile.txt"
            profiles_list_linker_node.append(ext.run_extended_analysis_node(results_dir, linker_node_profile, top_thresholds_list, "linkers"))

            # Run the Extended Analysis for the LINKER EDGE PROFILE

            linkers_subnetwork = results_dir+"/subnetwork_linkers.sif"
            profiles_list_linker_edge.append(ext.run_extended_analysis_edge(results_dir, linkers_subnetwork, top_thresholds_list, "linkers"))

            # Run the Extended Analysis for the LINKER FUNCTIONAL PROFILE

            linkers_functional_profile = results_dir+"/linker_functional_profile.txt"
            profiles_list_linker_functional.append(ext.run_extended_analysis_functional(results_dir, linkers_functional_profile, top_thresholds_list, "linkers"))


	#----------------------------#
	#   COMPARISON OF PROFILES   #
	#----------------------------#

    print("\n  DIANA INFO:\t#### STARTING COMPARISON OF PROFILES ####\n")


    print("  DIANA INFO:\tStarting comparison between node profiles.\n")

    if options.sif is None:
        results_drug1 = options.drug1+"/guild_results"
        results_drug2 = options.drug2+"/guild_results"
    else:
        results_drug1 = options.drug1+"/guild_results_using_sif"
        results_drug2 = options.drug2+"/guild_results_using_sif"

    if options.speed == 1:
        pvalue_file_drug1 = results_drug1+"/interactions_geneid.nodes"
        pvalue_file_drug2 = results_drug2+"/interactions_geneid.nodes"
    else:
        pvalue_file_drug1 = results_drug1+"/output_scores.sif.netcombo.pval"
        pvalue_file_drug2 = results_drug2+"/output_scores.sif.netcombo.pval"


    node_to_vals_drug1 = GU.get_values_from_pvalue_file(pvalue_file_drug1)
    node_to_vals_drug2 = GU.get_values_from_pvalue_file(pvalue_file_drug2)

    summary_node = comp.calculate_comparison(node_to_vals_drug1, node_to_vals_drug2)

    summary = {}
    # Keeping the Spearman results inside this dict
    summary["node"] = summary_node

    # If speed is 1 or 2, the analysis of the Edges and Functional profiles will be performed
    if options.speed is 1 or options.speed is 2:

        print("\n  DIANA INFO:\tStarting comparison between edge profiles.\n")

        edge_profile_drug1 = results_drug1+"/edge_profile.sif"
        edge_profile_drug2 = results_drug2+"/edge_profile.sif"

        # From the edge profile, obtain a dict with edges as keys and a set of Score and useless string as values
        edges_dict_drug1 = comp.parse_edge_profile(edge_profile_drug1)
        edges_dict_drug2 = comp.parse_edge_profile(edge_profile_drug2)

        summary_edges = comp.calculate_comparison(edges_dict_drug1, edges_dict_drug2)
        summary["edges"] = summary_edges

        print("\n  DIANA INFO:\tStarting comparison between functional enrichment profiles\n")

        enrichment_drug1 = results_drug1+"/enrichment.txt"
        enrichment_drug2 = results_drug2+"/enrichment.txt"

        # From the enrichment file, obtain a dict with GO terms as keys and a set of Log of odds and Adjusted pval as values
        enrichment_dict_drug1 = comp.parse_enrichment_file(enrichment_drug1)
        enrichment_dict_drug2 = comp.parse_enrichment_file(enrichment_drug2)

        summary_functional = comp.calculate_comparison(enrichment_dict_drug1, enrichment_dict_drug2)
        summary["functional"] = summary_functional


	#------------------------------------#
	#   DRUG-DRUG INTERACTIONS CHECKER   #
	#------------------------------------#

    # Check if drug-drug interactions has to be skipped or not
    if skip_ddi == False:
        # Check that the input drugs are not the names used by default
        if options.drug1 != 'drug1' and options.drug2 != 'drug2':

            print("\n  DIANA INFO:\t#### DRUG-DRUG INTERACTIONS CHECKER ####\n")

            print("  DIANA INFO:\tSearching if there is an experimental reported DDI in DrugBunk.\n")

            name1 = check_two_words_drug(options.drug1)
            name2 = check_two_words_drug(options.drug2)

            # Retrieve the description of the interaction between the two drugs if there is interaction
            # The description could be different depending on the order of the drugs, so the procedure is 
            # performed twice changing the order of the drugs
            description1 = dbs.search_DDI_in_drugbank(name1, name2)
            description2 = dbs.search_DDI_in_drugbank(name2, name1)

            if description1 == [] and description2 ==[]:
                print("  DIANA INFO:\tNo experimentally reported DDI found in DrugBunk between the two drugs.\n")
            else:
                print("  DIANA INFO:\tThe DDI between the two drugs is experimentally reported in DrugBank \n")
                if description1 == description2:
                    print("\t\tDescription of the DDI:\n")
                    print("\t\t%s\n" %(description1[0]))
                else:
                    print("\t\tDescription of %s vs. %s:\n\n\t\t%s\n" %(name1, name2, description1))
                    print("\t\tDescription of %s vs. %s:\n\n\t\t%s\n" %(name2, name1, description2))


    # Create a RESULTS FOLDER
    results_folder = 'results_'+options.drug1+'_'+options.drug2
    try:
        os.stat(results_folder)
        print("  DIANA INFO:\tIt already exists a directory for the results of %s and %s.\n" %(options.drug1, options.drug2))
    except:
        print("  DIANA INFO:\tCreating directory for %s and %s results.\n" %(options.drug1, options.drug2))
        os.mkdir(results_folder)


	#------------------------------------------------------#
	#   CREATE RESULTS DOCUMENT OF THE EXTENDED ANALYSIS   #
	#------------------------------------------------------#

    # EXTENDED ANALYSIS PROCEDURE: Output the comparison between the profiles of the different drugs

    output_file_general = results_folder+'/general_extended_analysis_%s_%s' %(options.drug1, options.drug2)
    output_file_linker = results_folder+'/linker_extended_analysis_%s_%s' %(options.drug1, options.drug2)
    if options.extended != None:
        ext.extended_analysis_comparison(profiles_list_node, profiles_list_edge, profiles_list_functional, top_thresholds_list, output_file_general)
        ext.extended_analysis_comparison(profiles_list_linker_node, profiles_list_linker_edge, profiles_list_linker_functional, top_thresholds_list, output_file_linker)

	#-----------------------------------------------------#
	#   CREATE RESULTS DOCUMENT OF THE GENERAL ANALYSIS   #
	#-----------------------------------------------------#

    # RESULTS: Output a document with the results of the general analysis
    if options.speed == 1:
        results_doc = results_folder+"/results_%s_%s.txt" %(options.drug1, options.drug2)
    else:
        results_doc = results_folder+"/results_SPEED%d_%s_%s.txt" %(options.speed, options.drug1, options.drug2)
    fr = open(results_doc,"w")
    fr.write("#### INPUTS ####\n")
    fr.write("Seeds for %s:\t%s\n" %(options.drug1, seeds1))
    fr.write("Seeds for %s:\t%s\n" %(options.drug2, seeds2))
    if options.sif is None:
        fr.write("SIF file not used\n")
    else:
        fr.write("SIF file used:\t%s" %(options.sif))
    fr.write("SPEED OF THE ANALYSIS: %s\n" %(options.speed))
    fr.write("%% of top used: %s. P-value used: %0.3f\n" %(options.top_threshold, float(options.pval_threshold)))
    fr.write("%% of top used (for edges): %s. P-value used (for edges): %0.3f\n" %(options.top_threshold_edges, float(options.pval_threshold_edges)))
    fr.write("%% of top used (for functional): %s. P-value used (for functional): %0.3f\n" %(options.top_threshold_functional, float(options.pval_threshold_functional)))
    fr.write("Radius of expansion network: %s. TaxID restriction: %s. Experiment restriction: %s\n" %(options.radius, options.taxid, options.restriction))

    fr.write("\n#### RESULTS ####\n")
    fr.write("## Node profile comparison ##\n")
    fr.write("Spearman: %0.3f (P-value: %0.3f)\n" %(summary["node"][1],summary["node"][2]))
    fr.write("Dot product: %0.3f\n" %(summary["node"][0]))

    if options.speed == 1 or options.speed == 2:
        fr.write("## Edge profile comparison ##\n")
        fr.write("Spearman: %0.3f (P-value: %0.3f)\n" %(summary["edges"][1],summary["edges"][2]))
        fr.write("Dot product: %0.3f\n" %(summary["edges"][0]))
        fr.write("## Functional profile comparison ##\n")
        fr.write("Spearman: %0.3f (P-value: %0.3f)\n" %(summary["functional"][1],summary["functional"][2]))
        fr.write("Dot product: %0.3f\n" %(summary["functional"][0]))

    if options.speed == 1:

        fr.write("\n#### NODE PROFILE common nodes ####\n")

        node_profile1 = results_drug1+'/node_profile.txt'
        node_profile2 = results_drug2+'/node_profile.txt'
        common_nodes = comp.compare_profiles('node', node_profile1, node_profile2)
        for node in common_nodes:
            fr.write("%s " %(node))

        fr.write("\n\n#### EDGE PROFILE common edges ####\n")

        edge_profile1 = results_drug1+'/edge_profile.sif'
        edge_profile2 = results_drug2+'/edge_profile.sif'
        common_edges = comp.compare_profiles('edge', edge_profile1, edge_profile2)
        for edge in common_edges:
            fr.write("%s " %(edge))

        fr.write("\n\n#### FUNCTIONAL PROFILE common GOs ####\n")

        functional_profile1 = results_drug1+'/enrichment.txt'
        functional_profile2 = results_drug2+'/enrichment.txt'
        common_GO = comp.compare_profiles('functional', functional_profile1, functional_profile2, pval = options.pval_threshold)
        for GO in common_GO:
            fr.write("%s " %(GO))
        fr.write("\n")

    # Check if drug-drug interactions has to be skipped or not
    if skip_ddi == False:
        # Check that the input drugs are not the names used by default
        if options.drug1 != 'drug1' and options.drug2 != 'drug2':
            fr.write("\n#### DRUG-DRUG INTERACTIONS CHECKER ####\n")
            if description1 == [] and description2 ==[]:
                fr.write("No experimentally reported DDI found in DrugBunk between the two drugs.\n")
            else:
                fr.write("The DDI between the two drugs is experimentally reported in DrugBank.\n")
                if description1 == description2:
                    fr.write("Description of the DDI:\n%s\n" %(description1[0]))
                else:
                    fr.write("Description of %s vs. %s:\n%s\n" %(option.drug1, option.drug2, description1))
                    fr.write("Description of %s vs. %s:\n%s\n" %(option.drug2, option.drug1, description2))


    # End marker for time
    end = time.time()
    fr.write("\n#### TIME OF EXECUTION ####\n")
    fr.write("The time of execution of the analysis is: %0.3f sec. or %0.3f minutes.\n" %(end - start, (end - start) / 60))

    print("\n  DIANA INFO:\t#### TIME OF EXECUTION ####\n")
    print ("  DIANA INFO:\tThe time of execution of the analysis is: %0.3f sec. or %0.3f minutes.\n" %(end - start, (end - start) / 60))


    print("\t\t--------------------------------------------\n")
    print("\t\tAnalysis finished. Thank you for using DIANA\n")
    print("\t\t--------------------------------------------\n")


if  __name__ == "__main__":
    main()

