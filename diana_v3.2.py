import sys
import os
import argparse
import re
import indigo_package.indigo as indigo
import toolbox
import toolbox.comparison as comp
import toolbox.database_searcher as dbs
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
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-d1','--drug1_name',dest='drug1',action = 'store',
                        help = """Name of the drug number 1. If you do not provide seeds for this drug or the number of seeds is not large enough,
                        the program will use this name to search for seeds in BIANA database. If seeds are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between 
                        double quotes""")
    parser.add_argument('-d2','--drug2_name',dest='drug2',action = 'store',
                        help = """Name of the drug number 2. If you do not provide seeds for this drug or the number of seeds is not large enough,
                        the program will use this name to search for seeds in BIANA database. If seeds are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between 
                        double quotes""")
    parser.add_argument('-s1','--seeds_drug1',dest='seeds1',action = 'store',default='drug1.seeds',
                        help = 'Input file with seeds of drug 1. Each seed must be separated by a newline character (default is drug1.seeds)')
    parser.add_argument('-s2','--seeds_drug2',dest='seeds2',action = 'store',default='drug2.seeds',
                        help = 'Input file with seeds of drug 2. Each seed must be separated by a newline character (default is drug2.seeds)')
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = 'Input file with a protein-protein interaction network in SIF format (if not introduced, networks of expansion are used)')
    parser.add_argument('-tis','--type_if_sif',dest='type_if_sif',action = 'store',
                        help = 'If sif file introduced, state the annotation type of the elements in the network (for example: uniprotentry, geneid, biana_codes)')
    parser.add_argument('-tra','--trans_if_sif',dest='trans_if_sif',action = 'store',
                        help = 'If sif file introduced, but the type of notation is not "geneid", introduce a translation file to GeneID if possible, so that it is easier to translate when necessary (optional)')
    parser.add_argument('-t','--top',dest='top_threshold',action = 'store',default='10',
                        help = 'Define the percentage of top scoring nodes to create the node / edge / functional profile (default is 10)')
    parser.add_argument('-p','--pvalue',dest='pval_threshold',action = 'store',default='1',
                        help = 'Define the pvalue threshold to filter the top nodes and create the node / edge / functional profile (default is 1)')
    parser.add_argument('-rad','--radius',dest='radius',action = 'store',default='3',
                        help = 'Define the radius of expansion for the creation of the network from seeds (default is 3)')
    parser.add_argument('-tax','--taxid',dest='taxid',action = 'store',default='9606',
                        help = 'Define the restriction of species for the creation of the network from seeds using a Taxonomy ID (default is 9606 (human))')
    parser.add_argument('-res','--restriction',dest='restriction',action = 'store',default='Y2H',
                        help = """Define an experiment restriction for the creation of the network from seeds (default is Y2H). 
                        Options:
                        - AFF: Use interactions at least described by affinity methods (i.e. Tandem Affinity Purification)
                        - Y2H: Use interactions at least described by yeast two hybrid methods (Y2H)
                        - eAFF: Use all interactions except those described by affinity methods (i.e. Tandem Affinity Purification)
                        - eY2H: Use all interactions except those described by yeast two hybrid methods (Y2H)
                        - None: Not use experiment restrictions
                        """)
    parser.add_argument('-ex','--extended',dest='extended',action = 'store', default=None,
                        help = """Extended analysis. It performs a comparison of the Spearman coefficients using several Top/P-val cut-offs. The input consists in:
                        - A file with different top threshold values separated by spaces. For example, a file called "top_thresholds.list" containing:
                        0.5 1 2 3 4 5 10 15 20 25 30 40 50
                        Default: None
                        """)
    parser.add_argument('-db','--database',dest='database',action = 'store',default='BIANA_JAN_2017',
                        help = """Define the database to use for the generation of the network of expansion. There are three options: 
                        - BIANA_MARCH_2013
                        - BIANA_2016
                        - BIANA_JAN_2017
                        (default is BIANA_JAN_2017)""")


    options=parser.parse_args()

    return options

#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def check_special_drug_name(drug):
    """
    If a drug contains a single quote (') --> it is replaced by "_squote_".
    If it is not composed of two words, returns the initial drug name
    """
    # If it contains single quote, it is replaced by _squote_
    drug = drug.replace("'", "_squote_")
    # If it contains parentheses, they are replaced by _openpar_ or _closepar_
    drug = drug.replace("(", "_openpar_")
    drug = drug.replace(")", "_closepar_")

    # If it a name made of more than one word, the space is replaced by _
    processed_name = drug.split(" ")

    if len(processed_name) > 1:
        return "_".join(processed_name)
    else:
        return drug

def detect_special_drug_name(drug):
    """
    Detect if a drug has a special name, and return the real name
    """
    # Detect single quotes
    drug = drug.replace("_squote_", "'")
    # If it contains parentheses, they are replaced by _openpar_ or _closepar_
    drug = drug.replace("_openpar_", "(")
    drug = drug.replace("_closepar_", ")")

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
        seed = str(seed)
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
    command = "cp "+new_seeds_file+" "+"data/"+drug+"/guild_results/seeds.txt"
    os.system(command)

    return new_seeds

def find_seeds_in_sif_file(drug, sif_file_items, initial_seeds, new_seeds_file):
    """Finds the seeds that are in a sif file and creates a new seeds file. Returns the seeds in a list"""

    fs = open(new_seeds_file,"w")
    new_seeds = []

    for seed in initial_seeds:
        if str(seed) in sif_file_items:
            new_seeds.append(seed)
            fs.write("%s 1.00000\n" %(seed))

    fs.close()

    # Produce a copy of the seeds file in the results directory
    command = "cp "+new_seeds_file+" "+"data/"+drug+"/guild_results_using_sif/seeds.txt"
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

def get_smiles_similarity(smiles1, smiles2, fp_type = "sub", metric = "tanimoto"):
    """
    fp_type: sim | sub
      - "sim" = similarity fingerprints: shorter
      - "sub" = substructure fingerprints: more descriptive
    metric: tanimoto | tversky
    """
    if len(smiles1) == 0 or len(smiles2) == 0:
        return None
    ind = indigo.Indigo()
    m = ind.loadMolecule(smiles1)
    m.aromatize()
    fp = m.fingerprint(fp_type)
    m2 = ind.loadMolecule(smiles2)
    m2.aromatize() # Aromatize molecules in case they are not in aromatic form
    fp2 = m2.fingerprint(fp_type) # Calculate similarity between "similarity" fingerprints
    d = ind.similarity(fp, fp2, metric)
    return d

def obtain_all_geneids_in_sif_file(sif_file):
    """Creates a set with all the geneids in a sif file"""

    fn = open(sif_file,"r")
    sif_file_geneids = set()

    for line in fn:
        words = line.split()
        sif_file_geneids.add(words[0])
        sif_file_geneids.add(words[2])

    fn.close()
    return sif_file_geneids

def process_seeds(seeds_file):
    """Takes the seeds from the input seeds file and stores them into a list"""

    fd = open(seeds_file,"r")
    seeds = []

    for line in fd:
        seeds.append(line.strip())

    fd.close()

    return seeds

def run_generate_netscore_files(drug, seeds_file, radius, taxid, restriction, seeds_type, translation_type, database):
    """Runs the program generate_netscore_files.py that generates a network of expansion from a seeds file"""

    nodes = "data/"+drug+"/human_nodes_BIANA_"+drug+".txt"
    loc = "data/"+drug+"/noLoc_bPPI_"+drug+".txt"
    edges = "data/"+drug+"/human_edges_BIANA_"+drug+".txt"
    trans = "data/"+drug+"/translation_BIANA_"+drug+".txt"
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
    if database == "BIANA_MARCH_2013":
        db = ""
    elif database == "BIANA_2016":
        db = "_v2016"
    elif database == "BIANA_JAN_2017":
        db = "_v2017"
    if seeds_type == "geneid":
        addgeneid = "_geneID"
    else:
        addgeneid = ""

    command = "/soft/devel/python-2.7/bin/python2.7 toolbox/generate_netscore_files%s%s.py -iseed %s -radius %d -stype %s -ttype %s -score 0.0 -taxid %d -node %s -loc %s -edge %s -trans %s -format guild %s " % (addgeneid, db, seeds_file, radius, seeds_type, translation_type, taxid, nodes, loc, edges, trans, restriction)
    os.system(command)

    return

def run_score_network(network_file, pvalue_file, seed_file, top_threshold, pval_threshold, data_dir):
    """
       Runs the program top_network.py which creates:
       - A sif file which contains the edges of the filtered nodes scored
    """

    top_threshold = float(top_threshold)
    pval_threshold = float(pval_threshold)

    command = "/soft/devel/python-2.7/bin/python toolbox/score_network.py {} {} {} {:f} {:f} {}".format( network_file, pvalue_file, seed_file, top_threshold, pval_threshold, data_dir )
    os.system(command)

    return

def run_top_scoring(data_dir, pvalue_file, top_threshold, pval_threshold, analysis_type):
    """
       Runs the program top_scoring_combination.py that has the following characteristics:
       - Takes the entries which are at a given percentage of the top scoring ones
       - Then, filters the top scoring entries by a given threshold of p-value
       - If analysis_type = "node": Creates a node profile containing Id/Score/P-value from a pvalue file in BIANA codes
       - If entry_type = "functional": Creates a functional enrichment profile from the top scoring nodes of a translated pvalue file in GeneID
    """

    top_threshold = float(top_threshold)
    pval_threshold = float(pval_threshold)

    command = "/soft/devel/python-2.7/bin/python  toolbox/top_scoring_combination.py %s %s %s %f %f %s " % (data_dir, pvalue_file, data_dir, top_threshold, pval_threshold, analysis_type)
    os.system(command)

    return

def run_transform_network(network_input_file, nodes_input_file, translation_file, output_file):
    """Runs the program TransformNetwork.py that translates a file in BIANA codes to a desired notation, using the previously obtained translation file from generate_netscore_files"""

    output_edges_file = output_file + ".edges"
    output_nodes_file = output_file + ".nodes"

    command = "/soft/devel/python-2.7/bin/python toolbox/TransformNetwork.py -i %s -n %s -trans %s -oe %s -on %s" % (network_input_file, nodes_input_file, translation_file, output_edges_file, output_nodes_file)
    os.system(command)

    return

def run_translate_network(network_input_file, nodes_input_file, taxid, translation_type, output_file, database):
    """Runs the program TranslateNetwork.py that translates a file in BIANA codes to a desired notation"""

    taxid = int(taxid)

    if database == "BIANA_MARCH_2013":
        db = ""
    elif database == "BIANA_2016":
        db = "_v2016"
    elif database == "BIANA_JAN_2017":
        db = "_v2017"

    command = "/soft/devel/python-2.7/bin/python toolbox/TranslateNetwork%s.py -i %s -n %s -x -taxid %d -ttype %s -o %s" % (db, network_input_file, nodes_input_file, taxid, translation_type, output_file)
    os.system(command)

    return

def translate_seeds(seeds_list, input_type, translation_type, taxid, drug, database):
    """
    Uses TranslateNetwork.py to translate a list of seeds from an input_type to a translation_type (uniprotentry to geneid or geneid to uniprotentry).
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
    output_file = "seeds/"+drug+".seeds."+input_type
    taxid = int(taxid)
    if database == "BIANA_MARCH_2013":
        db = ""
    elif database == "BIANA_2016":
        db = "_v2016"
    elif database == "BIANA_JAN_2017":
        db = "_v2017"

    if input_type == "uniprotentry":
        command = "/soft/devel/python-2.7/bin/python toolbox/TranslateNetwork%s.py -i %s -n %s -iformat %s -x -taxid %d -ttype %s -o %s" % (db, seeds_file, seeds_file, input_type, taxid, translation_type, output_file)
    elif input_type == "geneid":
        command = "/soft/devel/python-2.7/bin/python toolbox/TranslateNetwork%s.py -i %s -n %s -iformat %s -xf -taxid %d -ttype %s -o %s" % (db, seeds_file, seeds_file, input_type, taxid, translation_type, output_file)
    os.system(command)
    # Remove the edge file created, which is not useful
    command = "rm %s.edges" %(output_file)
    os.system(command)

    # Get the seeds from the node file created
    new_seeds = set()
    created_file = output_file+".nodes"
    f = open(created_file, 'r')
    for line in f:
        words = line.split()
        new_seed = words[0]
        new_seeds.add(new_seed)
    f.close()

    # Write a final seeds file with the seeds correctly written
    final_seeds_file = "seeds/"+drug+".seeds."+input_type+".trans"
    fs = open(final_seeds_file, 'w')
    for item in new_seeds:
        fs.write("%s\n" %(item))
    fs.close()

    # Remove the foo and the node file created, which is not useful anymore
    command = "rm %s" %(seeds_file)
    os.system(command)
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
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Thank you for your trust in our software\n")
    print("\t\t------------------------------------------------------------------------------------------------------------------------\n")

    # Start marker for time measure
    start = time.time()

    # Inform if EXTENDED analysis selected
    if options.extended != None:
        print("  DIANA INFO:\tEXTENDED ANALYSIS ACTIVATED\n")

    # Inform if SIF file provided
    if options.sif is not None:
        print("  DIANA INFO:\tSIF file provided.\n")
        # If the type of elements in the SIF file is not provided, the analysis is stopped
        if options.type_if_sif is None:
            print("  DIANA INFO:\tPlease, introduce the type of annotation of the elements in the SIF file in '--type_if_sif' (i.e. uniprotentry, geneid, biana_codes).\n")
            sys.exit(10)

    # If no names for drug1 or drug2 are provided, the analysis will be performed without looking for drug-drug interactions. The analysis will just consider two groups of seeds
    if options.drug1 == None and options.drug2 == None:
        print("  DIANA INFO:\tAs there are no drug names, the analysis will be performed without looking for a drug-drug interaction.\n\t\tDIANA will only compare and look for a relation between the two groups of seeds provided.\n")
        options.drug1 = "seeds1"
        options.drug2 = "seeds2"

    # If just one name for a drug is provided, the other will be called drug1 or drug2
    if options.drug1 == None and options.drug2 != None:
        print(" DIANA INFO:\tNo name has been provided for drug 1, so it will be called 'drug1'\n")
        options.drug1 = "drug1"

    if options.drug1 != None and options.drug2 == None:
        print(" DIANA INFO:\tNo name has been provided for drug 2, so it will be called 'drug2'\n")
        options.drug2 = "drug2"

    # Check if the names of the drugs provided have special characters (parentheses, single quotes) or contain more than one word
    if options.drug1 != None:
        options.drug1 = check_special_drug_name(options.drug1)
    if options.drug2 != None:
        options.drug2 = check_special_drug_name(options.drug2)


    # Create a directory for the data
    try:
        os.stat("data")
    except:
        os.mkdir("data")

    # Create a directory for the random networks
    try:
        os.stat("data/random_networks")
    except:
        os.mkdir("data/random_networks")

    # Create a directory for the results
    try:
        os.stat("results")
    except:
        os.mkdir("results")


    # Taking drugs and seeds in a list
    list_seeds_drugs = [ (options.drug1, options.seeds1), (options.drug2, options.seeds2) ]

    # Looping through each drug
    for (drug, seeds) in list_seeds_drugs:

        # Create a directory for each drug
        try:
            os.stat("data/"+drug)
            print("  DIANA INFO:\tIt already exists a directory for %s.\n" %(drug))
        except:
            print("  DIANA INFO:\tCreating directory for %s.\n" %(drug))
            os.mkdir("data/"+drug)

	    #----------------------#
	    #   SEEDS CONTROLLER   #
	    #----------------------#

	    # SEEDS CONTROLLER: Checks the seeds provided by the user. If necessary, performs a search 
	    # in BIANA database to obtain more seeds

        # Check if the seeds files are provided
        if not fileExist(seeds):

            # If the seed file is not provided but the name of the drug is known, seeds will be searched in BIANA database
            if drug != 'drug1' and drug != 'drug2' and drug != None:
                print("  DIANA INFO:\tFile with input seeds for %s is missing. Seeds will be searched in BIANA database.\n" % (drug))
                name = detect_special_drug_name(drug)
                seeds = dbs.search_seeds_in_file(name, provided_seeds = None)
                #seeds = dbs.search_seeds_in_drugbank(name, provided_seeds = None)
                seeds_type = check_type_of_seed(seeds)
                print("  DIANA INFO:\tThe seeds found in DrugBank for %s are:\n\t\t%s\n" % (name, seeds))

                # If no targets found, the program is stopped
                if len(seeds) < 1:
                    print("  DIANA INFO:\tNo targets are found for this drug.\n\t\t\tTry to find a drug with more reported targets.\n\t\t\tSorry for the inconvenience.\n")
                    sys.exit(10)
                # If the number of seeds is 1 to 4, the analysis will not be considered highly reliable but will carry on
                elif len(seeds) > 0 and len(seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 1 and 4.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
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
            if seeds_type == 'other':
                print("  DIANA INFO:\tFile with input seeds for %s contains an inappropriate type of seed.\n\t\tPlease, provide a file containing GeneID or UniprotEntry seeds.\n" % (drug))
                sys.exit(10)
            elif seeds_type == 'different':
                print("  DIANA INFO:\tFile with input seeds for %s contains different types of seed in the same file.\n\t\tPlease, provide a file containing only GeneID or UniprotEntry seeds.\n" % (drug))
                sys.exit(10)

            # If there are no seeds, new seeds will be searched in the database
            if len(seeds) < 1:
                print("  DIANA INFO:\tNo targets are found for this drug.\n\t\tWe will try to find targets in the database.\n")
                name = detect_special_drug_name(drug)
                seeds = dbs.search_seeds_in_drugbank(name, provided_seeds = seeds)
                print("  DIANA INFO:\tThe new list of seeds for %s are:\n\t\t%s\n" % (name, seeds))

                # If the number of seeds is still small, the program is stopped
                if len(seeds) < 1:
                    print("  DIANA INFO:\tNo targets are found for this drug.\n\t\t\tTry to find a drug with more reported targets.\n\t\t\tSorry for the inconvenience.\n")
                    sys.exit(10)
                # If the number of seeds is 1 to 4, the analysis will not be considered highly reliable but will carry on
                elif len(seeds) > 0 and len(seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 1 and 4.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
                # If the number of seeds is above 4, the analysis will be considered highly reliable
                else:
                    print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")

            # If the number of seeds is 1 to 4, the analysis will not be highly reliable, so more seeds will be searched in DrugBank
            elif len(seeds) > 0 and len(seeds) < 5:
                print("  DIANA INFO:\tThe number of seeds in the network is between 1 and 4.\n\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n\t\tMore seeds will be searched in DrugBank\n")
                name = detect_special_drug_name(drug)
                seeds_type = check_type_of_seed(seeds)
                seeds = dbs.search_seeds_in_drugbank(name, provided_seeds = seeds)
                print("  DIANA INFO:\tThe new list of seeds for %s are:\n\t\t%s\n" % (name, seeds))

                # If the number of seeds is 1 to 4, the analysis will not be considered highly reliable but will carry on
                if len(seeds) > 0 and len(seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 1 and 4.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
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
        seeds_file = 'seeds/'+drug+'.newseeds.'+seeds_type
        f = open(seeds_file,"w")
        for seed in seeds:
            f.write("%s\n" %(seed))
        f.close()

        # Assign the new seeds file to its corresponding variable
        if drug == list_seeds_drugs[0][0]:
            options.seeds1 = seeds_file
        else:
            options.seeds2 = seeds_file


	#----------------------#
	#   SEEDS COMPARISON   #
	#----------------------#

    print("  DIANA INFO:\t#### STARTING SEEDS COMPARISON ####\n")
    print("  DIANA INFO:\tStarting comparison between seed node profiles.\n")
    # Comparison using only the seeds and the Spearman / Dot product comparison
    print("  DIANA INFO:\t## Spearman's coefficient / Dot product comparison between targets ##\n")
    seeds1_dict = generate_seeds_dict_for_comparison(seeds1)
    seeds2_dict = generate_seeds_dict_for_comparison(seeds2)
    summary_seed_node = comp.calculate_comparison(seeds1_dict, seeds2_dict)
    seeds_jaccard = comp.calculate_jaccard_index(set(seeds1), set(seeds2))
    summary_seed_node.append(seeds_jaccard)
    print("  DIANA INFO:\tJaccard index for the comparison of targets: {}\n".format(seeds_jaccard))
    # Comparison using PFAM families of the seeds and the Spearman / Dot product comparison
    print("  DIANA INFO:\t## Spearman's coefficient / Dot product comparison between PFAMs ##\n")
    pfams1 = dbs.search_pfams_for_seeds(seeds1)
    pfams2 = dbs.search_pfams_for_seeds(seeds2)
    pfams1_dict = generate_seeds_dict_for_comparison(pfams1)
    pfams2_dict = generate_seeds_dict_for_comparison(pfams2)
    summary_seed_pfam = comp.calculate_comparison(pfams1_dict, pfams2_dict)
    pfams_jaccard = comp.calculate_jaccard_index(set(pfams1), set(pfams2))
    summary_seed_pfam.append(pfams_jaccard)
    print("  DIANA INFO:\tJaccard index for the comparison of PFAMs: {}\n".format(pfams_jaccard))

    print("  DIANA INFO:\tStarting comparison between seed functional profiles.\n")

    # Taking drugs and seeds in a list
    list_seeds_drugs = [ (options.drug1, options.seeds1), (options.drug2, options.seeds2) ]

    for (drug, seeds) in list_seeds_drugs:
        if drug == list_seeds_drugs[0][0]:
            seeds1 = process_seeds(seeds)
            seeds = seeds1
        else:
            seeds2 = process_seeds(seeds)
            seeds = seeds2

        # Check again seeds type (there could be changes)
        seeds_type = check_type_of_seed(seeds)

        # Perform a functional enrichment analysis using the seeds
        seed_functional_profile = "data/"+drug+"/seed_functional_profile.txt"

        if seeds_type == "geneid":
            seeds_func_analysis = seeds
        elif seeds_type == "uniprotentry":
            # If the seeds are in UniprotEntry, they are translated to GeneID using the function translate_seeds()
            # The new seeds are reassigned in a variable called "seeds_geneid" so that we do not loose the uniprotentry seeds
            (seeds_geneid, seeds_file_geneid) = translate_seeds(seeds, "uniprotentry", "geneid", options.taxid, drug, options.database)
            seeds_func_analysis = list(seeds_geneid)

            ## This is for the analysis when using a SIF file
            ## If the sif is in geneid and the seeds are in uniprotentry, the geneid translated seeds are reassigned as main seeds
            if options.sif != None and options.type_if_sif == 'geneid':
                if drug == list_seeds_drugs[0][0]:
                    options.seeds1 = seeds_file_geneid
                else:
                    options.seeds2 = seeds_file_geneid
                seeds_type = 'geneid'

        if not fileExist(seed_functional_profile):
            sf = open(seed_functional_profile, 'w')
            functional_enrichment.check_functional_enrichment(seeds_func_analysis, seeds_func_analysis, "geneid",
                                                                  sf.write, species = "Homo sapiens")
            sf.close()
        else:
            print("  DIANA INFO:\tFunctional profile was already created for %s, so the enrichment analysis has been skipped\n" %(drug))

    enrichment_drug1 = "data/"+options.drug1+"/seed_functional_profile.txt"
    enrichment_drug2 = "data/"+options.drug2+"/seed_functional_profile.txt"
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
                seeds = options.seeds1
            else:
                seeds = options.seeds2
            seeds_list = process_seeds(seeds)

            nodes_file = "data/"+drug+"/human_nodes_BIANA_"+drug+".txt"
            edges_file = "data/"+drug+"/human_edges_BIANA_"+drug+".txt"

            # If there are no nodes or edges file, the network will be generated. Else, this part will be skipped
            if not fileExist(nodes_file) and not fileExist(edges_file):

                # Generate network from seeds provided running the script "generate_netscore_files.py"
                print("  DIANA INFO:\tGenerating network for %s. This can take a few minutes...\n" %(drug))
                # The seeds type will indicate if the program uses the general script generate_netscore_files or the specific one for GENEID 
                # The translation type is important to be included in generate_netscore_files because it creates a translation
                # file that can be used later to translate the network with "TransformNetwork.py"
                # If the type of protein is uniprotentry, we will use geneid as translation type in order to translate to geneid
                # and compute the functional enrichment analysis (it cannot be done with uniprotentries!!)
                translation_type = "geneid"
                run_generate_netscore_files(drug, seeds, options.radius, options.taxid, options.restriction, seeds_type, translation_type, options.database)
                print("  DIANA INFO:\tNetwork finished for %s.\n" %(drug))

                # Change the name and folder of the translation of seeds
                translation = "data/"+drug+'/translation_'+drug+'_seeds2BIANAcodes.txt'
                command = 'mv translation_seeds_to_BIANA_codes.txt '+translation
                os.system(command)

                # Create the folder where the results will be stored
                try:
                    os.stat("data/"+drug+"/guild_results")
                except:
                    print("  DIANA INFO:\tCreating directory for GUILD results of %s.\n" %(drug))
                    os.mkdir("data/"+drug+"/guild_results")

                # Obtain the seeds which have been included inside the network
                # There could be seeds with no interactions which do not appear in the network
                # If there are no seeds in the network, the analysis is stopped
                print("  DIANA INFO:\tFinding the seeds which appear in the network of %s.\n" %(drug))
                nodes_file = "data/"+drug+"/human_nodes_BIANA_"+drug+".txt"
                new_seeds_file = "data/"+drug+"/"+drug+"_seeds_guild.txt"
                new_seeds = find_seeds_in_nodes_file(drug, nodes_file, new_seeds_file)
                print("  DIANA INFO:\tThere are %d seeds in the network from the %d initial ones.\n" %(len(new_seeds), len(seeds_list)) )
                if len(new_seeds) < 1:
                    print("  DIANA INFO:\tThere are no seeds in the network.\n\t\t\tTry to find a drug with more reported targets.\n\t\t\tSorry for the inconvenience.\n")
                    sys.exit(10)
                elif len(new_seeds) > 0 and len(new_seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 1 and 4.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
                else:
                    print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")


                # Transform the edges file to SIF format (network format)
                print("  DIANA INFO:\tTransforming the edges file to SIF format.\n")
                edges_file = "data/"+drug+"/human_edges_BIANA_"+drug+".txt"
                sif_file = "data/"+drug+"/network_"+drug+"_guild.sif"
                command = 'cat '+edges_file+'|cut -f 1,2|sed -e "s/\t//" > '+sif_file
                os.system(command)

                #--------------------------------#
                #  SCORING OF NETWORKS (GUILD)   #
                #--------------------------------#

                # Run GUILD
                print("  DIANA INFO:\tRunning GUILD (network scoring program). This can take a few minutes...\n")
                data_dir = "data/"+drug
                results_dir = "data/"+drug+'/guild_results'
                random_networks_dir = 'data/random_networks'

                guild_command = "/soft/devel/python-2.7/bin/python toolbox/run_guild.py %s %s %s %s %s" % (data_dir, new_seeds_file, sif_file, results_dir, random_networks_dir)
                os.system(guild_command)
                print("\n  DIANA INFO:\tGUILD has finished.\n")

            else:
                print("  DIANA INFO:\tThe network of expansion for %s was already done and it has been skipped.\n" %(drug))



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
            else:
                seeds = options.seeds2
            seeds_list = process_seeds(seeds)

            # Create the folder where the results will be stored
            results_dir = "data/"+drug+"/guild_results_using_sif"
            try:
                os.stat(results_dir)
            except:
                print("  DIANA INFO:\tCreating directory for GUILD results of %s.\n" %(drug))
                os.mkdir(results_dir)


            pval_file = results_dir+"/output_scores.sif.netcombo.pval"

            # If the p-value does not file exist, it means that the scoring network procedure with GUILD needs to be done
            # Else, it will be skipped
            if not fileExist(pval_file):

                # Obtain the seeds which have been included inside the network
                # There could be seeds with no interactions which do not appear in the network
                print("  DIANA INFO:\tFinding the seeds which appear in the network %s, from the seeds file %s.\n" %(options.sif, seeds))

                new_seeds_file = "data/"+drug+"/"+drug+"_seeds_guild_using_sif.txt"
                sif_file_geneids = obtain_all_geneids_in_sif_file(options.sif)
                new_seeds = find_seeds_in_sif_file(drug, sif_file_geneids, seeds_list, new_seeds_file)

                print("  DIANA INFO:\tThere are %d seeds in the network from the %d initial ones.\n" %(len(new_seeds), len(seeds_list)) )

                if len(new_seeds) < 1:
                    print("  DIANA INFO:\tThere are no seeds in the network.\n\t\t\tTry to find a drug with more reported targets.\n\t\t\tSorry for the inconvenience.\n")
                    sys.exit(10)
                elif len(new_seeds) > 0 and len(new_seeds) < 5:
                    print("  DIANA INFO:\tThe number of seeds in the network is between 1 and 4.\n\t\t\tThe analysis with GUILD can still be reliable, but it could improve with more seeds.\n")
                else:
                    print("  DIANA INFO:\tThere is a large enough number of seeds. The analysis will be highly reliable\n")

                #--------------------------------#
                #  SCORING OF NETWORKS (GUILD)   #
                #--------------------------------#

                # Run GUILD
                print("  DIANA INFO:\tRunning GUILD (network scoring program).\n")
                data_dir = "data/"+drug
                random_networks_dir = 'data/random_networks'

                guild_command = "/soft/devel/python-2.7/bin/python toolbox/run_guild.py %s %s %s %s %s" % (data_dir, new_seeds_file, options.sif, results_dir, random_networks_dir)
                os.system(guild_command)
                print("  DIANA INFO:\tGUILD has finished.\n")

            else:
                print("  DIANA INFO:\tThe scoring of the network with GUILD for %s was already done and it has been skipped.\n" %(drug))


    profiles_list_node = []
    profiles_list_functional = []
    profiles_list_edge = []
    pvalue_files = []
    scored_networks = []

    for drug in (options.drug1, options.drug2):

        if options.sif is None:
            results_dir = "data/"+drug+"/guild_results"
        else:
            results_dir = "data/"+drug+"/guild_results_using_sif"

		#-------------------------#
		#   TOP SCORING PROCESS   #
		#-------------------------#

        print("  DIANA INFO:\t#### STARTING TOP SCORING PROCESS FOR %s ####\n" %(drug.upper()))

        # TOP SCORE procedure: 
        #     - In a first run, obtain the top scoring nodes of the network in BIANA codes
        #     - Translate the node profile to GENEID, so that the functional enrichment analysis can be performed
        #     - In a second run, obtain a functional enrichment analysis from the top scoring nodes of a translated network in GeneID
        #     - If the type used is UniprotEntry, the node profile will be translated into UniprotEntry

        # Generate the NODE profile from the top % scoring nodes
        analysis_type = "node"
        pvalue_file = results_dir+"/output_scores.sif.netcombo.pval"
        top_threshold = options.top_threshold
        pval_threshold = options.pval_threshold
        run_top_scoring(results_dir, pvalue_file, top_threshold, pval_threshold, analysis_type)

        print("  DIANA INFO:\tNODE profile obtained from the top scoring nodes of the network.\n")

        # If we do not provide SIF OR we provide it but in a notation which is not GeneID but with a translation to GeneID file... we can transform the network
        if options.sif == None or (options.sif != None and options.trans_if_sif != None and options.type_if_sif != 'geneid'):
            ## Translation of the NODE PROFILE to GENEID
            translation_type = "geneid"
            network_input_file = results_dir+"/edge_scores.sif"
            nodes_input_file = results_dir+"/node_profile.txt"
            translation_file = "data/"+drug+"/translation_BIANA_%s.txt" %(drug)
            output_file = results_dir+"/node_profile_"+translation_type
            if options.sif != None and options.trans_if_sif != None:
                translation_file = options.trans_if_sif

            # If the output file does not exist, the translation will be performed. Else, it will be skipped
            if not fileExist(output_file+".nodes") and not fileExist(output_file+".edges"):
                print("  DIANA INFO:\tTranslating the NODE profile obtained to GENEID to generate a FUNCTIONAL profile.\n")
                run_transform_network(network_input_file, nodes_input_file, translation_file, output_file)
            else:
                print("  DIANA INFO:\tThe translation of the network to %s was already done and it has been skipped.\n" %(translation_type))

            ## Translation of the SEEDS in the network to the translation type
            seeds_input_file = results_dir+"/seeds.txt"
            output_file = results_dir+"/seeds_"+translation_type
            run_transform_network(network_input_file, seeds_input_file, translation_file, output_file)

            # Remove the edges file created, which is not useful anymore
            created_file = output_file+".edges"
            command = "rm %s" %(created_file)
            os.system(command)

            print("  DIANA INFO:\tTranslation process finished.\n")

        else:
            # If we provide a SIF file but in a notation which is not GeneID and without translation file... we need to translate the network
            if options.sif != None and options.trans_if_sif == None and options.type_if_sif != 'geneid':
                nodes_file = results_dir+"/node_profile.txt"
                nodes_input_file = results_dir+"/node_profile_mod.txt"
                translation_type = "geneid"
                output_file = results_dir+"/node_profile_"+translation_type

                # We create a nodes file without the header to use it as a network as well, and make faster the translation
                nf = open(nodes_file, 'r')
                nif = open(nodes_input_file, 'w')
                for line in nf:
                    if line[0] != "#":
                        nif.write(line)
                nf.close()
                nif.close()

                # If the output file does not exist, the translation will be performed. Else, it will be skipped
                if not fileExist(output_file+".nodes"):
                    print("  DIANA INFO:\tTranslating the NODE profile obtained for %s to UniprotEntry.\n" %(drug.upper()))
                    run_translate_network(nodes_input_file, nodes_input_file, options.taxid, translation_type, output_file, options.database)
                    print("  DIANA INFO:\tTranslation process finished.\n")
                else:
                    print("  DIANA INFO:\tThe translation of the NODE profile to %s was already done and it has been skipped.\n" %(translation_type))

                # Remove the node input file created, which is not useful anymore
                command = "rm %s" %(nodes_input_file)
                os.system(command)


        # Generate the FUNCTIONAL profile from the translated node profile
        node_profile_geneid = results_dir+"/node_profile_geneid.nodes"
        if options.sif != None and options.type_if_sif == 'geneid': # if we have used a geneid SIF, the node profile will be node_profile.txt
            node_profile_geneid = results_dir+"/node_profile.txt"
        enrichment_file = results_dir + "/enrichment.txt"
        # Get GUILD scores
        node_to_vals = GU.get_values_from_pvalue_file(node_profile_geneid)
        # Get all nodes
        all_nodes = node_to_vals.keys()
        # Calculate the enrichment analysis
        ef = open(enrichment_file, 'w')
        functional_enrichment.check_functional_enrichment(all_nodes, all_nodes, "geneid",
                                                              ef.write, species = "Homo sapiens")
        ef.close()

        print("  DIANA INFO:\tFUNCTIONAL profile obtained from the top scoring nodes of the network.\n")


        # ## Translation of the NODE PROFILE created into UniprotEntry
        # ## (only if there is no sif file and the initial seeds type was uniprotentry)
        # if seeds_type == "uniprotentry" and options.sif == None:
        #     nodes_file = results_dir+"/node_profile.txt"
        #     nodes_input_file = results_dir+"/node_profile_mod.txt"
        #     translation_type = "uniprotentry"
        #     output_file = results_dir+"/node_profile_"+translation_type

        #     # We create a nodes file without the header to use it as a network as well, and make faster the translation
        #     nf = open(nodes_file, 'r')
        #     nif = open(nodes_input_file, 'w')
        #     for line in nf:
        #         if line[0] != "#":
        #             nif.write(line)
        #     nf.close()
        #     nif.close()

        #     # If the output file does not exist, the translation will be performed. Else, it will be skipped
        #     if not fileExist(output_file+".nodes"):
        #         print("  DIANA INFO:\tTranslating the NODE profile obtained for %s to UniprotEntry.\n" %(drug.upper()))
        #         run_translate_network(nodes_input_file, nodes_input_file, options.taxid, translation_type, output_file, options.database)
        #         print("  DIANA INFO:\tTranslation process finished.\n")
        #     else:
        #         print("  DIANA INFO:\tThe translation of the NODE profile to %s was already done and it has been skipped.\n" %(translation_type))

        #     # Remove the node input file created, which is not useful anymore
        #     command = "rm %s" %(nodes_input_file)
        #     os.system(command)

        print("  DIANA INFO:\tTop scoring process finished for %s.\n" %(drug.upper()))


		#-----------------------------#
		#   NETWORK SCORING PROCESS   #
		#-----------------------------#

        print("  DIANA INFO:\t#### STARTING TOP NETWORK PROCESS FOR %s ####\n" %(drug.upper()))

        # TOP NETWORK procedure (edge profile creation): obtain a file of scored edges and a
        # subnetwork created using nodes ander a given threshold
        output_file = results_dir+"/edge_profile.sif"
        if not fileExist(output_file):
            pvalue_file = results_dir+"/output_scores.sif.netcombo.pval"
            network_file = results_dir+"/edge_scores.sif"
            seed_file = results_dir+"/seeds.txt"
            top_threshold = options.top_threshold
            pval_threshold = options.pval_threshold
            run_score_network(network_file, pvalue_file, seed_file, top_threshold, pval_threshold, results_dir)

            # Translation of LINKER NODE PROFILE
            # If we do not have SIF file, or we have SIF file and we have provided a translation file for geneids, we can use TransformNetwork
            if options.sif == None or (options.sif != None and options.trans_if_sif != None and options.type_if_sif != 'geneid'):
                translation_type = "geneid"
                network_input_file = results_dir+"/network_linkers.sif"
                nodes_input_file = results_dir+"/linker_node_profile.txt"
                translation_file = "data/"+drug+"/translation_BIANA_%s.txt" %(drug)
                if options.sif != None and options.trans_if_sif != None: # in the case of using a SIF file and a translation file, the translation will be the trans provided
                    translation_file = options.trans_if_sif
                translated_linker_node_profile = results_dir+"/linker_node_profile_"+translation_type

                run_transform_network(network_input_file, nodes_input_file, translation_file, translated_linker_node_profile)

            # In the case of SIF file, if we do not provide a translation file the translation must be done with TranslateNetwork
            elif options.sif != None and options.trans_if_sif == None and options.type_if_sif != 'geneid':
                if options.type_if_sif != "geneid":
                    translation_type = "geneid"
                    network_input_file = results_dir+"/network_linkers.sif"
                    nodes_input_file = results_dir+"/linker_node_profile.txt"
                    translated_linker_node_profile = results_dir+"/linker_node_profile_"+translation_type

                    # If the output file does not exist, the translation will be performed. Else, it will be skipped
                    if not fileExist(translated_linker_node_profile+".nodes"):
                        print("  DIANA INFO:\tTranslating the LINKER NODE profile obtained for %s to %s.\n" %(drug.upper(), translation_type))
                        run_translate_network(nodes_input_file, nodes_input_file, options.taxid, translation_type, translated_linker_node_profile, options.database)
                        print("  DIANA INFO:\tTranslation process finished.\n")
                    else:
                        print("  DIANA INFO:\tThe translation of the LINKER NODE profile to %s was already done and it has been skipped.\n" %(translation_type))


            # Parse the nodes of the linker node profile to create the LINKER FUNCTIONAL PROFILE
            if options.sif != None and options.type_if_sif == 'geneid': # if we use a geneid SIF network, the linker node profile is the initial one
                linker_node_profile = results_dir+"/linker_node_profile.txt"
                translation_type = 'geneid'
            else: # else, it is the translated one
                linker_node_profile = translated_linker_node_profile+'.nodes'
            linker_node_profile_list = []
            fl = open(linker_node_profile, 'r')
            for line in fl:
                linker_node = line.strip().split()[0]
                linker_node_profile_list.append(linker_node)
            fl.close()

            # Creation of the LINKER FUNCTIONAL PROFILE
            linker_functional_profile = results_dir+"/linker_functional_profile.txt"
            ff = open(linker_functional_profile, 'w')
            functional_enrichment.check_functional_enrichment(linker_node_profile_list, linker_node_profile_list, translation_type,
                                                                  ff.write, species = "Homo sapiens")
            ff.close()


            print("  DIANA INFO:\tTop network process finished for %s.\n" %(drug.upper()))
        else:
            print("  DIANA INFO:\tThe top network procedure was already done and it will be skipped.\n")


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

            pvalue_file = results_dir+"/output_scores.sif.netcombo.pval"
            profiles_list_node.append(ext.run_extended_analysis_node(results_dir, pvalue_file, top_thresholds_list))
            pvalue_files.append(pvalue_file)

            # Run the Extended Analysis for the FUNCTIONAL PROFILE

            profiles_list_functional.append(ext.run_extended_analysis_functional(results_dir, top_thresholds_list))
            
            # Run the Extended Analysis for the EDGE PROFILE

            network_scored = results_dir+"/network_scored.sif"
            profiles_list_edge.append(ext.run_extended_analysis_edge(results_dir, network_scored, pvalue_file, top_thresholds_list))
            scored_networks.append(network_scored)

            print("  DIANA INFO:\tExtended analysis finished for %s.\n" %(drug.upper()))


	#----------------------------#
	#   COMPARISON OF PROFILES   #
	#----------------------------#

    print("\n  DIANA INFO:\t#### STARTING COMPARISON OF PROFILES ####\n")

    # Comparison of the Node profile
    print("  DIANA INFO:\tStarting comparison between node profiles.\n")

    if options.sif is None:
        results_drug1 = "data/"+options.drug1+"/guild_results"
        results_drug2 = "data/"+options.drug2+"/guild_results"
    else:
        results_drug1 = "data/"+options.drug1+"/guild_results_using_sif"
        results_drug2 = "data/"+options.drug2+"/guild_results_using_sif"


    pvalue_file_drug1 = results_drug1+"/output_scores.sif.netcombo.pval"
    pvalue_file_drug2 = results_drug2+"/output_scores.sif.netcombo.pval"

    node_to_vals_drug1 = GU.get_values_from_pvalue_file(pvalue_file_drug1)
    node_to_vals_drug2 = GU.get_values_from_pvalue_file(pvalue_file_drug2)

    summary_node = comp.calculate_comparison(node_to_vals_drug1, node_to_vals_drug2)

    summary = {}
    # Keeping the Spearman results inside this dict
    summary["node"] = summary_node


    print("\n  DIANA INFO:\tStarting comparison between edge profiles.\n")

    network_scored_drug1 = results_drug1+"/network_scored.sif"
    network_scored_drug2 = results_drug2+"/network_scored.sif"

    # From the edge profile, obtain a dict with edges as keys and a set of Score and useless string as values
    edges_dict_drug1 = comp.parse_edge_profile(network_scored_drug1)
    edges_dict_drug2 = comp.parse_edge_profile(network_scored_drug2)

    summary_edges = comp.calculate_comparison(edges_dict_drug1, edges_dict_drug2)
    summary["edges"] = summary_edges

    # Comparison of the Fuctional profile
    print("\n  DIANA INFO:\tStarting comparison between functional enrichment profiles\n")

    enrichment_drug1 = results_drug1+"/enrichment.txt"
    enrichment_drug2 = results_drug2+"/enrichment.txt"

    # From the enrichment file, obtain a dict with GO terms as keys and a set of Log of odds and Adjusted pval as values
    enrichment_dict_drug1 = comp.parse_enrichment_file(enrichment_drug1)
    enrichment_dict_drug2 = comp.parse_enrichment_file(enrichment_drug2)

    # In this case we will use the comparison of seeds function, because we do not have a complete functions file or anything similar
    # We will compare the functions directly
    summary_functional = comp.calculate_comparison(enrichment_dict_drug1, enrichment_dict_drug2)
    summary["functional"] = summary_functional

    print("\n  DIANA INFO:\tStarting comparison between linker node profiles\n")

    # Comparison of the Linker Node profile
    #                   -----------
    linker_node_profile_drug1 = results_drug1+"/linker_node_profile.txt"
    linker_node_profile_drug2 = results_drug2+"/linker_node_profile.txt"

    linkers_to_vals_drug1 = GU.get_values_from_pvalue_file(linker_node_profile_drug1)
    linkers_to_vals_drug2 = GU.get_values_from_pvalue_file(linker_node_profile_drug2)

    summary_linker_node = comp.calculate_comparison(linkers_to_vals_drug1, linkers_to_vals_drug2)
    summary["linker_node"] = summary_linker_node

    print("\n  DIANA INFO:\tStarting comparison between linker edge profiles\n")

    # Comparison of the Linker Edge profile
    #                   -----------
    linker_edge_profile_drug1 = results_drug1+"/network_linkers.sif"
    linker_edge_profile_drug2 = results_drug2+"/network_linkers.sif"

    link_edges_dict_drug1 = comp.parse_edge_profile(linker_edge_profile_drug1)
    link_edges_dict_drug2 = comp.parse_edge_profile(linker_edge_profile_drug2)

    summary_linker_edge = comp.calculate_comparison(link_edges_dict_drug1, link_edges_dict_drug2)
    summary["linker_edge"] = summary_linker_edge

    print("\n  DIANA INFO:\tStarting comparison between linker functional profiles\n")

    # Comparison of the Linker Functional profile
    #                   -----------------
    linker_functional_profile_drug1 = results_drug1+"/linker_functional_profile.txt"
    linker_functional_profile_drug2 = results_drug2+"/linker_functional_profile.txt"

    link_functional_dict_drug1 = comp.parse_enrichment_file(linker_functional_profile_drug1)
    link_functional_dict_drug2 = comp.parse_enrichment_file(linker_functional_profile_drug2)

    summary_linker_functional = comp.calculate_comparison(link_functional_dict_drug1, link_functional_dict_drug2)
    summary["linker_functional"] = summary_linker_functional


    #-------------------------------------#
    #   CALCULATE STRUCTURAL SIMILARITY   #
    #-------------------------------------#

    smiles1 = dbs.get_smiles_from_dcdbid(options.drug1)
    smiles2 = dbs.get_smiles_from_dcdbid(options.drug2)

    if smiles1 != None and smiles2 != None:
        struc_similarity = get_smiles_similarity(smiles1, smiles2, fp_type = "sub", metric = "tanimoto")
        print("  DIANA INFO:\tStructural similarity between the two drugs: %0.3f\n" %(struc_similarity))
    else:
        struc_similarity = 'NA'
        print("  DIANA INFO:\tStructural similarity unavailable")


    # Create a RESULTS FOLDER
    results_folder = 'results/results_'+options.drug1+'_'+options.drug2
    results_folder2 = 'results/results_'+options.drug2+'_'+options.drug2
    try:
        os.stat(results_folder)
        print("  DIANA INFO:\tIt already exists a directory for the results of %s and %s.\n" %(options.drug1, options.drug2))
    except:
        try:
            os.stat(results_folder2)
            print("  DIANA INFO:\tIt already exists a directory for the results of %s and %s.\n" %(options.drug1, options.drug2))
            results_folder = results_folder2
        except:
            print("  DIANA INFO:\tCreating directory for %s and %s results.\n" %(options.drug1, options.drug2))
            os.mkdir(results_folder)


	#------------------------------------------------------#
	#   CREATE RESULTS DOCUMENT OF THE EXTENDED ANALYSIS   #
	#------------------------------------------------------#

    # EXTENDED ANALYSIS PROCEDURE: Output the comparison between the profiles of the different drugs

    output_file_general = results_folder+'/extended_analysis_%s_%s.txt' %(options.drug1, options.drug2)
    if options.extended != None:
        ext.extended_analysis_comparison(profiles_list_node, profiles_list_edge, profiles_list_functional, pvalue_files, scored_networks, top_thresholds_list, output_file_general)

	#-----------------------------------------------------#
	#   CREATE RESULTS DOCUMENT OF THE GENERAL ANALYSIS   #
	#-----------------------------------------------------#

    # RESULTS: Output a document with the results of the general analysis
    results_doc = results_folder+"/results_%s_%s.txt" %(options.drug1, options.drug2)
    fr = open(results_doc,"w")
    fr.write("#### INPUTS ####\n")
    fr.write("# Seeds for %s:\t%s\n" %(options.drug1, seeds1))
    fr.write("# Seeds for %s:\t%s\n" %(options.drug2, seeds2))
    if options.sif is None:
        fr.write("# SIF file not used\n")
    else:
        fr.write("# SIF file used:\t%s\n" %(options.sif))
    fr.write("# %% of top used: %s. P-value used: %0.3f\n" %(options.top_threshold, float(options.pval_threshold)))
    if options.sif is None:
        fr.write("# Radius of expansion network: %s. TaxID restriction: %s. Experiment restriction: %s\n" %(options.radius, options.taxid, options.restriction))

    fr.write("\n#### RESULTS ####\n")
    fr.write("# Spearman\tP-value\tDot Product\tJaccard index (for seeds)\n")
    fr.write("## Node profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\n" %(summary["node"][1],summary["node"][2], summary["node"][0]))

    fr.write("## Edge profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\n" %(summary["edges"][1],summary["edges"][2], summary["edges"][0]))

    fr.write("## Functional profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\n" %(summary["functional"][1],summary["functional"][2], summary["functional"][0]))

    fr.write("## Linker node profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\n" %(summary["linker_node"][1],summary["linker_node"][2], summary["linker_node"][0]))

    fr.write("## Linker edge profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\n" %(summary["linker_edge"][1],summary["linker_edge"][2], summary["linker_edge"][0]))

    fr.write("## Linker functional profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\n" %(summary["linker_functional"][1],summary["linker_functional"][2], summary["linker_functional"][0]))

    fr.write("## Seeds node profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\t%0.3f\n" %(summary_seed_node[1],summary_seed_node[2], summary_seed_node[0], summary_seed_node[3]))

    fr.write("## Seeds PFAM profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\t%0.3f\n" %(summary_seed_pfam[1],summary_seed_pfam[2], summary_seed_pfam[0], summary_seed_pfam[3]))

    fr.write("## Seeds functional profile comparison ##\n")
    fr.write("%0.3f\t%0.3f\t%0.3f\n" %(summary_seed_functional[1],summary_seed_functional[2], summary_seed_functional[0]))

    fr.write("## Structural similarity ##\n")
    fr.write("%s\n" %(struc_similarity))

    fr.write("\n#### NODE PROFILE common nodes ####\n")
    if options.sif == None:
        node_profile1 = results_drug1+'/node_profile_%s.nodes' %(seeds_type)
        node_profile2 = results_drug2+'/node_profile_%s.nodes' %(seeds_type)
    else: # if SIF file used, the initial node profile will be used to get the common nodes
        node_profile1 = results_drug1+'/node_profile.txt'
        node_profile2 = results_drug2+'/node_profile.txt'
    common_nodes = comp.compare_profiles('node', node_profile1, node_profile2)
    fr.write("# ")
    for node in common_nodes:
        fr.write("%s " %(node))

    fr.write("\n\n#### EDGE PROFILE common edges ####\n")
    edge_profile1 = results_drug1+'/edge_profile.sif'
    edge_profile2 = results_drug2+'/edge_profile.sif'
    common_edges = comp.compare_profiles('edge', edge_profile1, edge_profile2)
    fr.write("# ")
    for edge in common_edges:
        fr.write("%s " %(edge))

    fr.write("\n\n#### FUNCTIONAL PROFILE common GOs ####\n")
    functional_profile1 = results_drug1+'/enrichment.txt'
    functional_profile2 = results_drug2+'/enrichment.txt'
    common_GO = comp.compare_profiles('functional', functional_profile1, functional_profile2, pval = options.pval_threshold)
    fr.write("# ")
    for GO in common_GO:
        fr.write("%s " %(GO))
    fr.write("\n")


    # End marker for time
    end = time.time()
    fr.write("\n#### TIME OF EXECUTION ####\n")
    fr.write("# The time of execution of the analysis is: %0.3f sec. or %0.3f minutes.\n" %(end - start, (end - start) / 60))

    print("\n  DIANA INFO:\t#### TIME OF EXECUTION ####\n")
    print ("  DIANA INFO:\tThe time of execution of the analysis is: %0.3f sec. or %0.3f minutes.\n" %(end - start, (end - start) / 60))


    print("\t\t--------------------------------------------\n")
    print("\t\tAnalysis finished. Thank you for using DIANA\n")
    print("\t\t--------------------------------------------\n")


if  __name__ == "__main__":
    main()

