import argparse
import cPickle
import hashlib
import mysql.connector
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug


def main():

    options = parse_user_arguments()
    design_clinical_trials_standard(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of the input drug",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-dc','--drug_combinations_file',dest='drug_combinations_file',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace/drug_combination_pairs.tsv'),
                        help = """Define the file of drug combinations from ClinicalTrials.gov""")
    parser.add_argument('-cr','--crossings_file',dest='crossings_file',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace/clinicaltrials_crossings_file.txt'),
                        help = """Define the file where the drug crossings to be explored will be written""")
    parser.add_argument('-min','--minimum_targets',dest='minimum_targets',action = 'store',default=3,
                        help = """Define the minimum number of targets that the drugs need to have to be considered in the experiment""")
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory are stored.""")
    parser.add_argument('-db','--database',dest='database',action = 'store',default='BIANA_JUN_2017',
                        help = """Define the database to use:
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
    parser.add_argument('-up','--unification',dest='unification_protocol',action = 'store',default='geneid_seqtax_v2',
                        help = """Define the unification protocol used in BIANA database (default is BIANA_JUN_2017)""")

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def design_clinical_trials_standard(options):
    """
    Designs the drug crossings to be explored from the ClinicalTrials drug combinations.
    """

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')


    #-------------------------------------------------#
    #   GET THE LIST OF DRUGS AND THEIR DRUBANK IDs   #
    #-------------------------------------------------#

    drug_name_to_drugbank_ids = {}
    drug_combinations = []

    biana_cnx = mysql.connector.connect(user=options.db_user, password=options.db_pass,
                                        host=options.db_host,
                                        database=options.database)

    with open(options.drug_combinations_file, 'r') as clinical_trials_fd:

        first_line = clinical_trials_fd.readline()
        fields_dict = obtain_header_fields(first_line, separator='\t')
        #drug1   drug2   studies

        for line in clinical_trials_fd:

            fields = line.strip().split('\t')
            drug1 = fields[ fields_dict['drug1'] ]
            drug2 = fields[ fields_dict['drug2'] ]
            studies = fields[ fields_dict['studies'] ].split('|')
            drug_combinations.append(frozenset([drug1, drug2]))

            for drug_name in [drug1, drug2]:
                if drug_name not in drug_name_to_drugbank_ids:
                    drugbank_ids = diana_drug.find_drugbank_id_from_name(biana_cnx, options.unification_protocol, drug_name)
                    if drugbank_ids:
                        drug_name_to_drugbank_ids[drug_name] = drugbank_ids

    biana_cnx.close()

    #------------------------------------------#
    #   OBTAIN THE CROSSINGS IN DRUGBANK IDs   #
    #------------------------------------------#

    drugbank2targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    drugbank2targets = cPickle.load(open(drugbank2targets_file))

    # Check how many drugs have enough targets
    drugs_with_targets = set()
    for drug_name in drug_name_to_drugbank_ids:
        drugbank_ids = drug_name_to_drugbank_ids[drug_name]
        for db_id in drugbank_ids:
            if db_id in drugbank2targets:
                if len(drugbank2targets[db_id]) >= int(options.minimum_targets):
                    drugs_with_targets.add(db_id)
    print('\nThe number of drugs with more than {} targets in the network is: {}\n'.format(options.minimum_targets, len(drugs_with_targets)))


    # Check that all drugs have PFAM
    pfam_pickle_file = os.path.join(toolbox_dir, 'drugbank_target_to_pfams.pcl')
    geneid2pfams = cPickle.load(open(pfam_pickle_file))

    drugs_with_pfams = set()
    for drug in drugs_with_targets:
        for target in drugbank2targets[drug]:
            if target in geneid2pfams:
                drugs_with_pfams.add(drug)
    print('\nThe number of drugs with at least one PFAM is: {}\n'.format(len(drugs_with_pfams)))


    # Check how many drugs have SMILES
    smiles_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_smiles.pcl')
    drug2smiles = cPickle.load(open(smiles_pickle_file))

    drugs_considered = set(drug2smiles.keys()) & drugs_with_pfams
    print('\nThe number of drugs with at least one SMILES is: {}\n'.format(len(drugs_considered)))


    drug_combinations_considered = []
    for drug_pair in drug_combinations:
        drug1, drug2 = drug_pair
        if drug1 in drug_name_to_drugbank_ids and drug2 in  drug_name_to_drugbank_ids:
            drugbank_ids1 = drug_name_to_drugbank_ids[drug1]
            drugbank_ids2 = drug_name_to_drugbank_ids[drug2]
            for drugbank_id1 in drugbank_ids1:
                for drugbank_id2 in drugbank_ids2:
                    if drugbank_id1 in drugs_considered and drugbank_id2 in drugs_considered:
                        drug_combinations_considered.append(frozenset([drugbank_id1, drugbank_id2]))
    print('\nThe number of drug combinations in ClinicalTrials.gov considered is: {}\n'.format(len(drug_combinations_considered)))


    #-----------------------------------------------------------#
    #   GENERATE THE FILE WITH THE DRUG COMBINATION CROSSINGS   #
    #-----------------------------------------------------------#

    dc_not_done = 0
    comparisons_dir = os.path.join(options.workspace, 'comparisons')
    all_drugs = set()
    ct_crossings = set()

    with open(options.crossings_file,"w") as crossings_file_fd:
        for drug_pair in drug_combinations_considered:
            drug1, drug2 = drug_pair
            all_drugs.add(drug1)
            all_drugs.add(drug2)
            crossing = '{}---{}'.format(drug1, drug2)
            crossings_file_fd.write("{}\n".format(crossing))
            ct_crossings.add(crossing)


            # Check if already done
            drug_name1 = drug1.lower()
            targets1 = list(drugbank2targets[drug1.upper()])
            drug_id1 = generate_drug_id(drug_name1, targets1)

            drug_name2 = drug2.lower()
            targets2 = list(drugbank2targets[drug2.upper()])
            drug_id2 = generate_drug_id(drug_name2, targets2)

            # Check if the results table file is already created
            comp_dir = os.path.join(comparisons_dir, '{}---{}'.format(drug_id1, drug_id2))
            results_table = os.path.join(comp_dir, 'results_table.tsv')
            if not fileExist(results_table):
                dc_not_done+=1

    print('The number of drug combinations not already done is: {}\n'.format(dc_not_done))


    #-------------------------------------------------------#
    #   GENERATE THE FILE WITH THE ALL POSSIBLE CROSSINGS   #
    #-------------------------------------------------------#

    list_of_drugs = list(all_drugs)
    list_of_drugs2 = list(all_drugs)

    crossings = set()
    dc = 0
    non_dc = 0
    n = 0
    pair2comb = {}

    while (n < len(list_of_drugs)):
        i = 0
        while (i < len(list_of_drugs2)):
            drug1 = list_of_drugs[n]
            drug2 = list_of_drugs2[i]
            if drug1 == drug2:
                i+=1
                continue

            ddi_name1 = "%s---%s"%(drug1, drug2)
            ddi_name2 = "%s---%s"%(drug2, drug1)
            #print("%s vs. %s" %(drug1, drug2))

            # We check that none of the two possible names are in the crossings set, and we add it (this is not necessary, but it is added as security)
            if ddi_name1 in ct_crossings:
                crossings.add(ddi_name1)
                pair2comb[ddi_name1] = 1
            elif ddi_name2 in ct_crossings:
                crossings.add(ddi_name2)
                pair2comb[ddi_name2] = 1
            else:
                if ddi_name1 not in crossings and ddi_name2 not in crossings:
                    crossings.add(ddi_name1)
                    pair2comb[ddi_name1] = 0

            i+=1

        # We remove the first drug from the second list, so that we do not have to repeat pairings
        list_of_drugs2.remove(drug1)
        n+=1

    print('There are {} possible ClinicalTrials crossings\n'.format(len(crossings)))
    checking = len(list_of_drugs) * (len(list_of_drugs) - 1) / 2
    if len(crossings) != checking:
        print("THERE IS AN ERROR IN THE ANALYSIS. The number of crossings does not correspond to the theoretical number")
        sys.exit(10)


    # Save the dict containing if the pairings are drug combinations or not
    dump_file = os.path.join(toolbox_dir, 'pair2comb_clinicaltrials.pcl')
    cPickle.dump(pair2comb, open(dump_file, 'w')) 

    all_crossings_file = os.path.join(options.workspace, 'positives_and_negatives_clinical_trials.txt')
    with open(all_crossings_file,"w") as crossings_file_fd:
        for pair in crossings:
            crossings_file_fd.write("{}\n".format(pair))

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


def check_file(file):
    """
    Checks if a file exists and if not, raises FileNotFound exception
    """
    if not fileExist(file):
        raise FileNotFound(file)


class FileNotFound(Exception):
    """
    Exception raised when a file is not found.
    """
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return 'The file {} has not been found.\nTherefore, the comparison cannot be performed. Please, check that all the profiles have been correctly generated.\n'.format(self.file)

def obtain_header_fields(first_line, separator='\t'):
    """ 
    Obtain a dictionary: "field_name" => "position" 
    """
    fields_dict = {}
    header_fields = first_line.strip().split(separator)
    for x in xrange(0, len(header_fields)):
        fields_dict[header_fields[x].lower()] = x
    return fields_dict

def generate_drug_id(drug_name, targets):
    """
    Generates an id for the drug containing the drug name followed
    by an ID made from the processing of all the targets of the drug
    """
    m = hashlib.md5()
    targets = [str(x) for x in targets] # Transform the targets to strings
    targets_str = ''.join(sorted(targets)) # Join them in one string
    m.update(targets_str) # Introduce the string in the hashlib instance
    targets_id = m.hexdigest()[:12] # Obtain a unique ID from the targets string. Only get the first 12 characters
    drug_str = ''.join(drug_name.split('\s')) # Obtain the drug name, and remove any space in the name
    unique_id = drug_str+'_'+targets_id # Add the drug name with the targets ID creating a unique ID for the drug
    return unique_id



if  __name__ == "__main__":
    main()


