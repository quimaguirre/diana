import argparse
import ConfigParser
import cPickle
import csv
import mysql.connector
import time
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug


def main():

    options = parse_user_arguments()
    design_experiment_dcdb(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of the input drug",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-cr','--crossings_file',dest='crossings_file',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace/crossings_file.txt'),
                        help = """Define the file where the drug crossings to be explored will be written""")
    parser.add_argument('-min','--minimum_targets',dest='minimum_targets',action = 'store',default=1,
                        help = """Define the minimum number of targets that the drugs need to have to be considered in the experiment""")
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = """" Input file with the protein-protein interaction network in SIF format that will be used in the experiment. """)
    parser.add_argument('-se','--restrict_se',dest='restrict_se',action = 'store_true',
                        help = """" Restrict to drugs that have Side Effects / ATCs. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def design_experiment_dcdb(options):
    """
    Designs the drug crossings to be explored in the experiment of DCDB.
    """

    # Start marker for time measure
    start = time.time()

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')


    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)


    #-------------------------#
    #   CREATE PICKLE FILES   #
    #-------------------------#

    drugbank2targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    pfam_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_target_to_pfams.pcl')
    smiles_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_smiles.pcl')
    atc_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_atcs.pcl')
    sider_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_side_effects.pcl')
    int_to_drugs_file = os.path.join(toolbox_dir, 'drug_int_2_drugs.pcl')
    int_to_info_file = os.path.join(toolbox_dir, 'drug_int_2_info.pcl')
    dump_file = os.path.join(toolbox_dir, 'pair2comb.pcl')
    pubchem2drugbank_file = os.path.join(toolbox_dir, 'pubchem_to_drugbank.pcl')

    # Obtain the targets for all the DrugBank drugs
    if not fileExist(drugbank2targets_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( drugbank2targets_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        diana_drug.obtain_drugbank_to_targets(biana_cnx, config.get('BIANA', 'unification_protocol'), options.sif, drugbank2targets_file)
        biana_cnx.close()

    # Obtain all the PFAMs of the targets
    if not fileExist(pfam_drugbank_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( pfam_drugbank_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drugbank2targets = cPickle.load(open(drugbank2targets_file))
        all_targets = set()
        for drug in drugbank2targets:
            for target in drugbank2targets[drug]:
                all_targets.add(target)
        diana_drug.obtain_target_to_pfam(biana_cnx, config.get('BIANA', 'unification_protocol'), all_targets, pfam_drugbank_pickle_file)
        biana_cnx.close()

    # Obtain the SMILES of all the DrugBank drugs
    if not fileExist(smiles_drugbank_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( smiles_drugbank_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drug2targets = cPickle.load(open(drugbank2targets_file))
        all_drugs = set(drug2targets.keys())
        diana_drug.obtain_drug_to_smiles(biana_cnx, config.get('BIANA', 'unification_protocol'), all_drugs, 'drugbank', smiles_drugbank_pickle_file)
        biana_cnx.close()

    # Obtain the ATCs of all the DrugBank drugs
    if not fileExist(atc_drugbank_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( atc_drugbank_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drug2targets = cPickle.load(open(drugbank2targets_file))
        all_drugs = set(drug2targets.keys())
        diana_drug.obtain_drug_to_atcs(biana_cnx, config.get('BIANA', 'unification_protocol'), all_drugs, 'drugbank', atc_drugbank_pickle_file)
        biana_cnx.close()



    # Obtain the side effects of all the DrugBank drugs
    if not fileExist(sider_drugbank_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( sider_drugbank_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drug2targets = cPickle.load(open(drugbank2targets_file))
        all_drugs = set(drug2targets.keys())
        diana_drug.obtain_drug_to_side_effects(biana_cnx, config.get('BIANA', 'unification_protocol'), all_drugs, 'drugbank', sider_drugbank_pickle_file)
        biana_cnx.close()

    # Obtain the component drugs of the drug interactions in DCDB
    if not fileExist(int_to_drugs_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( int_to_drugs_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        diana_drug.obtain_drug_interaction_to_drugs(biana_cnx, int_to_drugs_file)
        biana_cnx.close()

    # Obtain the information of the drug interactions in DCDB
    if not fileExist(int_to_info_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( int_to_info_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        diana_drug.obtain_drug_interaction_to_info(biana_cnx, int_to_info_file)
        biana_cnx.close()


    #-------------------------------------------------------#
    #   OBTAIN THE NAMES OF THE DCDB DRUGS IN DRUGBANK ID   #
    #-------------------------------------------------------#

    dcdb_to_drugbank_file = os.path.join(toolbox_dir, 'dcdb_to_drugbank.pcl')
    check_file(dcdb_to_drugbank_file)
    dcdb_to_drugbank = cPickle.load(open(dcdb_to_drugbank_file))

    print('The number of DCDB IDs with DrugBank ID is: {}\n'.format(len(dcdb_to_drugbank.keys())))


    #------------------------------#
    #   GET THE DRUGS CONSIDERED   #
    #------------------------------#

    # Get the drugs with at least the minimum number of targets
    drugs_with_targets = set()
    drugbankIDs_with_targets = set()
    drugbank2targets = cPickle.load(open(drugbank2targets_file))
    for dcdbID in dcdb_to_drugbank:
        for drugbankID in dcdb_to_drugbank[dcdbID]:
            if drugbankID in drugbank2targets:
                if len(drugbank2targets[drugbankID]) >= int(options.minimum_targets):
                    drugs_with_targets.add(dcdbID)
                    drugbankIDs_with_targets.add(drugbankID)
    print('\nThe number of DCDB IDs with targets is: {}\n'.format(len(drugs_with_targets)))

    # Get the drugs that have at least one PFAM
    drugs_with_pfams = set()
    drugbankIDs_with_pfams = set()
    geneid2pfams = cPickle.load(open(pfam_drugbank_pickle_file))
    for dcdbID in drugs_with_targets:
        for drugbankID in dcdb_to_drugbank[dcdbID]:
            if drugbankID in drugbank2targets:
                for target in drugbank2targets[drugbankID]:
                    if target in geneid2pfams:
                        drugs_with_pfams.add(dcdbID)
                        drugbankIDs_with_pfams.add(drugbankID)
    print('The number of DCDB IDs with at least one PFAM is: {}\n'.format(len(drugs_with_pfams)))

    # Check how many drugs have side effects
    drugs_with_side_effects = set()
    drugbankIDs_with_side_effects = set()
    drug2side_effects = cPickle.load(open(sider_drugbank_pickle_file))
    for dcdbID in dcdb_to_drugbank:
        for drugbankID in dcdb_to_drugbank[dcdbID]:
            if drugbankID in drug2side_effects:
                drugs_with_side_effects.add(dcdbID)
                drugbankIDs_with_side_effects.add(drugbankID)
    print('The number of DCDB IDs with at least one SIDE EFFECT is: {}\n'.format(len(drugs_with_side_effects)))

    # Get the drugs that have at least one ATC
    drugs_with_atcs = set()
    drugbankIDs_with_atcs = set()
    drug2atcs = cPickle.load(open(atc_drugbank_pickle_file))
    for dcdbID in dcdb_to_drugbank:
        for drugbankID in dcdb_to_drugbank[dcdbID]:
            if drugbankID in drug2atcs:
                drugs_with_atcs.add(dcdbID)
                drugbankIDs_with_atcs.add(drugbankID)
    print('The number of DCDB IDs with at least one ATC is: {}\n'.format(len(drugs_with_atcs)))

    # Check how many drugs have SMILES
    drugs_with_smiles = set()
    drugbankIDs_with_smiles = set()
    drug2smiles = cPickle.load(open(smiles_drugbank_pickle_file))
    for dcdbID in dcdb_to_drugbank:
        for drugbankID in dcdb_to_drugbank[dcdbID]:
            if drugbankID in drug2smiles:
                drugs_with_smiles.add(dcdbID)
                drugbankIDs_with_smiles.add(drugbankID)
    print('The number of DCDB drugs with at least one SMILES is: {}\n'.format(len(drugs_with_smiles)))

    if options.restrict_se:

        # Get the drugs considered (with at least n targets, 1 PFAM, 1 SMILE, 1 ATC, 1 SE)
        drugs_considered = drugs_with_targets & drugs_with_pfams & drugs_with_smiles & drugs_with_atcs & drugs_with_side_effects
        print('\nThe number of DCDB IDs considered is: {}'.format(len(drugs_considered)))

        drugs_considered_drugbank = drugbankIDs_with_targets & drugbankIDs_with_pfams & drugbankIDs_with_smiles & drugbankIDs_with_atcs & drugbankIDs_with_side_effects
        print('\nThe number of DrugBank IDs considered is: {}\n'.format(len(drugs_considered_drugbank)))

    else:

        # Get the drugs considered (with at least n targets, 1 PFAM, 1 SMILE)
        drugs_considered = drugs_with_targets & drugs_with_pfams & drugs_with_smiles
        print('\nThe number of DCDB IDs considered is: {}'.format(len(drugs_considered)))

        drugs_considered_drugbank = drugbankIDs_with_targets & drugbankIDs_with_pfams & drugbankIDs_with_smiles
        print('\nThe number of DrugBank IDs considered is: {}\n'.format(len(drugs_considered_drugbank)))


    #-------------------------------------------------#
    #   DEFINE ALL POSSIBLE CROSSINGS BETWEEN PAIRS   #
    #-------------------------------------------------#

    # We need to copy the list using the list() method because if not, when you modify one of the lists, the other gets modified as well
    # This can also be done with copy.copy() method or copy.deepcopy() if the list contains objects and you want to copy them as well
    # More info: http://stackoverflow.com/questions/2612802/how-to-clone-or-copy-a-list 
    list_of_drugs = list(drugs_considered)
    list_of_drugs2 = list(drugs_considered)

    drug_int_2_drugs = cPickle.load(open(int_to_drugs_file))
    drug_int_2_info = cPickle.load(open(int_to_info_file))

    crossings = set()
    pair2comb = {}
    dc = 0
    non_dc = 0
    n = 0

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
            if ddi_name1 not in crossings and ddi_name2 not in crossings:
                crossings.add(ddi_name1)

            i+=1

        # We remove the first drug from the second list, so that we do not have to repeat pairings
        list_of_drugs2.remove(drug1)
        n+=1

    print('There are {} possible DCDB crossings\n'.format(len(crossings)))
    checking = len(list_of_drugs) * (len(list_of_drugs) - 1) / 2
    if len(crossings) != checking:
        print("THERE IS AN ERROR IN THE ANALYSIS. The number of crossings does not correspond to the theoretical number")
        sys.exit(10)
    #print(crossings)


    #--------------------------------#
    #   TRANSLATE DCDB TO DRUGBANK   #
    #--------------------------------#

    db_crossings = set()
    for crossing in crossings:
        drug1, drug2 = crossing.split('---')
        db_drugs1 = [ x for x in dcdb_to_drugbank[drug1] if x in drugbankIDs_with_targets ]
        db_drugs2 = [ x for x in dcdb_to_drugbank[drug2] if x in drugbankIDs_with_targets ]
        for db_drug1 in db_drugs1:
            for db_drug2 in db_drugs2:
                db_crossing1 = '{}---{}'.format(db_drug1, db_drug2)
                db_crossing2 = '{}---{}'.format(db_drug1, db_drug2)
                if db_crossing1 not in db_crossings and db_crossing2 not in db_crossings:
                    db_crossings.add(db_crossing1)

                    pair2comb[db_crossing1] = 0 # We will introduce 0 if it is not drug interaction
                    non_dc+=1
                    for drug_int in drug_int_2_drugs:
                        if drug1 in drug_int_2_drugs[drug_int] and drug2 in drug_int_2_drugs[drug_int]:
                            if drug_int_2_info[drug_int]['type'] == 'pharmacodynamical':
                                pair2comb[db_crossing1] = 1 # We will introduce 1 if it is a pharmacodynamical drug interaction
                                dc+=1
                                non_dc-=1
                            else:
                                pair2comb[db_crossing1] = 0 # We will introduce 0 if it is not a pharmacodynamical drug interaction
                            break

    print('There are {} possible DrugBank crossings\n'.format(len(db_crossings)))

    print('Number of pharmacodynamical drug interactions:\t{}\n'.format(dc))
    print('Number of non-pharmacodynamical drug interactions:\t{}\n'.format(non_dc))

    # Save the dict containing if the pairings are drug combinations or not
    cPickle.dump(pair2comb, open(dump_file, 'w')) 


    #------------------------------------------#
    #   GENERATE THE FILE WITH THE CROSSINGS   #
    #------------------------------------------#

    with open(options.crossings_file,"w") as crossings_file_fd:
        for pair in db_crossings:
            crossings_file_fd.write("{}\n".format(pair))


    #------------------------------------------------------------#
    #   GENERATE THE FILE WITH THE CROSSINGS OF THE COMPARISON   #
    #------------------------------------------------------------#

    if not fileExist(pubchem2drugbank_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( pubchem2drugbank_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        diana_drug.obtain_pubchem_to_drugbank(biana_cnx, config.get('BIANA', 'unification_protocol'), pubchem2drugbank_file)
        biana_cnx.close()

    pubchem_to_drugbank = cPickle.load(open(pubchem2drugbank_file))
    comparison_file = os.path.join(toolbox_dir, 'comparison_other_methods/DrugPairIndex_original.txt')
    ATC_file = os.path.join(toolbox_dir, 'comparison_other_methods/ATC_Profile.csv')
    DDI_file = os.path.join(toolbox_dir, 'comparison_other_methods/DDI_Profile.csv')
    DDI_interact_file = os.path.join(toolbox_dir, 'comparison_other_methods/DDI_Interact.csv')
    DTI_file = os.path.join(toolbox_dir, 'comparison_other_methods/DTI_Profile.csv')
    SE_file = os.path.join(toolbox_dir, 'comparison_other_methods/FeatureMat_SE.csv')
    comparison_output_file = os.path.join(toolbox_dir, 'comparison_other_methods/DrugPairIndex_drugbank.txt')
    comparison_crossings_file = os.path.join(toolbox_dir, 'comparison_other_methods/comparison_crossings.txt')

    ids_drugs_removed_ATC = read_matrix_of_drugs(ATC_file)
    ids_drugs_removed_DTI = read_matrix_of_drugs(DTI_file)
    ids_drugs_removed_DDI = read_matrix_of_drugs(DDI_file)
    ids_drugs_removed_DDI_Int = read_matrix_of_drugs(DDI_interact_file)
    ids_drugs_removed = ids_drugs_removed_ATC | ids_drugs_removed_DTI | ids_drugs_removed_DDI | ids_drugs_removed_DDI_Int
    print(ids_drugs_removed)
    print('Drugs removed in ATC: {}'.format(len(ids_drugs_removed_ATC)))
    print('Drugs removed in DTI: {}'.format(len(ids_drugs_removed_DTI)))
    print('Drugs removed in DDI: {}'.format(len(ids_drugs_removed_DDI)))
    print('Drugs removed in DDI (interact): {}'.format(len(ids_drugs_removed_DDI_Int)))
    print('Drugs removed in total: {}'.format(len(ids_drugs_removed)))

    original_crossings = set()
    filtered_crossings = set()
    new_crossings = set()
    drugs_comparison = set()
    drugs_not_found = set()
    with open(comparison_file, 'r') as comparison_fd, open(comparison_output_file, 'w') as comparison_output_fd:
        for line in comparison_fd:
            fields = line.strip().split('\t')
            pubchem1 = int(fields[0].split('CID')[1])
            pubchem2 = int(fields[1].split('CID')[1])
            drug_id1 = int(fields[2])
            drug_id2 = int(fields[3])
            original_crossings.add(frozenset([pubchem1, pubchem2]))
            if drug_id1 in ids_drugs_removed or drug_id2 in ids_drugs_removed:
                continue
            filtered_crossings.add(frozenset([pubchem1, pubchem2]))
            if pubchem1 in pubchem_to_drugbank and pubchem2 in pubchem_to_drugbank:
                for db1 in pubchem_to_drugbank[pubchem1]:
                    for db2 in pubchem_to_drugbank[pubchem2]:
                        if db1 in drugs_considered_drugbank and db2 in drugs_considered_drugbank:
                            new_crossings.add(frozenset([db1, db2]))
                            comparison_output_fd.write('{}\t{}\n'.format( db1, db2 ))
                            drugs_comparison.add(db1)
                            drugs_comparison.add(db2)
                        else:
                            if db1 not in drugs_considered_drugbank:
                                drugs_not_found.add(db1)
                            if db2 not in drugs_considered_drugbank:
                                drugs_not_found.add(db2)

    print('Drugs to compare: {}'.format(len(drugs_comparison)))

    possible_comparison_crossings = set()
    with open(comparison_crossings_file, 'w') as comparison_crossings_fd:
        for crossing in db_crossings:
            db1, db2 = crossing.split('---')
            if db1 in drugs_comparison and db2 in drugs_comparison:
                comparison_crossings_fd.write('{}---{}\n'.format( db1, db2 ))
                possible_comparison_crossings.add(crossing)

    # Checking that the number of crossings for the comparison with other methods is correct
    checking = len(drugs_comparison) * (len(drugs_comparison) - 1) / 2
    if len(possible_comparison_crossings) != checking:
        print(checking)
        print(len(possible_comparison_crossings))
        print("THERE IS AN ERROR IN THE ANALYSIS. The number of crossings does not correspond to the theoretical number")

    print('In the initial comparison there were: {} drug combinations'.format(len(original_crossings)))
    print('We filtered to: {} drug combinations'.format(len(filtered_crossings)))
    print('In the current comparison there are: {} drug combinations'.format(len(new_crossings)))
    print('And there are {} possible comparison crossings\n'.format(len(possible_comparison_crossings)))
    print('Drugs in BMC dataset not found:')
    drug2targets = cPickle.load(open(drugbank2targets_file))
    drug2smiles = cPickle.load(open(smiles_drugbank_pickle_file))
    drug2atcs = cPickle.load(open(atc_drugbank_pickle_file))
    drug2side_effects = cPickle.load(open(sider_drugbank_pickle_file))
    for drug in drugs_not_found:
        print(drug)
        if drug in drug2targets:
            print(drug2targets[drug])
        else:
            print('No targets')
        if drug in drug2smiles:
            print(drug2smiles[drug])
        else:
            print('No SMILES')
        if drug in drug2atcs:
            print(drug2atcs[drug])
        else:
            print('No ATCs')
        if drug in drug2side_effects:
            print(drug2side_effects[drug])
        else:
            print('No Side Effects')


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


def read_matrix_of_drugs(matrix_file):
    """
    Read a matrix file of the comparison.
    """
    drug_id = 1
    ids_drugs_removed = set()
    with open(matrix_file, 'r') as matrix_fd:
        data = csv.reader(matrix_fd, delimiter=',')
        for row in data:
            if not '1' in row:
                ids_drugs_removed.add(drug_id)
            drug_id+=1
    return ids_drugs_removed


if  __name__ == "__main__":
    main()


