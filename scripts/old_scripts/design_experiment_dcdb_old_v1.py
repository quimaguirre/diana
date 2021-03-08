import argparse
import ConfigParser
import cPickle
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
    parser.add_argument('-min','--minimum_targets',dest='minimum_targets',action = 'store',default=3,
                        help = """Define the minimum number of targets that the drugs need to have to be considered in the experiment""")
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = """" Input file with the protein-protein interaction network in SIF format that will be used in the experiment. """)

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

    dcdb2targets_file = os.path.join(toolbox_dir, 'dcdb2targets.pcl')
    drugbank2targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    pfam_dcdb_pickle_file = os.path.join(toolbox_dir, 'dcdb_target_to_pfams.pcl')
    smiles_dcdb_pickle_file = os.path.join(toolbox_dir, 'dcdb2smiles.pcl')
    atc_dcdb_pickle_file = os.path.join(toolbox_dir, 'dcdb2atcs.pcl')
    sider_dcdb_pickle_file = os.path.join(toolbox_dir, 'dcdb2side_effects.pcl')
    int_to_drugs_file = os.path.join(toolbox_dir, 'drug_int_2_drugs.pcl')
    int_to_info_file = os.path.join(toolbox_dir, 'drug_int_2_info.pcl')
    dump_file = os.path.join(toolbox_dir, 'pair2comb.pcl')
    pubchem2drugbank_file = os.path.join(toolbox_dir, 'pubchem_to_drugbank.pcl')

    # Obtain the targets for all the DCDB drugs
    if not fileExist(dcdb2targets_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( dcdb2targets_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        diana_drug.obtain_dcdb_to_targets(biana_cnx, config.get('BIANA', 'unification_protocol'), options.sif, dcdb2targets_file)
        biana_cnx.close()

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
    if not fileExist(pfam_dcdb_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( pfam_dcdb_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))
        dcdb2targets = cPickle.load(open(dcdb2targets_file))
        all_targets = set()
        for drug in dcdb2targets:
            for target in dcdb2targets[drug]:
                all_targets.add(target)
        diana_drug.obtain_target_to_pfam(biana_cnx, config.get('BIANA', 'unification_protocol'), all_targets, pfam_dcdb_pickle_file)
        biana_cnx.close()

    # Obtain the SMILES of all the DCDB drugs
    if not fileExist(smiles_dcdb_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( smiles_dcdb_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drug2targets = cPickle.load(open(dcdb2targets_file))
        all_drugs = set(drug2targets.keys())
        diana_drug.obtain_drug_to_smiles(biana_cnx, config.get('BIANA', 'unification_protocol'), all_drugs, 'dcdb', smiles_dcdb_pickle_file)
        biana_cnx.close()

    # Obtain the ATCs of all the DCDB drugs
    if not fileExist(atc_dcdb_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( atc_dcdb_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drug2targets = cPickle.load(open(dcdb2targets_file))
        all_drugs = set(drug2targets.keys())
        diana_drug.obtain_drug_to_atcs(biana_cnx, config.get('BIANA', 'unification_protocol'), all_drugs, 'dcdb', atc_dcdb_pickle_file)
        biana_cnx.close()

    # Obtain the side effects of all the DCDB drugs
    if not fileExist(sider_dcdb_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( sider_dcdb_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drug2targets = cPickle.load(open(dcdb2targets_file))
        all_drugs = set(drug2targets.keys())
        diana_drug.obtain_drug_to_side_effects(biana_cnx, config.get('BIANA', 'unification_protocol'), all_drugs, 'dcdb', sider_dcdb_pickle_file)
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


    #---------------------------#
    #   GET THE LIST OF DRUGS   #
    #---------------------------#

    dcdb2targets = cPickle.load(open(dcdb2targets_file))

    # Get the drugs with at least the minimum number of targets
    drugs_with_targets = set()
    for dcdb in dcdb2targets:
        if len(dcdb2targets[dcdb]) >= int(options.minimum_targets):
            drugs_with_targets.add(dcdb)

    print('\nThe number of drugs included is: {}\n'.format(len(drugs_with_targets)))

    # Get the drugs that have at least one PFAM
    drugs_with_pfams = set()
    geneid2pfams = cPickle.load(open(pfam_dcdb_pickle_file))
    for drug in drugs_with_targets:
        for target in dcdb2targets[drug]:
            if target in geneid2pfams:
                drugs_with_pfams.add(drug)
    print('The number of drugs with at least one PFAM is: {}\n'.format(len(drugs_with_pfams)))

    # Check how many drugs have side effects
    drug2side_effects = cPickle.load(open(sider_dcdb_pickle_file))
    drugs_with_side_effects = set(drug2side_effects.keys()) & drugs_with_pfams
    print('The number of drugs with at least one SIDE EFFECT is: {}\n'.format(len(drugs_with_side_effects)))

    drugs_no_se = [x for x in drugs_with_targets if x not in set(drug2side_effects.keys())]
    print(drugs_no_se)

    # Get the drugs that have at least one ATC
    dcdb2atcs = cPickle.load(open(atc_dcdb_pickle_file))
    drugs_with_atcs = set(dcdb2atcs.keys()) & drugs_with_side_effects
    print('The number of drugs with at least one ATC is: {}\n'.format(len(drugs_with_atcs)))

    # Check how many drugs have SMILES
    drug2smiles = cPickle.load(open(smiles_dcdb_pickle_file))
    drugs_with_smiles = set(drug2smiles.keys()) & drugs_with_atcs
    print('The number of DCDB drugs with at least one SMILES is: {}\n'.format(len(drugs_with_smiles)))



    #-------------------------------------------------------#
    #   OBTAIN THE NAMES OF THE DCDB DRUGS IN DRUGBANK ID   #
    #-------------------------------------------------------#

    dcdb_to_drugbank_file = os.path.join(toolbox_dir, 'dcdb_to_drugbank.pcl')
    check_file(dcdb_to_drugbank_file)
    dcdb_to_drugbank = cPickle.load(open(dcdb_to_drugbank_file))
    drugs_considered = set()
    drugs_considered_drugbank = set()
    for dcdbid in drugs_with_smiles:
        if dcdbid in dcdb_to_drugbank:
            drugs_considered.add(dcdbid)
            for db in dcdb_to_drugbank[dcdbid]:
                drugs_considered_drugbank.add(db)

    print('The number of DCDB drugs with DrugBank ID considered is: {}\n'.format(len(drugs_considered)))
    print('The number of DrugBank IDs considered is: {}\n'.format(len(drugs_considered_drugbank)))



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

    drugbank2targets = cPickle.load(open(drugbank2targets_file))
    db_crossings = set()
    for crossing in crossings:
        drug1, drug2 = crossing.split('---')
        db_drugs1 = dcdb_to_drugbank[drug1]
        db_drugs2 = dcdb_to_drugbank[drug2]
        if 'DB01258' in db_drugs1:
            print(drug1)
            print(dcdb2targets[drug1])
            print(drugbank2targets['DB01258'])
        if 'DB01258' in db_drugs2:
            print(drug2)
            print(dcdb2targets[drug2])
            print(drugbank2targets['DB01258'])
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

    print('NUMBER OF PHARMACODYNAMICAL DRUG INTERACTIONS:\t\t{}\n'.format(dc))
    print('NUMBER OF NON-PHARMACODYNAMICAL DRUG INTERACTIONS:\t{}\n'.format(non_dc))

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

    original_crossings = set()
    new_crossings = set()
    pubchem_to_drugbank = cPickle.load(open(pubchem2drugbank_file))
    comparison_file = os.path.join(toolbox_dir, 'DrugPairIndex_original.txt')
    comparison_output_file = os.path.join(toolbox_dir, 'DrugPairIndex_drugbank.txt')
    with open(comparison_file, 'r') as comparison_fd, open(comparison_output_file, 'w') as comparison_output_fd:
        for line in comparison_fd:
            fields = line.strip().split('\t')
            pubchem1 = int(fields[0].split('CID')[1])
            pubchem2 = int(fields[1].split('CID')[1])
            original_crossings.add(frozenset([pubchem1, pubchem2]))
            if pubchem1 in pubchem_to_drugbank and pubchem2 in pubchem_to_drugbank:
                for db1 in pubchem_to_drugbank[pubchem1]:
                    for db2 in pubchem_to_drugbank[pubchem2]:
                        if db1 in drugs_considered_drugbank and db2 in drugs_considered_drugbank:
                            new_crossings.add(frozenset([db1, db2]))
                            comparison_output_fd.write('{}\t{}\n'.format( db1, db2 ))

    print('In the initial comparison there were: {} crossings'.format(len(original_crossings)))
    print('In the current comparison there are: {} crossings'.format(len(new_crossings)))


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



if  __name__ == "__main__":
    main()


