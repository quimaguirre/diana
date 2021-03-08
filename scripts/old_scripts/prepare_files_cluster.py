import argparse
import cPickle
import ConfigParser
import mysql.connector
import time
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug



def main():

    options = parse_user_arguments()
    prepare_files_cluster(options)


def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Prepare the necessary files to run an analysis in the cluster",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = """" Input file with a protein-protein interaction network in SIF format.
                        If not introduced, the program will create a network of expansion using the targets as center and expanding as many neighbors
                        as specified in the parameter radius. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def prepare_files_cluster(options):
    """
    Creates the files necessary to run an analysis in the cluster.
    """

    # Start marker for time measure
    start = time.time()

    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = ConfigParser.ConfigParser()
    config.read(config_file)


    #----------------------#
    #   CREATE THE FILES   #
    #----------------------#

    dcdb2targets_file = os.path.join(toolbox_dir, 'dcdb2targets.pcl')
    drugbank2targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    pfam_dcdb_pickle_file = os.path.join(toolbox_dir, 'dcdb_target_to_pfams.pcl')
    pfam_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_target_to_pfams.pcl')
    smiles_dcdb_pickle_file = os.path.join(toolbox_dir, 'dcdb2smiles.pcl')
    smiles_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_smiles.pcl')


    if not fileExist(drugbank2targets_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( drugbank2targets_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        diana_drug.obtain_drugbank_to_targets(biana_cnx, config.get('BIANA', 'unification_protocol'), options.sif, drugbank2targets_file)
        biana_cnx.close()


    if not fileExist(dcdb2targets_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( dcdb2targets_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        diana_drug.obtain_dcdb_to_targets(biana_cnx, config.get('BIANA', 'unification_protocol'), options.sif, dcdb2targets_file)
        biana_cnx.close()


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


if  __name__ == "__main__":
    main()
