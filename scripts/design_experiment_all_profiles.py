import argparse
import ConfigParser
import cPickle
import mysql.connector
import ntpath
import time
import sys, os, re

from context import diana
import diana.classes.drug as diana_drug


def main():

    options = parse_user_arguments()
    design_experiment_all_profiles(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of the input drug",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace/drugbank_drugs_selected.txt'),
                        help = """File containing all the drugs in DrugBank that will be considered to generate all the profiles""")
    parser.add_argument('-on','--output_names_file',dest='output_names_file',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace/drugbank_names_selected.txt'),
                        help = """File containing all the drug names in DrugBank that will be considered to generate all the profiles""")
    parser.add_argument('-min','--minimum_targets',dest='minimum_targets',action = 'store',default=1,
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

def design_experiment_all_profiles(options):
    """
    Obtains the targets of the DrugBankIDs with more than the minimum number of specified targets in the network.
    Obtains also the PFAMs and SMILES.
    Outputs a file with all the drugs having at least the minimum number of targets, SMILES and PFAM.
    Prepares also the correspondency between drugbank drug and dcdb drug, and the correspondency between
    drug and its ID (created from the drug name and an ID obtained from the targets used).
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


    #----------------------#
    #   CREATE THE FILES   #
    #----------------------#

    drugbank2targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    pfam_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_target_to_pfams.pcl')
    smiles_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_smiles.pcl')
    atc_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_atcs.pcl')
    sider_drugbank_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_side_effects.pcl')
    bio_processes_file = os.path.join(toolbox_dir, 'target_to_bio_processes.pcl')
    pathways_file = os.path.join(toolbox_dir, 'target_to_pathways.pcl')
    drugbank_to_diana_id_file = os.path.join(toolbox_dir, 'drugbank_to_diana_id.pcl')
    diana_id_to_drugbank_file = os.path.join(toolbox_dir, 'diana_id_to_drugbank.pcl')
    dcdb_to_drugbank_file = os.path.join(toolbox_dir, 'dcdb_to_drugbank.pcl')
    drugbank_to_dcdb_file = os.path.join(toolbox_dir, 'drugbank_to_dcdb.pcl')
    drugbank_to_names_pickle_file = os.path.join(toolbox_dir, 'drugbank_to_names.pcl')


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



    # Obtain the side effects of all the DrugBank drugs
    if not fileExist(drugbank_to_names_pickle_file):

        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( drugbank_to_names_pickle_file ))

        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))

        drug2targets = cPickle.load(open(drugbank2targets_file))
        all_drugs = set(drug2targets.keys())
        diana_drug.obtain_drugbank_to_names(biana_cnx, all_drugs, drugbank_to_names_pickle_file)
        biana_cnx.close()



    # Obtain the biological processes for the targets of the drugs
    if not fileExist(bio_processes_file):
        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( bio_processes_file ))
        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))
        drugbank2targets = cPickle.load(open(drugbank2targets_file))
        all_targets = set()
        for drugbankid in drugbank2targets:
            for target in drugbank2targets[drugbankid]:
                all_targets.add(target)
        diana_drug.obtain_target_to_bio_processes(biana_cnx, config.get('BIANA', 'unification_protocol'), all_targets, bio_processes_file)
        biana_cnx.close()



    # Obtain the pathways for the targets of the drugs
    if not fileExist(pathways_file):
        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( pathways_file ))
        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))
        drugbank2targets = cPickle.load(open(drugbank2targets_file))
        all_targets = set()
        for drugbankid in drugbank2targets:
            for target in drugbank2targets[drugbankid]:
                all_targets.add(target)
        diana_drug.obtain_target_to_pathways(biana_cnx, config.get('BIANA', 'unification_protocol'), all_targets, pathways_file)
        biana_cnx.close()



    # Check the number of drugs considered
    drugbank2targets = cPickle.load(open(drugbank2targets_file))
    drugs_with_targets = set()
    for drugbankid in drugbank2targets:
        if len(drugbank2targets[drugbankid]) >= int(options.minimum_targets):
            drugs_with_targets.add(drugbankid)
    print('\nThe number of drugs with more than {} targets in the network is: {}\n'.format(options.minimum_targets, len(drugs_with_targets)))


    # Check that all drugs have PFAM
    drugs_with_pfams = set()
    geneid2pfams = cPickle.load(open(pfam_drugbank_pickle_file))
    for drug in drugs_with_targets:
        for target in drugbank2targets[drug]:
            if target in geneid2pfams:
                drugs_with_pfams.add(drug)
    print('\nThe number of drugs with at least one PFAM is: {}\n'.format(len(drugs_with_pfams)))


    # Check how many drugs have SMILES
    drug2smiles = cPickle.load(open(smiles_drugbank_pickle_file))
    drugs_with_smiles = set(drug2smiles.keys()) & drugs_with_targets
    print('\nThe number of drugs with at least one SMILES is: {}\n'.format(len(drugs_with_smiles)))


    # Check how many drugs have side effects
    drug2side_effects = cPickle.load(open(sider_drugbank_pickle_file))
    drugs_with_side_effects = set(drug2side_effects.keys()) & drugs_with_targets
    print('\nThe number of drugs with at least one SIDE EFFECT is: {}\n'.format(len(drugs_with_side_effects)))


    # Check how many drugs have ATCs
    drugbank2atcs = cPickle.load(open(atc_drugbank_pickle_file))
    drugs_with_atcs = set(drugbank2atcs.keys()) & drugs_with_targets
    print('\nThe number of drugs with at least one ATC is: {}\n'.format(len(drugs_with_atcs)))


    # Get the drugs considered (with at least n targets, 1 PFAM and 1 SMILE)
    drugs_considered = drugs_with_targets & drugs_with_pfams & drugs_with_smiles
    print('\nThe number of drugs considered is: {}\n'.format(len(drugs_considered)))


    # Output a file with all the drugs considered
    with open(options.output_file, 'w') as output_fd:
        for drug in drugs_considered:
            output_fd.write('{}\n'.format(drug))


    # Output also a file with all the drug names considered
    drugbank_to_names = cPickle.load(open(drugbank_to_names_pickle_file))
    with open(options.output_names_file, 'w') as output_fd:
        for drug in drugs_considered:
            if drug in drugbank_to_names:
                for name in drugbank_to_names[drug]:
                    output_fd.write('{}\n'.format(name))


    # Obtain the IDs of the drugs generated with the targets
    if not fileExist(drugbank_to_diana_id_file):
        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( drugbank_to_diana_id_file ))
        drugbank_to_diana_id = {}
        diana_id_to_drugbank = {}
        for drug in drugs_considered:
            targets = list(drugbank2targets[drug.upper()])
            network_filename = ntpath.basename(options.sif)
            drug_id = diana_drug.generate_drug_id(drug.lower(), targets, network_filename)
            drugbank_to_diana_id[drug] = drug_id
            diana_id_to_drugbank[drug_id] = drug

        cPickle.dump(drugbank_to_diana_id, open(drugbank_to_diana_id_file, 'w'))
        cPickle.dump(diana_id_to_drugbank, open(diana_id_to_drugbank_file, 'w'))


    # Obtain the correspondency of DCDB ID to DrugBank ID.
    if not fileExist(dcdb_to_drugbank_file):
        print( "  DIANA INFO:\tCreating pickle file {}.\n".format( dcdb_to_drugbank_file ))
        biana_cnx = mysql.connector.connect(user=config.get('BIANA', 'user'), 
                                            password=config.get('BIANA', 'password'),
                                            host=config.get('BIANA', 'host'),
                                            database=config.get('BIANA', 'database'))
        diana_drug.obtain_dcdb_to_drugbank(biana_cnx, config.get('BIANA', 'unification_protocol'), dcdb_to_drugbank_file)
        biana_cnx.close()

    # Create as well a drugbank_to_dcdb file!
    if not fileExist(drugbank_to_dcdb_file):
        dcdb_to_drugbank_file = os.path.join(toolbox_dir, 'dcdb_to_drugbank.pcl')
        dcdb_to_drugbank = cPickle.load(open(dcdb_to_drugbank_file))
        drugbank_to_dcdb = {}
        for dcdbid in dcdb_to_drugbank:
            for db in dcdb_to_drugbank[dcdbid]:
                drugbank_to_dcdb.setdefault(db, set())
                drugbank_to_dcdb[db].add(dcdbid)
        cPickle.dump(drugbank_to_dcdb, open(drugbank_to_dcdb_file, 'w'))
        print(drugbank_to_dcdb)


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

