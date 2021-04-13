import argparse
import configparser
import networkx as nx
import pandas as pd
import pymysql
import time
import sys, os, re


def main():

    options = parse_user_arguments()
    create_drug_mappings(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    python /home/quim/PHD/Projects/DIANA/diana/scripts/create_drug_mappings.py -n /home/quim/PHD/Projects/DIANA/diana/data/network_cheng.txt -d /home/quim/PHD/Projects/DIANA/diana/data -o /home/quim/PHD/Projects/DIANA/diana/outputs
    """

    parser = argparse.ArgumentParser(
        description = "Create drug mappings using BIANA database",
        epilog      = "@oliva's lab 2020")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store',
                        help = """" Input file with the protein-protein interaction network in SIF format that will be used in the experiment. """)
    parser.add_argument('-d','--data_dir',dest='data_dir',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'data'),
                        help = """Directory where the additional data is stored""")
    parser.add_argument('-o','--output_dir',dest='output_dir',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'mappings'),
                        help = """Directory where the output files will be stored""")

    options=parser.parse_args()

    return options



#################
#################
# MAIN FUNCTION #
#################
#################

def create_drug_mappings(options):
    """
    Creates the mappings of the drugs from the unification protocol of BIANA.
    """

    # Start marker for time measure
    start = time.time()

    # Get the inputs
    network_file = options.network_file
    data_dir = options.data_dir
    create_directory(data_dir)
    output_dir = options.output_dir
    create_directory(output_dir)

    # Get the outputs (all drugs / targets)
    output_drug_user_entities_file = os.path.join(output_dir, 'all_drug_user_entity_ids.txt')
    output_drug_mapping_file = os.path.join(output_dir, 'all_drug_mappings.txt')
    output_atc_file = os.path.join(output_dir, 'all_drug_atc.txt')
    output_smiles_file = os.path.join(output_dir, 'all_drug_smiles.txt')
    output_side_effects_file = os.path.join(output_dir, 'all_drug_side_effects.txt')
    output_target_mapping_file = os.path.join(output_dir, 'all_target_mappings.txt')
    output_drug_target_interactions_file = os.path.join(output_dir, 'all_drug_target_interactions.txt')
    output_drug_target_interactions_network_file = os.path.join(output_dir, 'all_drug_target_interactions_in_network.txt')
    output_drug_to_targets_file = os.path.join(output_dir, 'all_drug_to_targets.txt')

    # Get the outputs (filtered drugs / targets)
    filtered_drug_mapping_file = os.path.join(output_dir, 'filtered_drug_mappings.txt')
    filtered_atc_file = os.path.join(output_dir, 'filtered_drug_atc.txt')
    filtered_smiles_file = os.path.join(output_dir, 'filtered_drug_smiles.txt')
    filtered_side_effects_file = os.path.join(output_dir, 'filtered_drug_side_effects.txt')
    filtered_target_mapping_file = os.path.join(output_dir, 'filtered_target_mappings.txt')
    filtered_drug_target_interactions_file = os.path.join(output_dir, 'filtered_drug_target_interactions.txt')
    filtered_drug_target_interactions_network_file = os.path.join(output_dir, 'filtered_drug_target_interactions_in_network.txt')
    filtered_drug_to_targets_file = os.path.join(output_dir, 'filtered_drug_to_targets.txt')

    # Get the outputs (filtered with DrugbankID / GeneID)
    drugbank_drug_mapping_file = os.path.join(output_dir, 'drugbank_drug_mappings.txt')
    drugbank_atc_file = os.path.join(output_dir, 'drugbank_drug_atc.txt')
    drugbank_smiles_file = os.path.join(output_dir, 'drugbank_drug_smiles.txt')
    drugbank_side_effects_file = os.path.join(output_dir, 'drugbank_drug_side_effects.txt')
    geneid_target_mapping_file = os.path.join(output_dir, 'geneid_target_mappings.txt')
    drugbank_geneid_drug_target_interactions_file = os.path.join(output_dir, 'drugbank_geneid_drug_target_interactions.txt')
    drugbank_geneid_drug_target_interactions_network_file = os.path.join(output_dir, 'drugbank_geneid_drug_target_interactions_in_network.txt')
    drugbank_geneid_drug_to_targets_file = os.path.join(output_dir, 'drugbank_geneid_drug_to_targets.txt')
    diana_drugbank_benchmark_file = os.path.join(output_dir, 'DIANA_drugbank_benchmark.txt')
    diana_drugbank_benchmark_names_file = os.path.join(output_dir, 'DIANA_drugbank_benchmark_names.json')

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')

    # Check network
    if not fileExist(network_file):
        print('Network file does not exist: {}'.format(network_file))
        sys.exit(10)
    else:
        network_geneid = nx.Graph()
        with open(network_file, 'r') as net_fd:
            for line in net_fd:
                fields = line.strip().split('\t')
                if len(fields) == 2:
                    network_geneid.add_edge(fields[0], fields[1])
                elif len(fields) == 3:
                    network_geneid.add_edge(fields[0], fields[2])
                else:
                    print(fields)
                    print(len(fields))
                    print('Network in incorrect format!')
                    sys.exit(10)
        print('\nTHE INPUT NETWORK HAS {} EDGES AND {} NODES (GENEIDS)\n'.format(network_geneid.number_of_edges(), network_geneid.number_of_nodes()))


    #--------------------------------------#
    #   GET INFORMATION FROM CONFIG FILE   #
    #--------------------------------------#

    # Read the config file
    config_file = os.path.join(main_path, 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_file)


    #------------------------------#
    #   INICIATE BIANA CONNEXION   #
    #------------------------------#
    cnx = pymysql.connect(user=config.get('BIANA', 'user'), 
                          password=config.get('BIANA', 'password'),
                          host=config.get('BIANA', 'host'),
                          database=config.get('BIANA', 'database'))

    unification_table = return_unification_protocol_table(cnx, config.get('BIANA', 'unification_protocol').lower())

    cursor = cnx.cursor()


    #------------------------------#
    #   GET DRUG USER ENTITY IDS   #
    #------------------------------#

    if not fileExist(output_drug_user_entities_file):

        # Get drug user entity IDs
        user_entity_ids = get_drug_user_entity_ids(cursor, unification_table)

        # Output drug user entity IDs
        with open(output_drug_user_entities_file, 'w') as out_fd:
            for user_entity_id in user_entity_ids:
                out_fd.write('{}\n'.format(user_entity_id))
    
    else:

        drug_user_entities_df = pd.read_csv(output_drug_user_entities_file, sep='\t', index_col=None, header=None)
        user_entity_ids = set(drug_user_entities_df[0])


    #-----------------------#
    #   GET DRUG MAPPINGS   #
    #-----------------------#

    if not fileExist(output_drug_mapping_file):


        # Get drug names
        ueid_to_drugtype_to_drugnames = get_drug_names_of_user_entity_ids(cursor, unification_table, user_entity_ids)

        # Get drugbankids
        ueid_to_drugbankids = get_drugbankids_of_user_entity_ids(cursor, unification_table, user_entity_ids)

        # Get pubchembcompounds
        ueid_to_pubchemcompounds = get_pubchemcompounds_of_user_entity_ids(cursor, unification_table, user_entity_ids)

        # Get chemblids
        ueid_to_chemblids = get_chemblids_of_user_entity_ids(cursor, unification_table, user_entity_ids)

        # Output drug mappings
        output_drug_mappings(output_drug_mapping_file, user_entity_ids, ueid_to_drugtype_to_drugnames, ueid_to_drugbankids, ueid_to_pubchemcompounds, ueid_to_chemblids)


    #-------------------#
    #   GET DRUG ATCS   #
    #-------------------#

    if not fileExist(output_atc_file):

        # Get ATCs
        ueid_to_atc_to_databases = get_atcs_of_user_entity_ids(cursor, unification_table, user_entity_ids)

        # Output drug ATCs
        print('\nWRITING OUTPUT ATC FILE...\n')

        with open(output_atc_file, 'w') as out_fd:
            out_fd.write('#user_entity_id\tatc\tdatabases\n')
            for user_entity_id in sorted(user_entity_ids):
                if user_entity_id in ueid_to_atc_to_databases:
                    for atc in ueid_to_atc_to_databases[user_entity_id]:
                        databases = '; '.join(sorted(ueid_to_atc_to_databases[user_entity_id][atc]))
                        out_fd.write('{}\t{}\t{}\n'.format(user_entity_id, atc, databases))


    #---------------------#
    #   GET DRUG SMILES   #
    #---------------------#

    if not fileExist(output_smiles_file):

        # Get SMILES
        ueid_to_smiles_to_databases = get_smiles_of_user_entity_ids(cursor, unification_table, user_entity_ids)

        # Output SMILES
        print('\nWRITING OUTPUT SMILES FILE...\n')

        with open(output_smiles_file, 'w') as out_fd:
            out_fd.write('#user_entity_id\tsmiles\tdatabases\n')
            for user_entity_id in sorted(user_entity_ids):
                if user_entity_id in ueid_to_smiles_to_databases:
                    for smiles in ueid_to_smiles_to_databases[user_entity_id]:
                        databases = '; '.join(sorted(ueid_to_smiles_to_databases[user_entity_id][smiles]))
                        out_fd.write('{}\t{}\t{}\n'.format(user_entity_id, smiles, databases))


    #-----------------------------#
    #   GET SIDE EFFECTS SMILES   #
    #-----------------------------#

    if not fileExist(output_side_effects_file):

        # Get pubchem to side effects
        meddra_all_se_file = os.path.join(data_dir, 'meddra_all_se.tsv')
        pubchem_to_umls, umls_to_name = get_side_effects_from_sider(meddra_all_se_file)

        pubchemcompound_to_ueids = {}
        with open(output_drug_mapping_file, 'r') as inp_fd:
            first_line = inp_fd.readline()
            for line in inp_fd:
                ue_drug, type_identifier, identifier, type_name = line.strip().split('\t')
                if type_identifier == 'pubchemcompound':
                    pubchemcompound_to_ueids.setdefault(ue_drug, set()).add(identifier.upper())

        # Output side effects
        with open(output_side_effects_file, 'w') as out_fd:
            out_fd.write('#user_entity_id\tumls_id\tumls_name\n')
            for pubchemcompound in pubchemcompound_to_ueids:
                for user_entity_id in pubchemcompound_to_ueids[pubchemcompound]:
                    if pubchemcompound in pubchem_to_umls:
                        for umls_id in pubchem_to_umls[pubchemcompound]:
                            umls_name = umls_to_name[umls_id]
                            out_fd.write('{}\t{}\t{}\n'.format(user_entity_id, umls_id, umls_name))


    #----------------------------------#
    #   GET DRUG TARGET INTERACTIONS   #
    #----------------------------------#

    if not fileExist(output_drug_target_interactions_file):

        # Get drug targets from DrugBank
        drugbank_drug_target_interactions, drugbank_targets, drugbank_drug_to_targets, drugbank_drug_metab_interactions, drugbank_metab, drugbank_drug_to_metab = get_drug_targets_from_drugbank(cursor, unification_table, user_entity_ids)

        # Get drug targets from DrugCentral
        drugcentral_drug_target_interactions, drugcentral_targets, drugcentral_drug_to_targets = get_drug_targets_from_drugcentral(cursor, unification_table, user_entity_ids)

        # Get drug targets from DGIdb
        dgidb_drug_target_interactions, dgidb_targets, dgidb_drug_to_targets = get_drug_targets_from_dgidb(cursor, unification_table, user_entity_ids)

        # Get drug targets from ChEMBL
        chembl_drug_target_interactions, chembl_targets, chembl_drug_to_targets = get_drug_targets_from_chembl(cursor, unification_table, user_entity_ids)

        # Get drug targets from TTD
        ttd_drug_target_interactions, ttd_targets, ttd_drug_to_targets = get_drug_targets_from_ttd(cursor, unification_table, user_entity_ids)

        # Join interactions
        drug_target_interactions_with_metab = drugbank_drug_target_interactions | drugcentral_drug_target_interactions | dgidb_drug_target_interactions | chembl_drug_target_interactions | ttd_drug_target_interactions

        # Join targets
        targets_with_metab = drugbank_targets | drugcentral_targets | dgidb_targets | chembl_targets | ttd_targets

        print('NUMBER OF TARGET USER ENTITIES without removing metab: {}'.format(len(targets_with_metab)))
        print('NUMBER OF DRUG-TARGET INTERACTIONS without removing metab: {}'.format(len(drug_target_interactions_with_metab)))

        # Remove metab targets
        targets = targets_with_metab - drugbank_metab

        # Remove metab interactions
        drug_target_interactions = set()
        for interaction in drug_target_interactions_with_metab:
            (ueid_drug, ueid_target) = interaction
            if ueid_target not in drugbank_metab:
                drug_target_interactions.add(interaction)

        print('NUMBER OF TARGET USER ENTITIES after removing metab: {}'.format(len(targets)))
        print('NUMBER OF DRUG-TARGET INTERACTIONS after removing metab: {}'.format(len(drug_target_interactions)))

        # Output drug-target interactions file
        print('\nWRITING OUTPUT DRUG-TARGET INTERACTIONS FILE...\n')

        with open(output_drug_target_interactions_file, 'w') as out_fd:
            out_fd.write('#drug_user_entity_id\ttarget_user_entity_id\tdatabases\n')
            for interaction in drug_target_interactions:
                ue_drug, ue_target = interaction
                sources = []
                if interaction in drugbank_drug_target_interactions:
                    sources.append('drugbank')
                if interaction in dgidb_drug_target_interactions:
                    sources.append('dgidb')
                if interaction in drugcentral_drug_target_interactions:
                    sources.append('drugcentral')
                if interaction in chembl_drug_target_interactions:
                    sources.append('chembl')
                if interaction in ttd_drug_target_interactions:
                    sources.append('ttd')
                sources = '; '.join(sorted(sources))
                out_fd.write('{}\t{}\t{}\n'.format(ue_drug, ue_target, sources))

    # Read the drug target interactions
    drug_target_interactions = set()
    targets = set()
    drug_to_targets = {}
    with open(output_drug_target_interactions_file, 'r') as inp_fd:
        first_line = inp_fd.readline()
        for line in inp_fd:
            ue_drug, ue_target, sources = line.strip().split('\t')
            drug_to_targets.setdefault(ue_drug, set()).add(ue_target)
            interaction = (ue_drug, ue_target)
            drug_target_interactions.add(interaction)
            targets.add(ue_target)


    #-------------------------#
    #   GET TARGET MAPPINGS   #
    #-------------------------#

    if not fileExist(output_target_mapping_file):

        # Get gene symbols
        ueid_to_genesymbol_to_types = get_genesymbols_of_user_entity_ids(cursor, unification_table, targets)

        # Get gene IDs 
        ueid_to_geneid_to_types = get_geneids_of_user_entity_ids(cursor, unification_table, targets)

        # Get uniprot accessions
        ueid_to_uniprotaccession_to_types = get_uniprotaccessions_of_user_entity_ids(cursor, unification_table, targets)

        # Get uniprot entries
        ueid_to_uniprotentry_to_types = get_uniprotentries_of_user_entity_ids(cursor, unification_table, targets)

        # Get PFAM identifiers
        ueid_to_pfams = get_pfam_of_user_entity_ids(cursor, unification_table, targets)

        # Get taxonomy identifiers
        ueid_to_taxids = get_taxid_of_user_entity_ids(cursor, unification_table, targets)

        # Output target mappings
        output_target_mappings(output_target_mapping_file, targets, ueid_to_genesymbol_to_types, ueid_to_geneid_to_types, ueid_to_uniprotaccession_to_types, ueid_to_uniprotentry_to_types, ueid_to_pfams, ueid_to_taxids)


    #---------------------------------------------------#
    #   GET DRUG TO TARGETS & TARGETS IN NETWORK FILE   #
    #---------------------------------------------------#

    if not fileExist(output_drug_to_targets_file) or not fileExist(output_drug_target_interactions_network_file):

        # Output drug to targets file
        print('\nWRITING OUTPUT DRUG TO TARGETS FILE...\n')

        # Get Gene ID to UEID mapping
        geneid_to_ueids = {}
        with open(output_target_mapping_file, 'r') as inp_fd:
            for line in inp_fd:
                ue_target, type_identifier, identifier, type_name = line.strip().split('\t')
                if type_identifier == 'geneid':
                    geneid_to_ueids.setdefault(identifier, set()).add(ue_target)

        # Get network in UEID
        network_ueid = nx.Graph()
        for u,v in network_geneid.edges():
            if u in geneid_to_ueids and v in geneid_to_ueids:
                for u_ueid in geneid_to_ueids[u]:
                    for v_ueid in geneid_to_ueids[v]:
                        network_ueid.add_edge(str(u_ueid), str(v_ueid))

        with open(output_drug_to_targets_file, 'w') as out_fd:
            out_fd.write('#drug_user_entity_id\ttarget_user_entity_ids\ttargets_in_network\n')
            for user_entity_id in user_entity_ids:
                if (str(user_entity_id) in drug_to_targets):
                    user_entity_id_targets = [str(target) for target in drug_to_targets[str(user_entity_id)]]
                    ueid_targets_in_network = [str(target) for target in user_entity_id_targets if target in network_ueid.nodes()]
                    out_fd.write('{}\t{}\t{}\n'.format(str(user_entity_id), '; '.join(sorted(user_entity_id_targets)), '; '.join(sorted(ueid_targets_in_network))))

        with open(output_drug_target_interactions_file, 'r') as inp_fd, open(output_drug_target_interactions_network_file, 'w') as out_fd:
            out_fd.write('#drug_user_entity_id\ttarget_user_entity_id\tdatabases\n')
            for line in inp_fd:
                ue_drug, ue_target, sources = line.strip().split('\t')
                if ue_target in network_ueid.nodes():
                    out_fd.write('{}\t{}\t{}\n'.format(ue_drug, ue_target, sources))



    #------------------------------#
    #   FILTER DRUGS AND TARGETS   #
    #------------------------------#

    # We select the drugs that have a DrugBankID assigned
    # We select the targets that have an Entrez GeneID assigned and are human (taxid = 9606)

    print('\nFILTERING DRUG FILES...\n')

    # Filter drug mappings
    if not fileExist(filtered_drug_mapping_file):
        drug_mapping_df = pd.read_csv(output_drug_mapping_file, sep='\t', index_col=None)
        filtered_drugs = set(drug_mapping_df[drug_mapping_df['type_identifier'] == 'drugbankid']['#user_entity_id'])
        drug_mapping_filt_df = drug_mapping_df[drug_mapping_df['#user_entity_id'].isin(filtered_drugs)]
        drug_mapping_filt_df.to_csv(filtered_drug_mapping_file, sep='\t', index=False)

        # Check drugs with multiple drugbank IDs
        ueid_to_drugbankids = {}
        ueids_with_multiple_drugbankids = set()
        for index, row in drug_mapping_filt_df.iterrows():
            ueid = row['#user_entity_id']
            type_identifier = row['type_identifier']
            identifier = row['identifier']
            if type_identifier == 'drugbankid':
                ueid_to_drugbankids.setdefault(ueid, set()).add(identifier)
                if len(ueid_to_drugbankids[ueid]) > 1:
                    #print('Multiple DrugBankIDs for biana ID {}: {}'.format(ueid, ueid_to_drugbankids[ueid]))
                    ueids_with_multiple_drugbankids.add(ueid)
        print('DRUG USER ENTITIES WITH MULTIPLE DRUGBANK IDS: {}'.format(ueids_with_multiple_drugbankids))

    else:
        drug_mapping_filt_df = pd.read_csv(filtered_drug_mapping_file, sep='\t', index_col=None)
        filtered_drugs = set(drug_mapping_filt_df['#user_entity_id'])
    print('NUMBER OF DRUG USER ENTITIES ORIGINALLY: {}'.format(len(user_entity_ids)))
    print('NUMBER OF DRUG USER ENTITIES FILTERED: {}'.format(len(filtered_drugs)))

    # Filter drug ATCs
    if not fileExist(filtered_atc_file):
        drug_atc_df = pd.read_csv(output_atc_file, sep='\t', index_col=None)
        drug_atc_filt_df = drug_atc_df[drug_atc_df['#user_entity_id'].isin(filtered_drugs)]
        drug_atc_filt_df.to_csv(filtered_atc_file, sep='\t', index=False)

    # Filter drug SMILES
    if not fileExist(filtered_smiles_file):
        drug_smiles_df = pd.read_csv(output_smiles_file, sep='\t', index_col=None)
        drug_smiles_filt_df = drug_smiles_df[drug_smiles_df['#user_entity_id'].isin(filtered_drugs)]
        drug_smiles_filt_df.to_csv(filtered_smiles_file, sep='\t', index=False)

    # Filter drug side effects
    if not fileExist(filtered_side_effects_file):
        drug_side_effects_df = pd.read_csv(output_side_effects_file, sep='\t', index_col=None)
        drug_side_effects_filt_df = drug_side_effects_df[drug_side_effects_df['#user_entity_id'].isin(filtered_drugs)]
        drug_side_effects_filt_df.to_csv(filtered_side_effects_file, sep='\t', index=False)

    print('\nFILTERING TARGET FILES...\n')

    # Filter target mappings (by the ones that have at least one geneID and are from human)
    if not fileExist(filtered_target_mapping_file):
        target_mapping_df = pd.read_csv(output_target_mapping_file, sep='\t', index_col=None)
        original_targets = set(target_mapping_df['#user_entity_id'])
        targets_with_geneid = set(target_mapping_df[target_mapping_df['type_identifier'] == 'geneid']['#user_entity_id'])
        target_mapping_with_geneid_df = target_mapping_df[target_mapping_df['#user_entity_id'].isin(targets_with_geneid)]
        filtered_targets = set(target_mapping_with_geneid_df[ (target_mapping_with_geneid_df['type_identifier'] == 'taxid') & (target_mapping_with_geneid_df['identifier'] == '9606') ]['#user_entity_id'])
        target_mapping_filt_df = target_mapping_with_geneid_df[target_mapping_with_geneid_df['#user_entity_id'].isin(filtered_targets)]
        target_mapping_filt_df.to_csv(filtered_target_mapping_file, sep='\t', index=False)

        # Check targets with multiple geneids
        target_to_geneids = {}
        targets_with_multiple_geneids = set()
        for index, row in target_mapping_filt_df.iterrows():
            ueid = row['#user_entity_id']
            type_identifier = row['type_identifier']
            identifier = row['identifier']
            if type_identifier == 'geneid':
                target_to_geneids.setdefault(ueid, set()).add(identifier)
                if len(target_to_geneids[ueid]) > 1:
                    #print('Multiple Gene IDs for biana ID {}: {}'.format(ueid, target_to_geneids[ueid]))
                    targets_with_multiple_geneids.add(ueid)
        print('TARGET USER ENTITIES WITH MULTIPLE GENE IDS: {}'.format(targets_with_multiple_geneids))

    else:
        target_mapping_df = pd.read_csv(output_target_mapping_file, sep='\t', index_col=None)
        original_targets = set(target_mapping_df['#user_entity_id'])
        target_mapping_filt_df = pd.read_csv(filtered_target_mapping_file, sep='\t', index_col=None)
        filtered_targets = set(target_mapping_filt_df['#user_entity_id'])
    print('NUMBER OF TARGET USER ENTITIES ORIGINALLY: {}'.format(len(original_targets)))
    print('NUMBER OF TARGET USER ENTITIES FILTERED: {}'.format(len(filtered_targets)))

    print('\nFILTERING DRUG TARGET FILES...\n')

    # Filter drug-target interactions file
    if not fileExist(filtered_drug_target_interactions_file) or not fileExist(filtered_drug_target_interactions_network_file):
        # Drug-target interactions in general
        drug_target_interactions_df = pd.read_csv(output_drug_target_interactions_file, sep='\t', index_col=None)
        drug_target_interactions_filt_df = drug_target_interactions_df[ (drug_target_interactions_df['#drug_user_entity_id'].isin(filtered_drugs)) & (drug_target_interactions_df['target_user_entity_id'].isin(filtered_targets)) ]
        drug_target_interactions_filt_df.to_csv(filtered_drug_target_interactions_file, sep='\t', index=False)
        # Drug-target interactions in network
        drug_target_interactions_network_df = pd.read_csv(output_drug_target_interactions_network_file, sep='\t', index_col=None)
        drug_target_interactions_network_filt_df = drug_target_interactions_network_df[ (drug_target_interactions_network_df['#drug_user_entity_id'].isin(filtered_drugs)) & (drug_target_interactions_network_df['target_user_entity_id'].isin(filtered_targets)) ]
        drug_target_interactions_network_filt_df.to_csv(filtered_drug_target_interactions_network_file, sep='\t', index=False)
    else:
        drug_target_interactions_df = pd.read_csv(output_drug_target_interactions_file, sep='\t', index_col=None)
        drug_target_interactions_filt_df = pd.read_csv(filtered_drug_target_interactions_file, sep='\t', index_col=None)
        drug_target_interactions_network_df = pd.read_csv(output_drug_target_interactions_network_file, sep='\t', index_col=None)
        drug_target_interactions_network_filt_df = pd.read_csv(filtered_drug_target_interactions_network_file, sep='\t', index_col=None)
    print('NUMBER OF DRUG-TARGET INTERACTIONS ORIGINALLY: {}'.format(len(drug_target_interactions_df.index)))
    print('NUMBER OF DRUG-TARGET INTERACTIONS FILTERED: {}'.format(len(drug_target_interactions_filt_df.index)))
    print('NUMBER OF DRUG-TARGET INTERACTIONS IN NETWORK ORIGINALLY: {}'.format(len(drug_target_interactions_network_df.index)))
    print('NUMBER OF DRUG-TARGET INTERACTIONS IN NETWORK FILTERED: {}'.format(len(drug_target_interactions_network_filt_df.index)))

    # Filter drug to targets file
    if not fileExist(filtered_drug_to_targets_file):
        drug_to_targets_df = pd.read_csv(output_drug_to_targets_file, sep='\t', index_col=None)
        drug_to_targets_filt_df = pd.DataFrame(columns=['#drug_user_entity_id', 'target_user_entity_ids', 'targets_in_network'])
        filtered_targets = [str(target) for target in filtered_targets]
        filtered_drugs = [str(drug) for drug in filtered_drugs]
        for index, row in drug_to_targets_df.iterrows():
            drug_ueid = str(row['#drug_user_entity_id'])
            target_ueids = [ target for target in str(row['target_user_entity_ids']).split('; ') if target in filtered_targets ]
            targets_in_network = [ target for target in str(row['targets_in_network']).split('; ') if target in filtered_targets ]
            if( (drug_ueid in filtered_drugs) and (len(target_ueids) > 0) ):
                df2 = pd.DataFrame([[drug_ueid, '; '.join(sorted(target_ueids)), '; '.join(sorted(targets_in_network))]], columns=['#drug_user_entity_id', 'target_user_entity_ids', 'targets_in_network'])
                # Add the information to the main data frame
                drug_to_targets_filt_df = drug_to_targets_filt_df.append(df2)
        drug_to_targets_filt_df.to_csv(filtered_drug_to_targets_file, sep='\t', index=False)

    drug_to_targets_df = pd.read_csv(output_drug_to_targets_file, sep='\t', index_col=None)
    drug_to_targets_in_network = set( drug_to_targets_df[ drug_to_targets_df['targets_in_network'].notnull() ]['#drug_user_entity_id'] )
    drug_to_targets_filt_df = pd.read_csv(filtered_drug_to_targets_file, sep='\t', index_col=None)
    drug_to_targets_filt_in_network = set( drug_to_targets_filt_df[ drug_to_targets_filt_df['targets_in_network'].notnull() ]['#drug_user_entity_id'] )
    print('NUMBER OF DRUGS WITH AT LEAST ONE TARGET ORIGINALLY: {}'.format(len(drug_to_targets_df.index)))
    print('NUMBER OF DRUGS WITH AT LEAST ONE TARGET FILTERED: {}'.format(len(drug_to_targets_filt_df.index)))
    print('NUMBER OF DRUGS WITH AT LEAST ONE TARGET IN THE NETWORK ORIGINALLY: {}'.format(len(drug_to_targets_in_network)))
    print('NUMBER OF DRUGS WITH AT LEAST ONE TARGET IN THE NETWORK FILTERED: {}'.format(len(drug_to_targets_filt_in_network)))


    #---------------------------------------------------#
    #   MAKE MAPPING FILES WITH DRUGBANKID AND GENEID   #
    #---------------------------------------------------#

    # Get UEID to DrugBank mapping
    ueid_to_drugbankids = {}
    ueid_to_drugnames_unique = {}
    with open(filtered_drug_mapping_file, 'r') as inp_fd:
        first_line = inp_fd.readline()
        for line in inp_fd:
            ue_drug, type_identifier, identifier, type_name = line.strip().split('\t')
            if type_identifier == 'drugbankid':
                ueid_to_drugbankids.setdefault(ue_drug, set()).add(identifier.upper())
            if type_identifier == 'name' and type_name == 'unique':
                ueid_to_drugnames_unique.setdefault(ue_drug, set()).add(identifier.lower())
    with open(filtered_drug_mapping_file, 'r') as inp_fd, open(drugbank_drug_mapping_file, 'w') as out_fd:
        first_line = inp_fd.readline()
        out_fd.write('#drugbankid\ttype_identifier\tidentifier\ttype_name\n')
        for line in inp_fd:
            ue_drug, type_identifier, identifier, type_name = line.strip().split('\t')
            if type_identifier != 'drugbankid':
                if ue_drug in ueid_to_drugbankids:
                    for drugbankid in ueid_to_drugbankids[ue_drug]:
                        out_fd.write('{}\t{}\t{}\t{}\n'.format(drugbankid, type_identifier, identifier, type_name))

    # Mapping for drug ATCs
    with open(filtered_atc_file, 'r') as inp_fd, open(drugbank_atc_file, 'w') as out_fd:
        first_line = inp_fd.readline()
        out_fd.write('#drugbankid\tatc\tdatabases\n')
        for line in inp_fd:
            ue_drug, atc, databases = line.strip().split('\t')
            if ue_drug in ueid_to_drugbankids:
                for drugbankid in ueid_to_drugbankids[ue_drug]:
                    out_fd.write('{}\t{}\t{}\n'.format(drugbankid, atc, databases))

    # Mapping for drug SMILES
    with open(filtered_smiles_file, 'r') as inp_fd, open(drugbank_smiles_file, 'w') as out_fd:
        first_line = inp_fd.readline()
        out_fd.write('#drugbankid\tsmiles\tdatabases\n')
        for line in inp_fd:
            ue_drug, smiles, databases = line.strip().split('\t')
            if ue_drug in ueid_to_drugbankids:
                for drugbankid in ueid_to_drugbankids[ue_drug]:
                    out_fd.write('{}\t{}\t{}\n'.format(drugbankid, smiles, databases))

    # Mapping for drug side effects
    with open(filtered_side_effects_file, 'r') as inp_fd, open(drugbank_side_effects_file, 'w') as out_fd:
        first_line = inp_fd.readline()
        out_fd.write('#drugbankid\tumls_id\tumls_name\n')
        for line in inp_fd:
            ue_drug, umls_id, umls_name = line.strip().split('\t')
            if ue_drug in ueid_to_drugbankids:
                for drugbankid in ueid_to_drugbankids[ue_drug]:
                    out_fd.write('{}\t{}\t{}\n'.format(drugbankid, umls_id, umls_name))

    # Get UEID to GeneID mapping
    ueid_to_geneids = {}
    with open(filtered_target_mapping_file, 'r') as inp_fd:
        first_line = inp_fd.readline()
        for line in inp_fd:
            ue_target, type_identifier, identifier, type_name = line.strip().split('\t')
            if type_identifier == 'geneid':
                ueid_to_geneids.setdefault(ue_target, set()).add(identifier.upper())
    with open(filtered_target_mapping_file, 'r') as inp_fd, open(geneid_target_mapping_file, 'w') as out_fd:
        first_line = inp_fd.readline()
        out_fd.write('#geneid\ttype_identifier\tidentifier\ttype_name\n')
        for line in inp_fd:
            ue_target, type_identifier, identifier, type_name = line.strip().split('\t')
            if type_identifier != 'geneid':
                if ue_target in ueid_to_geneids:
                    for geneid in ueid_to_geneids[ue_target]:
                        out_fd.write('{}\t{}\t{}\t{}\n'.format(geneid, type_identifier, identifier, type_name))

    # Mapping for drug-target interactions
    with open(filtered_drug_target_interactions_file, 'r') as inp_fd, open(drugbank_geneid_drug_target_interactions_file, 'w') as out_fd:
        first_line = inp_fd.readline()
        out_fd.write('#drugbankid\tgeneid\tdatabases\n')
        for line in inp_fd:
            ue_drug, ue_target, databases = line.strip().split('\t')
            if ue_drug in ueid_to_drugbankids and ue_target in ueid_to_geneids:
                for drugbankid in ueid_to_drugbankids[ue_drug]:
                    for geneid in ueid_to_geneids[ue_target]:
                        out_fd.write('{}\t{}\t{}\n'.format(drugbankid, geneid, databases))

    # Mapping for drug-target interactions in network
    with open(filtered_drug_target_interactions_network_file, 'r') as inp_fd, open(drugbank_geneid_drug_target_interactions_network_file, 'w') as out_fd:
        first_line = inp_fd.readline()
        out_fd.write('#drugbankid\tgeneid\tdatabases\n')
        for line in inp_fd:
            ue_drug, ue_target, databases = line.strip().split('\t')
            if ue_drug in ueid_to_drugbankids and ue_target in ueid_to_geneids:
                for drugbankid in ueid_to_drugbankids[ue_drug]:
                    for geneid in ueid_to_geneids[ue_target]:
                        out_fd.write('{}\t{}\t{}\n'.format(drugbankid, geneid, databases))

    # Mapping for drug to targets
    drugbank_benchmark_names = set()
    with open(filtered_drug_to_targets_file, 'r') as inp_fd, open(drugbank_geneid_drug_to_targets_file, 'w') as out_fd, open(diana_drugbank_benchmark_file, 'w') as out2_fd:
        first_line = inp_fd.readline()
        out_fd.write('#drugbankid\tgeneids\tgeneids_in_network\n')
        for line in inp_fd:
            fields = line.strip().split('\t')
            if len(fields) == 3:
                ue_drug, ue_targets, ue_targets_network = fields
            elif len(fields) == 2:
                ue_drug, ue_targets = fields
                ue_targets_network = ''
            if ue_drug in ueid_to_drugbankids:
                drugbankids = ueid_to_drugbankids[ue_drug]
                if len(drugbankids) > 1:
                    print('Multiple DrugBank IDs for User Entity ID {}: {}'.format(ue_drug, drugbankids))
                    sys.exit(10)
                else:
                    drugbankid = list(drugbankids)[0]
                geneids = []
                geneids_in_network = []
                for ue_target in ue_targets.split('; '):
                    if ue_target in ueid_to_geneids:
                        for geneid in ueid_to_geneids[ue_target]:
                            geneids.append(geneid)
                            if ue_target in ue_targets_network.split('; '):
                                geneids_in_network.append(geneid)
                out_fd.write('{}\t{}\t{}\n'.format(drugbankid, '; '.join(geneids), '; '.join(geneids_in_network)))
                if len(geneids_in_network) > 0:
                    # Write the drugbank benchmark using drugbank ids to create the profiles
                    out2_fd.write('{}\n'.format(drugbankid))
                    # Keep the drug names for the web server autocomplete
                    if ue_drug in ueid_to_drugnames_unique:
                        for drugname in ueid_to_drugnames_unique[ue_drug]:
                            drugbank_benchmark_names.add(drugname)

    # Write the drugbank benchmark using drug names (in json for the web server autocomplete)
    drugbank_benchmark_names = sorted(list(drugbank_benchmark_names))
    with open(diana_drugbank_benchmark_names_file, 'w') as out_fd:
        out_fd.write('[\n')
        for x in range(len(drugbank_benchmark_names)):
            drugname = drugbank_benchmark_names[x]
            if x < (len(drugbank_benchmark_names) - 1):
                out_fd.write('"{}",\n'.format(drugname))
            else:
                out_fd.write('"{}"\n]'.format(drugname))


    #-----------------------------------------#
    #   COMPARE TARGETS WITH OTHER DATASETS   #
    #-----------------------------------------#

    venn_dir = os.path.join(output_dir, 'venn_drugtargets')
    venn_network_dir = os.path.join(output_dir, 'venn_drugtargets_network')
    create_directory(venn_dir)
    create_directory(venn_network_dir)
    cheng_drug_targets_file = os.path.join(data_dir, 'Cheng_NatCom19_DrugTargetInteractions.txt')
    cheng_output_file = os.path.join(venn_dir, 'cheng_drugtargets.txt')
    cheng_output_file_network = os.path.join(venn_network_dir, 'cheng_drugtargets_network.txt')
    biana_output_file = os.path.join(venn_dir, 'biana_drugtargets.txt')
    biana_output_file_network = os.path.join(venn_network_dir, 'biana_drugtargets_network.txt')
    ivenn_output_file = os.path.join(venn_dir, 'all_drugtargets.ivenn')
    ivenn_output_file_network = os.path.join(venn_network_dir, 'all_drugtargets_network.ivenn')
    if fileExist(cheng_drug_targets_file):

        print('\nCOMPARING DRUG-TARGET DATASETS...\n')

        # Read Cheng drug target interactions and write outputs
        cheng_drugtargets = set()
        cheng_drugtargets_network = set()
        with open(cheng_drug_targets_file, 'r') as inp_fd, open(cheng_output_file, 'w') as out1, open(cheng_output_file_network, 'w') as out2:
            first_line = inp_fd.readline()
            for line in inp_fd:
                drugbankid, geneid = line.strip().split('\t')
                interaction = '{}---{}'.format(drugbankid.upper(), geneid.upper())
                cheng_drugtargets.add(interaction)
                out1.write('{}\n'.format(interaction))
                if geneid in network_geneid.nodes():
                    cheng_drugtargets_network.add(interaction)
                    out2.write('{}\n'.format(interaction))
            print('NUMBER OF DRUG-TARGET INTERACTIONS IN CHENG: {}'.format(len(cheng_drugtargets)))
            print('NUMBER OF DRUG-TARGET INTERACTIONS IN NETWORK IN CHENG: {}'.format(len(cheng_drugtargets_network)))

        # Get UEID to GeneID mapping
        ueid_to_geneids = {}
        with open(output_target_mapping_file, 'r') as inp_fd:
            for line in inp_fd:
                ue_target, type_identifier, identifier, type_name = line.strip().split('\t')
                if type_identifier == 'geneid':
                    ueid_to_geneids.setdefault(ue_target, set()).add(identifier.upper())

        # Get UEID to DrugBankID mapping
        ueid_to_drugbankids = {}
        with open(output_drug_mapping_file, 'r') as inp_fd:
            for line in inp_fd:
                ue_drug, type_identifier, identifier, type_name = line.strip().split('\t')
                if type_identifier == 'drugbankid':
                    ueid_to_drugbankids.setdefault(ue_drug, set()).add(identifier.upper())

        # Read interactions file
        source_to_drug_to_targets = {}
        source_to_interactions = {}
        source_to_interactions_network = {}
        with open(output_drug_target_interactions_file, 'r') as inp_fd:
            for line in inp_fd:
                ue_drug, ue_target, sources = line.strip().split('\t')
                if ue_drug in ueid_to_drugbankids and ue_target in ueid_to_geneids:
                    for drugbankid in ueid_to_drugbankids[ue_drug]:
                        for geneid in ueid_to_geneids[ue_target]:
                            for source in sources.split('; '):
                                source_to_drug_to_targets.setdefault(source, {})
                                source_to_drug_to_targets[source].setdefault(drugbankid, set()).add(geneid)
                                interaction = '{}---{}'.format(drugbankid, geneid)
                                source_to_interactions.setdefault(source, set()).add(interaction)
                                if geneid in network_geneid.nodes():
                                    source_to_interactions_network.setdefault(source, set()).add(interaction)
        with open(biana_output_file, 'w') as out1, open(biana_output_file_network, 'w') as out2:
            for source in ['drugbank', 'drugcentral', 'dgidb', 'chembl', 'ttd']:
                drugtargets_output_file = os.path.join(venn_dir, '{}_drugtargets.txt'.format(source))
                drugtargets_network_output_file = os.path.join(venn_network_dir, '{}_drugtargets_network.txt'.format(source))
                with open(drugtargets_output_file, 'w') as out3, open(drugtargets_network_output_file, 'w') as out4:
                    for interaction in source_to_interactions[source]:
                        out1.write('{}\n'.format(interaction))
                        out3.write('{}\n'.format(interaction))
                    for interaction in source_to_interactions_network[source]:
                        out2.write('{}\n'.format(interaction))
                        out4.write('{}\n'.format(interaction))
        with open(ivenn_output_file, 'w') as out_ivenn, open(ivenn_output_file_network, 'w') as out_ivenn_net:
            for source in ['drugbank', 'drugcentral', 'dgidb', 'chembl', 'ttd']:
                out_ivenn.write('{}:{};\n'.format(source, ','.join(source_to_interactions[source])))
                out_ivenn_net.write('{}:{};\n'.format(source, ','.join(source_to_interactions_network[source])))
            out_ivenn.write('cheng:{};\n'.format(','.join(cheng_drugtargets)))
            out_ivenn_net.write('cheng:{};\n'.format(','.join(cheng_drugtargets_network)))

    # End marker for time
    end = time.time()
    print('\nDIANA INFO: TIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))

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


def return_unification_protocol_table(biana_cnx, unification_protocol):
    """
    Returns the table that contains the Unification Protocol
    introduced as query
    """

    query = (''' SELECT unificationProtocolID FROM userEntityUnificationProtocol
                 WHERE description = %s ''')

    cursor = biana_cnx.cursor() # Start cursor to MySQL
    cursor.execute(query, (unification_protocol,))
    up_ids = []
    for items in cursor:
        for up in items:
            up_ids.append(up)
    up_id = up_ids[0]
    up_table = 'userEntityUnification_protocol_'+str(up_id)
    cursor.close()

    return up_table


def get_drug_user_entity_ids(cursor, unification_table):
    """
    Get the BIANA user entity ids from drugs
    """
    query1 = ("""SELECT U.userEntityID
                 FROM externalEntity E, {} U 
                 WHERE U.externalEntityID = E.externalEntityID AND E.type = "drug"
              """.format(unification_table))

    print('\nRETRIEVING USER ENTITIES ASSOCIATED TO DRUGS...\n')

    user_entity_ids = set()

    cursor.execute(query1)
    for row in cursor:
        for ueid in row:
            #print(ueid)
            user_entity_ids.add(ueid)

    print('\nNUMBER OF USER ENTITIES ASSOCIATED WITH DRUGS: {}\n'.format(len(user_entity_ids)))

    return user_entity_ids


def get_drug_names_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get the drug names using their BIANA user entity ids
    """
    query_name = ("""SELECT N.value, N.type
                     FROM externalEntityName N, {} U 
                     WHERE U.externalEntityID = N.externalEntityID AND U.userEntityID = %s
                  """.format(unification_table))

    print('\nRETRIEVING DRUG NAMES ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_drugtype_to_drugnames = {}

    for ueid in user_entity_ids:
        cursor.execute(query_name, (ueid,))
        for row in cursor:
            drug_name, drug_type = row
            #print(ueid, drug_name, drug_type)
            ueid_to_drugtype_to_drugnames.setdefault(ueid, {})
            ueid_to_drugtype_to_drugnames[ueid].setdefault(drug_type.lower(), set()).add(drug_name.lower())

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH DRUG NAMES: {}'.format(len(ueid_to_drugtype_to_drugnames)))

    return ueid_to_drugtype_to_drugnames


def get_drugbankids_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get the drugbankids using their BIANA user entity ids
    """
    query_drugbankid = ("""SELECT DB.value
                           FROM externalEntityDrugBankID DB, {} U
                           WHERE U.externalEntityID = DB.externalEntityID AND U.userEntityID = %s
                        """.format(unification_table))

    print('\nRETRIEVING DRUGBANK IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_drugbankids = {}

    for ueid in user_entity_ids:
        cursor.execute(query_drugbankid, (ueid,))
        for row in cursor:
            for drugbankid in row:
                #print(ueid, drugbankid)
                ueid_to_drugbankids.setdefault(ueid, set()).add(drugbankid.upper())

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH DRUGBANK IDS: {}'.format(len(ueid_to_drugbankids)))

    return ueid_to_drugbankids


def get_chemblids_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get the chemblids using their BIANA user entity ids
    """
    query_chemblid = ("""SELECT CH.value
                           FROM externalEntityCHEMBL CH, {} U
                           WHERE U.externalEntityID = CH.externalEntityID AND U.userEntityID = %s
                        """.format(unification_table))

    print('\nRETRIEVING CHEMBL IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_chemblids = {}

    for ueid in user_entity_ids:
        cursor.execute(query_chemblid, (ueid,))
        for row in cursor:
            for chemblid in row:
                #print(ueid, chemblid)
                ueid_to_chemblids.setdefault(ueid, set()).add(chemblid.upper())

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH CHEMBL IDS: {}'.format(len(ueid_to_chemblids)))

    return ueid_to_chemblids


def get_pubchemcompounds_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get the pubchemcompounds using their BIANA user entity ids
    """
    query_pubchemcompound = ("""SELECT PC.value
                           FROM externalEntityPubChemCompound PC, {} U
                           WHERE U.externalEntityID = PC.externalEntityID AND U.userEntityID = %s
                        """.format(unification_table))

    print('\nRETRIEVING PUBCHEMCOMPOUND IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_pubchemcompounds = {}

    for ueid in user_entity_ids:
        cursor.execute(query_pubchemcompound, (ueid,))
        for row in cursor:
            for pubchemcompound in row:
                #print(ueid, pubchemcompound)
                ueid_to_pubchemcompounds.setdefault(ueid, set()).add(pubchemcompound)

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH PUBCHEMCOMPOUND IDS: {}'.format(len(ueid_to_pubchemcompounds)))

    return ueid_to_pubchemcompounds


def get_atcs_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get the atcs using their BIANA user entity ids
    """
    # query_atc = ("""SELECT ATC.value
    #                        FROM externalEntityATC ATC, {} U
    #                        WHERE U.externalEntityID = ATC.externalEntityID AND U.userEntityID = %s
    #                     """.format(unification_table))
    query_atc = ("""SELECT ATC.value, DB.databaseName
                    FROM externalEntityATC ATC, {} U, externalEntity E, externalDatabase DB
                    WHERE U.externalEntityID = ATC.externalEntityID AND ATC.externalEntityID = E.externalEntityID AND E.externalDatabaseID = DB.externalDatabaseID AND U.userEntityID = %s
                 """.format(unification_table))

    print('\nRETRIEVING ATC IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_atc_to_databases = {}

    for ueid in user_entity_ids:
        cursor.execute(query_atc, (ueid,))
        for row in cursor:
            atc, database = row
            #print(ueid, atc)
            ueid_to_atc_to_databases.setdefault(ueid, {})
            ueid_to_atc_to_databases[ueid].setdefault(atc, set()).add(database)

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH ATCS: {}'.format(len(ueid_to_atc_to_databases)))

    return ueid_to_atc_to_databases


def get_smiles_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get the smiles using their BIANA user entity ids
    """
    query_smiles = ("""SELECT SM.value, DB.databaseName
                       FROM externalEntitySMILES SM, {} U, externalEntity E, externalDatabase DB
                       WHERE U.externalEntityID = SM.externalEntityID AND SM.externalEntityID = E.externalEntityID AND E.externalDatabaseID = DB.externalDatabaseID AND U.userEntityID = %s
                    """.format(unification_table))

    print('\nRETRIEVING SMILES IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_smiles_to_databases = {}

    for ueid in user_entity_ids:
        cursor.execute(query_smiles, (ueid,))
        for row in cursor:
            smiles, database = row
            #print(ueid, smiles)
            ueid_to_smiles_to_databases.setdefault(ueid, {})
            ueid_to_smiles_to_databases[ueid].setdefault(smiles, set()).add(database)

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH SMILES IDS: {}'.format(len(ueid_to_smiles_to_databases)))

    return ueid_to_smiles_to_databases


def get_side_effects_from_sider(meddra_all_se_file):
    """
    Get the most frequent side effects from SIDER
    """
    pubchem_to_umls = {}
    umls_to_name = {}
    with open(meddra_all_se_file, 'r') as med_fd:
        for line in med_fd:
            fields = line.strip().split('\t')
            pubchem = str(int(fields[1].split('CID')[1]))
            concept_type = fields[3].upper()
            umls_id = fields[4].upper()
            umls_term = fields[5].lower()
            if concept_type == 'PT':
                pubchem_to_umls.setdefault(pubchem, set()).add(umls_id)
                umls_to_name[umls_id] = umls_term

    print('NUMBER OF PUBCHEM IDS ASSOCIATED WITH UMLS: {}'.format(len(pubchem_to_umls)))

    return pubchem_to_umls, umls_to_name


def get_side_effects_from_sider_by_frequency(meddra_freq_file, freq_cutoff=0.01):
    """
    Get the most frequent side effects from SIDER
    """
    pubchem_to_umls_to_freq = {}
    umls_to_name = {}
    with open(meddra_freq_file, 'r') as med_fd:
        for line in med_fd:
            fields = line.strip().split('\t')
            pubchem = str(int(fields[1].split('CID')[1]))
            freq_description = fields[4].lower()
            freq_low = fields[5].upper()
            freq_up = fields[6].upper()
            concept_type = fields[7].upper()
            umls_id = fields[8].upper()
            umls_term = fields[9].lower()
            if freq_low != '' and freq_up != '' and concept_type == 'PT':
                freq = ( float(freq_low) + float(freq_up) ) / 2
                if freq >= freq_cutoff:
                    pubchem_to_umls_to_freq.setdefault(pubchem, {})
                    pubchem_to_umls_to_freq[pubchem][umls_id] = freq
                    umls_to_name[umls_id] = umls_term

    print('NUMBER OF PUBCHEM IDS ASSOCIATED WITH UMLS: {}'.format(len(pubchem_to_umls_to_freq)))

    return pubchem_to_umls_to_freq, umls_to_name


def get_drug_targets_from_drugbank(cursor, unification_table, user_entity_ids):
    """
    Get the drug targets (therapeutic) of a list of drug user entity IDs from DrugBank
    """
    query_drugbank_targets = ("""SELECT U2.userEntityID, DT.value, DT.targetType 
                                 FROM {} U1, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, {} U2, externalEntityDrugBank_targetID DT
                                 WHERE U1.externalEntityID = R1.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityID != R2.externalEntityID AND R2.externalEntityID = U2.externalEntityID AND R2.externalEntityID = DT.externalEntityID AND U1.userEntityID = %s
              """.format(unification_table, unification_table))

    print('\nRETRIEVING DRUG TARGETS FROM DRUGBANK ASSOCIATED TO USER ENTITY IDS...\n')

    drugbank_drug_target_interactions = set()
    drugbank_drug_to_targets = {}
    drugbank_targets = set()

    drugbank_drug_metab_interactions = set()
    drugbank_drug_to_metab = {}
    drugbank_metab = set()

    for ueid1 in user_entity_ids:
        cursor.execute(query_drugbank_targets, (ueid1,))
        for row in cursor:
            ueid2, targetid, target_type = row
            if target_type.lower() == "therapeutic":
                #print(ueid1, ueid2, targetid, target_type)
                drugbank_targets.add(ueid2)
                interaction = (ueid1, ueid2)
                drugbank_drug_target_interactions.add(interaction)
                drugbank_drug_to_targets.setdefault(ueid1, set()).add(ueid2)
            else:
                drugbank_metab.add(ueid2)
                interaction = (ueid1, ueid2)
                drugbank_drug_metab_interactions.add(interaction)
                drugbank_drug_to_metab.setdefault(ueid1, set()).add(ueid2)

    print('NUMBER OF DRUG TARGET INTERACTIONS RETRIEVED FROM DRUGBANK: {}'.format(len(drugbank_drug_target_interactions)))
    print('NUMBER OF DRUG TARGETS RETRIEVED FROM DRUGBANK: {}'.format(len(drugbank_targets)))
    print('NUMBER OF DRUG METABOLIC TARGET INTERACTIONS RETRIEVED FROM DRUGBANK: {}'.format(len(drugbank_drug_metab_interactions)))
    print('NUMBER OF DRUG METABOLIC TARGETS RETRIEVED FROM DRUGBANK: {}'.format(len(drugbank_metab)))

    drug_target_interactions = drugbank_drug_target_interactions | drugbank_drug_metab_interactions
    targets = drugbank_targets | drugbank_metab
    print('NUMBER OF DRUG TARGET + METABOLIC INTERACTIONS RETRIEVED FROM DRUGBANK: {}'.format(len(drug_target_interactions)))
    print('NUMBER OF DRUG TARGETS + METABOLIC RETRIEVED FROM DRUGBANK: {}'.format(len(targets)))

    return drugbank_drug_target_interactions, drugbank_targets, drugbank_drug_to_targets, drugbank_drug_metab_interactions, drugbank_metab, drugbank_drug_to_metab


def get_drug_targets_from_drugcentral(cursor, unification_table, user_entity_ids):
    """
    Get the drug targets (tclin) of a list of drug user entity IDs from DrugCentral
    """
    query_drugcentral_targets = ("""SELECT U2.userEntityID, DCD.value 
                                    FROM {} U1, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, {} U2, externalEntityDrugCentral_druggability DCD
                                    WHERE U1.externalEntityID = R1.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityID != R2.externalEntityID AND R2.externalEntityID = U2.externalEntityID AND R2.externalEntityRelationID = DCD.externalEntityID AND DCD.value = 'tclin' AND U1.userEntityID = %s
                                 """.format(unification_table, unification_table))

    print('\nRETRIEVING DRUG TARGETS FROM DRUGCENTRAL ASSOCIATED TO USER ENTITY IDS...\n')

    drugcentral_drug_target_interactions = set()
    drugcentral_drug_to_targets = {}
    drugcentral_targets = set()

    for ueid1 in user_entity_ids:
        cursor.execute(query_drugcentral_targets, (ueid1,))
        for row in cursor:
            ueid2, druggability = row
            #print(ueid1, ueid2, druggability)
            drugcentral_targets.add(ueid2)
            interaction = (ueid1, ueid2)
            drugcentral_drug_target_interactions.add(interaction)
            drugcentral_drug_to_targets.setdefault(ueid1, set()).add(ueid2)

    print('NUMBER OF DRUG TARGET INTERACTIONS RETRIEVED FROM DRUGCENTRAL: {}'.format(len(drugcentral_drug_target_interactions)))
    print('NUMBER OF DRUG TARGETS RETRIEVED FROM DRUGCENTRAL: {}'.format(len(drugcentral_targets)))

    return drugcentral_drug_target_interactions, drugcentral_targets, drugcentral_drug_to_targets


def get_drug_targets_from_dgidb(cursor, unification_table, user_entity_ids):
    """
    Get the drug targets of a list of drug user entity IDs from DGIdb
    """
    query_dgidb_targets = ("""SELECT U2.userEntityID, DIS.value 
                              FROM {} U1, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, {} U2, externalEntityDGIDB_interaction_source DIS
                              WHERE U1.externalEntityID = R1.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityID != R2.externalEntityID AND R2.externalEntityID = U2.externalEntityID AND R2.externalEntityRelationID = DIS.externalEntityID AND U1.userEntityID = %s
                           """.format(unification_table, unification_table))

    print('\nRETRIEVING DRUG TARGETS FROM DGIDB ASSOCIATED TO USER ENTITY IDS...\n')

    dgidb_drug_target_interactions = set()
    dgidb_drug_to_targets = {}
    dgidb_targets = set()
    reliable_sources = ['chemblinteractions', 'guidetopharmacologyinteractions', 'tdgclinicaltrial', 'fda', 'tend', 'ttd']

    for ueid1 in user_entity_ids:
        cursor.execute(query_dgidb_targets, (ueid1,))
        for row in cursor:
            ueid2, source = row
            #print(ueid1, ueid2, source)
            if source in reliable_sources:
                dgidb_targets.add(ueid2)
                interaction = (ueid1, ueid2)
                dgidb_drug_target_interactions.add(interaction)
                dgidb_drug_to_targets.setdefault(ueid1, set()).add(ueid2)

    print('NUMBER OF DRUG TARGET INTERACTIONS RETRIEVED FROM DGIDB: {}'.format(len(dgidb_drug_target_interactions)))
    print('NUMBER OF DRUG TARGETS RETRIEVED FROM DGIDB: {}'.format(len(dgidb_targets)))

    return dgidb_drug_target_interactions, dgidb_targets, dgidb_drug_to_targets


def get_drug_targets_from_chembl(cursor, unification_table, user_entity_ids):
    """
    Get the drug targets of a list of drug user entity IDs from ChEMBL
    """
    query_chembl_targets = ("""SELECT U2.userEntityID, E2.type, DB.databaseName
                               FROM {} U1, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, {} U2, externalEntity E2, externalDatabase DB
                               WHERE U1.externalEntityID = R1.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityID != R2.externalEntityID AND R2.externalEntityID = U2.externalEntityID AND R2.externalEntityID = E2.externalEntityID AND E2.externalDatabaseID = DB.externalDatabaseID AND DB.databaseName = "chembl" AND E2.type = "protein" AND U1.userEntityID = %s
                            """.format(unification_table, unification_table))

    print('\nRETRIEVING DRUG TARGETS FROM CHEMBL ASSOCIATED TO USER ENTITY IDS...\n')

    chembl_drug_target_interactions = set()
    chembl_drug_to_targets = {}
    chembl_targets = set()

    for ueid1 in user_entity_ids:
        cursor.execute(query_chembl_targets, (ueid1,))
        for row in cursor:
            ueid2, ee_type, database = row
            #print(ueid1, ueid2, source)
            chembl_targets.add(ueid2)
            interaction = (ueid1, ueid2)
            chembl_drug_target_interactions.add(interaction)
            chembl_drug_to_targets.setdefault(ueid1, set()).add(ueid2)

    print('NUMBER OF DRUG TARGET INTERACTIONS RETRIEVED FROM CHEMBL: {}'.format(len(chembl_drug_target_interactions)))
    print('NUMBER OF DRUG TARGETS RETRIEVED FROM CHEMBL: {}'.format(len(chembl_targets)))

    return chembl_drug_target_interactions, chembl_targets, chembl_drug_to_targets


def get_drug_targets_from_ttd(cursor, unification_table, user_entity_ids):
    """
    Get the drug targets of a list of drug user entity IDs from TTD
    """
    query_ttd_targets = ("""SELECT U2.userEntityID, TT.value 
                            FROM {} U1, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, {} U2, externalEntityTTD_target_type TT
                            WHERE U1.externalEntityID = R1.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityID != R2.externalEntityID AND R2.externalEntityID = U2.externalEntityID AND R2.externalEntityID = TT.externalEntityID AND TT.value = "successful target" AND U1.userEntityID = %s
                         """.format(unification_table, unification_table))

    print('\nRETRIEVING DRUG TARGETS FROM TTD ASSOCIATED TO USER ENTITY IDS...\n')

    ttd_drug_target_interactions = set()
    ttd_drug_to_targets = {}
    ttd_targets = set()

    for ueid1 in user_entity_ids:
        cursor.execute(query_ttd_targets, (ueid1,))
        for row in cursor:
            ueid2, type_target = row
            #print(ueid1, ueid2, source)
            ttd_targets.add(ueid2)
            interaction = (ueid1, ueid2)
            ttd_drug_target_interactions.add(interaction)
            ttd_drug_to_targets.setdefault(ueid1, set()).add(ueid2)

    print('NUMBER OF DRUG TARGET INTERACTIONS RETRIEVED FROM TTD: {}'.format(len(ttd_drug_target_interactions)))
    print('NUMBER OF DRUG TARGETS RETRIEVED FROM TTD: {}'.format(len(ttd_targets)))

    return ttd_drug_target_interactions, ttd_targets, ttd_drug_to_targets


def get_geneids_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get the Entrez Gene IDs of targets using their BIANA user entity ids
    """
    query_geneid = ("""SELECT G.value, G.type
                       FROM externalEntityGeneID G, {} U 
                       WHERE U.externalEntityID = G.externalEntityID AND U.userEntityID = %s
                    """.format(unification_table))

    print('\nRETRIEVING GENE IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_geneid_to_types = {}

    for ueid in user_entity_ids:
        cursor.execute(query_geneid, (ueid,))
        for row in cursor:
            geneid, geneid_type = row
            #print(ueid, geneid, geneid_type)
            ueid_to_geneid_to_types.setdefault(ueid, {})
            ueid_to_geneid_to_types[ueid].setdefault(str(geneid), set()).add(geneid_type.lower())

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH GENE IDS: {}'.format(len(ueid_to_geneid_to_types)))

    return ueid_to_geneid_to_types


def get_uniprotaccessions_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get Uniprot Accessions of targets using their BIANA user entity ids
    """
    query_uniprotaccession = ("""SELECT UA.value, UA.type
                                 FROM externalEntityUniprotAccession UA, {} U 
                                 WHERE U.externalEntityID = UA.externalEntityID AND U.userEntityID = %s
                              """.format(unification_table))

    print('\nRETRIEVING UNIPROT ACCESSIONS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_uniprotaccession_to_types = {}

    for ueid in user_entity_ids:
        cursor.execute(query_uniprotaccession, (ueid,))
        for row in cursor:
            uniprotaccession, uniprotaccession_type = row
            #print(ueid, uniprotaccession, uniprotaccession_type)
            ueid_to_uniprotaccession_to_types.setdefault(ueid, {})
            ueid_to_uniprotaccession_to_types[ueid].setdefault(str(uniprotaccession), set()).add(uniprotaccession_type.lower())

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH UNIPROT ACCESSIONS: {}'.format(len(ueid_to_uniprotaccession_to_types)))

    return ueid_to_uniprotaccession_to_types


def get_uniprotentries_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get Uniprot Entries of targets using their BIANA user entity ids
    """
    query_uniprotentry = ("""SELECT UA.value, UA.type
                                 FROM externalEntityUniprotEntry UA, {} U 
                                 WHERE U.externalEntityID = UA.externalEntityID AND U.userEntityID = %s
                              """.format(unification_table))

    print('\nRETRIEVING UNIPROT ENTRIES ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_uniprotentry_to_types = {}

    for ueid in user_entity_ids:
        cursor.execute(query_uniprotentry, (ueid,))
        for row in cursor:
            uniprotentry, uniprotentry_type = row
            #print(ueid, uniprotentry, uniprotentry_type)
            ueid_to_uniprotentry_to_types.setdefault(ueid, {})
            ueid_to_uniprotentry_to_types[ueid].setdefault(str(uniprotentry), set()).add(uniprotentry_type.lower())

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH UNIPROT ENTRIES: {}'.format(len(ueid_to_uniprotentry_to_types)))

    return ueid_to_uniprotentry_to_types


def get_genesymbols_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get Gene Symbols of targets using their BIANA user entity ids
    """
    query_genesymbol = ("""SELECT UA.value, UA.type
                                 FROM externalEntityGeneSymbol UA, {} U 
                                 WHERE U.externalEntityID = UA.externalEntityID AND U.userEntityID = %s
                              """.format(unification_table))

    print('\nRETRIEVING GENE SYMBOLS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_genesymbol_to_types = {}

    for ueid in user_entity_ids:
        cursor.execute(query_genesymbol, (ueid,))
        for row in cursor:
            genesymbol, genesymbol_type = row
            #print(ueid, genesymbol, genesymbol_type)
            ueid_to_genesymbol_to_types.setdefault(ueid, {})
            ueid_to_genesymbol_to_types[ueid].setdefault(str(genesymbol), set()).add(genesymbol_type.lower())

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH GENE SYMBOLS: {}'.format(len(ueid_to_genesymbol_to_types)))

    return ueid_to_genesymbol_to_types


def get_taxid_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get Tax IDs of targets using their BIANA user entity ids
    """
    query_taxid = ("""SELECT TX.value, TX.type
                      FROM externalEntityTaxID TX, {} U 
                      WHERE U.externalEntityID = TX.externalEntityID AND U.userEntityID = %s
                   """.format(unification_table))

    print('\nRETRIEVING TAXONOMY IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_taxids = {}

    for ueid in user_entity_ids:
        cursor.execute(query_taxid, (ueid,))
        for row in cursor:
            taxid, taxid_type = row
            #print(ueid, taxid, taxid_type)
            ueid_to_taxids.setdefault(ueid, set()).add(taxid)

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH TAXONOMY IDS: {}'.format(len(ueid_to_taxids)))

    return ueid_to_taxids


def get_pfam_of_user_entity_ids(cursor, unification_table, user_entity_ids):
    """
    Get PFAM IDs of targets using their BIANA user entity ids
    """
    query_pfam = ("""SELECT PF.value, PF.type
                     FROM externalEntityPFAM PF, {} U 
                     WHERE U.externalEntityID = PF.externalEntityID AND U.userEntityID = %s
                  """.format(unification_table))

    print('\nRETRIEVING PFAM IDS ASSOCIATED TO USER ENTITY IDS...\n')

    ueid_to_pfams = {}

    for ueid in user_entity_ids:
        cursor.execute(query_pfam, (ueid,))
        for row in cursor:
            pfam, pfam_type = row
            #print(ueid, pfam, pfam_type)
            ueid_to_pfams.setdefault(ueid, set()).add(pfam)

    print('NUMBER OF USER ENTITIES ASSOCIATED WITH PFAMS: {}'.format(len(ueid_to_pfams)))

    return ueid_to_pfams


def output_drug_mappings(output_drug_mapping_file, user_entity_ids, ueid_to_drugtype_to_drugnames, ueid_to_drugbankids, ueid_to_pubchemcompounds, ueid_to_chemblids):
    """
    Output the drug mappings
    """

    print('\nWRITING OUTPUT DRUG MAPPING FILE...\n')

    with open(output_drug_mapping_file, 'w') as out_fd:
        out_fd.write('#user_entity_id\ttype_identifier\tidentifier\ttype_name\n')
        for user_entity_id in sorted(user_entity_ids):
            if user_entity_id in ueid_to_drugtype_to_drugnames:
                if 'unique' in ueid_to_drugtype_to_drugnames[user_entity_id]:
                    for drugname in ueid_to_drugtype_to_drugnames[user_entity_id]['unique']:
                        out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'name', drugname, 'unique'))
                for drugtype in ueid_to_drugtype_to_drugnames[user_entity_id]:
                    if drugtype != 'unique':
                        for drugname in ueid_to_drugtype_to_drugnames[user_entity_id][drugtype]:
                            out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'name', drugname, drugtype))
            if user_entity_id in ueid_to_drugbankids:
                for drugbankid in ueid_to_drugbankids[user_entity_id]:
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'drugbankid', drugbankid, '-'))
            if user_entity_id in ueid_to_pubchemcompounds:
                for pubchemcompound in ueid_to_pubchemcompounds[user_entity_id]:
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'pubchemcompound', pubchemcompound, '-'))
            if user_entity_id in ueid_to_chemblids:
                for chemblid in ueid_to_chemblids[user_entity_id]:
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'chemblid', chemblid, '-'))

    return


def output_target_mappings(output_target_mapping_file, targets, ueid_to_genesymbol_to_types, ueid_to_geneid_to_types, ueid_to_uniprotaccession_to_types, ueid_to_uniprotentry_to_types, ueid_to_pfams, ueid_to_taxids):
    """
    Output the target mappings
    """

    print('\nWRITING OUTPUT TARGET MAPPING FILE...\n')

    with open(output_target_mapping_file, 'w') as out_fd:
        out_fd.write('#user_entity_id\ttype_identifier\tidentifier\ttype_name\n')
        for user_entity_id in sorted(targets):
            if user_entity_id in ueid_to_genesymbol_to_types:
                for genesymbol in ueid_to_genesymbol_to_types[user_entity_id]:
                    types_genesymbol = '; '.join(sorted(ueid_to_genesymbol_to_types[user_entity_id][genesymbol]))
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'genesymbol', genesymbol, types_genesymbol))
            if user_entity_id in ueid_to_geneid_to_types:
                for geneid in ueid_to_geneid_to_types[user_entity_id]:
                    types_geneid = '; '.join(sorted(ueid_to_geneid_to_types[user_entity_id][geneid]))
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'geneid', geneid, types_geneid))
            if user_entity_id in ueid_to_uniprotaccession_to_types:
                for uniprotaccession in ueid_to_uniprotaccession_to_types[user_entity_id]:
                    types_uniprotaccession = '; '.join(sorted(ueid_to_uniprotaccession_to_types[user_entity_id][uniprotaccession]))
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'uniprotaccession', uniprotaccession, types_uniprotaccession))
            if user_entity_id in ueid_to_uniprotentry_to_types:
                for uniprotentry in ueid_to_uniprotentry_to_types[user_entity_id]:
                    types_uniprotentry = '; '.join(sorted(ueid_to_uniprotentry_to_types[user_entity_id][uniprotentry]))
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'uniprotentry', uniprotentry, types_uniprotentry))
            if user_entity_id in ueid_to_pfams:
                for pfam in ueid_to_pfams[user_entity_id]:
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'pfam', pfam, '-'))
            if user_entity_id in ueid_to_taxids:
                for taxid in ueid_to_taxids[user_entity_id]:
                    out_fd.write('{}\t{}\t{}\t{}\n'.format(user_entity_id, 'taxid', taxid, '-'))

    return


if  __name__ == "__main__":
    main()