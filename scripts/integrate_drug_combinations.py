import argparse
import pandas as pd
import time
import sys, os, re

import dcdbParser
import DrugCombDBParser
import FDAOrangeBookParser


def main():

    options = parse_user_arguments()
    integrate_drug_combinations(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    python /home/quim/PHD/Projects/DIANA/diana/scripts/integrate_drug_combinations.py -db /home/quim/Databases -m /home/quim/PHD/Projects/DIANA/diana/outputs -o /home/quim/PHD/Projects/DIANA/diana/outputs
    """

    parser = argparse.ArgumentParser(
        description = "Create drug mappings using BIANA database",
        epilog      = "@oliva's lab 2020")
    parser.add_argument('-db','--databases_dir',dest='databases_dir',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'databases'),
                        help = """Directory where the additional data is stored""")
    parser.add_argument('-m','--mappings_dir',dest='mappings_dir',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'mappings'),
                        help = """Directory where the mapping files of DrugBank to drug names, cross-references and drug targets are stored""")
    parser.add_argument('-o','--output_dir',dest='output_dir',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'mappings'),
                        help = """Directory where the output files will be stored""")

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def integrate_drug_combinations(options):
    """
    Creates a benchmark of drug combinations by integrating the information from different databases.
    """

    # Start marker for time measure
    start = time.time()

    # Get the inputs
    databases_dir = options.databases_dir
    create_directory(databases_dir)
    mappings_dir = options.mappings_dir
    output_dir = options.output_dir
    create_directory(output_dir)

    # Get database directories
    dcdb_dir = os.path.join(databases_dir, 'DCDB/DCDB2_plaintxt/dcdb')
    fdaorangebook_dir = os.path.join(databases_dir, 'FDAOrangeBook/EOBZIP_2020_04')
    drugcombdb_dir = os.path.join(databases_dir, 'DrugCombDB')
    cheng_dir = os.path.join(databases_dir, 'Cheng_NatCom19')

    # Check mapping files
    mappings_file = os.path.join(mappings_dir, 'drugbank_drug_mappings.txt')
    drug_to_targets_file = os.path.join(mappings_dir, 'drugbank_geneid_drug_to_targets.txt')
    if not fileExist(mappings_file):
        print('Mappings file does not exist: {}'.format(mappings_file))
        sys.exit(10)
    if not fileExist(drug_to_targets_file):
        print('Drug to targets file does not exist: {}'.format(drug_to_targets_file))
        sys.exit(10)


    #------------#
    # PARSE DCDB #
    #------------#

    dcdb_parser = dcdbParser.DCDB(dcdb_dir)
    dcdb_parser.parse()


    #-----------------------#
    # PARSE FDA ORANGE BOOK #
    #-----------------------#

    fda_parser = FDAOrangeBookParser.FDAOrangeBook(fdaorangebook_dir)
    fda_parser.parse_products_file()


    #------------------#
    # PARSE DRUGCOMBDB #
    #------------------#

    drugcombdb_parser = DrugCombDBParser.DrugCombDB(drugcombdb_dir)
    drugcombdb_parser.parse_drug_chemical_info()
    drugcombdb_parser.parse_asdcd_combinations(drugname_to_pubchem=drugcombdb_parser.drugname_to_pubchem)
    #drugcombdb_parser.parse_fda_combinations_all(drugname_to_pubchem=drugcombdb_parser.drugname_to_pubchem)
    drugcombdb_parser.parse_fda_combinations_2drugs(drugname_to_pubchem=drugcombdb_parser.drugname_to_pubchem)
    drugcombdb_parser.parse_textmining_combinations(drugname_to_pubchem=drugcombdb_parser.drugname_to_pubchem)


    #------------------------------------------------#
    # PARSE CHENG DRUG COMBINATIONS (BARABASI STUDY) #
    #------------------------------------------------#

    cheng_combinations_file = os.path.join(cheng_dir, 'Cheng_NatCom19_AllDrugCombinations.txt')
    cheng_raw_pairs = set()
    with open(cheng_combinations_file, 'r') as inp_fd:
        first_line = inp_fd.readline()
        for line in inp_fd:
            drug1, drug2 = line.strip().split('\t')
            combination = frozenset([drug1, drug2])
            cheng_raw_pairs.add(combination)


    #---------------------#
    # PARSE MAPPINGS FILE #
    #---------------------#

    print("\n.....PARSING MAPPINGS FILE.....\n")
    drugbankid_to_identifiers = {}
    drugbankid_to_typename_to_name = {}
    identifier_to_drugbankids = {}
    typename_to_name_to_drugbankids = {}
    with open(mappings_file, 'r') as inp_fd:
        #drugbankid	type_identifier	identifier	type_name
        first_line = inp_fd.readline()
        for line in inp_fd:
            drugbankid, type_identifier, identifier, type_name = line.strip().split('\t')
            drugbankid_to_identifiers.setdefault(drugbankid, {})
            drugbankid_to_identifiers[drugbankid].setdefault(type_identifier, set()).add(identifier)
            identifier_to_drugbankids.setdefault(type_identifier, {})
            identifier_to_drugbankids[type_identifier].setdefault(identifier, set()).add(drugbankid)
            if type_identifier == 'name':
                drugbankid_to_typename_to_name.setdefault(drugbankid, {})
                drugbankid_to_typename_to_name[drugbankid].setdefault(type_name, set()).add(identifier)
                typename_to_name_to_drugbankids.setdefault(type_name, {})
                typename_to_name_to_drugbankids[type_name].setdefault(identifier, set()).add(drugbankid)


    #----------------------------------#
    # PARSE MAPPINGS FILE WITH TARGETS #
    #----------------------------------#

    drugbankid_to_geneids_network = {}
    with open(drug_to_targets_file, 'r') as inp_fd:
        #drugbankid	type_identifier	identifier	type_name
        first_line = inp_fd.readline()
        for line in inp_fd:
            fields = line.strip().split('\t')
            if len(fields) == 3:
                drugbankid, geneids, geneids_network = fields
                geneids_network = geneids_network.split('; ')
                if len(geneids_network) > 0:
                    for geneid in geneids_network:
                        drugbankid_to_geneids_network.setdefault(drugbankid, set()).add(geneid)


    #------------------------#
    # INTEGRATE COMBINATIONS #
    #------------------------#

    print("\n.....INTEGRATING COMBINATIONS.....\n")
    combination_to_databases = {}
    combination_to_usage = {}
    multipledrugbankids_to_originalname = {}

    # DCDB
    # sources in DCDB: {'pubchem compound', 'pharmgkb', 'wikipedia', 'chebi', 'pdb', 'bindingdb', 'drugbank', 'kegg drug', 'pubchem substance', 'rxlist', 'kegg compound'}
    dcdb_combinations = set()
    dcdb_combinations_not_pairs = set()
    dcdb_drugs_found = set()
    dcdb_drugs_not_found = set()
    for combination in dcdb_parser.combination2component:
        skip = False
        combination_drugbankids = []
        for component in dcdb_parser.combination2component[combination]:
            # Get drug name and crossreferences for each component
            drugname = dcdb_parser.component2name[component]
            source_to_crossreferences = {}
            drugbankids = set()
            if component in dcdb_parser.component2outlink:
                for outlink in dcdb_parser.component2outlink[component]:
                    source = dcdb_parser.outlink2sourcelink[outlink]['source']
                    crossreference = dcdb_parser.outlink2sourcelink[outlink]['link']
                    source_to_crossreferences.setdefault(source, set()).add(crossreference)
            # Check if pubchem compounds are in our mapping file
            if len(drugbankids) == 0 and 'pubchem compound' in source_to_crossreferences:
                drugbankids = set()
                for pubchemcompound in source_to_crossreferences['pubchem compound']:
                    if str(pubchemcompound) in identifier_to_drugbankids['pubchemcompound']:
                        for drugbankid in identifier_to_drugbankids['pubchemcompound'][str(pubchemcompound)]:
                            drugbankids.add(drugbankid)
                if len(drugbankids) > 0:
                    combination_drugbankids.append('; '.join(sorted(drugbankids)))
                    multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(drugname)
            # Check if drugbank ids are in our mapping file
            if len(drugbankids) == 0 and 'drugbank' in source_to_crossreferences:
                drugbankids = set()
                for drugbankid in source_to_crossreferences['drugbank']:
                    if drugbankid in drugbankid_to_identifiers:
                        drugbankids.add(drugbankid)
                if len(drugbankids) > 0:
                    combination_drugbankids.append('; '.join(sorted(drugbankids)))
                    multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(drugname)
            # Check if name is in unique names of our mapping file
            if len(drugbankids) == 0 and drugname in typename_to_name_to_drugbankids['unique']:
                drugbankids = typename_to_name_to_drugbankids['unique'][drugname]
                combination_drugbankids.append('; '.join(sorted(drugbankids)))
                multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(drugname)
            # Check if name is in our mapping file
            if len(drugbankids) == 0 and drugname in identifier_to_drugbankids['name']:
                drugbankids = identifier_to_drugbankids['name'][drugname]
                combination_drugbankids.append('; '.join(sorted(drugbankids)))
                multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(drugname)
            # If not, annotate the drug as not found
            if len(drugbankids) == 0:
                skip = True
                dcdb_drugs_not_found.add(drugname)
        if not skip:
            if len(combination_drugbankids) == 2:
                # Add combinations
                dcdb_combinations.add(frozenset(combination_drugbankids))
                combination_to_databases.setdefault(frozenset(combination_drugbankids), set()).add('dcdb')
                for drug in combination_drugbankids:
                    dcdb_drugs_found.add(drug)
                # Get usage
                if combination in dcdb_parser.combination2usage:
                    for usage in dcdb_parser.combination2usage[combination]:
                        if usage in dcdb_parser.usage2icd10codes:
                            for icd10code in dcdb_parser.usage2icd10codes[usage]:
                                icd10name = dcdb_parser.icd10code2name[icd10code]
                                efficacy = '-'
                                effecttype = '-'
                                if usage in dcdb_parser.usage2efficacy:
                                    efficacy = dcdb_parser.usage2efficacy[usage]
                                if usage in dcdb_parser.usage2effecttype:
                                    effecttype = dcdb_parser.usage2effecttype[usage]
                                if efficacy != 'efficacious':
                                    continue
                                usage_complete = '{}|{}|{}|{}'.format(icd10code, icd10name, efficacy, effecttype)
                                combination_to_usage.setdefault(frozenset(combination_drugbankids), set()).add(usage_complete)
            else:
                dcdb_combinations_not_pairs.add(frozenset(combination_drugbankids))

    #print(dcdb_drugs_not_found)
    print('Number of drugs in DCDB not found: {}'.format(len(dcdb_drugs_not_found)))
    print('Number of drugs in DCDB found: {}'.format(len(dcdb_drugs_found)))
    #print(dcdb_combinations_not_pairs)
    print('Number of drug combinations in DCDB that are not pairs: {}'.format(len(dcdb_combinations_not_pairs)))
    print('Number of drug combinations in DCDB that are pairs: {}'.format(len(dcdb_combinations)))


    # FDAORANGEBOOK
    fda_combinations = set()
    fda_combinations_not_pairs = set()
    fda_drugs_found = set()
    fda_drugs_not_found = set()
    for appl_no in fda_parser.appl_numbers:
        skip = False
        combination_drugbankids = set()
        if appl_no in fda_parser.appl_no_to_ingredients:
            if len(fda_parser.appl_no_to_ingredients[appl_no]) > 1:
                for ingredient in fda_parser.appl_no_to_ingredients[appl_no]:
                    drugbankids = set()
                    # Check if name is in unique names of our mapping file
                    if len(drugbankids) == 0 and ingredient in typename_to_name_to_drugbankids['unique']:
                        drugbankids = typename_to_name_to_drugbankids['unique'][ingredient]
                        combination_drugbankids.add('; '.join(sorted(drugbankids)))
                        multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(ingredient)
                    # Check if name is in our mapping file
                    if len(drugbankids) == 0 and ingredient in identifier_to_drugbankids['name']:
                        drugbankids = identifier_to_drugbankids['name'][ingredient]
                        combination_drugbankids.add('; '.join(sorted(drugbankids)))
                        multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(ingredient)
                    # If not, annotate the drug as not found
                    if len(drugbankids) == 0:
                        skip = True
                        fda_drugs_not_found.add(ingredient)
                if not skip:
                    #print(len(combination_drugbankids), combination_drugbankids, len(fda_parser.appl_no_to_ingredients[appl_no]), fda_parser.appl_no_to_ingredients[appl_no])
                    if len(combination_drugbankids) == 2:
                        fda_combinations.add(frozenset(combination_drugbankids))
                        combination_to_databases.setdefault(frozenset(combination_drugbankids), set()).add('fdaorangebook')
                        for drug in combination_drugbankids:
                            fda_drugs_found.add(drug)
                    else:
                        fda_combinations_not_pairs.add(frozenset(combination_drugbankids))

    #print(fda_drugs_not_found)
    print('Number of drugs in FDAOrangeBook not found: {}'.format(len(fda_drugs_not_found)))
    print('Number of drugs in FDAOrangeBook found: {}'.format(len(fda_drugs_found)))
    #print(fda_combinations_not_pairs)
    print('Number of drug combinations in FDAOrangeBook that are not pairs: {}'.format(len(fda_combinations_not_pairs)))
    print('Number of drug combinations in FDAOrangeBook that are pairs: {}'.format(len(fda_combinations)))


    # DRUGCOMBDB
    drugcombdb_combinations = set()
    drugcombdb_combinations_not_pairs = set()
    drugcombdb_drugs_found = set()
    drugcombdb_drugs_not_found = set()
    for combination in drugcombdb_parser.combinations:
        skip = False
        combination_drugbankids = []
        drugcombdb_sources = drugcombdb_parser.combination_to_sources[combination]
        for drugname in combination:
            drugbankids = set()
            # Check if pubchem compounds are in our mapping file
            if len(drugbankids) == 0 and drugname in drugcombdb_parser.drugname_to_pubchem:
                drugbankids = set()
                for pubchemcompound in drugcombdb_parser.drugname_to_pubchem[drugname]:
                    if str(pubchemcompound) in identifier_to_drugbankids['pubchemcompound']:
                        for drugbankid in identifier_to_drugbankids['pubchemcompound'][str(pubchemcompound)]:
                            drugbankids.add(drugbankid)
                if len(drugbankids) > 0:
                    combination_drugbankids.append('; '.join(sorted(drugbankids)))
                    multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(drugname)
            # Check if name is in unique names of our mapping file
            if len(drugbankids) == 0 and drugname in typename_to_name_to_drugbankids['unique']:
                drugbankids = typename_to_name_to_drugbankids['unique'][drugname]
                combination_drugbankids.append('; '.join(sorted(drugbankids)))
                multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(drugname)
            # Check if name is in our mapping file
            if len(drugbankids) == 0 and drugname in identifier_to_drugbankids['name']:
                drugbankids = identifier_to_drugbankids['name'][drugname]
                combination_drugbankids.append('; '.join(sorted(drugbankids)))
                multipledrugbankids_to_originalname.setdefault('; '.join(sorted(drugbankids)), set()).add(drugname)
            # If not, annotate the drug as not found
            if len(drugbankids) == 0:
                skip = True
                drugcombdb_drugs_not_found.add(drugname)
        if not skip:
            if len(combination_drugbankids) == 2:
                drugcombdb_combinations.add(frozenset(combination_drugbankids))
                for source in drugcombdb_sources:
                    combination_to_databases.setdefault(frozenset(combination_drugbankids), set()).add('drugcombdb-'+source)
                for drug in combination_drugbankids:
                    drugcombdb_drugs_found.add(drug)
            else:
                drugcombdb_combinations_not_pairs.add(frozenset(combination_drugbankids))

    #print(drugcombdb_drugs_not_found)
    print('Number of drugs in DrugCombDB not found: {}'.format(len(drugcombdb_drugs_not_found)))
    print('Number of drugs in DrugCombDB found: {}'.format(len(drugcombdb_drugs_found)))
    #print(drugcombdb_combinations_not_pairs)
    print('Number of drug combinations in DrugCombDB that are not pairs: {}'.format(len(drugcombdb_combinations_not_pairs)))
    print('Number of drug combinations in DrugCombDB that are pairs: {}'.format(len(drugcombdb_combinations)))


    # CHENG DRUG COMBINATIONS
    cheng_combinations = set()
    cheng_combinations_not_pairs = set()
    cheng_drugs_found = set()
    cheng_drugs_not_found = set()
    for combination in cheng_raw_pairs:
        skip=False
        if len(combination) != 2:
            cheng_combinations_not_pairs.add(combination)
            continue
        for drugbankid in combination:
            if drugbankid in drugbankid_to_identifiers:
                cheng_drugs_found.add(drugbankid)
            else:
                cheng_drugs_not_found.add(drugbankid)
                skip=True
        if not skip:
            cheng_combinations.add(combination)
            combination_to_databases.setdefault(frozenset(combination), set()).add('cheng')

    #print(cheng_drugs_not_found)
    #print(cheng_combinations_not_pairs)
    print('Number of drugs in Cheng not found: {}'.format(len(cheng_drugs_not_found)))
    print('Number of drugs in Cheng found: {}'.format(len(cheng_drugs_found)))
    print('Number of raw drug combinations in Cheng: {}'.format(len(cheng_raw_pairs)))
    print('Number of drug combinations in Cheng that are pairs: {}'.format(len(cheng_combinations)))
    print('Number of drug combinations in Cheng that are not pairs: {}'.format(len(cheng_combinations_not_pairs)))


    # CHECK DRUGS WITH MULTIPLE DRUGBANK IDS
    multiple_drugbankids_to_databases = {}
    for combination in combination_to_databases:
        databases = combination_to_databases[combination]
        for drug in combination:
            if '; ' in drug:
                for database in databases:
                    multiple_drugbankids_to_databases.setdefault(drug, set()).add(database)
    output_file = os.path.join(output_dir, 'drug_combination_drugs_with_multiple_drugbankids.txt')
    with open(output_file, 'w') as out_fd:
        for drugbankids in multiple_drugbankids_to_databases:
            databases = '; '.join(sorted(multiple_drugbankids_to_databases[drugbankids]))
            drugnames = []
            for drugbankid in drugbankids.split('; '):
                drugnames.append(list(drugbankid_to_typename_to_name[drugbankid]['unique'])[0])
            originalnames = '; '.join(sorted(multipledrugbankids_to_originalname[drugbankids]))
            out_fd.write('{}\t{}\t{}\t{}\n'.format(drugbankids, '; '.join(drugnames), originalnames, databases))

    # Curated cases of multiple drugbank ids after revising the file generated above
    # There are still some cases that I consider that all drugbankIDs could be valid
    # In these cases, we will consider all of them as separate drugs
    multiple_drugbankids_to_unique = {
        'DB01373; DB13257' : 'DB13257',
        'DB00022; DB00041' : 'DB00022; DB00041',
        'DB06346; DB12186' : 'DB06346',
        'DB00316; DB00956; DB01050' : 'DB00956',
        'DB01914; DB09341' : 'DB01914',
        'DB00388; DB11086; DB11088' : 'DB00388',
        'DB00898; DB02325' : 'DB00898',
        'DB01694; DB09459' : 'DB01694; DB09459',
        'DB00286; DB09317' : 'DB00286',
        'DB00007; DB02325' : 'DB00007',
        'DB00005; DB06770; DB09145' : 'DB00005',
        'DB00882; DB06735' : 'DB00882',
        'DB14507; DB14509' : 'DB14507; DB14509',
        'DB00097; DB09126; DB09145' : 'DB00097; DB09126'
    }


    # CHECK DRUGS NOT FOUND
    drugs_not_found_to_databases = {}
    #drugs_not_found = dcdb_drugs_not_found | drugcombdb_drugs_not_found | fda_drugs_not_found | cheng_drugs_not_found
    for drugs_not_found, database in [(dcdb_drugs_not_found, 'dcdb'), (drugcombdb_drugs_not_found, 'drugcombdb'), (fda_drugs_not_found, 'fdaorangebook'), (cheng_drugs_not_found, 'cheng')]:
        for drug in drugs_not_found:
            drugs_not_found_to_databases.setdefault(drug, set()).add(database)
    output_file = os.path.join(output_dir, 'drug_combination_drugs_not_found.txt')
    with open(output_file, 'w') as out_fd:
        for drug in sorted(drugs_not_found_to_databases):
            databases = '; '.join(sorted(drugs_not_found_to_databases[drug]))
            out_fd.write('{}\t{}\n'.format(drug, databases))


    # WRITE FINAL FILE
    drugcombinations_file = os.path.join(output_dir, 'DIANA_drug_combinations.txt')
    with open(drugcombinations_file, 'w') as out_fd:
        out_fd.write('#DrugBankID1\tDrugBankID2\tDrugName1\tDrugName2\tDatabases\tICD10code|ICD10name|Efficacy|TypeEffect\n')
        for combination in combination_to_databases:
            databases = '; '.join(sorted(combination_to_databases[combination]))
            #print(combination, databases)
            usages = []
            if combination in combination_to_usage:
                usages = sorted(combination_to_usage[combination])
            drug1, drug2 = combination # Some drugs may have multiple drugbank IDs due to integration issues
            # Get unique IDs
            if drug1 in multiple_drugbankids_to_unique:
                drug1 = multiple_drugbankids_to_unique[drug1]
            if drug2 in multiple_drugbankids_to_unique:
                drug2 = multiple_drugbankids_to_unique[drug2]
            # Still there can be multiple drugbank ids. For these cases, we count them as separate drugs
            for drugbankid1 in drug1.split('; '):
                drugname1 = list(drugbankid_to_typename_to_name[drugbankid1]['unique'])[0]
                for drugbankid2 in drug2.split('; '):
                    drugname2 = list(drugbankid_to_typename_to_name[drugbankid2]['unique'])[0]
                    out_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(drugbankid1, drugbankid2, drugname1, drugname2, databases, '; '.join(usages)))


    # WRITE FINAL FILE WITH TARGETS
    drugcomb_with_targets_file = os.path.join(output_dir, 'DIANA_drug_combinations_with_targets_network.txt')
    drugs_with_targets = set(drugbankid_to_geneids_network.keys())
    drugcomb_df = pd.read_csv(drugcombinations_file, sep='\t', index_col=None)
    drugcomb_with_targets_df = drugcomb_df[(drugcomb_df['#DrugBankID1'].isin(drugs_with_targets)) & (drugcomb_df['DrugBankID2'].isin(drugs_with_targets))]
    drugcomb_with_targets_df.to_csv(drugcomb_with_targets_file, sep='\t', index=False)


    # CREATE IVENN FILE (TO CREATE A VENN DIAGRAM)
    database_to_combinations = {}
    for combination in combination_to_databases:
        drug1, drug2 = combination # Some drugs may have multiple drugbank IDs due to integration issues
        for drugbankid1 in drug1.split('; '):
            for drugbankid2 in drug2.split('; '):
                for database in combination_to_databases[combination]:
                    if database.startswith('drugcombdb'):
                        database = 'drugcombdb'
                    database_to_combinations.setdefault(database, set()).add('{}---{}'.format(drugbankid1, drugbankid2))

    ivenn_output_file_network = os.path.join(output_dir, 'drugcombinations_with_targets_network.ivenn')
    with open(ivenn_output_file_network, 'w') as out_ivenn_net:
        for database in ['dcdb', 'drugcombdb', 'fdaorangebook', 'cheng']:
            out_ivenn_net.write('{}:{};\n'.format(database, ','.join(database_to_combinations[database])))



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


if  __name__ == "__main__":
    main()