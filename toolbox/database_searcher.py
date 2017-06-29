import cPickle
import re
import sqlite3
import sys

def get_smiles_from_dcdbid(dcdbid):
    """
    Searches the SMILES of a DCDB drug using BIANA
    """

    dump_file = "toolbox/"+"dcdb2smiles.pcl"
    dcdb2smiles = cPickle.load(open(dump_file))

    try:
        return dcdb2smiles[dcdbid]
    except:
        return None

def result_to_list(result):
    """
    Transforms the result of the search to drugbank from a list of tuples to a sigle list
    """
    final_list = []
    for item in result:
        final_list.append(item[0])
    return final_list

def search_DDI_in_drugbank(drug1, drug2):
    """
    Performs a search of drug-drug interactions in DrugBank
    """

    # Connection to SQLite db
    conn = sqlite3.connect('drugbank.db')

    # Create a Cursor object in order to call its execute() method to perform SQL commands
    conn.text_factory = str
    c = conn.cursor()

    if process_drug_name(drug2) == 'name':
        c.execute('''SELECT drugbank_id
                      FROM Drugs
                              WHERE name = (?)
        ''', (drug2.lower(),))
        result = result_to_list(c.fetchall())
        result = result[0]
        c.execute('''SELECT i.description
                      FROM Drugs d
                        INNER JOIN Interactions i
                          ON i.Drugs_drugbank_id = d.drugbank_id
                              WHERE d.name = ? AND i.Drugs_drugbank_id1 = ?
        ''', (drug1.lower(), result,))
        command = c.fetchall()

    elif process_drug_name(drug2) == 'drugbank_id':
        c.execute('''SELECT i.description
                      FROM Drugs d
                        INNER JOIN Interactions i
                          ON i.Drugs_drugbank_id = d.drugbank_id
                              WHERE d.drugbank_id = (?) AND i.Drugs_drugbank_id1 = (?)
        ''', (drug1.lower(), result,))
        command = c.fetchall()


    # Close the connection
    conn.close()

    return command

def search_seeds_in_drugbank(drug, provided_seeds = None, protein_type = 'uniprotentry'):
    """
    Performs a search of targets in DrugBank using the name of the drug provided.
    If seeds are provided, it adds the seeds to the targets found.
    It returns the targets found including the seeds provided
    """

    # Connection to SQLite db
    conn = sqlite3.connect('drugbank.db')

    # Create a Cursor object in order to call its execute() method to perform SQL commands
    c = conn.cursor()

    # Depending on the type of protein chosen, the sql query will be different
    if protein_type == 'uniprotentry':
        protein_type = 'p.uniprot_accession'
    elif protein_type == 'geneid':
        protein_type = 'p.gene_id'

    # Process the drug name provided to see if it is a conventional drug name or a DrugBank ID
    # If it is a conventional name, the the query will be performed using the name
    if process_drug_name(drug) == 'name':
        drugbank_query = 'd.name'
    # If a DrugBank ID is provided, the query will be performed using DrugBank ID
    elif process_drug_name(drug) == 'drugbank_id':
        drugbank_query = 'd.drugbank_id'
    # If not, the program is terminated with an error
    else:
        if provided_seeds == None:
            print("  DIANA INFO:\tThe introduced drug %s is not found in DrugBank.\n\t\tTry to introduce a valid name for the drug.\n" % (drug))
        else:
            print("  DIANA INFO:\tThe introduced drug %s is not found in DrugBank.\n\t\tTry to introduce a valid name for the drug, or a higher number of seeds.\n" % (drug))
        sys.exit(10)

    # Define the sql query inside a variable, introducing the protein_type and drugbank_query variables
    sql_query = '''SELECT %s 
                      FROM Drugs d
                        INNER JOIN Drugs_has_Targets dht
                          ON dht.Drugs_drugbank_id = d.drugbank_id
                        INNER JOIN Targets t
                          ON  t.idTargets = dht.Targets_idTargets
                        INNER JOIN Polypeptides p
                          ON p.Targets_idTargets = t.idTargets
                              WHERE %s = (?)
        ''' %(protein_type, drugbank_query)

    c.execute(sql_query, (drug.lower(),))
    list_of_targets = result_to_list(c.fetchall())    

    # If a list of seeds is provided, add the seeds inside the list of targets
    if provided_seeds != None:
        for seed in provided_seeds:
            if seed not in list_of_targets:
                list_of_targets.append(seed)

    # Close the connection
    conn.close()

    return list_of_targets

def search_seeds_in_file(dcdbid, provided_seeds = None):
    """
    Performs a search of targets in the file extracted from the database
    If seeds are provided, it adds the seeds to the targets found.
    It returns the targets found including the seeds provided
    """

    dump_file = "toolbox/"+"dcdb2targets.pcl"
    dcdb2targets = cPickle.load(open(dump_file))
    list_of_targets = dcdb2targets[dcdbid]

    # If a list of seeds is provided, add the seeds inside the list of targets
    if provided_seeds != None:
        for seed in provided_seeds:
            if seed not in list_of_targets:
                list_of_targets.add(seed)

    return list_of_targets


def search_pfams_for_seeds(seeds):
    """
    Performs a search in the database to find the PFAMs of the seeds.
    Returns a set of PFAMs
    """

    dump_file = "toolbox/"+"geneid2pfam.pcl"
    geneid2pfam = cPickle.load(open(dump_file))
    
    all_pfams = set()

    for seed in seeds:
        pfams = geneid2pfam[seed]
        for pfam in pfams:
            all_pfams.add(pfam)

    return all_pfams


def process_drug_name(drug):
    """
    Processes the drug name provided in order to know if it is a conventional drug name or a
    DrugBank ID. Returns the strings 'drugbank_id' or 'name' according to the result
    """

    # DrugBank ID regular expression pattern
    p = re.compile('^DB\d{5}$')

    # If the drug name provided matches the pattern, it returns the string 'drugbank_id'
    if p.match(drug) != None:
        return 'drugbank_id'
    # If it does not match the pattern, it returns the string 'name'
    else:
        return 'name'