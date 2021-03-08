import sys, os, re
import csv


class DrugCombDB(object):

    def __init__(self, input_path):

        self.input_path = input_path
        self.drug_chemical_info_file = os.path.join(self.input_path, 'drug_chemical_info.csv')
        self.asdcd_file = os.path.join(self.input_path, 'SynDrugComb_external_ASDCDsynergism.txt')
        self.textmining_file = os.path.join(self.input_path, 'SynDrugComb_textmining_2drugs.txt')
        self.fda_file = os.path.join(self.input_path, 'SynDrugComb_fda_2drugs.txt')

        self.pubchems = set()
        self.pubchem_to_drugname = {}
        self.pubchem_to_drugnameofficial = {}
        self.drugname_to_pubchem = {}

        self.combinations = set()
        self.combination_to_sources = {}
        self.combination_to_dcdb = {}
        self.pubchem_to_dcdb = {}
        self.combination_to_pubchems = {}
        self.combination_to_pubmeds = {}

        return


    def parse_drug_chemical_info(self):
        """
        Parsing of drug chemical info file
        """

        print("\n.....PARSING DRUGCOMBDB DRUG CHEMICAL INFO FILE.....\n")

        num_line = 0

        with open(self.drug_chemical_info_file, 'r', encoding = "ISO-8859-1") as drug_chemical_info_fd:

            csvreader = csv.reader(drug_chemical_info_fd, delimiter=',')

            for fields in csvreader:

                num_line += 1

                # Obtain a dictionary: "field_name" : "position"
                if num_line == 1:
                    #drugName,cIds,drugNameOfficial,molecularWeight,smilesString
                    fields_dict = self.obtain_header_fields( ','.join(fields), separator=',')
                    continue

                # Get useful fields
                drugname = fields[ fields_dict['drugName'] ].lower().lstrip().rstrip()
                drugnameofficial = fields[ fields_dict['drugNameOfficial'] ].lower()
                pubchem_cid = fields[ fields_dict['cIds'] ].upper() # One per drug (e.g. CIDs00060750)

                # Check if application numbers, ingredients and trade names are available
                if drugname == '' or drugname == '#n/a':
                    print('Drug name unknown for PubChem {}'.format(pubchem_cid))
                    sys.exit(10)
                if pubchem_cid == '' or pubchem_cid == '#N/A':
                    print('Pubchem unknown for drug {}'.format(drugname))
                    sys.exit(10)

                # Add pubchem
                if pubchem_cid.startswith('CIDS'):
                    pubchem_cid = int(pubchem_cid.split('CIDS')[1])
                else:
                    print('Unknown format for PubChem CID format for drug {}: {}'.format(drugname, pubchem_cid))
                    continue
                    #sys.exit(10)
                self.pubchems.add(pubchem_cid)
                self.pubchem_to_drugname.setdefault(pubchem_cid, set()).add(drugname)
                self.drugname_to_pubchem.setdefault(drugname, set()).add(pubchem_cid)

                if drugnameofficial != '' and drugnameofficial != '#n/a':
                    self.pubchem_to_drugnameofficial.setdefault(pubchem_cid, set()).add(drugnameofficial)
                    self.drugname_to_pubchem.setdefault(drugnameofficial, set()).add(pubchem_cid)

        #print(self.pubchem_to_drugnameofficial)

        return


    def parse_asdcd_combinations(self, drugname_to_pubchem):
        """
        Parsing of ASDCD combinations file
        """

        print("\n.....PARSING DRUGCOMBDB ASDCD DRUG COMBINATIONS FILE.....\n")

        drugs_found = set()
        drugs_not_found = set()

        with open(self.asdcd_file, 'r', encoding = "ISO-8859-1") as asdcd_fd:

            csvreader = csv.reader(asdcd_fd, delimiter='\t')

            for fields in csvreader:

                drugname_A = fields[0].lower().lstrip().rstrip()
                drugname_B = fields[1].lower().lstrip().rstrip()
                pubmed = fields[2] # If there is no pubmed, the column has a number "1"
                pubmed_link = fields[3] # If there is no pubmed link, the column has a number "1"
                article = fields[4] # If there is no article, the column has a number "1"
                combination = frozenset([drugname_A, drugname_B])
                self.combinations.add(combination)
                self.combination_to_sources.setdefault(combination, set()).add('asdcd')
                self.combination_to_pubmeds.setdefault(combination, set()).add(pubmed)

                if drugname_A in drugname_to_pubchem:
                    drugs_found.add(drugname_A)
                else:
                    drugs_not_found.add(drugname_A)
                    #print('Drug without pubchem: {}'.format(drugname_A))
                if drugname_B in drugname_to_pubchem:
                    drugs_found.add(drugname_B)
                else:
                    drugs_not_found.add(drugname_B)
                    #print('Drug without pubchem: {}'.format(drugname_B))

                #drug_combination = frozenset([drugname_A, drugname_B])

        print('Number of drugs with PubChem found: {}'.format(len(drugs_found)))
        print('Number of drugs without PubChem found: {}'.format(len(drugs_not_found)))

        return


    def parse_fda_combinations_all(self, drugname_to_pubchem):
        """
        Parsing of FDA combinations file.
        These are the same drug combinations as in DCDB.
        """

        print("\n.....PARSING DRUGCOMBDB FDA DRUG COMBINATIONS FILE (ALL).....\n")

        drugs_found = set()
        drugs_not_found = set()
        num_line = 0

        with open(self.fda_file, 'r') as fda_fd:

            csvreader = csv.reader(fda_fd, delimiter='\t')

            for fields in csvreader:

                num_line += 1

                # Obtain a dictionary: "field_name" : "position"
                if num_line == 1:
                    #ID	PubChem_ID	Single PubChem_ID	DC_ID 	COMPONENT	MACHENISM	DCC_ID
                    fields_dict = self.obtain_header_fields( '\t'.join(fields), separator='\t')
                    continue

                comb_id = fields[ fields_dict['ID'] ]
                comb_pubchem_id = fields[ fields_dict['PubChem_ID'] ].upper()
                pubchem_ids = fields[ fields_dict['Single PubChem_ID'] ].upper().split('/')
                comb_dcdb_id = fields[ fields_dict['DC_ID'] ].upper()
                drug_names = fields[ fields_dict['COMPONENT'] ].lower()
                dcdb_ids = fields[ fields_dict['DCC_ID'] ].upper().split('/')

                if 'NULL' in pubchem_ids:
                    print('One of the drugs of the combination {} without pubchem ID: {}'.format(comb_id, pubchem_ids))
                else:
                    for pubchem_id in pubchem_ids:
                        combination.append(pubchem_id)


        return


    def parse_fda_combinations_2drugs(self, drugname_to_pubchem):
        """
        Parsing of FDA combinations file.
        These are the same drug combinations as in DCDB.
        """

        print("\n.....PARSING DRUGCOMBDB FDA DRUG COMBINATIONS FILE (2 DRUGS).....\n")

        drugs_found = set()
        drugs_not_found = set()
        num_line = 0

        with open(self.fda_file, 'r') as fda_fd:

            csvreader = csv.reader(fda_fd, delimiter='\t')

            for fields in csvreader:

                num_line += 1

                # Obtain a dictionary: "field_name" : "position"
                if num_line == 1:
                    #ID	Drug1	Drug2	Mechanism	Source
                    fields_dict = self.obtain_header_fields( '\t'.join(fields), separator='\t')
                    continue

                comb_id = fields[ fields_dict['ID'] ]
                drugname_A = fields[ fields_dict['Drug1'] ].lower().lstrip().rstrip()
                drugname_B = fields[ fields_dict['Drug2'] ].lower().lstrip().rstrip()
                mechanism = fields[ fields_dict['Mechanism'] ]
                pubmeds = fields[ fields_dict['Source'] ].upper().split('/') # e.g. 5282054/6047

                if drugname_A == '' or drugname_A == 'null':
                    print('Drug name not available!')
                    sys.exit(10)
                if drugname_B == '' or drugname_B == 'null':
                    print('Drug name not available!')
                    sys.exit(10)

                combination = frozenset([drugname_A, drugname_B])
                self.combinations.add(combination)
                self.combination_to_sources.setdefault(combination, set()).add('dcdb')
                for pubmed in pubmeds:
                    if pubmed != '' and pubmed != 'NULL':
                        self.combination_to_pubmeds.setdefault(combination, set()).add(pubmed)

                if drugname_A in drugname_to_pubchem:
                    drugs_found.add(drugname_A)
                else:
                    drugs_not_found.add(drugname_A)
                    #print('Drug without pubchem: {}'.format(drugname_A))
                if drugname_B in drugname_to_pubchem:
                    drugs_found.add(drugname_B)
                else:
                    drugs_not_found.add(drugname_B)
                    #print('Drug without pubchem: {}'.format(drugname_B))

        print('Number of drugs with PubChem found: {}'.format(len(drugs_found)))
        print('Number of drugs without PubChem found: {}'.format(len(drugs_not_found)))

        return


    def parse_textmining_combinations(self, drugname_to_pubchem):
        """
        Parsing of text mining combinations file
        """

        print("\n.....PARSING DRUGCOMBDB TEXT MINING DRUG COMBINATIONS FILE.....\n")

        drugs_found = set()
        drugs_not_found = set()
        num_line = 0

        with open(self.textmining_file, 'r') as textmining_fd:

            csvreader = csv.reader(textmining_fd, delimiter='\t')

            for fields in csvreader:

                num_line += 1

                # Obtain a dictionary: "field_name" : "position"
                if num_line == 1:
                    #drug1 	drug2	target	Source_Pubchem_
                    #curcumin	genistein	rhodesain of Trypanosoma brucei rhodesiense	29897253
                    fields_dict = self.obtain_header_fields( '\t'.join(fields), separator='\t')
                    continue

                drugname_A = fields[ fields_dict['drug1'] ].lower().lstrip().rstrip()
                drugname_B = fields[ fields_dict['drug2'] ].lower().lstrip().rstrip()
                comb_pubchem_id = fields[ fields_dict['Source_Pubchem_'] ].lower()

                if drugname_A == '' or drugname_A == 'null':
                    print('Drug name not available!')
                    sys.exit(10)
                if drugname_B == '' or drugname_B == 'null':
                    print('Drug name not available!')
                    sys.exit(10)
                if comb_pubchem_id == '' or comb_pubchem_id == 'null':
                    print('Pubchem not available for combination: {} - {}'.format(drugname_A, drugname_B))
                    sys.exit(10)

                combination = frozenset([drugname_A, drugname_B])
                self.combinations.add(combination)
                self.combination_to_sources.setdefault(combination, set()).add('textmining')
                self.combination_to_pubchems.setdefault(combination, set()).add(comb_pubchem_id)

                if drugname_A in drugname_to_pubchem:
                    drugs_found.add(drugname_A)
                else:
                    drugs_not_found.add(drugname_A)
                    #print('Drug without pubchem: {}'.format(drugname_A))
                if drugname_B in drugname_to_pubchem:
                    drugs_found.add(drugname_B)
                else:
                    drugs_not_found.add(drugname_B)
                    #print('Drug without pubchem: {}'.format(drugname_B))

        print('Number of drugs with PubChem found: {}'.format(len(drugs_found)))
        print('Number of drugs without PubChem found: {}'.format(len(drugs_not_found)))

        return


    def obtain_header_fields(self, first_line, separator='\t'):
        """ 
        Obtain a dictionary: "field_name" => "position" 
        """
        fields_dict = {}

        header_fields = first_line.strip().split(separator)
        for x in range(0, len(header_fields)):
            fields_dict[header_fields[x]] = x

        return fields_dict


if __name__ == "__main__":
    main()