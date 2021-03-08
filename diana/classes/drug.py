import os, sys, re
import pickle
import pandas as pd
import hashlib


class Drug(object):
    """
    Class defining a Drug object
    """

    def __init__(self, drug_name):
        """
        @param:    drug_name
        @pdef:     Name of the drug
        @ptype:    {String}

        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        self.drug_name = drug_name.lower()
        self.type_name = self.recognize_name(drug_name.lower())
        self.targets = []
        self.targets_in_network = []
        self.pfams = []
        self.smiles = []
        self.ATCs = []
        self.SEs = []
        self.target_type_id = None
        self.target_type_id_to_table = {
            'geneid' : 'externalEntityGeneID',
            'genesymbol' : 'externalEntityGeneSymbol',
            'uniprotentry' : 'externalEntityUniprotEntry',
            'uniprotaccession' : 'externalEntityUniprotAccession',
        }
        self.type_name_to_table = {
            'name' : 'externalEntityName',
            'drugbankid' : 'externalEntityDrugBankID',
            'dcdb' : 'externalEntityDCDB_drugID',
            'chemblid' : 'externalEntityCHEMBL',
            'pubchemcompound' : 'externalEntityPubChemCompound',
        }

    ###########
    # METHODS #
    ###########

    def obtain_targets_from_file(self, targets_file, target_type_id):
        """
        Obtains the targets from an input file and stores them into a list.
        The file must contain the names of the targets separated by new lines.
        The type of ID of the targets must be specified.
        """
        self.target_type_id = target_type_id.lower() # Annotate the type of ID of the targets
        with open(targets_file, 'r') as targets_file_fd:
            for line in targets_file_fd:
                self.targets.append(line.strip())

        # Check if the number of targets provided is sufficient for the analysis
        if len(self.targets) < 1:
            raise InsufficientTargets(self.targets)
        return

    def obtain_targets_from_pickle(self, drug2targets_file, target_type_id):
        """
        Obtains the targets from an input pickle file and stores them into a list.
        """
        self.target_type_id = target_type_id.lower() # Annotate the type of ID of the targets

        drug2targets = pickle.load(open(drug2targets_file))
        drug_id = self.drug_name.upper()

        if drug_id in drug2targets:
            self.targets = list(drug2targets[drug_id])
        else:
            raise InsufficientTargets(self.targets)

        # Check if the number of targets provided is sufficient for the analysis
        if len(self.targets) < 1:
            raise InsufficientTargets(self.targets)
        return

    def obtain_drugbankids_from_table(self, drug_mapping_file):
        """
        Obtains the drugbankids of a drug from an input table and stores them into a list.
        Usually, there is only one drugbankid, but there could be multiple ones in some occasions.
        """
        # Get DrugBankID for the input drug
        if self.type_name != 'drugbankid':
            drug_mapping_df = pd.read_csv(drug_mapping_file, sep='\t', index_col=None)
            if self.type_name == 'name':
                # Select drugbank ids with the input name
                drugnames_df = drug_mapping_df[(drug_mapping_df['type_identifier'] == 'name') & (drug_mapping_df['identifier'] == self.drug_name)]
                drugbankids = set(drugnames_df['#drugbankid'].tolist())
                if len(drugbankids) == 0:
                    raise DrugNameNotFound(self.drug_name, self.type_name)
                elif len(drugbankids) > 1:
                    # Check if the input name is unique
                    if 'unique' in drugnames_df['type_name'].tolist():
                        drugbankids = set(drugnames_df.loc[drugnames_df['type_name'] == 'unique', '#drugbankid'].tolist())
                        if len(drugbankids) == 0:
                            drugbankids = set(drugnames_df['#drugbankid'].tolist())
            else:
                drugbankids = set(drug_mapping_df.loc[(drug_mapping_df['type_identifier'] == self.type_name) & (drug_mapping_df['identifier'] == self.drug_name), '#drugbankid'].tolist())
                if len(drugbankids) == 0:
                    raise DrugNameNotFound(self.drug_name, self.type_name)
        else:
            drugbankids = [self.drug_name.upper()]
        return drugbankids

    def obtain_targets_from_table(self, drugbankids, drug_to_targets_file, target_type_id='geneid'):
        """
        Obtains the targets from an input table and stores them into a list.
        """
        self.target_type_id = target_type_id.lower() # Annotate the type of ID of the targets
        # Get targets
        targets = set()
        drug_to_targets_df = pd.read_csv(drug_to_targets_file, sep='\t', index_col=None)
        for group_targets in drug_to_targets_df.loc[drug_to_targets_df['#drugbankid'].isin(drugbankids), 'geneids'].tolist():
            targets = targets | set(group_targets.split('; '))

        # Check if the number of targets provided is sufficient for the analysis
        if len(targets) < 1:
            raise InsufficientTargets(self.targets)
        else:
            self.targets = targets

        return

    def obtain_targets_from_BIANA(self, biana_cnx, target_type_id, unification_protocol):
        """
        Obtains the targets from BIANA database using as query the drug name.
        The type of ID of the targets must be specified.
        "biana_cnx" parameter stands for the variable containing the connexion to the MySQL database.
        """
        self.target_type_id = target_type_id.lower() # Annotate the type of ID of the targets

        target_type_id_table = self.return_targets_biana_table(self.target_type_id) # Obtain the table containing the type of ID introduced
        type_name_table = self.return_drug_biana_table(self.type_name) # Obtain the table containing the type of name introduced

        up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

        cursor = biana_cnx.cursor() # Start cursor to MySQL

        # Select the external entity ID of the DCDB drug
        query1 = (''' SELECT externalEntityID FROM {} WHERE value = %s
                 '''.format(type_name_table))

        # Select the geneID targets of the drug, only therapeutic ones, and only from DrugBank database!
        query2 = (''' SELECT G.value FROM externalEntity E1, {} U1, {} U2, externalEntity E2, externalEntityRelationParticipant R2, externalEntityRelationParticipant R3, externalEntityDrugBank_targetID T, {} U3, {} U4, externalEntityGeneID G
                     WHERE E1.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = E2.externalEntityID AND E2.type = 'drug'
                     AND E2.externalEntityID = R2.externalEntityID AND R2.externalEntityRelationID = R3.externalEntityRelationID AND R3.externalEntityID = T.externalEntityID AND T.targetType = "therapeutic"
                     AND R3.externalEntityID = U3.externalEntityID AND U3.userEntityID = U4.userEntityID AND U4.externalEntityID = G.externalEntityID AND E1.externalEntityID = %s
                 '''.format(up_table, up_table, up_table, up_table, target_type_id_table))

        cursor.execute(query1, (self.drug_name,))

        external_entities = set()
        geneids = set()

        # Search for the external entities corresponding to the name of the drug
        for items in cursor:
            for ee in items:
                external_entities.add(ee)
        # Search for the geneIDs interacting with the drug
        if len(external_entities) > 0:
            for ee in external_entities:
                cursor.execute(query2, (ee,))
                for items in cursor:
                    for geneid in items:
                        geneids.add(geneid)
        else:
            raise DrugNameNotFound(self.drug_name, self.type_name)
        # Why in two steps?
        # Because as the table "externalEntityName" is too large, do one complex command can be very time-consuming
        # It is better to split the search in two commands

        cursor.close()

        self.targets = list(geneids)

        # Check if the number of targets provided is sufficient for the analysis
        if len(self.targets) < 1:
            raise InsufficientTargets(self.targets)

        return

    def recognize_name(self, drug_name):
        """
        Recognizes the type of name of the drug
        (dcdb, drugbank or name)
        """

        dcdb_pattern = re.compile('^dcc[0-9]{4}$')
        drugbank_pattern = re.compile('^db[0-9]{5}$')
        chembl_pattern = re.compile('^chembl[0-9]+$')
        pubchem_pattern = re.compile('^[0-9]+$')
        diana_pattern = re.compile('^diana_.*$')

        if dcdb_pattern.match(drug_name):
            self.drug_name = drug_name.upper()
            return 'dcdb'
        elif drugbank_pattern.match(drug_name):
            self.drug_name = drug_name.upper()
            return 'drugbankid'
        elif chembl_pattern.match(drug_name):
            self.drug_name = drug_name.upper()
            return 'chemblid'
        elif pubchem_pattern.match(drug_name):
            return 'pubchemcompound'
        elif diana_pattern.match(drug_name):
            return 'diana'
        else:
            return 'name'

    def return_drug_biana_table(self, type_name):
        """
        Returns the table in BIANA where the type of drug name
        introduced is stored.
        """
        if type_name in self.type_name_to_table:
            return self.type_name_to_table[type_name]

    def return_targets_biana_table(self, target_type_id):
        """
        Returns the table in BIANA where the annotations of the type of ID
        introduced are stored.
        """
        if target_type_id in self.target_type_id_to_table:
            return self.target_type_id_to_table[target_type_id]
        else:
            raise IncorrectTypeID(target_type_id, self.target_type_id_to_table)

    def obtain_pfams_from_file(self, pfam_file):
        """
        Obtains the pfams from an input file and stores them into a list.
        The file must contain the names of the pfams separated by new lines.
        """
        with open(pfam_file, 'r') as pfam_file_fd:
            for line in pfam_file_fd:
                self.pfams.append(line.strip())
        return

    def obtain_pfams_from_pickle(self, pfam_pickle_file, output_file):
        """
        Obtains the pfams from an input pickle file and stores them into a list.
        """
        geneid2pfam = pickle.load(open(pfam_pickle_file))

        all_pfams = set()

        for target in self.targets:
            if target in geneid2pfam:
                pfams = geneid2pfam[target]
                for pfam in pfams:
                    all_pfams.add(pfam)

        if len(all_pfams) > 0:
            self.pfams = list(all_pfams)
            with open(output_file, 'w') as pfam_fd:
                for pfam in self.pfams:
                    pfam_fd.write('{}\n'.format(pfam))
        else:
            print('No PFAMS found for the targets introduced: {}.\n'.format(', '.join(self.targets)))

        return

    def obtain_pfams_from_geneid_target_table(self, geneids, geneid_target_mapping_file):
        """
        Obtains the pfams of a list of targets (in gene ID) from an input table and stores them into a list.
        """
        # Get pfams
        geneid_mappings_df = pd.read_csv(geneid_target_mapping_file, sep='\t', index_col=None)
        pfams_df = geneid_mappings_df[(geneid_mappings_df['#geneid'].isin(geneids)) & (geneid_mappings_df['type_identifier'] == 'pfam')]
        self.pfams = set([pfam.upper() for pfam in pfams_df['identifier'].tolist()])
        return

    def obtain_pfams_from_targets(self, biana_cnx, output_file, unification_protocol):
        """
        Obtains the pfams from BIANA database using as query the targets.
        "biana_cnx" parameter stands for the variable containing the connexion to the MySQL database.
        Stores the PFAMs found in an output file
        """

        target_type_id_table = self.return_targets_biana_table(self.target_type_id) # Obtain the table containing the type of ID introduced

        up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

        query = (''' SELECT P.value FROM {} G, {} U1, {} U2, externalEntityPFAM P
                     WHERE G.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = P.externalEntityID AND G.value = %s
                 '''.format(target_type_id_table, up_table, up_table))

        if len(self.targets) > 0:

            cursor = biana_cnx.cursor() # Start cursor to MySQL

            for target in self.targets:
                cursor.execute(query, (target,))

                pfams = set()
                for items in cursor:
                    for pfam in items:
                        pfams.add(pfam.upper())

            cursor.close()

        else:
            print('There are no targets, so it is impossible to get the PFAMs!\n')
            sys.exit(10)

        if len(pfams) > 0:
            self.pfams = list(pfams)
            with open(output_file, 'w') as pfam_fd:
                for pfam in self.pfams:
                    pfam_fd.write('{}\n'.format(pfam))
        else:
            print('No PFAMS found for the targets introduced: {}.\n'.format(', '.join(self.targets)))

        return

    def obtain_SMILES_from_file(self, smiles_file):
        """
        Obtains the SMILES from an input file and stores them into a list.
        The file must contain the SMILES separated by new lines.
        """
        with open(smiles_file, 'r') as smiles_file_fd:
            for line in smiles_file_fd:
                self.smiles.append(line.strip())
        return

    def obtain_SMILES_from_table(self, drugbankids, drugbank_smiles_file):
        """
        Obtains the SMILES of a drug from an input table and stores them into a list.
        """
        drugbank_smiles_df = pd.read_csv(drugbank_smiles_file, sep='\t', index_col=None)
        self.smiles = set(drugbank_smiles_df.loc[drugbank_smiles_df['#drugbankid'].isin(drugbankids), 'smiles'].tolist())
        return

    def obtain_SMILES_from_pickle(self, smiles_pickle_file, output_file):
        """
        Obtains the SMILES from an input pickle file and stores them into a list.
        """
        drug2smiles = pickle.load(open(smiles_pickle_file))
        drug_id = self.drug_name.upper()
        if drug_id in drug2smiles:
            self.smiles = drug2smiles[drug_id]
        else:
            self.smiles = None

        if len(self.smiles) > 0:
            with open(output_file, 'w') as smiles_fd:
                for result in self.smiles:
                    smiles_fd.write('{}\n'.format(result))
        else:
            print('No SMILES found for the drug {}.\n'.format(self.drug_name))

        return

    def obtain_SMILES_from_BIANA(self, biana_cnx, output_file, unification_protocol):
        """
        Obtains the SMILES from BIANA database using as query the name of the drug.
        "biana_cnx" parameter stands for the variable containing the connexion to the MySQL database.
        Stores the SMILES in an output file.
        If there is more than one different SMILES, they are printed separated by new lines.
        """

        up_table = return_unification_protocol_table(biana_cnx, unification_protocol)
        type_name_table = self.return_drug_biana_table(self.type_name) # Obtain the table containing the type of name introduced

        query = (''' SELECT S.value FROM {} N, {} U1, {} U2, externalEntitySMILES S
                     WHERE N.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = S.externalEntityID AND N.value = %s
                 '''.format(type_name_table, up_table, up_table))

        cursor = biana_cnx.cursor() # Start cursor to MySQL

        cursor.execute(query, (self.drug_name,))

        smiles = set()
        for items in cursor:
            for result in items:
                smiles.add(result)

        cursor.close()

        if len(smiles) > 0:
            self.smiles = list(smiles)
            with open(output_file, 'w') as smiles_fd:
                for result in self.smiles:
                    smiles_fd.write('{}\n'.format(result))
        else:
            print('No SMILES found for the drug {}.\n'.format(self.drug_name))

        return

    def obtain_ATCs_from_file(self, ATCs_file):
        """
        Obtains the ATCs from an input file and stores them into a list.
        The file must contain the names of the ATCs separated by new lines.
        """
        with open(ATCs_file, 'r') as ATC_file_fd:
            for line in ATC_file_fd:
                self.ATCs.append(line.strip())
        return

    def obtain_ATCs_from_table(self, drugbankids, drugbank_atc_file):
        """
        Obtains the pfams of a list of targets (in gene ID) from an input table and stores them into a list.
        """
        drugbank_atc_df = pd.read_csv(drugbank_atc_file, sep='\t', index_col=None)
        atcs = drugbank_atc_df.loc[drugbank_atc_df['#drugbankid'].isin(drugbankids), 'atc'].tolist()
        self.ATCs = set([atc.upper() for atc in atcs])
        return

    def obtain_ATCs_from_BIANA(self, biana_cnx, output_file, unification_protocol):
        """
        Obtains the ATCs from BIANA database using as query the targets.
        "biana_cnx" parameter stands for the variable containing the connexion to the MySQL database.
        Stores the ATCs found in an output file.
        """

        up_table = return_unification_protocol_table(biana_cnx, unification_protocol)
        type_name_table = self.return_drug_biana_table(self.type_name) # Obtain the table containing the type of name introduced

        query = (''' SELECT A.value FROM {} N, {} U1, {} U2, externalEntityATC A
                     WHERE N.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = A.externalEntityID AND N.value = %s
                 '''.format(type_name_table, up_table, up_table))

        cursor = biana_cnx.cursor()
        cursor.execute(query, (self.drug_name,))

        ATCs = set()
        for items in cursor:
            for ATC in items:
                ATCs.add(ATC.upper())
        cursor.close()

        if len(ATCs) > 0:
            self.ATCs = list(ATCs)
            with open(output_file, 'w') as ATCs_fd:
                for ATCs in self.ATCs:
                    ATCs_fd.write('{}\n'.format(ATCs))
        else:
            print('  DIANA INFO:\tNo ATCs for the drug introduced: {}.\n'.format(self.drug_name))

        return

    def obtain_ATCs_from_pickle(self, atc_pickle_file, output_file):
        """
        Obtains the ATCs from an input pickle file and stores them into a list.
        """
        drug2atcs = pickle.load(open(atc_pickle_file))
        drug_id = self.drug_name.upper()
        if drug_id in drug2atcs:
            self.ATCs = drug2atcs[drug_id]
        else:
            self.ATCs = set()

        if len(self.ATCs) > 0:
            with open(output_file, 'w') as atc_fd:
                for result in self.ATCs:
                    atc_fd.write('{}\n'.format(result))
        else:
            print('No ATCs found for the drug {}.\n'.format(self.drug_name))

        return

    def obtain_SE_from_file(self, SE_file):
        """
        Obtains the SE from an input file and stores them into a list.
        The file must contain the names of the SE separated by new lines.
        """
        with open(SE_file, 'r') as SE_file_fd:
            for line in SE_file_fd:
                self.SEs.append(line.strip())
        return

    def obtain_SE_from_table(self, drugbankids, drugbank_side_effects_file):
        """
        Obtains the side effects of a drug from an input table and stores them into a list.
        """
        drugbank_side_effects_df = pd.read_csv(drugbank_side_effects_file, sep='\t', index_col=None)
        self.SEs = set(drugbank_side_effects_df.loc[drugbank_side_effects_df['#drugbankid'].isin(drugbankids), 'umls_id'].tolist())
        return

    def obtain_SE_from_BIANA(self, biana_cnx, output_file, unification_protocol):
        """
        Obtains the SE from BIANA database using as query the targets.
        "biana_cnx" parameter stands for the variable containing the connexion to the MySQL database.
        Stores the SE found in an output file.
        """

        target_type_id_table = self.return_targets_biana_table(self.target_type_id)

        up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

        query = (''' SELECT UD.value FROM {} D, {} U1, {} U2, externalEntityPubChemCompound PC, externalEntityRelationParticipant P1, externalEntityRelationParticipant P2, externalEntityRelation R, externalEntityUMLS_diseaseID UD
                     WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = PC.externalEntityID AND PC.externalEntityID = P1.externalEntityID AND D.value = %s
                     AND P1.externalEntityRelationID = R.externalEntityRelationID AND P2.externalEntityRelationID = R.externalEntityRelationID AND P1.externalEntityID != P2.externalEntityID AND R.type = "drug_phenotype_association"
                     AND P2.externalEntityID = UD.externalEntityID
                 '''.format(target_type_id_table, up_table, up_table))

        cursor = biana_cnx.cursor()
        cursor.execute(query, (self.drug_name,))

        SEs = set()
        for items in cursor:
            for SEs in items:
                SEs.add(SEs.upper())
        cursor.close()

        if len(SEs) > 0:
            self.SEs = list(SEs)
            with open(output_file, 'w') and SEs_fd:
                for SEs in self.SEs:
                    SEs_fd.write('{}\n'.format(SEs))
        else:
            print('  DIANA INFO:\tNo Side Effects for the drug introduced: {}.\n'.format(self.drug_name))

        return

    def obtain_SE_from_pickle(self, se_pickle_file, output_file):
        """
        Obtains the side effects from an input pickle file and stores them into a list.
        """
        drug2side_effects = pickle.load(open(se_pickle_file))
        drug_id = self.drug_name.upper()
        if drug_id in drug2side_effects:
            self.SEs = drug2side_effects[drug_id]
        else:
            self.SEs = set()

        if len(self.SEs) > 0:
            with open(output_file, 'w') as se_fd:
                for result in self.SEs:
                    se_fd.write('{}\n'.format(result))
        else:
            print('No SEs found for the drug {}.\n'.format(self.drug_name))

        return





class InsufficientTargets(Exception):
    """
    Exception raised when the number of targets is below 3.
    This exception is raised because the analyses of GUILD with less than 3
    targets are not reliable.
    """
    def __init__(self, targets, limit_targets=1):
        self.targets = targets
        self.limit_targets = limit_targets

    def __str__(self):
        return 'The number of targets provided ({}) is insufficient.\nGUILD must have at least 1 target to run a reliable analysis.\n'.format(len(self.targets), self.limit_targets)

class DrugNameNotFound(Exception):
    """
    Exception raised when the drug name is not found in BIANA.
    """
    def __init__(self, drug_name, type_name):
        self.drug_name = drug_name
        self.type_name = type_name

    def __str__(self):
        return 'The drug {} {} has not been found in BIANA.\nTherefore, any target could be found. Please, introduce another name or the targets of the drug.\n'.format(self.type_name, self.drug_name)

class IncorrectTypeID(Exception):
    """
    Exception that raises when a type of IDs of the proteins is not admitted for
    the program.
    """
    def __init__(self, target_type_id, target_type_id_to_table):
        self.target_type_id = target_type_id
        self.target_type_id_to_table = target_type_id_to_table

    def __str__(self):
        return 'The initial type of IDs of the proteins ({}) is not admitted.\nThe types of ID admitted in DIANA are: {}\n'.format(self.target_type_id, ', '.join(self.target_type_id_to_table.keys()))




def generate_diana_id(drug_name, targets, network_name):
    """
    Generates an ID for the drug using the drug name, the sorted targets and
    the network name.
    """
    id_str = 'diana_'
    drug_str = ''.join(drug_name.split('\s')) # Obtain the drug name, and remove any space in the name
    id_str += drug_str.lower() # Add the drug name in the ID string
    targets = [str(x) for x in targets] # Transform the targets to strings
    targets_str = ''.join(sorted(targets)) # Join the targets in one string
    id_str += targets_str.lower() # Add the targets in the ID string
    id_str += network_name.lower() # Add the network name in the ID string
    id_str=id_str.encode('utf-8') # Encode it by utf-8
    m = hashlib.md5()
    m.update(id_str) # Introduce the string in the hashlib instance
    unique_id = m.hexdigest()[:12] # Obtain a unique ID from the string. Only get the first 12 characters
    return unique_id

def old_generate_drug_id(drug_name, targets, network_name):
    """
    Generates an ID for the drug using the drug name, the sorted targets and
    the network name.
    """
    id_str = ''
    drug_str = ''.join(drug_name.split('\s')) # Obtain the drug name, and remove any space in the name
    id_str += drug_str.lower() # Add the drug name in the ID string
    targets = [str(x) for x in targets] # Transform the targets to strings
    targets_str = ''.join(sorted(targets)) # Join the targets in one string
    id_str += targets_str.lower() # Add the targets in the ID string
    id_str += network_name.lower() # Add the network name in the ID string
    id_str=id_str.encode('utf-8') # Encode it by utf-8
    m = hashlib.md5()
    m.update(id_str) # Introduce the string in the hashlib instance
    unique_id = m.hexdigest()[:12] # Obtain a unique ID from the string. Only get the first 12 characters
    return unique_id

def create_targets_file(targets, file_name):
    """
    Creates a targets file, containing the targets separated by new line characters
    """
    with open(file_name, 'w') as fw:
        for target in targets:
            fw.write('{}\n'.format(target))
    return

def create_number_of_targets_file(drugbank_geneid_mapping_file, number_of_targets_file):
    """
    Creates a file with the total number of different targets.
    """
    drugbank_geneid_mappings_df = pd.read_csv(drugbank_geneid_mapping_file, sep='\t', index_col=None)
    targets = set(drugbank_geneid_mappings_df['geneid'].tolist())
    with open(number_of_targets_file, 'w') as number_of_targets_fd:
        number_of_targets_fd.write('{}\n'.format(len(targets)))
    return len(targets)

def create_number_of_pfams_file(geneid_target_mapping_file, number_of_pfams_file):
    """
    Creates a file with the total number of different PFAMs associated to targets.
    """
    geneid_mappings_df = pd.read_csv(geneid_target_mapping_file, sep='\t', index_col=None)
    pfams_df = geneid_mappings_df[geneid_mappings_df['type_identifier'] == 'pfam']
    #pfams_df = geneid_mappings_df[(geneid_mappings_df['#geneid'].isin(network_nodes)) & (geneid_mappings_df['type_identifier'] == 'pfam')]
    pfams = set([pfam.upper() for pfam in pfams_df['identifier'].tolist()])
    with open(number_of_pfams_file, 'w') as number_of_pfams_fd:
        number_of_pfams_fd.write('{}\n'.format(len(pfams)))
    return len(pfams)

def get_all_targets_from_mappings(drugbank_geneid_mapping_file):
    """
    Get all targets from the drugbank geneid mapping file.
    """
    drugbank_geneid_mappings_df = pd.read_csv(drugbank_geneid_mapping_file, sep='\t', index_col=None)
    targets = set(map(str, drugbank_geneid_mappings_df['geneid'].tolist()))
    return targets

def get_all_pfams_from_mappings(geneid_target_mapping_file):
    """
    Get all PFAMs from the geneid target mapping file
    """
    geneid_mappings_df = pd.read_csv(geneid_target_mapping_file, sep='\t', index_col=None)
    pfams_df = geneid_mappings_df[geneid_mappings_df['type_identifier'] == 'pfam']
    pfams = set([pfam.upper() for pfam in pfams_df['identifier'].tolist()])
    return pfams

def get_all_atcs_from_mappings(drugbank_atc_file):
    """
    Get all ATCs from the drugbank ATC mapping file
    """
    drugbank_atc_df = pd.read_csv(drugbank_atc_file, sep='\t', index_col=None)
    atcs = set(drugbank_atc_df['atc'].tolist())
    return atcs

def get_all_ses_from_mappings(drugbank_side_effects_file):
    """
    Get all side effects from the drugbank side effects mapping file
    """
    drugbank_side_effects_df = pd.read_csv(drugbank_side_effects_file, sep='\t', index_col=None)
    ses = set(drugbank_side_effects_df['umls_id'].tolist())
    return ses

def read_number_file(number_file):
    """
    Reads the file with the total number of targets/PFAMs.
    """
    with open(number_file, 'r') as number_fd:
        number_of_entities = int(number_fd.readline().strip("\n"))
    return number_of_entities

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

def obtain_drugbank_to_targets(biana_cnx, unification_protocol, sif_file, output_pickle_file):
    """
    Obtains a file containing the targets of every DrugBank drug
    """

    # Get all the nodes in the network
    all_nodes = set()
    with open(sif_file, 'r') as sif_file_fd:
        for line in sif_file_fd:
            node1, score, node2 = line.strip().split('\t')
            all_nodes.add(int(node1))
            all_nodes.add(int(node2))

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query1 = (''' SELECT value FROM externalEntityDrugBankID ''')
    query2 = (''' SELECT G.value FROM externalEntityDrugBankID D, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, externalEntityDrugBank_targetID T, {} U2, {} U3, externalEntityGeneID G
                 WHERE D.externalEntityID = R1.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND  R2.externalEntityID = T.externalEntityID AND T.targetType = "therapeutic"
                 AND R2.externalEntityID = U2.externalEntityID AND U2.userEntityID = U3.userEntityID AND U3.externalEntityID = G.externalEntityID AND D.value = %s
             '''.format(up_table, up_table))

    cursor.execute(query1)

    drugbank2targets = {}
    drugbank_ids = set()

    # Search for all the DrugbankIDs
    for items in cursor:
        for drugbankid in items:
            drugbank_ids.add(drugbankid)
    # Search for the geneIDs interacting with the drugs
    for drugbankid in drugbank_ids:
        geneids = set()
        cursor.execute(query2, (drugbankid,))
        for items in cursor:
            for geneid in items:
                geneids.add(geneid)
        if len(geneids) > 0:
            for geneid in geneids:
                if geneid in all_nodes:
                    drugbank2targets.setdefault(drugbankid, set())
                    drugbank2targets[drugbankid].add(geneid)

    cursor.close()

    print(drugbank2targets)

    pickle.dump(drugbank2targets, open(output_pickle_file, 'wb'))

    return drugbank2targets

def obtain_dcdb_to_targets(biana_cnx, unification_protocol, sif_file, output_pickle_file):
    """
    Obtains a file containing the targets of every DCDB drug
    """

    # Get all the nodes in the network
    all_nodes = set()
    with open(sif_file, 'r') as sif_file_fd:
        for line in sif_file_fd:
            node1, score, node2 = line.strip().split('\t')
            all_nodes.add(int(node1))
            all_nodes.add(int(node2))

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query1 = (''' SELECT value FROM externalEntityDCDB_drugID ''')
    query2 = (''' SELECT G.value FROM externalEntityDCDB_drugID D, {} U1, {} U2, externalEntity E2, externalEntityRelationParticipant R2, externalEntityRelationParticipant R3, externalEntityDrugBank_targetID T, {} U3, {} U4, externalEntityGeneID G
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = E2.externalEntityID AND E2.type = 'drug'
                 AND E2.externalEntityID = R2.externalEntityID AND R2.externalEntityRelationID = R3.externalEntityRelationID AND R3.externalEntityID = T.externalEntityID AND T.targetType = "therapeutic"
                 AND R3.externalEntityID = U3.externalEntityID AND U3.userEntityID = U4.userEntityID AND U4.externalEntityID = G.externalEntityID AND D.value = %s
             '''.format(up_table, up_table, up_table, up_table))

    cursor.execute(query1)

    dcdb2targets = {}
    dcdb_ids = set()

    # Search for all the DCDB drug IDs
    for items in cursor:
        for dcdbid in items:
            dcdb_ids.add(dcdbid)
    # Search for the geneIDs interacting with the drugs
    for dcdbid in dcdb_ids:
        geneids = set()
        cursor.execute(query2, (dcdbid,))
        for items in cursor:
            for geneid in items:
                geneids.add(geneid)
        if len(geneids) > 0:
            for geneid in geneids:
                if geneid in all_nodes:
                    dcdb2targets.setdefault(dcdbid, set())
                    dcdb2targets[dcdbid].add(geneid)

    cursor.close()

    print(dcdb2targets)

    pickle.dump(dcdb2targets, open(output_pickle_file, 'wb'))

    return dcdb2targets

def obtain_target_to_pfam(biana_cnx, unification_protocol, all_targets, output_pickle_file):
    """
    Obtains a file containing the PFAMs of every target
    """

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query = (''' SELECT P.value FROM externalEntityGeneID G, {} U1, {} U2, externalEntityPFAM P
                 WHERE G.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = P.externalEntityID AND G.value = %s
             '''.format(up_table, up_table))

    geneid2pfam = {}

    # Search the PFAMS of the targets
    for target in all_targets:
        pfams = set()
        cursor.execute(query, (target,))
        for items in cursor:
            for pfam in items:
                pfams.add(pfam)
        if len(pfams) > 0:
            geneid2pfam[target] = pfams

    cursor.close()

    print(geneid2pfam)

    pickle.dump(geneid2pfam, open(output_pickle_file, 'wb'))

    return geneid2pfam

def obtain_drug_to_smiles(biana_cnx, unification_protocol, all_drugs, type_drug_name, output_pickle_file):
    """
    Obtains a file containing the SMILES of every drug.
    The type of drug name must be indicated (drugbank, dcdb, name)
    """

    type_name_to_table = {
        'name' : 'externalEntityName',
        'drugbank' : 'externalEntityDrugBankID',
        'dcdb' : 'externalEntityDCDB_drugID',
    }

    type_drug_table = type_name_to_table[type_drug_name] # Obtain the table containing the type of name introduced
    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    # Obtain the SMILES for the drugs
    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query = (''' SELECT S.value FROM {} D, {} U1, {} U2, externalEntitySMILES S
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = S.externalEntityID AND D.value = %s
             '''.format(type_drug_table, up_table, up_table))

    drug2smiles = {}

    for drug in all_drugs:
        smiles = set()
        cursor.execute(query, (drug,))
        for items in cursor:
            for result in items:
                smiles.add(result)
        if len(smiles) > 0:
            drug2smiles[drug] = smiles

    cursor.close()


    # Obtain the PubChem IDs for the drugs
    cursor = biana_cnx.cursor()

    query = (''' SELECT D.value, P.value FROM {} D, {} U1, {} U2, externalEntityPubChemCompound P
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = P.externalEntityID
             '''.format(type_drug_table, up_table, up_table))

    cursor.execute(query)

    drug2pubchems = {}

    for items in cursor:
        drug2pubchems.setdefault(items[0], set())
        drug2pubchems[items[0]].add(items[1])

    cursor.close()


    for drug in all_drugs:
        if drug not in drug2smiles:
            if drug in drug2pubchems:
                for pubchem in drug2pubchems[drug]:
                    command = 'wget https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/CanonicalSMILES/XML -O out.xml'.format(pubchem)
                    os.system(command)
                    smiles_regex = re.compile("<CanonicalSMILES>(.+)</CanonicalSMILES>")
                    with open('out.xml', 'r') as f:
                        for line in f:
                            m = smiles_regex.search(line)
                            if m:
                                smiles = m.group(1)
                                print(smiles)
                                drug2smiles.setdefault(drug, set())
                                drug2smiles[drug].add(smiles)
                    os.system('rm out.xml')


    print(drug2smiles)

    pickle.dump(drug2smiles, open(output_pickle_file, 'wb'))

    return drug2smiles

def obtain_drug_to_atcs(biana_cnx, unification_protocol, all_drugs, type_drug_name, output_pickle_file):
    """
    Obtains a file containing the ATCs of every drug.
    The type of drug name must be indicated (drugbank, dcdb, name)
    """

    type_name_to_table = {
        'name' : 'externalEntityName',
        'drugbank' : 'externalEntityDrugBankID',
        'dcdb' : 'externalEntityDCDB_drugID',
    }

    type_drug_table = type_name_to_table[type_drug_name] # Obtain the table containing the type of name introduced
    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    # Obtain the ATCs for the drugs
    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query = (''' SELECT A.value FROM {} D, {} U1, {} U2, externalEntityATC A
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = A.externalEntityID AND D.value = %s
             '''.format(type_drug_table, up_table, up_table))

    drug2atcs = {}

    for drug in all_drugs:
        ATCs = set()
        cursor.execute(query, (drug,))
        for items in cursor:
            for result in items:
                ATCs.add(result)
        if len(ATCs) > 0:
            drug2atcs[drug] = ATCs

    print(drug2atcs)
    pickle.dump(drug2atcs, open(output_pickle_file, 'wb'))

    cursor.close()

def obtain_drug_to_side_effects(biana_cnx, unification_protocol, all_drugs, type_drug_name, output_pickle_file):
    """
    Obtains a file containing the side effects associated to every drug.
    The type of drug name must be indicated (drugbank, dcdb, name).
    """

    type_name_to_table = {
        'name' : 'externalEntityName',
        'drugbank' : 'externalEntityDrugBankID',
        'dcdb' : 'externalEntityDCDB_drugID',
    }

    type_drug_table = type_name_to_table[type_drug_name] # Obtain the table containing the type of name introduced
    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    # Obtain the ATCs for the drugs
    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query = (''' SELECT UD.value FROM {} D, {} U1, {} U2, externalEntityPubChemCompound PC, externalEntityRelationParticipant P1, externalEntityRelationParticipant P2, externalEntityRelation R, externalEntityUMLS_diseaseID UD
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = PC.externalEntityID AND PC.externalEntityID = P1.externalEntityID AND D.value = %s
                 AND P1.externalEntityRelationID = R.externalEntityRelationID AND P2.externalEntityRelationID = R.externalEntityRelationID AND P1.externalEntityID != P2.externalEntityID AND R.type = "drug_phenotype_association"
                 AND P2.externalEntityID = UD.externalEntityID
             '''.format(type_drug_table, up_table, up_table))

    drug2side_effects = {}

    for drug in all_drugs:
        side_effects = set()
        cursor.execute(query, (drug,))
        for items in cursor:
            for result in items:
                side_effects.add(result)
        if len(side_effects) > 0:
            drug2side_effects[drug] = side_effects

    print(drug2side_effects)
    pickle.dump(drug2side_effects, open(output_pickle_file, 'wb'))

    cursor.close()

def obtain_drugbank_to_names(biana_cnx, drugs, output_pickle_file):
    """
    Obtains the names associated to DrugbankIDs.
    """

    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query = (''' SELECT N.value 
                 FROM externalEntityDrugBankID D, externalEntityName N 
                 WHERE D.externalEntityID = N.externalEntityID AND (N.type = "unique" OR N.type = "brand") AND D.value = %s
             ''')

    drugbank_to_names = {}
    for drug in drugs:
        cursor.execute(query, (drug,))
        names = set()
        for items in cursor:
            for name in items:
                names.add(name)
        if len(names)>0:
            for name in names:
                drugbank_to_names.setdefault(drug, set())
                drugbank_to_names[drug].add(name.lower())

    cursor.close()

    print(drugbank_to_names)
    pickle.dump(drugbank_to_names, open(output_pickle_file, 'wb'))

    return drugbank_to_names

def obtain_drug_interaction_to_drugs(biana_cnx, output_pickle_file):
    """
    Obtain a dictionary drug_interaction : [drugs]
    """

    cursor = biana_cnx.cursor()

    query = ('''SELECT I.value, D.value from externalEntityRelation R, externalEntityRelationParticipant P, externalEntityDCDB_drugID D, externalEntityDCDB_druginteractionID I
                WHERE R.externalEntityRelationID = P.externalEntityRelationID AND R.type = "interaction" AND P.externalEntityID = D.externalEntityID AND R.externalEntityRelationID = I.externalEntityID
             ''')

    cursor.execute(query)

    drug_int_2_drugs = {}

    for items in cursor:
        drug_int = items[0]
        drug = items[1]
        drug_int_2_drugs.setdefault(drug_int, [])
        drug_int_2_drugs[drug_int].append(drug)

    cursor.close()

    print(drug_int_2_drugs)

    pickle.dump(drug_int_2_drugs, open(output_pickle_file, 'wb'))

    return drug_int_2_drugs

def obtain_drug_interaction_to_info(biana_cnx, output_pickle_file):
    """
    Obtain a dictionary drug_interaction : { 'type' : ... , 'classification' : ... }
    """

    cursor = biana_cnx.cursor()

    query = (''' SELECT value, interactionType, classification from externalEntityDCDB_druginteractionID ''')

    cursor.execute(query)

    drug_int_2_info = {}

    for items in cursor:
        drug_int = items[0]
        type_int = items[1]
        class_int = items[2]
        drug_int_2_info.setdefault(drug_int, {})
        drug_int_2_info[drug_int]['type'] = type_int
        drug_int_2_info[drug_int]['classification'] = class_int

    cursor.close()

    print(drug_int_2_info)

    pickle.dump(drug_int_2_info, open(output_pickle_file, 'wb'))

    return drug_int_2_info

def obtain_dcdb_to_drugbank(biana_cnx, unification_protocol, output_pickle_file):
    """
    Obtain a dictionary {dcdb : drugbank}
    """

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    query = ('''SELECT DC.value, DB.value FROM externalEntityDCDB_drugID DC, {} U1, {} U2, externalEntityDrugBankID DB
                WHERE DC.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = DB.externalEntityID
             '''.format(up_table, up_table))

    cursor = biana_cnx.cursor()
    cursor.execute(query)

    dcdb_to_drugbank = {}

    for items in cursor:
        dcdb = items[0]
        drugbank = items[1]
        dcdb_to_drugbank.setdefault(dcdb, set())
        dcdb_to_drugbank[dcdb].add(drugbank)

    cursor.close()

    print(dcdb_to_drugbank)

    pickle.dump(dcdb_to_drugbank, open(output_pickle_file, 'wb'))

    return dcdb_to_drugbank

def obtain_pubchem_to_drugbank(biana_cnx, unification_protocol, output_pickle_file):
    """
    Obtain a dictionary {pubchem : drugbank}
    """

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    query = ('''SELECT PC.value, DB.value FROM externalEntityPubChemCompound PC, {} U1, {} U2, externalEntityDrugBankID DB
                WHERE PC.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = DB.externalEntityID
             '''.format(up_table, up_table))

    cursor = biana_cnx.cursor()
    cursor.execute(query)

    pubchem_to_drugbank = {}

    for items in cursor:
        pubchem = items[0]
        drugbank = items[1].upper()
        pubchem_to_drugbank.setdefault(pubchem, set())
        pubchem_to_drugbank[pubchem].add(drugbank)

    cursor.close()

    print(pubchem_to_drugbank)

    pickle.dump(pubchem_to_drugbank, open(output_pickle_file, 'wb'))

    return pubchem_to_drugbank

def obtain_target_to_bio_processes(biana_cnx, unification_protocol, all_targets, output_pickle_file):
    """
    Obtains a dictionary containing the targets (in Entrez GeneID) and their corresponding biological processes.
    Record a pickle file.
    """

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query1 = (''' SELECT GO.value FROM externalEntityGeneID G, {} U1, {} U2, externalEntityGO GO
                  WHERE G.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = GO.externalEntityID AND G.value = %s
             '''.format(up_table, up_table))

    query2 = (''' SELECT GO.value FROM externalEntityGO GO, externalEntityGO_type T
                  WHERE GO.externalEntityID = T.externalEntityID AND T.value = 'biological_process' AND GO.value = %s
             ''')

    target_to_bio_processes = {}
    go_to_bio = {}

    num_tar = len(all_targets)
    for target in all_targets:
        print(target)
        GOs = set()
        cursor.execute(query1, (target,))
        for items in cursor:
            for go in items:
                GOs.add(go)
        bio_proc = set()
        for go in GOs:
            if go in go_to_bio:
                for bp in go_to_bio[go]:
                    bio_proc.add(bp)
            else:
                cursor.execute(query2, (go,))
                for items in cursor:
                    for bp in items:
                        bio_proc.add(bp)
                        go_to_bio.setdefault(go, set())
                        go_to_bio[go].add(bp)
        if len(bio_proc) > 0:
            target_to_bio_processes[target] = bio_proc
            print(bio_proc)
        num_tar -= 1
        print('{} left'.format(num_tar))

    cursor.close()

    print(target_to_bio_processes)

    pickle.dump(target_to_bio_processes, open(output_pickle_file, 'wb'))

    return target_to_bio_processes

def obtain_target_to_pathways(biana_cnx, unification_protocol, all_targets, output_pickle_file):
    """
    Obtains a dictionary containing the targets (in Entrez GeneID) and their corresponding pathways.
    Record a pickle file.
    """

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query = (''' SELECT Re.value FROM externalEntityGeneID G, {} U1, {} U2, externalEntityRelationParticipant P1, externalEntityRelationParticipant P2, externalEntityReactome Re
                  WHERE G.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID
                  AND U2.externalEntityID = P1.externalEntityID AND P1.externalEntityRelationID = P2.externalEntityRelationID AND P2.externalEntityID = Re.externalEntityID AND G.value = %s
             '''.format(up_table, up_table))

    target_to_pathways = {}

    num_tar = len(all_targets)
    for target in all_targets:
        print(target)
        pathways = set()
        cursor.execute(query, (target,))
        for items in cursor:
            for path in items:
                pathways.add(path)
        if len(pathways) > 0:
            print(pathways)
            target_to_pathways[target] = pathways
        num_tar -= 1
        print('{} left'.format(num_tar))

    cursor.close()

    print(target_to_pathways)

    pickle.dump(target_to_pathways, open(output_pickle_file, 'wb'))

    return target_to_pathways

def find_drugbank_id_from_name(biana_cnx, unification_protocol, drug_name):
    """
    Obtains the DrugBank ID of the drug from its name, if it is in the database.
    If it is not in the Database, it returns None.
    """

    up_table = return_unification_protocol_table(biana_cnx, unification_protocol)

    cursor = biana_cnx.cursor() # Start cursor to MySQL

    query = (''' SELECT D.value FROM externalEntityName N, {} U1, {} U2, externalEntityDrugBankID D
                 WHERE N.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = D.externalEntityID AND N.value = %s
             '''.format(up_table, up_table))

    cursor.execute(query, (drug_name,))
    drugbank_ids = set()
    for items in cursor:
        for db in items:
            drugbank_ids.add(db)
    cursor.close()

    if len(drugbank_ids) == 0:
        return None
    else:
        return drugbank_ids


