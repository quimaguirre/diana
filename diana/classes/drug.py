import os, sys, re
import hashlib
import network_analysis


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
        self.targets = []
        self.pfams = []
        self.type_id = None
        self.type_id_to_table = {
            'geneid' : 'externalEntityGeneID',
            'uniprotentry' : 'externalEntityUniprotEntry',
            'uniprotaccession' : 'externalEntityUniprotAccession',
        }

    ###########
    # METHODS #
    ###########

    def obtain_targets_from_file(self, targets_file, type_id):
        """
        Obtains the targets from an input file and stores them into a list.
        The file must contain the names of the targets separated by new lines.
        The type of ID of the targets must be specified.
        """
        self.type_id = type_id.lower() # Annotate the type of ID of the targets
        with open(targets_file, 'r') as targets_file_fd:
            for line in targets_file_fd:
                self.targets.append(line.strip())

        # Check if the number of targets provided is sufficient for the analysis
        if len(self.targets) < 3:
            raise InsufficientTargets(self.targets)
        return

    def obtain_targets_from_BIANA(self, biana_cnx, type_id, unification_protocol):
        """
        Obtains the targets from BIANA database using as query the drug name.
        The type of ID of the targets must be specified.
        "biana_cnx" parameter stands for the variable containing the connexion to the MySQL database.
        """
        self.type_id = type_id.lower() # Annotate the type of ID of the targets

        type_id_table = self.return_targets_biana_table(self.type_id) # Obtain the table containing the type of ID introduced

        up_table = self.return_unification_protocol_table(biana_cnx, unification_protocol)

        cursor = biana_cnx.cursor() # Start cursor to MySQL

        query1 = (''' SELECT externalEntityID FROM externalEntityName D WHERE value = %s
                 ''')
        query2 = (''' SELECT G.value FROM externalEntity E1, {} U1, {} U2, externalEntity E2, externalEntityRelationParticipant R2, externalEntityRelationParticipant R3, {} U3, {} U4, {} G
                     WHERE E1.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = E2.externalEntityID AND E2.type = 'drug' AND E2.externalEntityID = R2.externalEntityID AND R2.externalEntityRelationID = R3.externalEntityRelationID
                     AND R3.externalEntityID = U3.externalEntityID AND U3.userEntityID = U4.userEntityID AND U4.externalEntityID = G.externalEntityID AND E1.externalEntityID = %s
                 '''.format(up_table, up_table, up_table, up_table, type_id_table))

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
            raise DrugNameNotFound(self.drug_name)
        # Why in two steps? 
        # Because as the table "externalEntityName" is too large, do one complex command can be very time-consuming
        # It is better to split the search in two commands

        cursor.close()

        self.targets = list(geneids)

        # Check if the number of targets provided is sufficient for the analysis
        if len(self.targets) < 3:
            raise InsufficientTargets(self.targets)

        return

    def return_targets_biana_table(self, type_id):
        """
        Returns the table in BIANA where the annotations of the type of ID 
        introduced are stored.
        """
        if type_id in self.type_id_to_table:
            return self.type_id_to_table[type_id]
        else:
            raise IncorrectTypeID(type_id, self.type_id_to_table)

    def return_unification_protocol_table(self, biana_cnx, unification_protocol):
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


    def obtain_pfams_from_file(self, pfam_file):
        """
        Obtains the pfams from an input file and stores them into a list.
        The file must contain the names of the pfams separated by new lines.
        """
        with open(pfam_file, 'r') as pfam_file_fd:
            for line in pfam_file_fd:
                self.pfams.append(line.strip())
        return

    def obtain_pfams_from_targets(self, biana_cnx, output_file, unification_protocol):
        """
        Obtains the pfams from BIANA database using as query the targets.
        "biana_cnx" parameter stands for the variable containing the connexion to the MySQL database.
        Stores the PFAMs found in an output file
        """

        type_id_table = self.return_targets_biana_table(self.type_id) # Obtain the table containing the type of ID introduced

        up_table = self.return_unification_protocol_table(biana_cnx, unification_protocol)

        query = (''' SELECT P.value FROM {} G, {} U1, {} U2, externalEntityPFAM P 
                     WHERE G.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = P.externalEntityID AND G.value = %s
                 '''.format(type_id_table, up_table, up_table))

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




class InsufficientTargets(Exception):
    """
    Exception raised when the number of targets is below 3.
    This exception is raised because the analyses of GUILD with less than 3
    targets are not reliable.
    """
    def __init__(self, targets):
        self.targets = targets

    def __str__(self):
        return 'The number of targets provided ({}) is insufficient.\nGUILD must have at least 3 targets to run a reliable analysis.\n'.format(len(self.targets))

class DrugNameNotFound(Exception):
    """
    Exception raised when the drug name is not found in BIANA.
    """
    def __init__(self, drug_name):
        self.drug_name = drug_name

    def __str__(self):
        return 'The drug name ({}) has not been found in BIANA.\nTherefore, any target could be found. Please, introduce another name or the targets of the drug.\n'.format(self.drug_name)

class IncorrectTypeID(Exception):
    """
    Exception that raises when a type of IDs of the proteins is not admitted for
    the program.
    """
    def __init__(self, type_id, type_id_to_table):
        self.type_id = type_id
        self.type_id_to_table = type_id_to_table

    def __str__(self):
        return 'The initial type of IDs of the proteins ({}) is not admitted.\nThe types of ID admitted in DIANA are: {}\n'.format(self.type_id, ', '.join(self.type_id_to_table.keys()))




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

def create_targets_file(targets, file_name):
    """
    Creates a targets file, containing the targets separated by new line characters
    """
    with open(file_name, 'w') as fw:
        for target in targets:
            fw.write('{}\n'.format(target))
    return



