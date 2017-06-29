import xml.etree.ElementTree as ET
import sqlite3
import sys
import time


# Start marker for time measure
start = time.time()

# Connection to SQLite db
conn = sqlite3.connect('drugbank.db')

# Create a Cursor object in order to call its execute() method to perform SQL commands
c = conn.cursor()

# Activate the foreign key restriction
c.execute("PRAGMA foreign_keys = ON")

# Create tables
c.execute('''CREATE TABLE Drugs(
             drugbank_id TEXT PRIMARY KEY, 
             name TEXT NOT NULL UNIQUE
             )''')

c.execute('''CREATE TABLE Drugs_has_PubChem(
             Drugs_drugbank_id TEXT NOT NULL, 
             pubchem_compound TEXT NOT NULL, 
             FOREIGN KEY(Drugs_drugbank_id) REFERENCES Drugs(drugbank_id), 
             PRIMARY KEY(pubchem_compound)
             )''')

c.execute('''CREATE TABLE Targets(
             idTargets TEXT PRIMARY KEY, 
             name TEXT NOT NULL,
             type TEXT NOT NULL
             )''')

c.execute('''CREATE TABLE Polypeptides(
             Targets_idTargets TEXT NOT NULL,
             uniprot_accession TEXT NOT NULL,
             FOREIGN KEY(Targets_idTargets) REFERENCES Targets(idTargets),
             PRIMARY KEY(Targets_idTargets, uniprot_accession)
             )''')

c.execute('''CREATE TABLE Drugs_has_Targets(
             Drugs_drugbank_id TEXT NOT NULL, 
             Targets_idTargets TEXT NOT NULL, 
             FOREIGN KEY(Drugs_drugbank_id) REFERENCES Drugs(drugbank_id), 
             FOREIGN KEY(Targets_idTargets) REFERENCES Targets(idTargets),
             PRIMARY KEY(Drugs_drugbank_id, Targets_idTargets)
             )''')

c.execute('''CREATE TABLE Interactions(
             description TEXT NOT NULL, 
             Drugs_drugbank_id TEXT NOT NULL, 
             Drugs_drugbank_id1 TEXT NOT NULL, 
             FOREIGN KEY(Drugs_drugbank_id) REFERENCES Drugs(drugbank_id), 
             FOREIGN KEY(Drugs_drugbank_id1) REFERENCES Drugs(drugbank_id),
             PRIMARY KEY(description, Drugs_drugbank_id, Drugs_drugbank_id1)
             )''')



# DrugBank directory
drugbank = "/home/quim/Databases/DrugBank/full_database_old_parser.xml"
drugbank_old = "/home/quim/project/DrugBank/old/full_database_2016_07.xml"
proof = "/home/quim/project/DrugBank/drugbank_proof.xml"
proof2 = "/home/quim/project/DrugBank/drugbank_proof2.xml"

# ------------------------------------------------------------------------------
# NOTE: REMBEMBER TO MODIFY THE ROOT ELEMENT SO THAT IT DOES NOT HAVE ATTRIBUTES
# Example: 
# <drugbank xmlns="http://www.drugbank.ca" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.drugbank.ca http://www.drugbank.ca/docs/drugbank.xsd" version="5.0" exported-on="2016-10-01">
# It must be changed to:
# <drugbank>
# ------------------------------------------------------------------------------

# Parsing of drugbank
tree = ET.parse(drugbank)
root = tree.getroot()
#drug_num = 1


for drug in root.findall('drug'):
    # Control variables
    #drug_num_saved = None
    #save = False
    # Inserting drugbank id and drug name into the table Drugs
    drugbank_id = drug.find('drugbank-id').text
    name = drug.find('name').text.lower()
    drug_set = (drugbank_id, name)
    print("MAIN DRUG: ID: %s NAME: %s" %(drugbank_id, name))

    try:
        c.execute('INSERT INTO Drugs VALUES (?,?)', drug_set)
    except sqlite3.IntegrityError:
        print("The DRUG name %s with id %s is already in the table" %(name, drugbank_id))
        pass


    drug_external_identifiers = drug.find('external-identifiers')

    for drug_external_id in drug_external_identifiers:
        resource = drug_external_id.findall('resource')
        if resource != None:
            if resource[0].text == "PubChem Compound":
                pubchem = drug_external_id.find('identifier').text.upper()

                print("DrugBankID: %s Name: %s PUBCHEMCOMPOUND: %s" %(drugbank_id, name, pubchem))
                drug_pubchem_set = (drugbank_id, pubchem)

                try:
                    c.execute('INSERT INTO Drugs_has_PubChem VALUES (?,?)', drug_pubchem_set)
                except sqlite3.IntegrityError:
                    print("The DRUG name %s with id %s and pubchem compound %s is already in the table" %(name, drugbank_id, pubchem))
                    pass

    # c.execute('SELECT idDrugs FROM Drugs WHERE drugbank_id = (?)', (drugbank_id,))
    # command = c.fetchall()
    # if command == []:
    #     c.execute('INSERT INTO Drugs VALUES (?,?,?)', drug_set)
    #     print("MAIN DRUG: NUM: %s NAME: %s" %(drug_num, name))
    # else:
    #     command = command[0][0]
    #     drug_num_saved = drug_num
    #     drug_num = command
    #     print("MAIN DRUG: NUM: %s NAME: %s" %(command, name))

    targets = drug.find('targets')

    # Search for targets id and name
    for target in targets:
        target_id = target.find('id').text.upper()
        target_name = target.find('name').text.lower()
        target_type = 'target'
        #print("TARGET NAME: %s. ID: %s"%(target_name, target_id))
        target_set = (target_id, target_name, target_type)
        print("TARGET ID: %s NAME: %s TYPE: %s" %(target_id, target_name, target_type))
        drug_target_set = (drugbank_id, target_id)
        try:
            c.execute('INSERT INTO Targets VALUES (?,?,?)', target_set)
        except sqlite3.IntegrityError:
            #print("The TARGET was already in the table")
            pass
        try:
            c.execute('INSERT INTO Drugs_has_Targets VALUES (?,?)', drug_target_set)
        except sqlite3.IntegrityError:
            pass

        # Search for the polypeptides of the target and store their uniprot accessions
        for polypeptides in target.findall('polypeptide'):
            for polypeptide in polypeptides:
                for element in polypeptide.iter('external-identifiers'):
                    for external_id in element.findall('external-identifier'):
                        resource = external_id.findall('resource')
                        if resource != None:
                            if resource[0].text == "UniProt Accession":
                                uniprot = external_id.find('identifier').text.upper()
                                print("POLYPEPTIDE UNIPROT: %s" %(uniprot))
                    polypeptides_set = (target_id, uniprot)
                    try:
                        c.execute('INSERT INTO Polypeptides VALUES (?,?)', polypeptides_set)
                    except sqlite3.IntegrityError:
                        pass

    enzymes = drug.find('enzymes')

    # Search for targets id and name
    for enzyme in enzymes:
        target_id = enzyme.find('id').text.upper()
        target_name = enzyme.find('name').text.lower()
        target_type = 'enzyme'
        #print("TARGET NAME: %s. ID: %s"%(target_name, target_id))
        target_set = (target_id, target_name, target_type)
        print("TARGET ID: %s NAME: %s TYPE: %s" %(target_id, target_name, target_type))
        drug_target_set = (drugbank_id, target_id)
        try:
            c.execute('INSERT INTO Targets VALUES (?,?,?)', target_set)
        except sqlite3.IntegrityError:
            #print("The TARGET was already in the table")
            pass
        try:
            c.execute('INSERT INTO Drugs_has_Targets VALUES (?,?)', drug_target_set)
        except sqlite3.IntegrityError:
            pass

        # Search for the polypeptides of the enzyme and store their uniprot accessions
        for polypeptides in enzyme.findall('polypeptide'):
            for polypeptide in polypeptides:
                for element in polypeptide.iter('external-identifiers'):
                    for external_id in element.findall('external-identifier'):
                        resource = external_id.findall('resource')
                        if resource != None:
                            if resource[0].text == "UniProt Accession":
                                uniprot = external_id.find('identifier').text.upper()
                                print("POLYPEPTIDE UNIPROT: %s" %(uniprot))
                    polypeptides_set = (target_id, uniprot)
                    try:
                        c.execute('INSERT INTO Polypeptides VALUES (?,?)', polypeptides_set)
                    except sqlite3.IntegrityError:
                        pass

    carriers = drug.find('carriers')

    # Search for targets id and name
    for carrier in carriers:
        target_id = carrier.find('id').text.upper()
        target_name = carrier.find('name').text.lower()
        target_type = 'carrier'
        #print("TARGET NAME: %s. ID: %s"%(target_name, target_id))
        target_set = (target_id, target_name, target_type)
        print("TARGET ID: %s NAME: %s TYPE: %s" %(target_id, target_name, target_type))
        drug_target_set = (drugbank_id, target_id)
        try:
            c.execute('INSERT INTO Targets VALUES (?,?,?)', target_set)
        except sqlite3.IntegrityError:
            #print("The TARGET was already in the table")
            pass
        try:
            c.execute('INSERT INTO Drugs_has_Targets VALUES (?,?)', drug_target_set)
        except sqlite3.IntegrityError:
            pass

        # Search for the polypeptides of the carrier and store their uniprot accessions
        for polypeptides in carrier.findall('polypeptide'):
            for polypeptide in polypeptides:
                for element in polypeptide.iter('external-identifiers'):
                    for external_id in element.findall('external-identifier'):
                        resource = external_id.findall('resource')
                        if resource != None:
                            if resource[0].text == "UniProt Accession":
                                uniprot = external_id.find('identifier').text.upper()
                                print("POLYPEPTIDE UNIPROT: %s" %(uniprot))
                    polypeptides_set = (target_id, uniprot)
                    try:
                        c.execute('INSERT INTO Polypeptides VALUES (?,?)', polypeptides_set)
                    except sqlite3.IntegrityError:
                        pass


    transporters = drug.find('transporters')

    # Search for targets id and name
    for transporter in transporters:
        target_id = transporter.find('id').text.upper()
        target_name = transporter.find('name').text.lower()
        target_type = 'transporter'
        #print("TARGET NAME: %s. ID: %s"%(target_name, target_id))
        target_set = (target_id, target_name, target_type)
        print("TARGET ID: %s NAME: %s TYPE: %s" %(target_id, target_name, target_type))
        drug_target_set = (drugbank_id, target_id)
        try:
            c.execute('INSERT INTO Targets VALUES (?,?,?)', target_set)
        except sqlite3.IntegrityError:
            #print("The TARGET was already in the table")
            pass
        try:
            c.execute('INSERT INTO Drugs_has_Targets VALUES (?,?)', drug_target_set)
        except sqlite3.IntegrityError:
            pass

        # Search for the polypeptides of the transporter and store their uniprot accessions
        for polypeptides in transporter.findall('polypeptide'):
            for polypeptide in polypeptides:
                for element in polypeptide.iter('external-identifiers'):
                    for external_id in element.findall('external-identifier'):
                        resource = external_id.findall('resource')
                        if resource != None:
                            if resource[0].text == "UniProt Accession":
                                uniprot = external_id.find('identifier').text.upper()
                                print("POLYPEPTIDE UNIPROT: %s" %(uniprot))
                    polypeptides_set = (target_id, uniprot)
                    try:
                        c.execute('INSERT INTO Polypeptides VALUES (?,?)', polypeptides_set)
                    except sqlite3.IntegrityError:
                        pass

    # i = 0
    # for element in targets.iter('external-identifiers'):
    #     for external_id in element.findall('external-identifier'):
    #         resource = external_id.findall('resource')
    #         if resource != None:
    #             if resource[0].text == "UniProt Accession":
    #                 identifier = external_id.find('identifier').text.upper()
    #                 try:
    #                     target_set = (target_num, targets_list[i], identifier)
    #                     drug_target_set = (drug_num, target_num)
    #                     c.execute('INSERT INTO Targets VALUES (?,?,?)', target_set)
    #                     c.execute('INSERT INTO Drugs_has_Targets VALUES (?,?)', drug_target_set)
    #                     i += 1
    #                     target_num += 1
    #                 except:
    #                     print("There is a problem with: %s, in drug %s" %(identifier, drug))
    #                     continue

    # if drug_num_saved == None:
    #     interacting_drug_num = drug_num
    # else:
    #     interacting_drug_num = drug_num_saved

    for element in drug.iter('drug-interactions'):
        for interaction in element.findall('drug-interaction'):
            interaction_drugbank_id = interaction.find('drugbank-id').text
            interaction_name = interaction.find('name').text.lower()
            interaction_description = interaction.find('description').text
            try:
                c.execute('INSERT INTO Drugs VALUES (?,?)', (interaction_drugbank_id, interaction_name))
                #print("Interaction DRUG: ID: %s NAME: %s" %(interaction_drugbank_id, interaction_name))
            except sqlite3.IntegrityError:
                #print("This drug_interaction was already in DRUGS")
                pass
            try:
                c.execute('INSERT INTO Interactions VALUES (?,?,?)', (interaction_description, drugbank_id, interaction_drugbank_id))
            except:
                print("ERROR WHEN INSERTING THE INTERACTION %s" %(interaction_name))
                sys.exit(10)

            # c.execute('SELECT idDrugs FROM Drugs WHERE drugbank_id = (?)', (interaction_drugbank_id,))
            # command = c.fetchall()
            # # If drug is already in Drugs table, just insert the interaction
            # if command != []:
            #     command = command[0][0]
            #     c.execute('INSERT INTO Interactions VALUES (?,?,?)', (interaction_description, drug_num, command))
            # # If drug is not in Drugs table, add the drug in Drugs table and insert the interaction
            # else:
            #     if drug_num_saved == None:
            #         interacting_drug_num += 1
            #         c.execute('INSERT INTO Drugs VALUES (?,?,?)', (interaction_drugbank_id, interaction_name))
            #         c.execute('INSERT INTO Interactions VALUES (?,?,?)', (interaction_description, drugbank_id, interaction_drugbank_id))
            #         #print("NUM: %s NAME: %s" %(interacting_drug_num, interaction_name))
            #     else:
            #         c.execute('INSERT INTO Drugs VALUES (?,?,?)', (interacting_drug_num, interaction_name, interaction_drugbank_id))
            #         c.execute('INSERT INTO Interactions VALUES (?,?,?)', (interaction_description, drug_num, interacting_drug_num))
            #         #print("DRUG SAVED: NUM: %s NAME: %s" %(interacting_drug_num, interaction_name))
            #         interacting_drug_num += 1
            #         save = True
            #     #print(interacting_drug_num)
    # if save == False:
    #     drug_num = interacting_drug_num + 1
    # else:
    #     drug_num = interacting_drug_num



# Save (commit) the changes
conn.commit()

# We can also close the connection if we are done with it.
# Just be sure any changes have been committed or they will be lost.
conn.close()

# End marker for time
end = time.time()
print ("  DIANA INFO:\tThe time of execution of the parsing is: %0.3f sec. or %0.3f minutes.\n" %(end - start, (end - start) / 60))
