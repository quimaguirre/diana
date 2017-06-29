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

# Create tables
c.execute('''CREATE TABLE Drugs(
             drugbank_id TEXT PRIMARY KEY, 
             name TEXT NOT NULL
             )''')

c.execute('''CREATE TABLE Targets(
             idTargets TEXT PRIMARY KEY, 
             name TEXT NOT NULL)''')

c.execute('''CREATE TABLE Polypeptides(
             Targets_idTargets TEXT NOT NULL,
             uniprot_accession TEXT NOT NULL,
             gene_id TEXT NOT NULL,
             FOREIGN KEY(Targets_idTargets) REFERENCES Targets(idTargets),
             PRIMARY KEY(Targets_idTargets, uniprot_accession, gene_id)
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
             PRIMARY KEY(Drugs_drugbank_id, Drugs_drugbank_id1)
             )''')



# DrugBank directory
drugbank = "/home/quim/project/DrugBank/full_database.xml"
proof = "drugbank_proof.xml"

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
        #print("The DRUG ID is already in the table")
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
        #print("TARGET NAME: %s. ID: %s"%(target_name, target_id))
        target_set = (target_id, target_name)
        print("TARGET ID: %s NAME: %s" %(target_id, target_name))
        drug_target_set = (drugbank_id, target_id)
        try:
            c.execute('INSERT INTO Targets VALUES (?,?)', target_set)
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
                            if resource[0].text == "GenBank Protein Database":
                                gene_id = external_id.find('identifier').text.upper()
                                print("POLYPEPTIDE GENEID: %s" %(gene_id))
                            elif resource[0].text == "UniProt Accession":
                                uniprot = external_id.find('identifier').text.upper()
                                print("POLYPEPTIDE UNIPROT: %s" %(uniprot))
                    polypeptides_set = (target_id, uniprot, gene_id)
                    try:
                        c.execute('INSERT INTO Polypeptides VALUES (?,?,?)', polypeptides_set)
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
