import mysql.connector
import cPickle
import os
import re

def obtain_all_geneids_in_sif_file(sif_file):
    """Creates a list with all the geneids in a sif file"""

    fn = open(sif_file,"r")
    sif_file_geneids = set()

    for line in fn:
        words = line.split()
        sif_file_geneids.add(words[0])
        sif_file_geneids.add(words[2])

    fn.close()
    return sif_file_geneids


def obtain_atc_for_dcdb_drugs(cnx):
    """
    Obtain the ATC code for each DCDB drug
    Returns ATC lvl 1 and lvl2
    Example: R03AB02 --> lvl 1: R, lvl 2: R03
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, A.value FROM externalEntityDCDB_drugID D, userEntityUnification_protocol_2 U1, userEntityUnification_protocol_2 U2, externalEntityATC A
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = A.externalEntityID
             ''')

    cursor.execute(query)

    dcdb2atclvl1 = {}
    dcdb2atclvl2 = {}

    for items in cursor:
        dcdb = items[0]
        atc = items[1]
        atc_lvl1 = atc[:1]
        atc_lvl2 = atc[:3]
        dcdb2atclvl1.setdefault(dcdb, set())
        dcdb2atclvl1[dcdb].add(atc_lvl1)
        dcdb2atclvl2.setdefault(dcdb, set())
        dcdb2atclvl2[dcdb].add(atc_lvl2)

    cursor.close()

    return dcdb2atclvl1, dcdb2atclvl2


def obtain_DCDB_drugs(cnx):
    """
    Obtain the ID, external entity and the user entity in which belongs for every DCDB drug
    Store it in a dictionary "external entity" => [ "user entity", "dcdb drug ID" ]
    """

    cursor = cnx.cursor()

    query = (''' SELECT E.externalEntityID, U.userEntityID, D.value FROM externalEntity E, userEntityUnification_protocol_2 U, externalEntityDCDB_drugID D
                 WHERE E.externalEntityID = U.externalEntityID AND E.externalEntityID = D.externalEntityID AND E.externalDatabaseID = 17 AND E.type = 'drug'
             ''')

    cursor.execute(query)

    ee2values = {}

    for items in cursor:
        ee2values.setdefault(items[0], [])
        ee2values[items[0]] = [items[1], items[2]]

    cursor.close()

    return ee2values


def obtain_targets(cnx, user_entity):
    """
    Obtain the user entities of all the targets from a user entity that corresponds to a drug
    """

    cursor = cnx.cursor()

    query = ('''SELECT U2.userEntityID FROM userEntityUnification_protocol_2 U, userEntityUnification_protocol_2 U2, externalEntity E1, externalEntity E2, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2
             WHERE U.externalEntityID = E1.externalEntityID AND E1.externalEntityID = R1.externalEntityID AND R1.externalEntityRelationID = R2.externalEntityRelationID AND R2.externalEntityID = E2.externalEntityID AND E2.type = "protein" AND E2.externalEntityID = U2.externalEntityID AND U.userEntityID = %s
             ''')

    cursor.execute(query, (user_entity,))

    user_entities = set()

    for items in cursor:
        ue = items[0]
        user_entities.add(ue)

    cursor.close()

    return user_entities


def obtain_geneids(cnx, user_entity):
    """
    Obtain the geneids of a user entity corresponding to a protein
    """

    cursor = cnx.cursor()

    query = (''' SELECT G.value FROM userEntityUnification_protocol_2 U, externalEntity E, externalEntityGeneID G 
             WHERE U.externalEntityID = E.externalEntityID AND E.externalEntityID = G.externalEntityID AND U.userEntityID = %s
             ''')

    cursor.execute(query, (user_entity,))

    geneids = set()

    for items in cursor:
        for geneid in items:
            geneids.add(geneid)

    cursor.close()

    return geneids


def obtain_geneids_from_dcdbid(cnx, dcdbid):
    """
    Obtain the geneids of a user entity corresponding to a protein
    """

    cursor = cnx.cursor()

    query = (''' SELECT G.value FROM externalEntityDCDB_drugID D, externalEntityRelationParticipant P1, externalEntityRelationParticipant P2, externalEntity E2, userEntityUnification_protocol_2 U1, userEntityUnification_protocol_2 U2, externalEntityGeneID G 
             WHERE D.externalEntityID = P1.externalEntityID AND P1.externalEntityRelationID = P2.externalEntityRelationID AND P2.externalEntityID = E2.externalEntityID AND E2.type = "protein" AND E2.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = G.externalEntityID AND D.value = %s
             ''')

    cursor.execute(query, (dcdbid,))

    geneids = set()

    for items in cursor:
        for geneid in items:
            geneids.add(geneid)

    cursor.close()

    return geneids


def obtain_pfams(cnx, geneid):
    """
    Obtain the geneids of a user entity corresponding to a protein
    """

    cursor = cnx.cursor()

    query = (''' SELECT P.value FROM externalEntityGeneID G, userEntityUnification_protocol_2 U1, userEntityUnification_protocol_2 U2, externalEntityPFAM P 
             WHERE G.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = P.externalEntityID AND G.value = %s
             ''')

    cursor.execute(query, (geneid,))

    pfams = set()

    for items in cursor:
        for pfam in items:
            pfams.add(pfam.upper())

    cursor.close()

    return pfams

def obtain_pubchem(cnx):
    """
    Obtain the SMILES code for each DCDB drug
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, P.value FROM externalEntityDCDB_drugID D, userEntityUnification_protocol_2 U1, userEntityUnification_protocol_2 U2, externalEntityPubChemCompound P
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = P.externalEntityID
             ''')

    cursor.execute(query)

    dcdb2pubchem = {}

    for items in cursor:
        dcdb2pubchem.setdefault(items[0], set())
        dcdb2pubchem[items[0]].add(items[1])

    cursor.close()

    return dcdb2pubchem

def obtain_smiles(cnx):
    """
    Obtain the SMILES code for each DCDB drug
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, S.value FROM externalEntityDCDB_drugID D, userEntityUnification_protocol_2 U1, userEntityUnification_protocol_2 U2, externalEntitySMILES S
                 WHERE D.externalEntityID = U1.externalEntityID AND U1.userEntityID = U2.userEntityID AND U2.externalEntityID = S.externalEntityID
             ''')

    cursor.execute(query)

    dcdbid2smiles = {}

    for items in cursor:
        dcdbid2smiles[items[0]] = items[1]

    cursor.close()

    return dcdbid2smiles

def obtain_drug_combinations(cnx):
    """
    Obtain the drug combinations in DCDB
    """

    cursor = cnx.cursor()

    query = ('''SELECT C.value, D.value FROM externalEntityRelation R, externalEntityRelationParticipant P, externalEntityDCDB_drugID D, externalEntityDCDB_combinationID C 
                WHERE R.externalEntityRelationID = P.externalEntityRelationID AND R.type = "drug_combination" AND P.externalEntityID = D.externalEntityID AND R.externalEntityRelationID = C.externalEntityID
             ''')

    cursor.execute(query)

    drug_comb_2_drugs = {}

    for items in cursor:
        drug_comb = items[0]
        drug = items[1]
        drug_comb_2_drugs.setdefault(drug_comb, [])
        drug_comb_2_drugs[drug_comb].append(drug)

    cursor.close()

    return drug_comb_2_drugs

def obtain_drug_interaction_to_drugs(cnx):
    """
    Obtain a dictionary drug_interaction : [drugs]
    """

    cursor = cnx.cursor()

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

    return drug_int_2_drugs

def obtain_drug_interaction_to_info(cnx):
    """
    Obtain a dictionary drug_interaction : { 'type' : ... , 'classification' : ... }
    """

    cursor = cnx.cursor()

    query = ('''SELECT value, interactionType, classification from externalEntityDCDB_druginteractionID
             ''')

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

    return drug_int_2_info



cnx = mysql.connector.connect(user='quim', password="",
                              host='localhost',
                              database='BIANA_JAN_2017')

sif_file = "sif/human_eAFF_geneid_2017.sif"
sif_file_geneids = obtain_all_geneids_in_sif_file(sif_file)
ee2values = obtain_DCDB_drugs(cnx)
dcdb2targets = {}
targets = set()

#### Obtain the ids of the DCDB drugs and the GeneIDs of the targets ####

for ee in ee2values:
    ue_drug, dcdbid = ee2values[ee]
    print("External entity: {}    User entity: {}    Drug ID: {}".format(ee, ue_drug, dcdbid))

    geneids = obtain_geneids_from_dcdbid(cnx, dcdbid)
    geneids_in_sif = set()
    for geneid in geneids:
        if str(geneid) in sif_file_geneids:
            geneids_in_sif.add(geneid)
    if len(geneids_in_sif) > 0:
        dcdb2targets.setdefault(dcdbid.upper(), set())
        print("Geneids for {}: {}".format(dcdbid, geneids_in_sif))
        for geneid in geneids_in_sif:
            dcdb2targets[dcdbid].add(geneid)
            targets.add(geneid)
    else:
        print("No geneids found for {}".format(dcdbid))

print("Number of DCDB ids with at least 1 GeneID: {}".format(len(dcdb2targets)))

#### Dump the data structure in a pickle file ####

dump_file = "toolbox/"+"dcdb2targets.pcl"

cPickle.dump(dcdb2targets, open(dump_file, 'w')) 

# dcdb2targets = {}
# print(dcdb2targets)

# dcdb2targets = cPickle.load(open(dump_file))

# print(dcdb2targets)


#### Obtain the PFAM families of the targets found ####

geneid2pfam = {}

for target in targets:
    pfams = obtain_pfams(cnx, target)
    geneid2pfam.setdefault(target, set())
    geneid2pfam[target] = pfams
    #print("PFAM for {}: {}".format(target, geneid2pfam[target]))

#print(geneid2pfam)

pfam_file = "toolbox/"+"geneid2pfam.pcl"

cPickle.dump(geneid2pfam, open(pfam_file, 'w')) 


#### Obtain all the drug combinations ####

drug_comb_2_drugs = obtain_drug_combinations(cnx)

#print(drug_comb_2_drugs)

drug_comb_file = "toolbox/"+"drug_comb_2_drugs.pcl"

cPickle.dump(drug_comb_2_drugs, open(drug_comb_file, 'w'))


#### Obtain all the drug interactions ####

# drug_int_2_drugs = obtain_drug_interaction_to_drugs(cnx)

# drug_int_2_info = obtain_drug_interaction_to_info(cnx)

# print(drug_int_2_drugs)

# #print(drug_int_2_info)

# drugint_drug_file = "toolbox/"+"drug_int_2_drugs.pcl"
# cPickle.dump(drug_int_2_drugs, open(drugint_drug_file, 'w'))

# drugint_info_file = "toolbox/"+"drug_int_2_info.pcl"
# cPickle.dump(drug_int_2_info, open(drugint_info_file, 'w'))


cnx.close()


cnx2 = mysql.connector.connect(user='quim', password="",
                              host='localhost',
                              database='test_2017')


#### Obtain ATC groups for drugs ####

dcdb2atc_lvl1, dcdb2atc_lvl2 = obtain_atc_for_dcdb_drugs(cnx2)
#print(dcdb2atc_lvl1)
#print(dcdb2atc_lvl2)
atc1_file = "toolbox/"+"dcdb2atc_lvl1.pcl"
cPickle.dump(dcdb2atc_lvl1, open(atc1_file, 'w')) 
atc2_file = "toolbox/"+"dcdb2atc_lvl2.pcl"
cPickle.dump(dcdb2atc_lvl2, open(atc2_file, 'w')) 


#### Obtain the SMILES of the drugs ####

dcdb2smiles = obtain_smiles(cnx2)
#print(dcdb2smiles)

dcdb2pubchem = obtain_pubchem(cnx2)
#print(dcdb2pubchem)

cnx2.close()

for dcdb in dcdb2targets:
    if dcdb not in dcdb2smiles:
        if dcdb in dcdb2pubchem:
            for pubchem in dcdb2pubchem[dcdb]:
                command = 'wget https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/CanonicalSMILES/XML -O out.xml'.format(pubchem)
                os.system(command)
                smiles_regex = re.compile("<CanonicalSMILES>(.+)</CanonicalSMILES>")
                f = open('out.xml', 'r')
                for line in f:
                    m = smiles_regex.search(line)
                    if m:
                        smiles = m.group(1)
                        print(smiles)
                        dcdb2smiles[dcdb] = smiles
                f.close()
                os.system('rm out.xml')

print(dcdb2smiles)
smiles_file = "toolbox/"+"dcdb2smiles.pcl"
cPickle.dump(dcdb2smiles, open(smiles_file, 'w')) 

