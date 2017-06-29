import matplotlib.pyplot as plt
import mysql.connector


def obtain_database_id(cnx, database_name):

    cursor = cnx.cursor()

    query = ('''SELECT externalDatabaseID FROM externalDatabase
                WHERE DatabaseName = %s
             ''')

    cursor.execute(query, (database_name,))

    result = []

    for db_id in cursor:
        result.append(db_id[0])
    cursor.close()
    print result
    return result[0]


def obtain_DrugCombinationComponentID_list(cnx, database_id):

    cursor = cnx.cursor()

    query = ('''SELECT D.value FROM externalEntityDrugCombinationComponentID D, externalEntity E
                WHERE D.externalEntityID = E.externalEntityID AND E.externalDatabaseID = %s
             ''')

    cursor.execute(query, (database_id,))

    dcc_set = set()

    for dcc_id in cursor:
        dcc_set.add(dcc_id[0])

    cursor.close()

    return dcc_set


def obtain_DrugBankID_list(cnx, database_id):

    cursor = cnx.cursor()

    query = ('''SELECT D.value FROM externalEntityDrugBankID D, externalEntity E
                WHERE D.externalEntityID = E.externalEntityID AND E.externalDatabaseID = %s
             ''')

    cursor.execute(query, (database_id,))

    drugbankid_set = set()

    for drugbankid in cursor:
        drugbankid_set.add(drugbankid[0])

    cursor.close()

    return drugbankid_set


def obtain_drug_external_entity(cnx, drugbankid):

    cursor = cnx.cursor()

    query = ('''SELECT externalEntityID FROM externalEntityDrugBankID
                WHERE value = %s
             ''')

    cursor.execute(query, (drugbankid,))

    result = []

    for ee in cursor:
        result.append(ee[0])
    cursor.close()

    return result[0]


def obtain_num_uniprotentries_per_drug(cnx, external_entity_id, database_id):

    cursor = cnx.cursor()

    query = ('''SELECT COUNT(U.externalEntityID) FROM externalEntityRelationParticipant P1, externalEntityRelationParticipant P2, externalEntityUniprotEntry U, externalEntityDrugBankID D, externalEntity E 
                WHERE P1.externalEntityRelationID = P2.externalEntityRelationID AND P1.externalEntityID != P2.externalEntityID 
                AND P1.externalEntityID = D.externalEntityID AND P2.externalEntityID = U.externalEntityID AND P1.externalEntityID = E.externalEntityID 
                AND E.externalDatabaseID = %s AND D.externalEntityID = %s
             ''')

    cursor.execute(query, (database_id, external_entity_id,))

    result = []

    for num in cursor:
        result.append(num[0])
    cursor.close()

    return result[0]


def obtain_drug_target_interactions(cnx, database_id):

    cursor = cnx.cursor()

    query = ("""SELECT D.value, U.value FROM externalEntity E1, externalEntity E2, externalEntityRelationParticipant P1, externalEntityRelationParticipant P2, externalEntityDrugBankID D, externalEntityUniprotEntry U
       WHERE E1.externalEntityID = P1.externalEntityID AND E2.externalEntityID = P2.externalEntityID AND P1.externalEntityRelationID = P2.externalEntityRelationID AND P1.externalEntityID != P2.externalEntityID  
       AND E1.externalEntityID = D.externalEntityID AND E2.externalEntityID = U.externalEntityID
       AND E1.type = "drug" AND E2.type = "protein" AND E1.externalDatabaseID = %s
    """)

    cursor.execute(query, (database_id,))

    drug2targets = {}

    for result in cursor:
        (drug, target) = result
        print("DRUG {}, TARGET {}".format(drug, target))
        drug2targets.setdefault(drug, [])
        drug2targets[drug].append(target)
    cursor.close()
    print("QUERY FINISHED")
    return drug2targets


def obtain_drug_target_interactions_dcdb(cnx, database_id):

    cursor = cnx.cursor()

    query = ("""SELECT D.value, T.value FROM externalEntity E1, externalEntity E2, externalEntityRelationParticipant P1, externalEntityRelationParticipant P2, externalEntityDrugCombinationComponentID D, externalEntityDrugCombinationTargetID T
       WHERE E1.externalEntityID = P1.externalEntityID AND E2.externalEntityID = P2.externalEntityID AND P1.externalEntityRelationID = P2.externalEntityRelationID AND P1.externalEntityID != P2.externalEntityID  
       AND E1.externalEntityID = D.externalEntityID AND E2.externalEntityID = T.externalEntityID
       AND E1.type = "drug" AND E2.type = "protein" AND E1.externalDatabaseID = %s
    """)

    cursor.execute(query, (database_id,))

    drug2targets = {}

    for result in cursor:
        (drug, target) = result
        print("DRUG {}, TARGET {}".format(drug, target))
        drug2targets.setdefault(drug, [])
        drug2targets[drug].append(target)
    cursor.close()
    print("QUERY FINISHED")

    return drug2targets



cnx = mysql.connector.connect(user='quim', password="",
                              host='localhost',
                              database='test_2017')

############################
### ANALYSIS OF DRUGBANK ###
############################

database_name = 'drugbank'
print("DATABASE: {}".format(database_name))

database_id = obtain_database_id(cnx, database_name)
drugbankid_set = obtain_DrugBankID_list(cnx, database_id)

print("NUMBER OF DRUGS: {}".format(len(drugbankid_set)))

num_targets_per_drug = []

# num = 0
# for drugbankid in drugbankid_set:
#     print("DRUGBANKID {}".format(drugbankid))
#     ee_id = obtain_drug_external_entity(cnx, drugbankid)
#     num_uniprotentries = obtain_num_uniprotentries_per_drug(cnx, ee_id, database_id)
#     print("NUM TARGETS {}".format(num_uniprotentries))
#     num_targets_per_drug.append(num_uniprotentries)
#     num += 1
#     if num % 20 == 0:
#         print("We are at {:.3f}% to finish".format(float(num) / float(len(drugbankid_set)) * 100))
# cnx.close()

drug2targets = obtain_drug_target_interactions(cnx, database_id)
drugs_with_target = set()

for drug in drug2targets:
    num_targets_per_drug.append(len(drug2targets[drug]))
    drugs_with_target.add(drug)

drugs_without_target = drugbankid_set - drugs_with_target
print("NUMBER OF DRUGS WITH TARGET: {}, NUMBER OF DRUGS WITHOUT TARGET: {}".format(drugs_with_target, drugs_without_target))

cnx.close()

for drug in drugs_without_target:
    num_targets_per_drug.append(0)

# Histogram 
n, bins, patches = plt.hist(num_targets_per_drug, 
                            histtype='bar')
plt.title(r'$\mathrm{Number\ of\ targets\ per\ drug\ DrugBank}$')
plt.xlabel('Number of targets')
plt.ylabel('Frequency of drugs')
plt.savefig('num_targets_per_drug_drugbank.png')
#plt.show()

plt.clf()

############################
### ANALYSIS OF DCDB ###
############################

# database_name = 'dcdb'
# print("DATABASE: {}".format(database_name))

# database_id = obtain_database_id(cnx, database_name)
# dcc_set = obtain_DrugCombinationComponentID_list(cnx, database_id)

# print("NUMBER OF DRUGS: {}".format(len(dcc_set)))

# num_targets_per_drug = []

# drug2targets = obtain_drug_target_interactions_dcdb(cnx, database_id)
# drugs_with_target = set()

# for drug in drug2targets:
#     num_targets_per_drug.append(len(drug2targets[drug]))
#     drugs_with_target.add(drug)

# drugs_without_target = dcc_set - drugs_with_target
# print("NUMBER OF DRUGS WITH TARGET: {}, NUMBER OF DRUGS WITHOUT TARGET: {}".format(drugs_with_target, drugs_without_target))

# cnx.close()

# for drug in drugs_without_target:
#     num_targets_per_drug.append(0)

# # Histogram 
# n, bins, patches = plt.hist(num_targets_per_drug, 
#                             histtype='bar')
# plt.title(r'$\mathrm{Number\ of\ targets\ per\ drug\ DrugBank}$')
# plt.xlabel('Number of targets')
# plt.ylabel('Frequency of drugs')
# plt.savefig('num_targets_per_drug_drugbank.png')
# #plt.show()

# plt.clf()

