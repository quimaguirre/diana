import sqlite3
import random
import diana
import toolbox.drugbank_searcher as dbs


def process_drug_name_for_design_experiment(drug):
    """
    Processes the drug name provided in order to know if it is a conventional drug name or a
    DrugBank ID. Returns the strings 'drugbank_id' or 'name' according to the result
    """

    if "'" in drug or "(" in drug or ")" in drug:
        drug = '"%s"' %(drug)
        return drug

    processed_name = drug.split(" ")

    if len(processed_name) > 1:
        drug = '"%s"' %(drug)
        return drug
    else:
        return drug


#------------#
#   INPUTS   #
#------------#

# The minimum number of targets that we want our drugs to have
num_target_threshold = 15

# The name of the final bash file to run the experiment
bash_file = "diana_100_analyses.sh"

# The database you want to use for the analysis
database = "BIANA_MARCH_2013"

# The SIF file needed
sif_file = "sif/human_Y2H_geneid.sif"

# The translation file needed
translation_file = "translation_geneid_Y2H_2013.trans"

#-----------------------------------------------------------#
#   FILTER THE DRUGS OF THE DATABASE BY NUMBER OF TARGETS   #
#-----------------------------------------------------------#

# Connection to SQLite db
conn = sqlite3.connect('drugbank.db')

# Create a Cursor object in order to call its execute() method to perform SQL commands
c = conn.cursor()

# Define the sql query inside a variable, introducing the minimum number of targets as variable
sql_query = '''SELECT d.name, COUNT(dht.Drugs_drugbank_id) as num_targets
                  FROM Drugs d
                    INNER JOIN Drugs_has_Targets dht
                      ON dht.Drugs_drugbank_id = d.drugbank_id
                    INNER JOIN Targets t
                      ON  t.idTargets = dht.Targets_idTargets
                  GROUP BY dht.Drugs_drugbank_id
                  HAVING num_targets >= %d
        ''' %(num_target_threshold)

c.execute(sql_query)

list_of_drugs = dbs.result_to_list(c.fetchall())

# Close the connection
conn.close()

print ("The number of drugs with more than %d targets is: %d" %(num_target_threshold-1, len(list_of_drugs)))


#-----------------------------------------------#
#   CHECK IF THERE IS A MINIMUM NUMBER OF DDI   #
#-----------------------------------------------#

list_of_drugs2 = list(list_of_drugs)

ddi = set()
no_ddi = set()
ddi_pairs = set()
no_ddi_pairs = set()

num_ddi = 0
num_no_ddi = 0

n = 0
while (n < len(list_of_drugs)):
    i = 0
    while (i < len(list_of_drugs2)):
        drug1 = list_of_drugs[n]
        drug2 = list_of_drugs2[i]
        if drug1 == drug2:
            i+=1
            continue


        ddi_name1 = "%s---%s"%(drug1, drug2)
        ddi_name2 = "%s---%s"%(drug2, drug1)
        print("COMPARING: %s vs. %s" %(drug1, drug2))

        # Retrieve the description of the interaction between the two drugs if there is interaction
        # The description could be different depending on the order of the drugs, so the procedure is 
        # performed twice changing the order of the drugs
        description1 = dbs.search_DDI_in_drugbank(drug1, drug2)
        description2 = dbs.search_DDI_in_drugbank(drug2, drug1)

        # If the descriptions are empty, there is NO DDI
        if description1 == [] and description2 ==[]:
            print("\t\tNo DDI")
            no_ddi.add(drug1)
            no_ddi.add(drug2)
            num_no_ddi+=1

            if ddi_name1 not in no_ddi_pairs and ddi_name2 not in no_ddi_pairs:
                no_ddi_pairs.add(ddi_name1)

        # If the descriptions are full, there is DDI
        else:
            print(description1)
            ddi.add(drug1)
            ddi.add(drug2)
            num_ddi+=1

            if ddi_name1 not in ddi_pairs and ddi_name2 not in ddi_pairs:
                ddi_pairs.add(ddi_name1)

        i+=1
    list_of_drugs2.remove(drug1)
    n+=1

print("There are %d DDI and %d no reported DDI" %(num_ddi, num_no_ddi))
print("There are %d drugs participating in DDI and %d drugs not participating in DDI" %(len(ddi), len(no_ddi)))


#------------------------------------------------#
#   CHOOSE RANDOMLY THE COMPARISONS TO ANALYZE   #
#------------------------------------------------#

new_ddi_pairs = random.sample(list(ddi_pairs), 50)
new_no_ddi_pairs = random.sample(list(no_ddi_pairs), 50)

#----------------------------------------------#
#   GENERATE THE BASH FILE OF THE EXPERIMENT   #
#----------------------------------------------#

f = open(bash_file,"w")

# Add a comment to the bash file
f.write("#!/bin/bash\n")
f.write("\n# DDI EXPERIMENTS\n\n")

for ddi_pair in new_ddi_pairs:

    # Obtain the two drugs participating in the DDI
    (drug1, drug2) = ddi_pair.split('---')

    # Check the type of drug name and modify it if necessary
    drug1 = process_drug_name_for_design_experiment(drug1)
    drug2 = process_drug_name_for_design_experiment(drug2)

    # Add a comment to the bash file
    f.write("# %s vs. %s\n" %(drug1, drug2))

    # Add the command of the comparison
    #command = "/soft/devel/python-2.7/bin/python diana_vbiana.py -d1 %s -d2 %s -sp 1 -ex top_thresholds.list -db %s\n" %(drug1, drug2, database)
    command = "/soft/devel/python-2.7/bin/python diana_vbiana.py -d1 %s -d2 %s -sif %s -tis geneid -tra %s -sp 1 -ex top_thresholds.list -db %s\n" %(drug1, drug2, sif_file, translation_file, database)
    f.write(command)

# Add a comment to the bash file
f.write("\n# NO DDI EXPERIMENTS\n\n")

for no_ddi_pair in new_no_ddi_pairs:

    # Obtain the two drugs participating in the DDI
    (drug1, drug2) = no_ddi_pair.split('---')

    # Check the type of drug name and modify it if necessary
    drug1 = process_drug_name_for_design_experiment(drug1)
    drug2 = process_drug_name_for_design_experiment(drug2)

    # Add a comment to the bash file
    f.write("# %s vs. %s\n" %(drug1, drug2))

    # Add the command of the comparison
    #command = "/soft/devel/python-2.7/bin/python diana_vbiana.py -d1 %s -d2 %s -sp 1 -ex top_thresholds.list -db %s\n" %(drug1, drug2, database)
    command = "/soft/devel/python-2.7/bin/python diana_vbiana.py -d1 %s -d2 %s -sif %s -tis geneid -tra %s -sp 1 -ex top_thresholds.list -db %s\n" %(drug1, drug2, sif_file, translation_file, database)
    f.write(command)

