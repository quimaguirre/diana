import random
import sqlite3
import sys
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

# The number of analysis that we want to do for DDI and for non-DDI
num_analyses = 50

# The name of the final bash file to run the experiment
bash_file = "diana_100_analyses_31_10_16.sh"

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

#-------------------------------------------------#
#   DEFINE ALL POSSIBLE CROSSINGS BETWEEN PAIRS   #
#-------------------------------------------------#

# We need to copy the list using the list() method because if not, when you modify one of the lists, the other gets modified as well
# This can also be done with copy.copy() method or copy.deepcopy() if the list contains objects and you want to copy them as well
# More info: http://stackoverflow.com/questions/2612802/how-to-clone-or-copy-a-list 
list_of_drugs2 = list(list_of_drugs)

crossings = set()

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
        #print("%s vs. %s" %(drug1, drug2))

        # We check that none of the two possible names are in the crossings set, and we add it (this is not necessary, but it is added as security)
        if ddi_name1 not in crossings and ddi_name2 not in crossings:
            crossings.add(ddi_name1)

        i+=1

    # We remove the first drug from the second list, so that we do not have to repeat pairings
    list_of_drugs2.remove(drug1)
    n+=1

print("There are %d possible crossings" %(len(crossings)))
checking = len(list_of_drugs) * (len(list_of_drugs) - 1) / 2
if len(crossings) != checking:
    print("THERE IS AN ERROR IN THE ANALYSIS. The number of crossings does not correspond to the theoretical number")
    sys.exit(10)
#print(crossings)


#------------------------------------------------#
#   CHOOSE RANDOMLY THE COMPARISONS TO ANALYZE   #
#------------------------------------------------#

# Define the sets where DDIs and non-DDIs will be stored
ddi_pairs = set()
no_ddi_pairs = set()

# We will run this loop until the DDIs set and non-DDIs set are filled with the defined number of analyses
while ( len(ddi_pairs) != num_analyses or len(no_ddi_pairs) != num_analyses ):

    # We choose one random crossing
    choice = random.choice(list(crossings))
    # We remove the selected item from the crossings set
    crossings.remove(choice)

    # Obtain the two drugs participating in the crossing
    (drug1, drug2) = choice.split('---')
    print("COMPARING: %s vs. %s" %(drug1, drug2))

    # Retrieve the description of the interaction between the two drugs if there is interaction
    # The description could be different depending on the order of the drugs, so the procedure is 
    # performed twice changing the order of the drugs
    description1 = dbs.search_DDI_in_drugbank(drug1, drug2)
    description2 = dbs.search_DDI_in_drugbank(drug2, drug1)

    # If the descriptions are empty, there is NO DDI
    if description1 == [] and description2 ==[]:

        # If the set of non-DDIs is not full, we add the crossing
        if len(no_ddi_pairs) < num_analyses:
            print("\t\tNo DDI")
            no_ddi_pairs.add(choice)
            print("%d NON-DDIs added\n" %(len(no_ddi_pairs)))

    # If the descriptions are full, there is DDI
    else:

        if len(ddi_pairs) < num_analyses:
            print(description1)
            ddi_pairs.add(choice)
            print("%d DDIs added\n" %(len(ddi_pairs)))


#----------------------------------------------#
#   GENERATE THE BASH FILE OF THE EXPERIMENT   #
#----------------------------------------------#

f = open(bash_file,"w")

# Add a comment to the bash file
f.write("#!/bin/bash\n")
f.write("\n# DDI EXPERIMENTS\n\n")

for ddi_pair in ddi_pairs:

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

for no_ddi_pair in no_ddi_pairs:

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

