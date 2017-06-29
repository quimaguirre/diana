import cPickle
import random
import sqlite3
import sys
import toolbox.database_searcher as dbs


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
num_target_threshold = 1

# The number of analysis that we want to do for DDI and for non-DDI
num_analyses = 50

# The name of the final bash file to run the experiment
bash_file = "diana_crossings_test_2.txt"

# The database you want to use for the analysis
database = "BIANA_JAN_2017"

# The SIF file needed
sif_file = "sif/human_eAFF_geneid_2017.sif"

# The translation file needed
translation_file = "translation_geneid_eAFF_2017.txt"

#-----------------------------------------------------------#
#   FILTER THE DRUGS OF THE DATABASE BY NUMBER OF TARGETS   #
#-----------------------------------------------------------#

dump_file = "toolbox/"+"dcdb2targets.pcl"
dcdb2targets = cPickle.load(open(dump_file))
list_of_drugs = dcdb2targets.keys()

print ("The number of drugs with at least %d targets is: %d" %(num_target_threshold, len(list_of_drugs)))

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
dump_file = "toolbox/"+"drug_comb_2_drugs.pcl"
drug_comb_2_drugs = cPickle.load(open(dump_file))
ddi_pairs = set()
no_ddi_pairs = set()

# We will run this loop until the DDIs set and non-DDIs set are filled with the defined number of analyses
while ( len(ddi_pairs) < num_analyses or len(no_ddi_pairs) < num_analyses ):

    ddi = False
    # We choose one random crossing
    choice = random.choice(list(crossings))
    # We remove the selected item from the crossings set
    crossings.remove(choice)

    # Obtain the two drugs participating in the crossing
    (drug1, drug2) = choice.split('---')
    #print("COMPARING: %s vs. %s" %(drug1, drug2))

    # Retrieve the description of the interaction between the two drugs if there is interaction
    # The description could be different depending on the order of the drugs, so the procedure is 
    # performed twice changing the order of the drugs
    for drug_comb in drug_comb_2_drugs:
        if drug1 in drug_comb_2_drugs[drug_comb] and drug2 in drug_comb_2_drugs[drug_comb]:
            if len(ddi_pairs) < num_analyses:
                ddi = True
                ddi_pairs.add(choice)
                print("%d drug comb added\n" %(len(ddi_pairs)))
                break
    if ddi == False:
        if len(no_ddi_pairs) < num_analyses:
            no_ddi_pairs.add(choice)
            print("%d non-drug comb added\n" %(len(no_ddi_pairs)))



#----------------------------------------------#
#   GENERATE THE BASH FILE OF THE EXPERIMENT   #
#----------------------------------------------#

f = open(bash_file,"w")

for ddi_pair in ddi_pairs:

    f.write("{}\n".format(ddi_pair))

for no_ddi_pair in no_ddi_pairs:

    f.write("{}\n".format(no_ddi_pair))


f.close()
