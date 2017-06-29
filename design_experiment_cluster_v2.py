import cPickle
import random
import sqlite3
import sys
import toolbox.database_searcher as dbs


#------------#
#   INPUTS   #
#------------#

# The name of the final bash file to run the experiment
crossings_file = "diana_crossings_cluster_3targets.txt"

# Minimum number of targets of the drugs
min_num_targets = 3

#---------------------------#
#   GET THE LIST OF DRUGS   #
#---------------------------#

dump_file = "toolbox/"+"dcdb2targets.pcl"
dcdb2targets = cPickle.load(open(dump_file))

list_of_drugs = []
for dcdb in dcdb2targets:
    if len(dcdb2targets[dcdb]) >= min_num_targets:
        list_of_drugs.append(dcdb)

print ("The number of drugs included is: %d" %(len(list_of_drugs)))

#-------------------------------------------------#
#   DEFINE ALL POSSIBLE CROSSINGS BETWEEN PAIRS   #
#-------------------------------------------------#

# We need to copy the list using the list() method because if not, when you modify one of the lists, the other gets modified as well
# This can also be done with copy.copy() method or copy.deepcopy() if the list contains objects and you want to copy them as well
# More info: http://stackoverflow.com/questions/2612802/how-to-clone-or-copy-a-list 
list_of_drugs2 = list(list_of_drugs)

dump_file = "toolbox/drug_int_2_drugs.pcl"
drug_int_2_drugs = cPickle.load(open(dump_file))
dump_file = "toolbox/drug_int_2_info.pcl"
drug_int_2_info = cPickle.load(open(dump_file))

crossings = set()
pair2comb = {}
dc = 0
non_dc = 0
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

            pair2comb[ddi_name1] = 0 # We will introduce 0 if it is not drug interaction
            non_dc+=1
            for drug_int in drug_int_2_drugs:
                if drug1 in drug_int_2_drugs[drug_int] and drug2 in drug_int_2_drugs[drug_int]:
                    if drug_int_2_info[drug_int]['type'] == 'pharmacodynamical':
                        pair2comb[ddi_name1] = 1 # We will introduce 1 if it is a pharmacodynamical drug interaction
                        dc+=1
                        non_dc-=1
                    else:
                        pair2comb[ddi_name1] = 0 # We will introduce 0 if it is not a pharmacodynamical drug interaction
                    break

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

print('NUMBER OF PHARMACODYNAMICAL DRUG INTERACTIONS:\t\t{}'.format(dc))
print('NUMBER OF NON-PHARMACODYNAMICAL DRUG INTERACTIONS:\t{}'.format(non_dc))

# Save the dict containing if the pairings are drug combinations or not
dump_file = "toolbox/"+"pair2comb.pcl"
cPickle.dump(pair2comb, open(dump_file, 'w')) 


#------------------------------------------#
#   GENERATE THE FILE WITH THE CROSSINGS   #
#------------------------------------------#

f = open(crossings_file,"w")

for pair in crossings:

    f.write("{}\n".format(pair))