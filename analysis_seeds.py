from os import listdir
from os.path import isfile, isdir, join
import re
import diana_vbiana as diana
import toolbox
import toolbox.drugbank_searcher as dbs
import matplotlib.pyplot as plt


analysis_file = "100_analyses/diana_100_analyses_27_10_16.sh"

# Parse the analysis file in order to obtain the names of the DDI pairs and non-DDI pairs
fa = open(analysis_file,"r")
ddi = False
non_ddi = False
pair_regex = re.compile('# "*(.+?)"* vs\. "*(.+?)"*\n')


drugs = {}

for line in fa:
    if line == '# DDI EXPERIMENTS\n':
        ddi = True
        non_ddi = False
    elif line == '# NO DDI EXPERIMENTS\n':
        ddi = False
        non_ddi = True

    if ddi == True:
        m = pair_regex.search(line)
        if m != None:
            drug_ori1 = m.group(1)
            drug_ori2 = m.group(2)
            drug1 = diana.check_special_drug_name(drug_ori1)
            drug2 = diana.check_special_drug_name(drug_ori2)
            drugs[drug_ori1] = drug1
            drugs[drug_ori2] = drug2
    elif non_ddi == True:
        m = pair_regex.search(line)
        if m != None:
            drug_ori1 = m.group(1)
            drug_ori2 = m.group(2)
            drug1 = diana.check_special_drug_name(drug_ori1)
            drug2 = diana.check_special_drug_name(drug_ori2)
            drugs[drug_ori1] = drug1
            drugs[drug_ori2] = drug2

fa.close()

main_directory = "100_analyses/data"

diff_seeds = []

for drug in drugs:
    directory = main_directory+'/'+drugs[drug]

    ori_seeds = dbs.search_seeds_in_drugbank(drug, provided_seeds = None)

    seeds_file = directory+'/guild_results_using_sif/seeds.txt'

    final_seeds = []
    fs = open(seeds_file, 'r')
    for seed in fs:
        final_seeds.append(seed)
    fs.close()

    print("%d    %d"%(len(ori_seeds), len(final_seeds)))
    diff = len(ori_seeds) - len(final_seeds)
    diff_seeds.append(diff)

print(diff_seeds)


# Histogram 
n, bins, patches = plt.hist(diff_seeds, 30, 
                            histtype='bar')
plt.title(r'$\mathrm{Loss\ of\ seeds}$')
plt.xlabel('Number of seeds not found in the network')
plt.ylabel('Frequency')
plt.savefig('100_analyses/loss_of_seeds.png')
#plt.show()


