from diana.goatools.obo_parser import GODag
from diana.goatools.associations import read_ncbi_gene2go
from diana.goatools.go_enrichment import GOEnrichmentStudy
import diana.toolbox.guild_utilities as GU
#from toolbox import functional_enrichment
import math
import os


def main():

    # results_dir = "data/DCC0303/guild_results_using_sif"

    # pvalue_file = results_dir+"/output_scores.sif.netcombo.pval"
    # node_profile = results_dir+"/node_profile.txt"

    # node_to_vals = GU.get_values_from_pvalue_file(pvalue_file)
    # all_nodes = node_to_vals.keys()
    # all_nodes = [ int(x) for x in all_nodes ]

    # top_node_to_vals = GU.get_values_from_pvalue_file(node_profile)
    # top_nodes = top_node_to_vals.keys()
    # top_nodes = [ int(x) for x in top_nodes ]

    # obodag = GODag("go-basic.obo")
    # geneid2gos_human = read_ncbi_gene2go("gene2go", taxids=[9606])

    # temp_file = results_dir+"/temp_enrichment_goatools.txt"
    # output_file = results_dir+"/enrichment_goatools.txt"

    # calculate_functional_enrichment_profile(obodag, geneid2gos_human, top_nodes, all_nodes, temp_file, output_file)

    return


def calculate_log_of_odds_ratio(q, k, m, t):
    """Calculates the log of the odds ratio"""
    # print("q: {}".format(q))
    # print("k: {}".format(k))
    # print("m: {}".format(m))
    # print("t: {}".format(t))
    odds_ratio = ( (float(q)/float(k)) / (float(m)/float(t)) )
    #odds_ratio = ( (float(q)/(float(k)-float(q))) / (float(m)/((float(t)-float(m)))) )
    if odds_ratio == 0:
        return -float('inf')
    else:
        return float(math.log(odds_ratio, 2))


def calculate_functional_enrichment_profile(obodag, geneid2go, top_nodes, all_nodes, temp_file, output_file):

    print("{N:,} annotated human genes".format(N=len(geneid2go)))

    # Create a dictionary of GOs to their geneids (it will be used for the Log of Odds calculus)
    go2geneids = {}
    for geneid in geneid2go:
        for function in geneid2go[geneid]:
            go2geneids.setdefault(function, set())
            go2geneids[function].add(geneid)

    # Define the GOEnrichmentStudy object
    goeaobj = GOEnrichmentStudy(
            all_nodes, # List of background genes
            geneid2go, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.9999, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method

    # Writing the temorary file (without the Log of odds ratio)
    goea_results_all = goeaobj.run_study(top_nodes)
    goeaobj.wr_txt(temp_file, goea_results_all)

    results = []
    ft = open(temp_file, 'r')
    k = len(top_nodes) # num genes in top
    t = len(all_nodes) # num genes in network

    # Calculating the Log of odds ratio
    for line in ft:
        words = line.strip().split()
        go = words[0]
        type_func = words[1]
        pvalue = words[2]
        name = ' '.join(words[4:])
        q = int(words[3]) # num genes assoc with function in top
        m = len(go2geneids[go]) # num genes assoc with func in network
        log_of_odds = calculate_log_of_odds_ratio(q, k, m, t)
        results.append([q, m, log_of_odds, '-', pvalue, go, name, type_func])

    ft.close()

    # Writing the output file
    fo = open(output_file, 'w')
    fo.write('# num of genes\tnum of total genes\tLog of odds ratio\tP-value\tAdjusted p-value\tGO term ID\tGO term name\tGO type\n')

    for line in sorted(results, key=lambda x: x[2], reverse=True):
        if line[0] <= 0: # Skip the functions that do not have at leaast 1 top gene associated (the log of odds ratio is -infinite!!)
            continue
        if line[2] < 0: # Skip the log of odds ratio that are negative (we will only use functions with log of odds ratio from 0 to inf)
            continue
        if line[7].upper() == 'BP' or line[7].upper() == 'MF':
            new_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7])
            fo.write(new_line)

    fo.close()

    # Removing the temporary file
    command = "rm {}".format(temp_file)
    os.system(command)

    return


if  __name__ == "__main__":
    main()

