from collections import defaultdict
from scipy.stats import fisher_exact, hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import sys, os, re


def calculate_functional_enrichment_profile(top_geneids, type_correction, associations_file, output_file=None):
    """
    Calculate the functional enrichment of the seeds, obtaining GO terms significantly enriched.
    Calculate how many of the nodes in the ranking are part of the enriched GO terms.
    Parameters:
        @top_geneids:           Set of Entrez gene IDs to check the enrichment
        @type_correction:       The multiple test p-value correction (fdr_bh / bonferroni)
        @output_file:           Resulting file which will contain the functions enriched. If None, the file is not created.
        @associations_file:     File containing the function-gene associations
    """
    # Read associations file
    dicbp, dicbp2, term2name = load_functional_terms(associations_file)
    N = len(dicbp.keys())
    NB = len(dicbp2.keys())

    # Calculate the functional enrichment
    test_passed_f1, terms_l1, term_to_values = functional_enrichment(dicbp2, N, list(set(top_geneids)), type_correction)

    # Write the output file
    if output_file:
        if len(top_geneids) > 0:
            write_enrichment_file_Carlota(term_to_values, term2name, output_file)
        else:
            with open(output_file, 'w') as out_fd:
                out_fd.write('# Term ID\tTerm name\tnum of genes\tnum of total genes\tP-value\tP-value corrected')
    return term_to_values


def functional_enrichment(dicbp2, N, geneset_toenrich, type_correction, pvalue_threshold=0.05):
    """
    Calculate the functional enrichment using the method described by Carlota.
    Parameters:
        @dicbp2:                A dictionary containing the associations of functions to their corresponding genes
        @N:                     The total number of genes
        @geneset_toenrich:      The gene set to check the enrichment
        @type_correction:       The multiple test p-value correction (fdr_bh / bonferroni)
        @pvalue_threshold:      P-value threshold of the multiple test p-value correction
    """
    pvals = {}
    genes_enriched = {}
    go_enriched = False
    k = len(geneset_toenrich)
    terms_l = []
    test_passed_f = False
    term_to_values = {}
    for term in dicbp2.keys():
        m = len(dicbp2[term])
        xl = [y for y in dicbp2[term] if y in geneset_toenrich]
        x = len(xl)
        if x != 0:
            go_enriched = True
            xlist = []
            for i in range(x, m + 1):
                xlist.append(i)
            # calculation of the hypervalue
            dhypervalue = hypergeom.pmf(xlist, N, m, k)
            # threshold of enrichment
            pvals[term] = sum(dhypervalue)
            genes_enriched[term] = [x, m] # Quim: addition to get genes enriched
    if go_enriched:
        pvals_values = list(pvals.values())
        terms = list(pvals.keys())
        pvals_corrected = multipletests(pvals_values, alpha=0.05, method=type_correction, is_sorted=False,
                                        returnsorted=False)
        for i in range(0, len(terms)):
            if list(pvals_corrected[1])[i] < pvalue_threshold:
                # Quim: addition to get the p-values and genes enriched
                pval = pvals_values[i]
                pval_corrected = list(pvals_corrected[1])[i]
                term = terms[i]
                x, m = genes_enriched[term]
                term_to_values[term] = [pval, pval_corrected, x, m]
                ####################################
                test_passed_f = True
                terms_l.append(terms[i])
    return test_passed_f, terms_l, term_to_values


def load_functional_terms(associations_file):
    """
    Loads the gene - function associations from the associations file.
    """
    dicbp = defaultdict(list)
    dicbp2 = defaultdict(list)
    term2name = {}
    with open(associations_file, 'r') as associations_fd:
        for line in associations_fd.readlines():
            fields = line.rstrip().split("\t")
            if len(fields) == 2:
                gobp = fields[0]
                geneid = fields[1]
            elif len(fields) == 3:
                gobp = fields[0]
                name = fields[1]
                geneid = fields[2]
                term2name[gobp] = name

            if gobp not in dicbp[geneid]:
                dicbp[geneid].append(gobp)

            if geneid not in dicbp2[gobp]:
                dicbp2[gobp].append(geneid)
    return dicbp, dicbp2, term2name


def write_enrichment_file_Carlota(term_to_values, term2name, enrichment_file):
    with open(enrichment_file, 'w') as out_fd:
        out_fd.write('# Term ID\tTerm name\tnum of genes\tnum of total genes\tP-value\tP-value corrected\n')
        for term, values in sorted(term_to_values.items(),key=lambda x: x[1][1],reverse=False):
            pval, pval_corrected, x, m = values
            if term in term2name:
                term_name = term2name[term]
            else:
                term_name = '-'
            out_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(term, term_name, x, m, pval, pval_corrected))
    return


def create_association_file(all_geneids, type_function, taxID, output_file, functions_data_dir):
    """
    Create a file containing gene-function associations:
        function ID <tab> function name <tab> gene
    all_geneids = the geneids considered to create the associations file
    type_function = reactome OR gobp OR gomf
    Evidence codes in GO: 'EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'ISS', 'ISA', 'ISM', 'ISO'
    Evidence codes in Reactome: 'IEA', 'TAS'
    """

    # Parse the raw association files
    go_dir = os.path.join(functions_data_dir, 'GO')
    reactome_dir = os.path.join(functions_data_dir, 'Reactome')
    evidence_codes_reactome = ['IEA', 'TAS']
    evidence_codes_go = ['EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'ISS', 'ISA', 'ISM', 'ISO']
    if type_function.lower() == 'reactome':
        reactome_file = os.path.join(reactome_dir, 'NCBI2Reactome.txt')
        term2genes, term2name = parse_reactome(reactome_file, evidence_codes_reactome, taxID)
    elif type_function.lower() == 'gobp':
        term2genes, term2name = parse_gene2go(go_dir, 'Process', evidence_codes_go, taxID)
    elif type_function.lower() == 'gomf':
        term2genes, term2name = parse_gene2go(go_dir, 'Function', evidence_codes_go, taxID)
    else:
        print('Incorrect type of function introduced: {}'.format(type_function))
        sys.exit(10)

    # Output the processed association file
    with open(output_file, 'w') as out_fd:
        for termID in term2genes:
            term_name = term2name[termID]
            for geneID in term2genes[termID]:
                if geneID in all_geneids:
                    out_fd.write('{}\t{}\t{}\n'.format(termID, term_name, geneID))
    return


def parse_gene2go(go_dir, type_go, evidence_codes_used, taxID):
    """
    Parse associations file from GO.
    """
    # In a normal species we use the file gene2go
    # Downloaded from: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
    if int(taxID) != 4932:
        go2name = {}
        go2geneids = {}
        gene2go_file = os.path.join(go_dir, 'gene2go')
        with open(gene2go_file, 'r') as inp_fd:
            first_line = inp_fd.readline()
            for line in inp_fd:
                #tax_id GeneID  GO_ID   Evidence    Qualifier   GO_term PubMed  Category
                tax_id, geneID, goID, evidence, qualifier, go_name, pubmed, category  = line.strip().split('\t')
                if category != type_go or tax_id != str(taxID):
                    continue
                if evidence in evidence_codes_used:
                    go2geneids.setdefault(goID, set()).add(geneID)
                    go2name[goID] = go_name
    # In the exceptional case of Saccharomyces cerevisiae (Taxonomy ID 4932) we use the file sgd.gaf
    # Downloaded from Gene Ontology
    # We parse it using GOAtools
    else:
        taxID_to_gaf_file = {
            9606 : 'goa_human.gaf',
            10090 : 'mgi.gaf',
            10116 : 'rgd.gaf',
            6239 : 'wb.gaf',
            7227 : 'fb.gaf',
            3702 : 'tair.gaf',
            4932 : 'sgd.gaf',
        } 
        taxID_to_geneinfo_file = {
            9606 : 'Homo_sapiens.gene_info',
            10090 : 'Mus_musculus.gene_info',
            10116 : 'Rattus_norvegicus.gene_info',
            6239 : 'Caenorhabditis_elegans.gene_info',
            7227 : 'Drosophila_melanogaster.gene_info',
            3702 : 'Arabidopsis_thaliana.gene_info',
            4932 : 'Saccharomyces_cerevisiae.gene_info',
        } 
        from diana.classes.goatools.obo_parser import GODag
        from diana.classes.goatools.associations import read_gaf
        go_obo_file = os.path.join(go_dir, 'go-basic.obo')
        obodag = GODag(go_obo_file)
        gaf_file = os.path.join(go_dir, 'goa_files/{}'.format(taxID_to_gaf_file[int(taxID)]))
        geneinfo_file = os.path.join(go_dir, 'goa_files/{}'.format(taxID_to_geneinfo_file[int(taxID)]))
        gene2gos = read_gaf(gaf_file, taxids=[int(taxID)], evidence_set=evidence_codes_used)
        gos = set()
        for gene in gene2gos:
            for go in gene2gos[gene]:
                gos.add(go)
        print('Genes: {}, GOs: {}'.format(len(gene2gos), len(gos)))

        geneid2gos = {}
        go2name = {}
        with open(geneinfo_file, 'r') as inp_fd:
            first_line = inp_fd.readline()
            for line in inp_fd:
                fields = line.strip().split('\t')
                geneID = fields[1]
                gene = fields[5]
                #print(geneID, gene)
                if ':' in gene:
                    gene = gene.split(':')[1]
                if gene != '-' and geneID != '-':
                    if gene in gene2gos:
                        for go in gene2gos[gene]:
                            geneid2gos.setdefault(geneID, set()).add(go)
                            go_name = obodag[go].name
                            go2name[go] = go_name
                    else:
                        #print('Gene {} not found'.format(gene))
                        pass
        gos = set()
        for gene in geneid2gos:
            for go in geneid2gos[gene]:
                gos.add(go)
        print('GeneIDs: {}, GOs: {}'.format(len(geneid2gos), len(gos)))
        go2geneids = associations.get_b2aset(geneid2gos)
    return go2geneids, go2name


def parse_reactome(reactome_file, evidence_codes, taxID):
    """
    Parse Reactome to NCBI gene file (Lowest level pathway diagram / Subset of the pathway).
    Downloaded from: https://reactome.org/download/current/NCBI2Reactome.txt
    """
    taxID2species = {
        9606 : 'Homo sapiens',
        10090 : 'Mus musculus',
        10116 : 'Rattus norvegicus',
        6239 : 'Caenorhabditis elegans',
        7227 : 'Drosophila melanogaster',
        3702 : 'Arabidopsis thaliana',
        4932 : 'Saccharomyces cerevisiae'
    }
    selected_species = taxID2species[int(taxID)]
    reactome2genes = {}
    reactome2name = {}
    with open(reactome_file, 'r') as inp_fd:
        for line in inp_fd:
            fields = line.strip().split('\t')
            # 1 R-HSA-114608    https://reactome.org/PathwayBrowser/#/R-HSA-114608  Platelet degranulation  TAS Homo sapiens
            geneID, reactomeID, reactome_url, reactome_name, evidence_code, species_name = fields
            if evidence_code not in evidence_codes:
                continue
            if selected_species == species_name:
                reactome2genes.setdefault(reactomeID, set()).add(geneID)
                reactome2name[reactomeID] = reactome_name
    return reactome2genes, reactome2name


def create_number_of_functions_file(associations_file, number_of_functions_file):
    """
    Create a file with the number of functions associated to the genes of a network.
    This will be used in the calculation of overlap of functions.
    """
    dicbp, dicbp2, term2name = load_functional_terms(associations_file)
    functions = dicbp2.keys()
    with open(number_of_functions_file, 'w') as number_of_functions_fd:
        number_of_functions_fd.write('{}\n'.format(len(functions)))
    return


def get_functions_from_associations_file(associations_file):
    """
    Get the functions from the associations file.
    """
    dicbp, dicbp2, term2name = load_functional_terms(associations_file)
    functions = dicbp2.keys()
    return functions


def read_number_of_functions_file(number_of_functions_file):
    """
    Reads the file with the total number of functions.
    """
    with open(number_of_functions_file, 'r') as num_functions_fd:
        number_of_functions = int(num_functions_fd.readline().strip("\n"))
    return number_of_functions


def calculate_functions_threshold(seed_geneids, geneid_to_score, type_correction, associations_file, output_sliding_window_file, output_seeds_enrichment_file=None, seed_functional_enrichment=False):
    """
    Calculate the threshold of a GUILD scored list of genes by their functional similarity with their seeds.
    Parameters:
        @seed_geneids:                  Set of seed Entrez gene IDs to check the enrichment
        @type_correction:               The multiple test p-value correction (fdr_bh / bonferroni)
        @associations_file:             File containing the function-gene associations
        @output_sliding_window_file     File containing the sliding window information.
        @output_seeds_enrichment_file:  Resulting file which will contain the functions enriched for the seeds. If None, the file is not created
        @seed_functional_enrichment:    If true, a functional enrichment of the seeds is calculated.
                                        If false, all the functions associated to the seeds are considered.
    """
    # Read associations file
    dicbp, dicbp2, term2name = load_functional_terms(associations_file)
    N = len(dicbp.keys())
    NB = len(dicbp2.keys())

    # Calculate the terms associated with the seeds
    terms_enriched = set()
    genes_in_terms_enriched = set()
    if seed_functional_enrichment == True:
        # Calculate functional enrichment
        test_passed_f1, terms_l1, term_to_values = functional_enrichment(dicbp2, N, list(set(seed_geneids)), type_correction)
        # Write the output seeds file
        if output_seeds_enrichment_file:
            if len(seed_geneids) > 0:
                write_enrichment_file_Carlota(term_to_values, term2name, output_seeds_enrichment_file)
            else:
                with open(output_seeds_enrichment_file, 'w') as out_fd:
                    out_fd.write('# Term ID\tTerm name\tnum of genes\tnum of total genes\tP-value\tP-value corrected')
        # Get the genes of the functions enriched
        for term in term_to_values:
            pval, pval_corrected, x, m = term_to_values[term]
            if float(pval_corrected) < 0.05:
                terms_enriched.add(term)
                if term in dicbp2:
                    for gene in dicbp2[term]:
                        genes_in_terms_enriched.add(gene)
    else:
        for gene in seed_geneids:
            if gene in dicbp:
                genes_in_terms_enriched.add(gene)
                for term in dicbp[gene]:
                    terms_enriched.add(term)

    # Order genes by highest score
    ranking_geneids = sorted(geneid_to_score, key=geneid_to_score.get, reverse=True)

    # Get the ranking positions of the genes that are associated with significantly enriched GOs
    enriched_positions = []
    non_enriched_positions = []
    seed_positions = []
    position = 0
    for geneid in ranking_geneids:
        position += 1
        if geneid in genes_in_terms_enriched:
            enriched_positions.append(position)
        else:
            non_enriched_positions.append(position)
        if geneid in seed_geneids:
            seed_positions.append(position)

    # Get the maximum number of sliding windows to check
    # We put it as 10 times the length of the sliding window (or the number of seeds)
    #maximum_rank = len(seed_positions) * 10
    maximum_rank = 500
    # Calculate enrichment in sliding window
    cutoff_central_position, cutoff_right_interval = calculate_enrichment_in_all_sliding_window_positions(output_file=output_sliding_window_file, enriched_positions=enriched_positions, non_enriched_positions=non_enriched_positions, seed_positions=seed_positions, maximum_rank=maximum_rank)

    return cutoff_central_position, cutoff_right_interval


def calculate_enrichment_in_all_sliding_window_positions(output_file, enriched_positions, non_enriched_positions, seed_positions, maximum_rank=500):
    """
    Calculate the enrichment in all the positions of the sliding window.
    """

    num_seeds = len(seed_positions)

    # Calculate all the central positions of the sliding window
    all_positions = range(calculate_initial_position_of_sliding_window(num_seeds), calculate_final_position_of_sliding_window(num_seeds, maximum_rank)+1)

    # When below_cutoff = 0, it has not found any p-value under 0.05
    # When below_cutoff = 1, it is finding for the first time values under 0.05
    # When below_cutoff = 2, it has gone over 0.05
    below_cutoff = 0
    cutoff_central_position = 0
    cutoff_right_interval = 0

    with open(output_file, 'w') as output_fd:
        output_fd.write('#Position\tOddsRatio\tP-value\tType\n')
        for position in all_positions:
            left_interval = position - int( float(num_seeds) / float(2) ) # Define left interval position
            right_interval = position + int( float(num_seeds) / float(2) ) # Define right interval position
            window = set(range(left_interval, right_interval+1)) # Define window positions
            non_window = set(all_positions) - window # Define non-window positions
            enriched_positions_window = set(enriched_positions) & window
            non_enriched_positions_window = set(non_enriched_positions) & window
            enriched_positions_non_window = set(enriched_positions) & non_window
            non_enriched_positions_non_window = set(non_enriched_positions) & non_window
            # Calculate if the enrichment is significant
            contingency_table = [[len(enriched_positions_window), len(non_enriched_positions_window)], [len(enriched_positions_non_window), len(non_enriched_positions_non_window)]]
            oddsratio, pvalue = fisher_exact(contingency_table)
            if position in seed_positions:
                cutoff_central_position = position # Get the last position where the cut-off is below 0.05
                cutoff_right_interval = right_interval # Get the last node which is included in the interval
                output_fd.write('{}\t{}\t{}\t{}\n'.format(position, oddsratio, pvalue, 'seed'))
                continue
            else:
                output_fd.write('{}\t{}\t{}\t{}\n'.format(position, oddsratio, pvalue, 'non-seed'))
            if below_cutoff == 0 and pvalue <= 0.05:
                below_cutoff = 1
            elif below_cutoff == 0 and pvalue > 0.05: # If in the first non-seed position, the value is already above 0.05, we stop getting the cut-off
                cutoff_central_position = position # Get the last position where the cut-off is below 0.05
                cutoff_right_interval = right_interval # Get the last node which is included in the interval
                below_cutoff = 2 # When pvalue is above 0.05, we stop getting the cutoff
                #print('Stop at position {}. central: {}. right: {}'.format(position, cutoff_central_position, cutoff_right_interval))
            if below_cutoff == 1 and pvalue <= 0.05:
                cutoff_central_position = position # Get the last position where the cut-off is below 0.05
                cutoff_right_interval = right_interval # Get the last node which is included in the interval
            elif below_cutoff == 1 and pvalue > 0.05:
                below_cutoff = 2 # When pvalue is above 0.05, we stop getting the cutoff

    return cutoff_central_position, cutoff_right_interval


def calculate_initial_position_of_sliding_window(num_seeds):
    """
    Calculate the initial position where the sliding position
    will start its calculations.
    """
    # Calculate the initial position of the sliding window:
    # We know that the left part of the first interval in the sliding window is:
    # i - num_seeds / 2 = 1
    # So, if we want to know i (the initial position) we have to calculate:
    # i = 1 + num_seeds / 2
    initial_position = 1 + int( float(num_seeds) / float(2) )
    #print('Initial position: {}'.format(initial_position))
    return initial_position


def calculate_final_position_of_sliding_window(num_seeds, maximum_rank=500):
    """
    Calculate the final position where the sliding position
    will start its calculations.
    maximum_rank = The maximum number of positions that we want to include
                   in the calculations
    """
    # Calculate the final position of the sliding window:
    # We know that the right part of the last interval will be:
    # i + num_seeds / 2 = last_rank
    # So the final position will be:
    # i = last_rank - num_seeds / 2
    last_rank = maximum_rank
    final_position = last_rank - int( float(num_seeds) / float(2) )
    #print('Final position: {}'.format(final_position))
    return final_position


def plot_sliding_window(sliding_window_file, output_plot_file, maximum_rank=500):
    """
    Plot the sliding window in matplotlib
    """
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import pylab
    
    # Read sliding window file
    sliding_window_df = pd.read_csv(sliding_window_file, sep='\t', index_col=None)
    sliding_window_df['Cut-off'] = 0.05

    # Plot LINE chart
    fig = pylab.figure(dpi=300)
    ax = pylab.axes()
    ax.plot(sliding_window_df['#Position'], sliding_window_df['P-value'])
    ax.plot(sliding_window_df['#Position'], sliding_window_df['Cut-off'])
    ax.grid(True)
    ax.set_xlabel('Ranking')
    ax.set_ylabel('P-values')
    ax.set_xlim(0, maximum_rank)
    plt.setp(ax)
    pylab.savefig(output_plot_file, format='png')
    return


def read_sliding_window_file(output_sliding_window_file, num_seeds):
    """
    Get the central  and right position of the sliding window from the sliding window file.
    """
    # When below_cutoff = 0, it has not found any p-value under 0.05
    # When below_cutoff = 1, it is finding for the first time values under 0.05
    # When below_cutoff = 2, it has gone over 0.05
    below_cutoff = 0
    cutoff_central_position = 0
    cutoff_right_interval = 0

    with open(output_sliding_window_file, 'r') as slide_fd:
        for line in slide_fd:
            if line[0] == '#':
                continue
            fields = line.strip().split('\t')
            position = int(fields[0])
            pvalue = float(fields[2])
            type_protein = fields[3]
            left_interval = position - int( float(num_seeds) / float(2) ) # Define left interval position
            right_interval = position + int( float(num_seeds) / float(2) ) # Define right interval position
            if type_protein == 'seed':
                cutoff_central_position = position # Get the last position where the cut-off is below 0.05
                cutoff_right_interval = right_interval # Get the last node which is included in the interval
                continue
            if below_cutoff == 0 and pvalue <= 0.05:
                below_cutoff = 1
            elif below_cutoff == 0 and pvalue > 0.05: # If in the first non-seed position, the value is already above 0.05, we stop getting the cut-off
                cutoff_central_position = position # Get the last position where the cut-off is below 0.05
                cutoff_right_interval = right_interval # Get the last node which is included in the interval
                below_cutoff = 2 # When pvalue is above 0.05, we stop getting the cutoff
                #print('Stop at position {}. central: {}. right: {}'.format(position, cutoff_central_position, cutoff_right_interval))
            if below_cutoff == 1 and pvalue <= 0.05:
                cutoff_central_position = position # Get the last position where the cut-off is below 0.05
                cutoff_right_interval = right_interval # Get the last node which is included in the interval
            elif below_cutoff == 1 and pvalue > 0.05:
                below_cutoff = 2 # When pvalue is above 0.05, we stop getting the cutoff
                break

    return cutoff_central_position, cutoff_right_interval


