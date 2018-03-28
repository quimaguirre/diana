import argparse
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
import scipy
import time
import sys, os, re

from context import diana
import diana.classes.comparison as diana_comparison
import diana.classes.analysis as diana_analysis




def main():

    options = parse_user_arguments()
    analysis_results(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of the input drug",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-th','--threshold_list',dest='threshold_list',action = 'store',
                        help = """List of percentages that will be used as cut-offs to define the profiles of the drugs. It has to be a file containing:
                        - Different numbers that will be the threshold values separated by newline characters. 
                        For example, a file called "top_threshold.list" containing:
                        0.1
                        0.5
                        1
                        5
                        10
                        """)
    parser.add_argument('-f','--formula',dest='formula',action = 'store',default='simpson',
                        help = """Define the formula used to classify. It can be: simpson, jaccard""")
    parser.add_argument('-se','--consider_se',dest='consider_se',action = 'store_true',
                        help = """" Consider Side Effects / ATCs. """)
    parser.add_argument('-ws','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory will be created""")

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def analysis_results(options):
    """
    Analyzes the results of the comparisons
    """

    # Start marker for time measure
    start = time.time()

    print("\n\t\t-------------------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Analysis of results: Classify drug combinations\n")
    print("\t\t-------------------------------------------------------------------------------------------------------------------------------\n")

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')

    # Check the directory of the profiles, comparisons and analysis
    data_dir = os.path.join(options.workspace, "profiles")
    check_directory(data_dir)

    results_dir = os.path.join(options.workspace, "comparisons")
    check_directory(results_dir)

    analysis_dir = os.path.join(options.workspace, "analysis")
    check_directory(analysis_dir)

    # Get the list of thresholds to create the profiles
    if options.threshold_list and fileExist(options.threshold_list):
        threshold_list = get_values_from_threshold_file(options.threshold_list)
    else:
        threshold_list = [1, 5, 10, 20, 50]

    # Do we consider Side Effects/ATC?
    if options.consider_se:
        consider_se = True
    else:
        consider_se = False

    # Get the names of the columns
    columns = diana_analysis.obtain_columns(threshold_list, ATC_SE=consider_se)



    #-----------------------------------------------------#
    #   PARSE THE RESULTS AND CREATE A PANDAS DATAFRAME   #
    #-----------------------------------------------------#

    pair2comb_file = os.path.join(toolbox_dir, 'pair2comb.pcl')
    pair2comb = cPickle.load(open(pair2comb_file))

    diana_id_to_drugbank_file = os.path.join(toolbox_dir, 'diana_id_to_drugbank.pcl')
    diana_id_to_drugbank = cPickle.load(open(diana_id_to_drugbank_file))

    ddi = sum(1 for x in pair2comb.values() if x == 1)
    non_ddi = sum(1 for x in pair2comb.values() if x == 0)

    print('NUMBER OF DRUG COMBINATIONS:\t\t{}\n'.format(ddi))
    print('NUMBER OF NON-DRUG COMBINATIONS:\t{}\n'.format(non_ddi))

    output_dataframe = os.path.join(analysis_dir, 'dcdb_comparisons.csv')

    if not fileExist(output_dataframe):

        # Create a data frame to store the results
        df = pd.DataFrame(columns=columns)


        # Obtain all the results subfolders of the results main folder
        results_dir_list = [f for f in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, f))]

        for comparison in results_dir_list:

            drug_id1, drug_id2 = comparison.split('---')
            comparison_dir = os.path.join(results_dir, comparison)
            results_table = os.path.join(comparison_dir, 'results_table.tsv')

            # Add the Comb field (if it is drug combination or not)
            drug1 = diana_id_to_drugbank[drug_id1].upper()
            drug2 = diana_id_to_drugbank[drug_id2].upper()
            comparison_without_id = '{}---{}'.format(drug1, drug2)
            if comparison_without_id in pair2comb:
                combination_field = pair2comb[comparison_without_id]
            else:
                print('The comparison {} is not in the pair2comb dictionary!\n'.format(comparison_without_id))
                print(pair2comb)
                sys.exit(10)

            if not fileExist(results_table):
                print('The comparison {} has not been executed properly!\n'.format(comparison))
                sys.exit(10)

            results = diana_analysis.get_results_from_table(results_table, columns, combination_field)

            df2 = pd.DataFrame([results], columns=columns, index=[comparison])
            # Add the information to the main data frame
            df = df.append(df2)

        # Output the Pandas dataframe in a CSV file
        df.to_csv(output_dataframe)

    else:
        df = pd.read_csv(output_dataframe, index_col=0)



    #---------------------------#
    #   REMOVE MISSING VALUES   #
    #---------------------------#

    # Replace the None values in dcstructure by nan
    if 'None' in df['dcstructure']:
        df = df.replace(to_replace={'dcstructure':{'None':np.nan}})

    # Remove the nan values in dcstructure
    df = df.dropna()

    # Count the number of drug combinations / non-drug combinations
    dc_data = df[df['combination'] == 1]
    ndc_data = df[df['combination'] == 0]
    num_dc = len(dc_data.index)
    num_ndc = len(ndc_data.index)
    print('Number of drug combinations after removing missing values:\t{}\n'.format(num_dc))
    print('Number of non-drug combinations after removing missing values:\t{}\n'.format(num_ndc))



    #---------------------------#
    #   IDENTIFY ME-TOO DRUGS   #
    #---------------------------#

    me_too_dir = os.path.join(analysis_dir, 'me_too_drugs')
    create_directory(me_too_dir)
    me_too_drugs_table = os.path.join(me_too_dir, 'me_too_drugs.tsv')
    me_too_drug_combs_table = os.path.join(me_too_dir, 'me_too_drug_combinations.tsv')

    me_too_drug_pairs_file = os.path.join(me_too_dir, 'me_too_drug_pairs.pcl')
    me_too_drug_comb_pairs_file = os.path.join(me_too_dir, 'me_too_drug_comb_pairs.pcl')

    if not fileExist(me_too_drug_pairs_file) or not fileExist(me_too_drug_comb_pairs_file):

        df_struc = df[['dcstructure']]
        df_struc = df_struc.astype(float)
        me_too_drug_pairs, me_too_drug_comb_pairs = diana_analysis.obtain_me_too_drugs_and_combinations(df_struc, columns, me_too_drugs_table, me_too_drug_combs_table)
        cPickle.dump(me_too_drug_pairs, open(me_too_drug_pairs_file, 'w'))
        cPickle.dump(me_too_drug_comb_pairs, open(me_too_drug_comb_pairs_file, 'w'))

    else:

        me_too_drug_pairs = cPickle.load(open(me_too_drug_pairs_file))
        me_too_drug_comb_pairs = cPickle.load(open(me_too_drug_comb_pairs_file))

    # Process me-too drug combination pairs
    me_too_drug_combinations = set()
    drug_pair_to_me_too_times = {}
    for pair in me_too_drug_comb_pairs:
        drug_comb1, drug_comb2 = pair.split('___')
        me_too_drug_combinations.add(frozenset([drug_comb1, drug_comb2]))
        drug_pair_to_me_too_times.setdefault(drug_comb1, 0)
        drug_pair_to_me_too_times.setdefault(drug_comb2, 0)
        drug_pair_to_me_too_times[drug_comb1] += 1
        drug_pair_to_me_too_times[drug_comb2] += 1
    removed_drug_pairs = set()
    for pair in me_too_drug_comb_pairs:
        drug_comb1, drug_comb2 = pair.split('___')
        if drug_comb1 in removed_drug_pairs or drug_comb2 in removed_drug_pairs:
            continue
        if drug_pair_to_me_too_times[drug_comb1] > drug_pair_to_me_too_times[drug_comb2]:
            removed_drug_pairs.add(drug_comb1)
        else:
            removed_drug_pairs.add(drug_comb2)

    # Remove the drug pairs which appear in me-too pairs of drug pairs more times
    df = df.loc[~df.index.isin(list(removed_drug_pairs))]

    # Count the number of drug combinations / non-drug combinations
    dc_data = df[df['combination'] == 1]
    ndc_data = df[df['combination'] == 0]
    num_dc = len(dc_data.index)
    num_ndc = len(ndc_data.index)
    print('Number of drug combinations after removing me-too conflictive drug pairs:\t{}\n'.format(num_dc))
    print('Number of non-drug combinations after removing me-too conflictive drug pairs:\t{}\n'.format(num_ndc))



    img_dir = os.path.join(analysis_dir, 'figures')
    create_directory(img_dir)
    fig_format = 'png'

    #-----------------------------------------------------#
    #   PLOT DISTRIBUTION OF NUMBER OF TARGETS PER DRUG   #
    #-----------------------------------------------------#

    # Plot distribution of comparisons of targets
    drugbank2targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    drugbank_to_targets = cPickle.load(open(drugbank2targets_file))
    plot_distribution_targets = os.path.join(img_dir, 'distribution_number_targets.{}'.format(fig_format))
    targets = [len(x) for x in drugbank_to_targets.values()]
    n, bins, patches = plt.hist(np.array(targets), bins=50, weights=np.zeros_like(np.array(targets)) + 1. / np.array(targets).size, facecolor='r')
    plt.xlabel('Number of targets per drug')
    plt.ylabel('Relative frequency')
    plt.title('Distribution of the number of targets per drug')
    plt.savefig(plot_distribution_targets, format=fig_format, dpi=300)
    plt.clf()

    #----------------------------------------------------------------------------------------------#
    #   EVALUATE OVERLAP BETWEEN TARGETS, BIOLOGICAL PROCESSES AND PATHWAYS IN DRUG COMBINATIONS   #
    #----------------------------------------------------------------------------------------------#

    tables_dir = os.path.join(analysis_dir, 'tables')
    create_directory(tables_dir)

    if options.formula != 'jaccard' and options.formula != 'simpson':
        print('Please, introduce a correct formula to classify drug combinations: jaccard or simpson!\n')
        sys.exit(10)

    # Plot of distribution of comparisons of Targets
    plot_ji_targets = os.path.join(img_dir, 'distribution_{}_index_targets.{}'.format(options.formula, fig_format))

    # Plot of distribution of comparisons of Biological Processes
    plot_ji_bp = os.path.join(img_dir, 'distribution_{}_index_biological_processes.{}'.format(options.formula, fig_format))

    # Plot of distribution of comparisons of Pathways
    plot_ji_pathways = os.path.join(img_dir, 'distribution_{}_index_pathways.{}'.format(options.formula, fig_format))

    # Output pickle file of the classification
    classification_targets_bp_file = os.path.join(toolbox_dir, 'classification_targets_bp.pcl')
    classification_targets_pathways_file = os.path.join(toolbox_dir, 'classification_targets_pathways.pcl')


    # Get the classification files
    drug_int_2_drugs_file = os.path.join(toolbox_dir, 'drug_int_2_drugs.pcl')
    drug_int_2_drugs = cPickle.load(open(drug_int_2_drugs_file))
    drug_int_2_info_file = os.path.join(toolbox_dir, 'drug_int_2_info.pcl')
    drug_int_2_info = cPickle.load(open(drug_int_2_info_file))
    drugbank_to_dcdb_file = os.path.join(toolbox_dir, 'drugbank_to_dcdb.pcl')
    drugbank_to_dcdb = cPickle.load(open(drugbank_to_dcdb_file))

    bio_processes_file = os.path.join(toolbox_dir, 'target_to_bio_processes.pcl')
    target_to_bio_processes = cPickle.load(open(bio_processes_file))
    pathways_file = os.path.join(toolbox_dir, 'target_to_pathways.pcl')
    target_to_pathways = cPickle.load(open(pathways_file))


    target_comparisons = []
    bp_comparisons = []
    pathway_comparisons = []

    dc_to_target_ji = {}
    dc_to_bp_ji = {}
    dc_to_pathway_ji = {}
    all_drugs = set()

    for index, row in dc_data.iterrows():

        (drug_id1, drug_id2) = index.split('---')
        drug1 = diana_id_to_drugbank[drug_id1].upper()
        drug2 = diana_id_to_drugbank[drug_id2].upper()
        all_drugs.add(drug1)
        all_drugs.add(drug2)


        if drug1 in drugbank_to_targets and drug2 in drugbank_to_targets:
            targets1 = drugbank_to_targets[drug1]
            targets2 = drugbank_to_targets[drug2]
            if options.formula == 'jaccard':
                result_targets = diana_comparison.calculate_jaccard_index(targets1, targets2)
            elif options.formula == 'simpson':
                result_targets = diana_comparison.calculate_simpson_index(targets1, targets2)
            target_comparisons.append(result_targets)
            dc_to_target_ji[index] = result_targets

            bio_proc1 = get_results_from_dict_of_sets(targets1, target_to_bio_processes)
            bio_proc2 = get_results_from_dict_of_sets(targets2, target_to_bio_processes)
            if options.formula == 'jaccard':
                result_bp = diana_comparison.calculate_jaccard_index(bio_proc1, bio_proc2)
            elif options.formula == 'simpson':
                result_bp = diana_comparison.calculate_simpson_index(bio_proc1, bio_proc2)
            bp_comparisons.append(result_bp)
            dc_to_bp_ji[index] = result_bp

            pathways1 = get_results_from_dict_of_sets(targets1, target_to_pathways)
            pathways2 = get_results_from_dict_of_sets(targets2, target_to_pathways)
            if options.formula == 'jaccard':
                result_pathways = diana_comparison.calculate_jaccard_index(pathways1, pathways2)
            elif options.formula == 'simpson':
                result_pathways = diana_comparison.calculate_simpson_index(pathways1, pathways2)
            pathway_comparisons.append(result_pathways)
            dc_to_pathway_ji[index] = result_pathways

    # Plot distribution of comparisons of targets
    n, bins, patches = plt.hist(np.array(target_comparisons), bins=50, weights=np.zeros_like(np.array(target_comparisons)) + 1. / np.array(target_comparisons).size, facecolor='r')
    plt.xlabel('{} Index of Targets'.format(options.formula.capitalize()))
    plt.ylabel('Relative frequency')
    plt.title('Distribution of {} Index of Targets in drug combinations'.format(options.formula.capitalize()))
    plt.savefig(plot_ji_targets, format=fig_format, dpi=300)
    plt.clf()

    # Plot distribution of comparisons of biological processes
    n, bins, patches = plt.hist(np.array(bp_comparisons), bins=50, weights=np.zeros_like(np.array(bp_comparisons)) + 1. / np.array(bp_comparisons).size, facecolor='b')
    plt.xlabel('{} Index of Biological Processes'.format(options.formula.capitalize()))
    plt.ylabel('Relative frequency')
    plt.title('Distribution of {} Index of Biological Processes in drug combinations'.format(options.formula.capitalize()))
    plt.savefig(plot_ji_bp, format=fig_format, dpi=300)
    plt.clf()

    # Plot distribution of comparisons of pathways
    n, bins, patches = plt.hist(np.array(pathway_comparisons), bins=50, weights=np.zeros_like(np.array(pathway_comparisons)) + 1. / np.array(pathway_comparisons).size, facecolor='g')
    plt.xlabel('{} Index of Pathways'.format(options.formula.capitalize()))
    plt.ylabel('Relative frequency')
    plt.title('Distribution of {} Index of Pathways in drug combinations'.format(options.formula.capitalize()))
    plt.savefig(plot_ji_pathways, format=fig_format, dpi=300)
    plt.clf()


    #------------------------------------#
    #   CLASSIFY THE DRUG COMBINATIONS   #
    #------------------------------------#

    # Similar targets   --> ji >  0.25
    # Different targets --> ji <= 0.25
    target_cut_off = 0.5

    # Similar biological processes   --> ji >= 0.25
    # Different biological processes --> ji <  0.25
    bp_cut_off = 0.5

    # Similar pathways   --> ji >= 0.5
    # Different pathways --> ji < 0.5
    pathway_cut_off = 0.5

    classification_tar_bp = {}

    st = 0
    dt = 0
    st_sbp = 0
    st_dbp = 0
    dt_sbp = 0
    dt_dbp = 0
    for dc in dc_to_target_ji:

        # Classify by targets and biological processes
        if dc in dc_to_bp_ji:

            ji_tar = dc_to_target_ji[dc]
            ji_bp = dc_to_bp_ji[dc]

            if ji_tar > target_cut_off:
                classification_tar_bp[dc] = 'similar_targets'
                st += 1
                if ji_bp > bp_cut_off:
                    st_sbp += 1
                elif ji_bp <= bp_cut_off:
                    st_dbp += 1
            elif ji_tar <= target_cut_off:
                dt += 1
                if ji_bp > bp_cut_off:
                    dt_sbp += 1
                    classification_tar_bp[dc] = 'different_targets_similar_bp'
                elif ji_bp <= bp_cut_off:
                    dt_dbp += 1
                    classification_tar_bp[dc] = 'different_targets_different_bp'

    print('Similar targets {}: similar bp {}, diff bp {}\n'.format(st, st_sbp, st_dbp))
    print('Different targets {}: similar bp {}, diff bp {}\n'.format(dt, dt_sbp, dt_dbp))

    cPickle.dump(classification_tar_bp, open(classification_targets_bp_file, 'w'))

    classification_tar_pathway = {}

    st = 0
    dt = 0
    st_spath = 0
    st_dpath = 0
    dt_spath = 0
    dt_dpath = 0
    for dc in dc_to_target_ji:

        # Classify by targets and biological processes
        if dc in dc_to_pathway_ji:

            ji_tar = dc_to_target_ji[dc]
            ji_path = dc_to_pathway_ji[dc]

            if ji_tar > target_cut_off:
                classification_tar_pathway[dc] = 'similar_targets'
                st += 1
                if ji_path > pathway_cut_off:
                    st_spath += 1
                elif ji_path <= pathway_cut_off:
                    st_dpath += 1
            elif ji_tar <= target_cut_off:
                dt += 1
                if ji_path > pathway_cut_off:
                    dt_spath += 1
                    classification_tar_pathway[dc] = 'different_targets_similar_pathways'
                elif ji_path <= pathway_cut_off:
                    dt_dpath += 1
                    classification_tar_pathway[dc] = 'different_targets_different_pathways'

    print('Similar targets {}: similar pathways {}, diff pathways {}\n'.format(st, st_spath, st_dpath))
    print('Different targets {}: similar pathways {}, diff pathways {}\n'.format(dt, dt_spath, dt_dpath))

    cPickle.dump(classification_tar_pathway, open(classification_targets_pathways_file, 'w'))


    # Get number of drugs in drug combinations per number of targets
    targets = [len(drugbank_to_targets[drug]) for drug in drugbank_to_targets if drug in all_drugs]
    numtargets_to_numdrugs = {}
    for target in targets:
        numtargets_to_numdrugs.setdefault(target, 0)
        numtargets_to_numdrugs[target] += 1

    print('Number of drugs in drug combination: {}. Divided by four: {}'.format(len(all_drugs), len(all_drugs)/4))
    for numtar, numdrug in sorted(numtargets_to_numdrugs.iteritems(), key=lambda (x, y): x, reverse = True):
        print(numtar, numdrug)

    # End marker for time
    end = time.time()
    print('\n  DIANA INFO:\tTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))



    return




#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def check_file(file):
    """
    Checks if a file exists and if not, raises FileNotFound exception
    """
    if not fileExist(file):
        raise FileNotFound(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


def check_directory(directory):
    """
    Checks if a directory exists and if not, raises DirNotFound exception
    """
    try:
        os.stat(directory)
    except:
        raise DirNotFound(directory)


class FileNotFound(Exception):
    """
    Exception raised when a file is not found.
    """
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return 'The file {} has not been found.\nTherefore, the comparison cannot be performed. Please, check that all the profiles have been correctly generated.\n'.format(self.file)


class DirNotFound(Exception):
    """
    Exception raised when a directory is not found.
    """
    def __init__(self, directory):
        self.directory = directory

    def __str__(self):
        return 'The directory {} has not been found.\nTherefore, the comparison cannot be performed. Please, check that all the parameters have been correctly introduced and the profiles have been correctly generated.\n'.format(self.directory)


def get_results_from_dict_of_sets(list_of_elements, dict_of_sets):
    """
    We have a list of elements that are in a dict of elements, and every element have a set with results.
    We want to extract the results corresponding to our elements.
    """
    results = set()
    for element in list_of_elements:
        if element in dict_of_sets:
            for result in dict_of_sets[element]:
                results.add(result)
    return results



if  __name__ == "__main__":
    main()


