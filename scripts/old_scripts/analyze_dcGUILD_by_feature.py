import argparse
import copy
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
#from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, xlabel, ylabel
import scipy
import time
import sys, os, re

from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier

from context import diana
import diana.classes.drug as diana_drug
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

    print("\n\t\t------------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Analysis of results: Analysis by targets\n")
    print("\t\t------------------------------------------------------------------------------------------------------------------------\n")

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
            drug1 = drug_id1.split('_')[0].upper()
            drug2 = drug_id2.split('_')[0].upper()
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

            results = get_results_from_table(results_table, columns, combination_field)

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



    #-------------------------#
    #   EVALUATE PERFORMANCE  #
    #-------------------------#

    img_dir = os.path.join(analysis_dir, 'figures')
    create_directory(img_dir)
    fig_format = 'png'

    tables_dir = os.path.join(analysis_dir, 'tables')
    create_directory(tables_dir)

    # Machine learning parameters
    repetitions = 25 # Number of repetititons
    n_fold = 2     # Number of folds
    min_num_dc_group = 10
    greater_or_smaller = 'greater'
    classifier = 'SVC'
    classifiers = {
        'KNeighbors' : KNeighborsClassifier(3),
        'SVC' : SVC(probability=True),
        'SVC linear' : SVC(kernel="linear", C=0.025),
        'SVC rbf' : SVC(gamma=2, C=1),
        'DecisionTree' : DecisionTreeClassifier(max_depth=5),
        'RandomForest' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        'MLP' : MLPClassifier(alpha=1),
        'AdaBoost' : AdaBoostClassifier(),
        'GaussianNB' : GaussianNB(),
        'QuadraticDiscr.' : QuadraticDiscriminantAnalysis(),
        'SVC best 1' : SVC(kernel="linear", C=0.1, probability=True),
        'SVC best 2' : SVC(kernel="rbf", gamma=0.01, C=100.0, probability=True)
    }

    # Plot of distributions of AUC
    plot_name = os.path.join(img_dir, 'dcGUILD_1_threshold_auc.{}'.format(fig_format))

    # Get the targets file
    drugbank_to_targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    drugbank_to_targets = cPickle.load(open(drugbank_to_targets_file))

    # Get the DIANA IDs file
    diana_id_to_drugbank_file = os.path.join(toolbox_dir, 'diana_id_to_drugbank.pcl')
    diana_id_to_drugbank = cPickle.load(open(diana_id_to_drugbank_file))


    print('\nEVALUATION OF DCGUILD\n')
    repetitions = 25
    n_fold = 10
    analysis_results = {}

    # Obtain the different non-drug combination groups to repeat the analysis
    ndc_repetitions = diana_analysis.obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times


    # dcGUILD_features = [str(x) for x in threshold_list]
    # dcGUILD_feature_to_columns = {}
    # # Get dcGUILD columns
    # for top_threshold in threshold_list:
    #     for data_type in ['node', 'edge', 'function']:
    #         for scoring_function in ['dot_product', 'spearman', 'jaccard']:
    #             col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
    #             dcGUILD_feature_to_columns.setdefault(str(top_threshold), [])
    #             dcGUILD_feature_to_columns[str(top_threshold)].append(col)
    #     dcGUILD_feature_to_columns[str(top_threshold)].append('combination')

    dcGUILD_features = []
    dcGUILD_feature_to_columns = {}
    # Get dcGUILD columns
    for top_threshold in [1]:
        for data_type in ['node', 'edge', 'function']:
            for scoring_function in ['dot_product', 'spearman', 'jaccard']:
                col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
                dcGUILD_features.append(col)
                dcGUILD_feature_to_columns[col]=[col, 'combination']


    for feature in dcGUILD_features:

        df_method = df[dcGUILD_feature_to_columns[feature]]

        dc_data = df_method[df_method['combination'] == 1]
        ndc_data = df_method[df_method['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)

        print(feature)
        print('Building {} repetition groups of {} (same) DC and {} (different) non-DC'.format(repetitions,num_dc,num_dc))
        ndc_repetitions = diana_analysis.obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times

        mean_aucs = [] # Here we will store the means of AUCs from the cross-validations
        std_aucs = [] # Here we will store the standard deviations of the AUCs from the cross-validations
        all_aucs = [] # Here we will store ALL the AUCs
        all_probs = [] # Here we store all the probabilities and labels

        num_repetitions=0
        for ndc_data_equal in ndc_repetitions:

            num_repetitions+=1
            num_items_group = int( float(num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation
            if num_repetitions == 1:
                print('Building {} fold groups of {} DC and {} non-DC x {} repetitions'.format(n_fold,num_items_group,num_items_group, repetitions))

            dc_groups = diana_analysis.obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group, me_too_drug_combinations) # Defining the drug combination groups in each cross-validation step
            ndc_groups = diana_analysis.obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group, me_too_drug_combinations) # Defining the non-drug combination groups in each cross-validation step
            merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

            mean, var, std, list_auc, list_prob = diana_analysis.run_nfold_crossvalidation_scikit_with_prob(n_fold, merged_groups, classifiers[classifier])

            mean_aucs.append(mean)
            std_aucs.append(std)
            all_aucs = all_aucs + list_auc
            all_probs = all_probs + list_prob

        final_mean = np.mean(mean_aucs)
        mean_std = np.mean(std_aucs)
        std_means = np.std(mean_aucs)
        std = np.std(all_aucs)
        print('FINAL MEAN: {}'.format(final_mean))
        print('MEAN of STD: {}'.format(mean_std))
        print('STD: {}\n'.format(std))

        # Store the distribution of AUCs in the dictionary
        analysis_results[feature] = all_aucs



    #------------------------------#
    #   PLOT DISTRIBUTION OF AUC   #
    #------------------------------#

    fig = pylab.figure(dpi=300)
    ax = pylab.axes()
    #pylab.hold(True)
    pos = 1
    col_num = 0

    xticks = [] # Define the places in which the labels will be
    xlabels = [] # Define the labels (the names of the features)
    #colors = [ ['#9ed0ff, blue'], ['#32f232', 'green'], ['#fbc562', '#d48900'], ['#ff7373', '#b80000'], ['grey', 'black'] ]

    for feature in dcGUILD_features:

        positions = []
        positions.append(pos) # Define the positions of the boxplots
        pos+=2 # Add separation between boxplots
        xlabels.append(feature) # Add the feature used at the x axis

        # Boxplot group
        #bp = boxplot(data, positions = positions, widths = 0.6)
        bp = pylab.boxplot(analysis_results[feature], positions = positions, widths = 0.6, patch_artist=True)

        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # Set axes limits and labels
    pylab.xlim(0,pos-1)
    pylab.ylim(0,1)
    ax.set_xticklabels(xlabels)
    ax.set_xticks(xticks)
    pylab.xlabel('Features')
    pylab.ylabel('Distribution of AUC values')

    fig.autofmt_xdate()
    pylab.savefig(plot_name, format=fig_format)
    pylab.show()


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


if  __name__ == "__main__":
    main()
