import argparse
import copy
import cPickle
import matplotlib.pyplot as plt
import ntpath
import numpy as np
import pandas as pd
import pylab
#from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, xlabel, ylabel
import scipy
import time
import sys, os, re

from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.grid_search import GridSearchCV
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
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
    parser.add_argument('-cr','--crossings_file',dest='crossings_file',action = 'store',
                        help = """Define the file where the drug crossings to be explored have been written""")
    parser.add_argument('-sif','--sif_file',dest='sif',action = 'store',
                        help = """" Input file with the protein-protein interaction network in SIF format that will be used in the experiment. """)
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
    parser.add_argument('-pca','--pca',dest='pca',action = 'store_true',
                        help = """" Make a PCA to reduce dimensionality. """)
    parser.add_argument('-cp','--comparison',dest='comparison_other_methods',action = 'store_true',
                        help = """" If we are considering a dataset to compare with other methods. """)
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

    print("\n\t\t----------------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Analysis of results: Selection of classifier\n")
    print("\t\t----------------------------------------------------------------------------------------------------------------------------\n")

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')

    # Check the directory of the profiles and comparisons
    data_dir = os.path.join(options.workspace, "profiles")
    check_directory(data_dir)

    results_dir = os.path.join(options.workspace, "comparisons")
    check_directory(results_dir)

    # Create a directory for the analysis inside the workspace
    analysis_dir = os.path.join(options.workspace, "analysis")
    create_directory(analysis_dir)

    # Create a directory for the analysis of the comparison with other methods
    if options.comparison_other_methods:
        analysis_dir = os.path.join(options.workspace, "analysis_comparison")
        create_directory(analysis_dir)

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

    # Change the name of the output file if we are doing a comparison with other methods
    if options.comparison_other_methods:
        output_dataframe = os.path.join(analysis_dir, 'comparison_other_methods.csv')

    if not fileExist(output_dataframe):

        # Create a data frame to store the results
        df = pd.DataFrame(columns=columns)

        # Prepare files
        network_filename = ntpath.basename(options.sif)
        drugbank2targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
        drug2targets = cPickle.load(open(drugbank2targets_file))

        # Open the crossings file
        crossings_file = options.crossings_file
        with open(crossings_file, 'r') as crossings_file_fd:
            for line in crossings_file_fd:
                crossing = line.strip()
                drug1, drug2 = crossing.split('---')

                # Get drug IDs
                targets1 = list(drug2targets[drug1.upper()])
                drug_id1 = diana_drug.generate_drug_id(drug1, targets1, network_filename)
                targets2 = list(drug2targets[drug2.upper()])
                drug_id2 = diana_drug.generate_drug_id(drug2, targets2, network_filename)

                # Check results table
                comparison = '{}---{}'.format(drug_id1, drug_id2)
                comparison_dir = os.path.join(results_dir, comparison)
                results_table = os.path.join(comparison_dir, 'results_table.tsv')
                if not fileExist(results_table):
                    print('The comparison of {} ({}) and {} ({}) has not been executed properly!\n'.format(drug1, drug_id1, drug2, drug_id2))
                    sys.exit(10)

                if crossing in pair2comb:
                    combination_field = pair2comb[crossing]
                else:
                    print('The comparison {} is not in the pair2comb dictionary!\n'.format(crossing))
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


    #-----------------------------------------------------------#
    #   DIVIDE THE DATASET IN A TRAINING AND A VALIDATION SET   #
    #-----------------------------------------------------------#

    training_dataframe = os.path.join(analysis_dir, 'dcdb_comparisons_training.csv')
    validation_dataframe = os.path.join(analysis_dir, 'dcdb_comparisons_validation.csv')
    proportion_training = 0.8

    # Change the name of the output file if we are doing a comparison with other methods
    if options.comparison_other_methods:
        training_dataframe = os.path.join(analysis_dir, 'comparison_other_methods_training.csv')
        validation_dataframe = os.path.join(analysis_dir, 'comparison_other_methods_validation.csv')

    if not fileExist(training_dataframe) or not fileExist(validation_dataframe):

        num_dc_training = int(round(num_dc*proportion_training))
        num_ndc_training = int(round(num_ndc*proportion_training))
        print('Training set (positives): {} out of {} ({}%)\n'.format(num_dc_training, num_dc, proportion_training*100))
        print('Training set (negatives): {} out of {} ({}%)\n'.format(num_ndc_training, num_ndc, proportion_training*100))

        dc_data_training = dc_data.sample(n=num_dc_training) # Get a random sample
        ndc_data_training = ndc_data.sample(n=num_ndc_training)
        dc_data_validation = dc_data.loc[~dc_data.index.isin(dc_data_training.index)]  # Remove the sample that we have taken from the dataframe
        ndc_data_validation = ndc_data.loc[~ndc_data.index.isin(ndc_data_training.index)]

        df_training = pd.concat([dc_data_training, ndc_data_training])
        df_validation = pd.concat([dc_data_validation, ndc_data_validation])

        # Output the Pandas dataframes in a CSV file
        df_training.to_csv(training_dataframe)
        df_validation.to_csv(validation_dataframe)

        # Define the variables for the training dataset
        df = df_training
        dc_data = df[df['combination'] == 1]
        ndc_data = df[df['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)
        print('Number of drug combinations after getting the training dataset:\t{}\n'.format(num_dc))
        print('Number of non-drug combinations after getting the training dataset:\t{}\n'.format(num_ndc))

    else:
        df_training = pd.read_csv(training_dataframe, index_col=0)
        df_validation = pd.read_csv(validation_dataframe, index_col=0)

        # Define the variables for the training dataset
        df = df_training
        dc_data = df[df['combination'] == 1]
        ndc_data = df[df['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)
        print('Number of drug combinations after getting the training dataset:\t{}\n'.format(num_dc))
        print('Number of non-drug combinations after getting the training dataset:\t{}\n'.format(num_ndc))



    #------------------------------------------------------------------#
    #   SELECT RELEVANT FEATURES / REDUCE DIMENSIONALITY OF THE DATA   #
    #------------------------------------------------------------------#

    if options.pca:

        # Strategy:
        # We calculate the explained variance ratio for all the features.
        # We define a a cut-off threshold of the minimum variance ratio that we consider relevant.
        # We will count the number of features with explained variance higher than the cut-off defined.
        # Then, we will reduce the dimensionality to the number of features with variance higher than the cut-off.

        variance_cut_off = 0.01
        num_components = 0
        df_raw = df.drop('combination', axis=1)
        raw_columns = copy.copy(columns)
        raw_columns.remove('combination')
        pca = PCA(n_components=None)
        pca.fit(df_raw)
        values_trans = pca.transform(df_raw)
        explained_variance = pca.explained_variance_ratio_
        for var in explained_variance:
            if var > variance_cut_off:
                num_components += 1

        if num_components < len(raw_columns):

            print('Number of features:\t{}\n'.format(len(raw_columns)))
            print('Reduction to {} components\n'.format(num_components))

            pca = PCA(n_components=num_components)
            pca.fit(df_raw)
            values_trans = pca.transform(df_raw)
            indexes = df.index.values
            df_trans = pd.DataFrame.from_records(values_trans, index=indexes)
            df_comb = df[['combination']]
            df_new = pd.concat([df_trans, df_comb], axis=1)
            df = df_new
            dc_data = df[df['combination'] == 1]
            ndc_data = df[df['combination'] == 0]
            num_dc = len(dc_data.index)
            num_ndc = len(ndc_data.index)

    else:

        # Manually introduced features
        guild_thresholds = [1, 5]
        rank_scoring = ['spearman', 'dot_product']
        list_scoring = ['jaccard']
        selected_columns = diana_analysis.obtain_columns_best_features(guild_thresholds, rank_scoring, list_scoring, ATC_SE=consider_se)
        print('Selected columns: {}\n'.format(', '.join(selected_columns)))
        print('Number of selected features: {}\n'.format(len(selected_columns)-1)) # We take away the combinations column

        # Define the new table with the selected columns
        df = df[selected_columns]
        dc_data = df[df['combination'] == 1]
        ndc_data = df[df['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)
 

    #------------------------------------------#
    #   TUNE THE ALGORITHM OF THE CLASSIFIER   #
    #------------------------------------------#

    tables_dir = os.path.join(analysis_dir, 'tables')
    create_directory(tables_dir)
    results_table = os.path.join(tables_dir, 'tuning_results.tsv')
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
        'QuadraticDiscr.' : QuadraticDiscriminantAnalysis()
    }
    classifier = 'SVC'
    pipe_svc = Pipeline([('slc', StandardScaler()),
                         ('clf', SVC(random_state=1))])

    param_range_c = [1.0, 10.0, 100]
    param_range_gamma = [1e-4, 1e-3, 0.01, 0.1]
    param_grid = [{'clf__C': param_range_c,
                   'clf__kernel': ['linear']},
                  {'clf__C': param_range_c,
                   'clf__gamma': param_range_gamma,
                   'clf__kernel': ['rbf']}]


    print('TUNNING THE ALGORITHM OF {}\n'.format(classifier.upper()))
    rounds = 2
    repetitions = 25
    dict_results = {}

    for n_round in xrange(rounds):

        print('ROUND NUMBER {}\n'.format(n_round+1))

        # Obtain the different non-drug combination groups to repeat the analysis
        ndc_training_groups = diana_analysis.obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times

        for ndc_training_data in ndc_training_groups:

            merged_groups = pd.concat([dc_data, ndc_training_data])
            X_train, y_train = merged_groups.iloc[:, :-1], merged_groups.iloc[:, -1]
            grid_search = GridSearchCV(estimator=pipe_svc,
                          param_grid=param_grid,
                          scoring='accuracy',
                          cv=10,
                          n_jobs=-1)
            grid = grid_search.fit(X_train, y_train)
            print(grid)
            # summarize the results of the grid search
            print('Grid best score: {}'.format(grid.best_score_))
            result = str(grid.best_params_)
            print('Grid best parameters: {}\n'.format(result))
            dict_results.setdefault(result, 0)
            dict_results[result] += 1

    print('\nFINAL RESULT\n')

    with open(results_table, 'w') as results_table_fd:

        for param_comb in sorted(dict_results, reverse = True):

            print('{}\t{}\n'.format(param_comb, dict_results[param_comb]))
            results_table_fd.write('{}\t{}\n'.format(param_comb, dict_results[param_comb]))



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


