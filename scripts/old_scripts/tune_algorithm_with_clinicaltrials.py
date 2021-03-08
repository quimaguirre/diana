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
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.grid_search import GridSearchCV
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier

from context import diana
import diana.classes.drug as diana_drug
import diana.classes.network_analysis as network_analysis




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

    print("\n\t\t---------------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Analysis of results: Tunning the classifier\n")
    print("\t\t---------------------------------------------------------------------------------------------------------------------------\n")

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolbox_dir = os.path.join(main_path, 'diana/toolbox')

    # Check the directory of the profiles and comparisons
    data_dir = os.path.join(options.workspace, "profiles")
    check_directory(data_dir)

    results_dir = os.path.join(options.workspace, "comparisons")
    check_directory(results_dir)

    results_ct_dir = os.path.join(options.workspace, "comparisons_clinicaltrials")
    check_directory(results_ct_dir)

    # Get the list of thresholds to create the profiles
    if options.threshold_list and fileExist(options.threshold_list):
        threshold_list = get_values_from_threshold_file(options.threshold_list)
    else:
        threshold_list = [1, 5, 10, 20, 50]

    # Get the names of the columns
    columns = obtain_columns(threshold_list, clinicaltrials=True)



    #-----------------------------------------------------#
    #   PARSE THE RESULTS AND CREATE A PANDAS DATAFRAME   #
    #-----------------------------------------------------#

    pair2comb_file = os.path.join(toolbox_dir, 'pair2comb.pcl')
    pair2comb = cPickle.load(open(pair2comb_file))

    ddi = sum(1 for x in pair2comb.values() if x == 1)
    non_ddi = sum(1 for x in pair2comb.values() if x == 0)

    print('NUMBER OF DRUG COMBINATIONS:\t\t{}\n'.format(ddi))
    print('NUMBER OF NON-DRUG COMBINATIONS:\t{}\n'.format(non_ddi))

    output_dataframe = os.path.join(options.workspace, 'dcdb_comparisons.csv')

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

            results = get_results_from_table(results_table, columns, combination_field, clinicaltrials_field = 0)

            df2 = pd.DataFrame([results], columns=columns, index=[comparison])
            # Add the information to the main data frame
            df = df.append(df2)

        # Output the Pandas dataframe in a CSV file
        df.to_csv(output_dataframe)

    else:
        df = pd.read_csv(output_dataframe, index_col=0)



    #-----------------------------------------------------------------#
    #   PARSE CLINICAL TRIALS RESULTS AND CREATE A PANDAS DATAFRAME   #
    #-----------------------------------------------------------------#

    pair2comb_ct_file = os.path.join(toolbox_dir, 'pair2comb_clinicaltrials.pcl')
    pair2comb_ct = cPickle.load(open(pair2comb_ct_file))
    output_dataframe = os.path.join(options.workspace, 'clinicaltrials_comparisons.csv')

    if not fileExist(output_dataframe):


        # Obtain all the results subfolders of the results main folder
        results_dir_list = [f for f in os.listdir(results_ct_dir) if os.path.isdir(os.path.join(results_ct_dir, f))]

        for comparison in results_dir_list:

            drug_id1, drug_id2 = comparison.split('---')
            comparison_dir = os.path.join(results_ct_dir, comparison)
            results_table = os.path.join(comparison_dir, 'results_table.tsv')

            # Add the Comb field (if it is drug combination or not)
            drug1 = drug_id1.split('_')[0].upper()
            drug2 = drug_id2.split('_')[0].upper()
            comparison_without_id = '{}---{}'.format(drug1, drug2)
            if comparison_without_id in pair2comb_ct:
                combination_field = pair2comb_ct[comparison_without_id]
            else:
                print('The comparison {} is not in the pair2comb_ct dictionary!\n'.format(comparison_without_id))
                print(pair2comb_ct)
                sys.exit(10)

            # Obtain the comparisons that are in DCDB and change their value by 2.0
            # (indicating Clinical Trial drug combination)
            if comparison in df.index:
                #if df.loc[comparison]['combination'] == 0:
                df.set_value(comparison, 'combination', combination_field)
                df.set_value(comparison, 'clinicaltrials', 1)
                #print(df.loc[comparison]['combination'])
                continue

            if not fileExist(results_table):
                print('The comparison {} has not been executed properly!\n'.format(comparison))
                sys.exit(10)

            results = get_results_from_table(results_table, columns, combination_field, clinicaltrials_field = 1)

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
    df = df.replace(to_replace={'dcstructure':{'None':np.nan}})

    # Remove the nan values in dcstructure
    df = df.dropna()

    # Count the number of drug combinations / non-drug combinations
    dcdb_df = df[df['clinicaltrials'] == 0]
    ct_df = df[df['clinicaltrials'] == 1]
    dc_data = dcdb_df[dcdb_df['combination'] == 1]
    ndc_data = dcdb_df[dcdb_df['combination'] == 0]
    dc_ct_data = ct_df[ct_df['combination'] == 1] # Clinical Trials data
    ndc_ct_data = ct_df[ct_df['combination'] == 0] # Clinical Trials data
    num_dc = len(dc_data.index)
    num_ndc = len(ndc_data.index)
    num_ct_dc = len(dc_ct_data.index)
    num_ct_ndc = len(ndc_ct_data.index)
    print('Number of drug combinations after removing missing values:\t{}\n'.format(num_dc))
    print('Number of non-drug combinations after removing missing values:\t{}\n'.format(num_ndc))
    print('Number of drug combinations after removing missing values in Clinical Trials:\t{}\n'.format(num_ct_dc))
    print('Number of non-drug combinations after removing missing values in Clinical Trials:\t{}\n'.format(num_ct_ndc))



    #---------------------------#
    #   IDENTIFY ME-TOO DRUGS   #
    #---------------------------#

    me_too_dir = os.path.join(options.workspace, 'me_too_drugs')
    create_directory(me_too_dir)
    me_too_drugs_table = os.path.join(me_too_dir, 'me_too_drugs.tsv')
    me_too_drug_combs_table = os.path.join(me_too_dir, 'me_too_drug_combinations.tsv')

    me_too_drug_pairs_file = os.path.join(me_too_dir, 'me_too_drug_pairs.pcl')
    me_too_drug_comb_pairs_file = os.path.join(me_too_dir, 'me_too_drug_comb_pairs.pcl')

    if not fileExist(me_too_drug_pairs_file) or not fileExist(me_too_drug_comb_pairs_file):

        df_struc = df[['dcstructure']]
        df_struc = df_struc.astype(float)
        me_too_drug_pairs, me_too_drug_comb_pairs = obtain_me_too_drugs_and_combinations(df_struc, columns, me_too_drugs_table, me_too_drug_combs_table)
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
    dcdb_df = df[df['clinicaltrials'] == 0]
    ct_df = df[df['clinicaltrials'] == 1]
    dc_data = dcdb_df[dcdb_df['combination'] == 1]
    ndc_data = dcdb_df[dcdb_df['combination'] == 0]
    dc_ct_data = ct_df[ct_df['combination'] == 1] 
    ndc_ct_data = ct_df[ct_df['combination'] == 0]
    num_dc = len(dc_data.index)
    num_ndc = len(ndc_data.index)
    num_ct_dc = len(dc_ct_data.index)
    num_ct_ndc = len(ndc_ct_data.index)
    print('Number of drug combinations after removing me-too conflictive drug pairs:\t{}\n'.format(num_dc))
    print('Number of non-drug combinations after removing me-too conflictive drug pairs:\t{}\n'.format(num_ndc))
    print('Number of drug combinations after removing me-too conflictive drug pairs in Clinical Trials:\t{}\n'.format(num_ct_dc))
    print('Number of non-drug combinations after removing me-too conflictive drug pairs in Clinical Trials:\t{}\n'.format(num_ct_ndc))




    #------------------------------------------------------------------#
    #   SELECT RELEVANT FEATURES / REDUCE DIMENSIONALITY OF THE DATA   #
    #------------------------------------------------------------------#

    # Strategy:
    # We calculate the explained variance ratio for all the features.
    # We define a a cut-off threshold of the minimum variance ratio that we consider relevant.
    # We will count the number of features with explained variance higher than the cut-off defined.
    # Then, we will reduce the dimensionality to the number of features with variance higher than the cut-off.

    variance_cut_off = 0.01
    num_components = 0
    df_raw = df.drop('combination', axis=1)
    df_raw = df_raw.drop('clinicaltrials', axis=1)
    raw_columns = copy.copy(columns)
    raw_columns.remove('combination')
    raw_columns.remove('clinicaltrials')
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
        df_clin = df[['clinicaltrials']]
        df_new = pd.concat([df_trans, df_comb], axis=1)
        df_new = pd.concat([df_new, df_clin], axis=1)
        df = df_new
        dcdb_df = df[df['clinicaltrials'] == 0]
        dcdb_df = dcdb_df.drop('clinicaltrials', axis=1)
        ct_df = df[df['clinicaltrials'] == 1]
        ct_df = ct_df.drop('clinicaltrials', axis=1)
        dc_data = dcdb_df[dcdb_df['combination'] == 1]
        ndc_data = dcdb_df[dcdb_df['combination'] == 0]
        dc_ct_data = ct_df[ct_df['combination'] == 1] 
        ndc_ct_data = ct_df[ct_df['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)
        num_ct_dc = len(dc_ct_data.index)
        num_ct_ndc = len(ndc_ct_data.index)




    #------------------------------------------#
    #   TUNE THE ALGORITHM OF THE CLASSIFIER   #
    #------------------------------------------#

    tables_dir = os.path.join(options.workspace, 'tables')
    create_directory(tables_dir)
    results_table = os.path.join(tables_dir, 'tuning_results.tsv')
    classifier = 'SVC rbf'

    pipe_svc = Pipeline([('slc', StandardScaler()),
                         ('clf', SVC(random_state=1))])

    param_range = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]
    param_grid = [{'clf__C': param_range,
                   'clf__kernel': ['linear']},
                  {'clf__C': param_range,
                   'clf__gamma': param_range,
                   'clf__kernel': ['rbf']}]


    print('TUNNING THE ALGORITHM OF {}\n'.format(classifier.upper()))
    rounds = 10
    repetitions = 10
    dict_results = {}

    for n_round in xrange(rounds):

        print('ROUND NUMBER {}\n'.format(n_round+1))

        # Obtain the different non-drug combination groups to repeat the analysis
        ndc_training_groups = obtain_n_groups_of_k_length(ndc_ct_data, repetitions, num_ct_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times

        for ndc_training_data in ndc_training_groups:

            merged_groups = pd.concat([dc_ct_data, ndc_training_data])
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


def obtain_header_fields(first_line, separator='\t'):
    """ 
    Obtain a dictionary: "field_name" => "position" 
    """
    fields_dict = {}
    header_fields = first_line.strip().split(separator)
    for x in xrange(0, len(header_fields)):
        fields_dict[header_fields[x].lower()] = x
    return fields_dict


def obtain_columns(threshold_list, clinicaltrials=False):
    """ 
    Obtain the name of the column from the values of the method, data_type and threshold.
    """

    columns = []

    for data_type in ['target', 'pfam', 'function']:
        for scoring_function in ['dot_product', 'spearman', 'jaccard']:
            col = 'dct'+'_'+data_type+'_'+scoring_function
            columns.append(col)
    
    for top_threshold in threshold_list:
        for data_type in ['node', 'edge', 'function']:
            for scoring_function in ['dot_product', 'spearman', 'jaccard']:
                col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
                columns.append(col)

    columns.append('dcstructure')
    columns.append('combination')
    if clinicaltrials:
        columns.append('clinicaltrials')

    return columns


def get_results_from_table(results_table, columns, combination_field, clinicaltrials_field = 0):
    """ 
    Obtain the results from the results table of the comparison.
    """

    column_to_results = {}
    with open(results_table, 'r') as results_table_fd:

        # method    data_type   threshold   dot_product spearman    jaccard
        first_line = results_table_fd.readline()
        fields_dict = obtain_header_fields(first_line, separator='\t')

        for line in results_table_fd:

            fields = line.strip().split('\t')
            method = fields[ fields_dict['method'] ]
            data_type = fields[ fields_dict['data_type'] ]
            threshold = fields[ fields_dict['threshold'] ]
            dot_product = fields[ fields_dict['dot_product'] ]
            spearman = fields[ fields_dict['spearman'] ]
            jaccard = fields[ fields_dict['jaccard'] ]

            if method == 'dctargets':
                for scoring_function, result in [['dot_product', dot_product], ['spearman', spearman], ['jaccard', jaccard]]:
                    col = 'dct'+'_'+data_type+'_'+scoring_function
                    column_to_results[col] = result
            elif method == 'dcguild':
                for scoring_function, result in [['dot_product', dot_product], ['spearman', spearman], ['jaccard', jaccard]]:
                    col = 'dcg'+'_'+data_type+'_'+str(threshold)+'_'+scoring_function
                    column_to_results[col] = result
            elif method == 'dcstructure':
                column_to_results['dcstructure'] = dot_product

    results = []
    for column in columns:
        if column in column_to_results:
            results.append(column_to_results[column])
        elif column == 'combination':
            results.append(combination_field)
        elif column == 'clinicaltrials':
            results.append(clinicaltrials_field)
        else:
            print('The column {} is not among the result columns!'.format(column))
            print('Predefined columns: {}'.format(sorted(columns)))
            print('Result columns: {}\n'.format(sorted(column_to_results.keys())))
            sys.exit(10)

    return results


def obtain_me_too_drugs_and_combinations(df, columns, me_too_drugs_table, me_too_drug_combs_table):
    """ 
    Obtain me-too drugs and me-to drug combinations in the dataset.
    """

    df_me_too = pd.DataFrame(columns=columns)

    me_too_drug_pairs = set()
    me_too_drug_comb_pairs = set()

    me_too_drugs_dict = {}
    me_too_drug_combs_dict = {}

    num_metoo_dc = 0
    num_nonmetoo_dc = 0

    done_pairings = []

    for index, row in df.iterrows():

        score = row['dcstructure']

        if score >= 0.7:
            me_too_drug_pairs.add(index)
            me_too_drugs_dict[index] = score

            df2 = pd.DataFrame([row], columns=columns, index=[index])
            df_me_too = df_me_too.append(df2)


    for index1, row1 in df_me_too.iterrows():

        score = row1['dcstructure']


        for index2, row2 in df_me_too.iterrows():

            if index1 == index2:
                continue

            combpair1 = '___'.join([index1, index2])
            combpair2 = '___'.join([index2, index1])

            if combpair1 in done_pairings or combpair2 in done_pairings:
                continue

            done_pairings.append(combpair1)
            done_pairings.append(combpair2)

            (drug11, drug12) = index1.split('---')
            (drug21, drug22) = index2.split('---')

            pairing11_1 = '---'.join([drug11, drug21])
            pairing11_2 = '---'.join([drug21, drug11])
            pairing12_1 = '---'.join([drug12, drug22])
            pairing12_2 = '---'.join([drug22, drug12])

            pairing21_1 = '---'.join([drug11, drug22])
            pairing21_2 = '---'.join([drug22, drug11])
            pairing22_1 = '---'.join([drug12, drug21])
            pairing22_2 = '---'.join([drug21, drug12])

            group1 = []
            no_pairing = False
            for possib1, possib2 in [ (pairing11_1,pairing11_2), (pairing12_1, pairing12_2) ]:

                if possib1 in df_me_too.index:
                    pairing = df_me_too.loc[[possib1]]
                    group1.append(pairing)
                elif possib2 in df_me_too.index:
                    pairing = df_me_too.loc[[possib2]]
                    group1.append(pairing)
                else:
                    #print('No pairing found!')
                    num_nonmetoo_dc+=1
                    no_pairing = True

            group2 = []
            for possib1, possib2 in [ (pairing21_1,pairing21_2), (pairing22_1, pairing22_2) ]:

                if possib1 in df_me_too.index:
                    pairing = df_me_too.loc[[possib1]]
                    group2.append(pairing)
                elif possib2 in df_me_too.index:
                    pairing = df_me_too.loc[[possib2]]
                    group2.append(pairing)
                else:
                    #print('No pairing found!')
                    if no_pairing == False:
                        num_nonmetoo_dc+=1
                    no_pairing = True

            if no_pairing:
                continue

            score11 = group1[0].iloc[0]['dcstructure']
            score12 = group1[1].iloc[0]['dcstructure']

            score21 = group2[0].iloc[0]['dcstructure']
            score22 = group2[1].iloc[0]['dcstructure']


            if (score11 < 0.7 and score12 < 0.7) or (score21 < 0.7 and score22 < 0.7):
                num_nonmetoo_dc+=1
            else:
                num_metoo_dc+=1
                me_too_drug_comb_pairs.add(combpair1)

                if (score11 >= 0.7 and score12 >= 0.7):
                    me_too_drug_combs_dict.setdefault(combpair1, {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_1', {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_2', {})
                    me_too_drug_combs_dict[combpair1]['me_too_1'][group1[0].index[0]] = score11
                    me_too_drug_combs_dict[combpair1]['me_too_2'][group1[1].index[0]] = score12
                elif (score21 < 0.7 and score22 < 0.7):
                    me_too_drug_combs_dict.setdefault(combpair2, {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_1', {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_2', {})
                    me_too_drug_combs_dict[combpair1]['me_too_1'][group2[0].index[0]] = score21
                    me_too_drug_combs_dict[combpair1]['me_too_2'][group2[1].index[0]] = score22

    print('Number of me-too drug combinations:\t{}\n'.format(num_metoo_dc))
    print('Number of non me-too drug combinations:\t{}\n'.format(num_nonmetoo_dc))

    me_too_drugs_fd = open(me_too_drugs_table, 'w')
    me_too_drug_comb_fd = open(me_too_drug_combs_table, 'w')

    # Write the results of me-too drug pairs
    for drug_pair, score in sorted(me_too_drugs_dict.iteritems(), key=lambda (x, y): y, reverse=True):
        (drug1, drug2) = drug_pair.split('---')
        name1 = '-'
        name2 = '-'
        #name1 = dcdb2name[drug1]
        #name2 = dcdb2name[drug2]
        me_too_drugs_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(drug1, name1, drug2, name2, score))

    # Write the results of me-too drug combination pairs
    for drug_comb_pair in me_too_drug_combs_dict:
        (dc1, dc2) = drug_comb_pair.split('___')
        (drug1, drug2) = dc1.split('---')
        name1 = drug1
        name2 = drug2
        #name1 = dcdb2name[drug1]
        #name2 = dcdb2name[drug2]
        (drug3, drug4) = dc2.split('---')
        name3 = drug3
        name4 = drug4
        #name3 = dcdb2name[drug3]
        #name4 = dcdb2name[drug4]

        me_too_1 = me_too_drug_combs_dict[drug_comb_pair]['me_too_1'].keys()[0]
        score1 = me_too_drug_combs_dict[drug_comb_pair]['me_too_1'][me_too_1]
        (mtd1, mtd2) = me_too_1.split('---') # mtd = me too drug
        mtn1 = mtd1
        mtn2 = mtd2
        #mtn1 = dcdb2name[mtd1] # mtn = me too drug name
        #mtn2 = dcdb2name[mtd2]

        me_too_2 = me_too_drug_combs_dict[drug_comb_pair]['me_too_2'].keys()[0]
        score2 = me_too_drug_combs_dict[drug_comb_pair]['me_too_2'][me_too_2]
        (mtd3, mtd4) = me_too_2.split('---')
        mtn3 = mtd3
        mtn4 = mtd4
        #mtn4 = dcdb2name[mtd3]
        #mtn4 = dcdb2name[mtd4]

        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name1,name2,name3,name4,mtn1,mtn2,score1,mtn3,mtn4,score2))
        me_too_drug_comb_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name1,name2,name3,name4,mtn1,mtn2,score1,mtn3,mtn4,score2))

    me_too_drugs_fd.close()
    me_too_drug_comb_fd.close()

    return me_too_drug_pairs, me_too_drug_comb_pairs


def obtain_n_groups_of_k_length(my_df, n, k, me_too_drug_combinations=None):
    """
    Obtain n number of groups of length k.
    If me_too_drug_combinations contains a list of the me-too drug combination pairs,
    it will check that there are no me-too drug combinations in different groups.
    If me_too_drug_combinations = None, it will not do anything.
    """

    repeat = True

    while repeat:
        groups = []
        repeat = False
        curr_df = my_df.copy() # Copy the dataframe, so that if we repeat we have the initial one
        for y in xrange(n):
            new_df = curr_df.sample(n=k) # Get a random sample of length k from the main dataframe
            curr_df = curr_df.loc[~curr_df.index.isin(new_df.index)]  # Remove the sample that we have taken from the main dataframe
            groups.append(new_df) # Append the sample in the list gropus

        # Check if two me-too drug combinations are part of two different groups
        # If this happens, we will repeat the process (because they could be used in training / testing at the same time)
        if me_too_drug_combinations:
            for pair in me_too_drug_combinations:
                drug_comb1, drug_comb2 = pair
                comb1_group = None
                comb2_group = None
                for x in xrange(len(groups)):
                    indexes = groups[x].index.values
                    if drug_comb1 in groups[x].index.values:
                        comb1_group = x
                    if drug_comb2 in groups[x].index.values:
                        comb2_group = x
                if comb1_group and comb2_group and comb1_group != comb2_group:
                    repeat = True
                    break

    return groups


def run_nfold_crossvalidation_scikit(n, groups, classifier):
    """
    n = number of folds
    groups = list with the balanced groups in each fold of the cross-validation
    classifier = classifier used in the machine learning approach
    """

    all_auc = []
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    stdsc = StandardScaler()
    #pca = PCA()

    for x in xrange(n):

        test = groups[x]
        train_groups = [item for index,item in enumerate(groups) if index != x]
        train = pd.concat(train_groups)

        X_train, y_train = train.iloc[:, :-1], train.iloc[:, -1]
        X_test, y_test = test.iloc[:, :-1], test.iloc[:, -1]

        X_train = stdsc.fit_transform(X_train)
        X_test = stdsc.transform(X_test)

        clf = classifier.fit(X_train, y_train)
        y_pred = clf.predict(X_test)

        if hasattr(clf, "decision_function"):
            y_score = clf.decision_function(X_test)
        else:
            prob = clf.predict_proba(X_test)
            classes = clf.classes_ # This is the order of the classes. The probabilities are given in this order
            for index in xrange(len(classes)):
                if classes[index] == 1:
                    dc_index = index # Obtain in which position is located the probability of being drug combination
            y_score = []
            for p in xrange(len(prob)):
                dc_prob = prob[p][dc_index] # We use the index to obtain the probability of being drug combination
                y_score.append(dc_prob) # Append the array in all_prob

        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        auc = metrics.roc_auc_score(y_test, y_pred)
        #print('SCIKIT AUC: {}\n'.format(auc))
        all_auc.append(auc)

    mean_tpr /= n
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    mean_auc2 = np.mean(all_auc)
    #print('Mean AUC: {}'.format(mean_auc))
    #print('Mean AUC2: {}'.format(mean_auc2))

    var_auc = np.var(all_auc)
    std_auc = np.std(all_auc)
    #print('Var AUC: {}'.format(var_auc))
    #print('Std AUC: {}'.format(std_auc))

    return mean_auc, var_auc, std_auc, all_auc


def run_machine_learning_scikit(classifier, repetitions, positive_training_data, negative_training_groups, positive_testing_data, negative_testing_groups):
    """
    classifier = classifier used in the machine learning approach
    repetitions = number of repetitions
    positive_training_data = pandas table with the positive training data
    negative_training_groups = list of pandas tables with the negative training data
    positive_testing_data = pandas table with the positive testing data
    negative_testing_groups = list of pandas tables with the negative testing data
    """

    all_auc = []
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    stdsc = StandardScaler()
    #pca = PCA()

    for x in xrange(repetitions):

        negative_training_data = negative_training_groups[x]
        negative_testing_data = negative_testing_groups[x]
        train = pd.concat([positive_training_data,negative_training_data])
        test = pd.concat([positive_testing_data,negative_testing_data])

        X_train, y_train = train.iloc[:, :-1], train.iloc[:, -1]
        X_test, y_test = test.iloc[:, :-1], test.iloc[:, -1]

        X_train = stdsc.fit_transform(X_train)
        X_test = stdsc.transform(X_test)

        clf = classifier.fit(X_train, y_train)
        y_pred = clf.predict(X_test)

        if hasattr(clf, "decision_function"):
            y_score = clf.decision_function(X_test)
        else:
            prob = clf.predict_proba(X_test)
            classes = clf.classes_ # This is the order of the classes. The probabilities are given in this order
            for index in xrange(len(classes)):
                if classes[index] == 1:
                    dc_index = index # Obtain in which position is located the probability of being drug combination
            y_score = []
            for p in xrange(len(prob)):
                dc_prob = prob[p][dc_index] # We use the index to obtain the probability of being drug combination
                y_score.append(dc_prob) # Append the array in all_prob

        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        auc = metrics.roc_auc_score(y_test, y_pred)
        #print('SCIKIT AUC: {}\n'.format(auc))
        all_auc.append(auc)

    mean_tpr /= repetitions
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    mean_auc2 = np.mean(all_auc)
    #print('Mean AUC: {}'.format(mean_auc))
    #print('Mean AUC2: {}'.format(mean_auc2))

    var_auc = np.var(all_auc)
    std_auc = np.std(all_auc)
    #print('Var AUC: {}'.format(var_auc))
    #print('Std AUC: {}'.format(std_auc))

    return mean_auc, var_auc, std_auc, all_auc



if  __name__ == "__main__":
    main()


