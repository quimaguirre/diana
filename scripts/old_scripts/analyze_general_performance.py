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
                        help = """ Consider Side Effects / ATCs. """)
    parser.add_argument('-datc','--different_atc',dest='different_atc',action = 'store_true',
                        help = """ Consider only drug combinations with different ATCs. """)
    parser.add_argument('-woutrep','--without_repetition',dest='without_repetition',action = 'store_true',
                        help = """ Make the cross-validation without repetitions. """)
    parser.add_argument('-pca','--pca',dest='pca',action = 'store_true',
                        help = """" Make a PCA to reduce dimensionality. """)
    parser.add_argument('-cp','--comparison',dest='comparison_other_methods',action = 'store_true',
                        help = """ If we are considering a dataset to compare with other methods. """)
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

    # Make a cross-validation with the validation set (True)
    # or make a training with the training and a validation with the validation (False)
    cross_validation = False

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


    #----------------------------#
    #   GET THE VALIDATION SET   #
    #----------------------------#

    training_dataframe = os.path.join(analysis_dir, 'dcdb_comparisons_training.csv')
    validation_dataframe = os.path.join(analysis_dir, 'dcdb_comparisons_validation.csv')
    df_training = pd.read_csv(training_dataframe, index_col=0)
    df_validation = pd.read_csv(validation_dataframe, index_col=0)

    if cross_validation:
        df = df_validation
        dc_data = df[df['combination'] == 1]
        ndc_data = df[df['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)
        print('Number of drug combinations after getting the validation dataset:\t{}\n'.format(num_dc))
        print('Number of non-drug combinations after getting the validation dataset:\t{}\n'.format(num_ndc))
    else:
        # Define the variables for the training dataset
        df = df_training
        dc_data = df[df['combination'] == 1]
        ndc_data = df[df['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)
        print('Number of drug combinations after getting the training dataset:\t{}\n'.format(num_dc))
        print('Number of non-drug combinations after getting the training dataset:\t{}\n'.format(num_ndc))

        # Define the variables for the validation dataset
        dc_data_val = df_validation[df_validation['combination'] == 1]
        ndc_data_val = df_validation[df_validation['combination'] == 0]
        num_dc_val = len(dc_data_val.index)
        num_ndc_val = len(ndc_data_val.index)
        print('Number of drug combinations after getting the validation dataset:\t{}\n'.format(num_dc_val))
        print('Number of non-drug combinations after getting the validation dataset:\t{}\n'.format(num_ndc_val))




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
    n_fold = 10     # Number of folds
    min_num_dc_group = 10
    classifier = 'SVC best 1'
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
        'SVC best 1' : SVC(kernel="rbf", gamma=0.01, C=100, probability=True),
        'SVC best 2' : SVC(kernel="rbf", gamma=0.1, C=1.0, probability=True)
    }

    if options.pca:
        pca_str = '_withPCA'
    else:
        pca_str = '_withoutPCA'
    if options.different_atc:
        atc_str = '_diff_ATC'
    else:
        atc_str = ''

    # Plot of distributions of AUC
    plot_name = os.path.join(img_dir, 'general_performance_by_methods{}{}.{}'.format(atc_str, pca_str, fig_format))

    # Get the targets file
    drugbank_to_targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    drugbank_to_targets = cPickle.load(open(drugbank_to_targets_file))

    # Get the ATC file
    drugbank_to_atcs_file = os.path.join(toolbox_dir, 'drugbank_to_atcs.pcl')
    drugbank_to_atcs = cPickle.load(open(drugbank_to_atcs_file))

    # Get the DIANA IDs file
    diana_id_to_drugbank_file = os.path.join(toolbox_dir, 'diana_id_to_drugbank.pcl')
    diana_id_to_drugbank = cPickle.load(open(diana_id_to_drugbank_file))

    print('\nEVALUATION OF GENERAL PERFORMANCE\n')
    repetitions = 25
    n_fold = 10
    analysis_results = {}
    method_to_results = {}
    method_to_probs = {}

    # Get columns for each method
    if consider_se:
        dct_columns, dcg_columns, dcs_columns, dcatc_columns, dcse_columns = diana_analysis.obtain_method_to_columns(threshold_list, ATC_SE=consider_se)
    else:
        dct_columns, dcg_columns, dcs_columns = diana_analysis.obtain_method_to_columns(threshold_list, ATC_SE=consider_se)

    # Remove ATC columns if different ATC
    if options.different_atc:
        columns = [col for col in columns if col not in dcatc_columns or col == 'combination']

    if consider_se:
        if options.different_atc:
            list_methods = [ ['Combination', columns], ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['dcatc', dcatc_columns], ['dcse', dcse_columns], ['random', columns] ]
            methods_ordered = ['Combination', 'dctargets', 'dcguild', 'dcstructure', 'dcatc', 'dcse', 'random']
            method_to_label = {'Combination':'All', 'dctargets':'Target', 'dcguild':'PPI', 'dcstructure':'Structure', 'dcatc':'ATC', 'dcse':'Side Effects', 'random':'Random'}
            colors_ordered = [['yellow', 'black'], ['#ff7373', 'red'], ['#32f232', 'green'], ['#4f4f4f', 'black'], ['#e59600', '#966200'], ['#aeaeae', 'black']] # yellow, red, green, black, orange, grey            
        else:
            list_methods = [ ['Combination', columns], ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['dcatc', dcatc_columns], ['dcse', dcse_columns], ['random', columns] ]
            methods_ordered = ['Combination', 'dctargets', 'dcguild', 'dcstructure', 'dcatc', 'dcse', 'random']
            method_to_label = {'Combination':'All', 'dctargets':'Target', 'dcguild':'PPI', 'dcstructure':'Structure', 'dcatc':'ATC', 'dcse':'Side Effects', 'random':'Random'}
            colors_ordered = [['yellow', 'black'], ['#ff7373', 'red'], ['#32f232', 'green'], ['#4f4f4f', 'black'], ['#22a9bd', '#0049e5'], ['#e59600', '#966200'], ['#aeaeae', 'black']] # yellow, red, green, black, blue, orange, grey
    else:
        list_methods = [ ['Combination', columns], ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['random', columns] ]
        methods_ordered = ['Combination', 'dctargets', 'dcguild', 'dcstructure', 'random']
        method_to_label = {'Combination':'All', 'dctargets':'Target', 'dcguild':'PPI', 'dcstructure':'Structure', 'random':'Random'}
        colors_ordered = [['white', 'black'], ['#ff7373', 'red'], ['#32f232', 'green'], ['#4f4f4f', 'black'], ['#aeaeae', 'black']] # white, red, green, black, grey


    #-------------------------------------------------#
    #   SELECT DRUG COMBINATIONS WITH DIFFERENT ATC   #
    #-------------------------------------------------#

    if options.different_atc:

        selected_rows = []
        for index, row in df.iterrows():
            (drug_id1, drug_id2) = index.split('---')
            drug1 = diana_id_to_drugbank[drug_id1].upper()
            drug2 = diana_id_to_drugbank[drug_id2].upper()

            atcs_drug1 = set([ atc[0] for atc in drugbank_to_atcs[drug1] ])
            atcs_drug2 = set([ atc[0] for atc in drugbank_to_atcs[drug2] ])
            intersection = atcs_drug1 & atcs_drug2
            if len(intersection) == 0:
                selected_rows.append(index)

        df = df.ix[selected_rows]
        dc_data = df[df['combination'] == 1]
        ndc_data = df[df['combination'] == 0]
        num_dc = len(dc_data.index)
        num_ndc = len(ndc_data.index)
        print('Num drug combinations after removing the ones with same ATC in training: {}'.format(num_dc))
        print('Num non-drug combinations after removing the ones with same ATC in training: {}'.format(num_ndc))

        selected_rows = []
        for index, row in df_validation.iterrows():
            (drug_id1, drug_id2) = index.split('---')
            drug1 = diana_id_to_drugbank[drug_id1].upper()
            drug2 = diana_id_to_drugbank[drug_id2].upper()

            atcs_drug1 = set([ atc[0] for atc in drugbank_to_atcs[drug1] ])
            atcs_drug2 = set([ atc[0] for atc in drugbank_to_atcs[drug2] ])
            intersection = atcs_drug1 & atcs_drug2
            if len(intersection) == 0:
                selected_rows.append(index)

        df_validation = df_validation.ix[selected_rows]
        dc_data_val = df_validation[df_validation['combination'] == 1]
        ndc_data_val = df_validation[df_validation['combination'] == 0]
        num_dc_val = len(dc_data_val.index)
        num_ndc_val = len(ndc_data_val.index)
        print('Num drug combinations (in validation) after removing the ones with same ATC in training: {}'.format(num_dc_val))
        print('Num non-drug combinations (in validation) after removing the ones with same ATC in training: {}'.format(num_ndc_val))



    #--------------------------#
    #   EVALUATE EACH METHOD   #
    #--------------------------#

    for method, columns_method in list_methods:

        print('Evaluating method {}\n'.format(method))


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
            scoring_methods = ['spearman', 'dot_product', 'jaccard']
            df_method = df[columns_method]
            df_val = df_validation[columns_method]
            df_all = pd.concat([df_method,df_val])
            df_raw = df_all.drop('combination', axis=1)
            raw_columns = copy.copy(columns_method)
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
                indexes = df_all.index.values
                df_trans = pd.DataFrame.from_records(values_trans, index=indexes)
                df_comb = df_all[['combination']]
                df_pca = pd.concat([df_trans, df_comb], axis=1)
                train_indexes = df_method.index.values
                val_indexes = df_val.index.values

                df_method = df_pca.loc[train_indexes]
                dc_data = df_method[df_method['combination'] == 1]
                ndc_data = df_method[df_method['combination'] == 0]
                num_dc = len(dc_data.index)
                num_ndc = len(ndc_data.index)

                df_val = df_pca.loc[val_indexes]
                dc_data_val = df_val[df_val['combination'] == 1]
                ndc_data_val = df_val[df_val['combination'] == 0]
                num_dc_val = len(dc_data_val.index)
                num_ndc_val = len(ndc_data_val.index)

        else:

            # Manually introduced features
            guild_thresholds = [1, 5]
            rank_scoring = ['spearman', 'dot_product']
            list_scoring = ['jaccard']
            if method == 'Combination' or method == 'random':
                selected_columns = diana_analysis.obtain_columns_best_features(guild_thresholds, rank_scoring, list_scoring, ATC_SE=consider_se)
            else:
                selected_columns = diana_analysis.obtain_columns_best_features_for_specific_method(method, guild_thresholds, rank_scoring, list_scoring)

            # Remove ATC columns if different ATC
            if options.different_atc and consider_se:
                if method != 'dcatc':
                    selected_columns = [col for col in selected_columns if col not in dcatc_columns or col == 'combination']

            print('Selected columns: {}\n'.format(', '.join(selected_columns)))
            print('Number of selected features: {}\n'.format(len(selected_columns)-1)) # We take away the combinations column

            # Define the new table with the selected columns
            df_method = df[selected_columns]
            dc_data = df_method[df_method['combination'] == 1]
            ndc_data = df_method[df_method['combination'] == 0]
            num_dc = len(dc_data.index)
            num_ndc = len(ndc_data.index)

            # Define also the validation table with new columns
            df_val = df_validation[selected_columns]
            dc_data_val = df_val[df_val['combination'] == 1]
            ndc_data_val = df_val[df_val['combination'] == 0]
            num_dc_val = len(dc_data_val.index)
            num_ndc_val = len(ndc_data_val.index)


        #-------------------------#
        #   CLASSIFY DRUG PAIRS   #
        #-------------------------#

        if options.without_repetition:

            # from sklearn.model_selection import train_test_split
            # data, target = df_method.iloc[:, :-1], df_method.iloc[:, -1]
            # X_train, X_test, y_train, y_test = train_test_split(
            #     data, target, test_size=0.1, random_state=0)
            # clf = classifiers[classifier].fit(X_train, y_train)
            # y_pred = clf.predict(X_test)
            # fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
            # auc = metrics.roc_auc_score(y_test, y_pred)
            # analysis_results[method] = auc
            # print('Method: {}. AUC {}.'.format(method, auc))
            # print(fpr)
            # print(tpr)

            # Calculate the number of items per group
            num_items_group_dc = int( float(num_dc) / float(n_fold) )
            num_items_group_ndc = int( float(num_ndc) / float(n_fold) )
            print('Building {} groups of {} drug combinations'.format(n_fold,num_items_group_dc))
            dc_groups = diana_analysis.obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group_dc, me_too_drug_combinations)
            print('Building {} groups of {} non-drug combinations'.format(n_fold,num_items_group_ndc))
            ndc_groups = diana_analysis.obtain_n_groups_of_k_length(ndc_data, n_fold, num_items_group_ndc, me_too_drug_combinations)
            merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

            if method == 'random':
                mean, var, std, list_auc, list_prob = diana_analysis.run_nfold_crossvalidation_dummy(n_fold, merged_groups, classifiers[classifier])
            else:
                mean, var, std, list_auc, list_prob = diana_analysis.run_nfold_crossvalidation_scikit_with_prob(n_fold, merged_groups, classifiers[classifier])

            analysis_results[method] = mean
            method_to_results[method] = (mean, std)
            method_to_probs[method] = list_prob
            print('Method: {}. AUC mean {}. AUC results: {}'.format(method, mean, list_auc))

        else:

            print('Building {} repetition groups of {} (same) DC and {} (different) non-DC'.format(repetitions,num_dc,num_dc))
            ndc_repetitions = diana_analysis.obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times

            mean_aucs = [] # Here we will store the means of AUCs from the cross-validations
            std_aucs = [] # Here we will store the standard deviations of the AUCs from the cross-validations
            all_aucs = [] # Here we will store ALL the AUCs
            all_probs = [] # Here we store all the probabilities and labels

            if cross_validation:
                num_repetitions=0
                for ndc_data_equal in ndc_repetitions:

                    num_repetitions+=1
                    num_items_group = int( float(num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation
                    if num_repetitions == 1:
                        print('Building {} fold groups of {} DC and {} non-DC x {} repetitions'.format(n_fold,num_items_group,num_items_group, repetitions))

                    dc_groups = diana_analysis.obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group, me_too_drug_combinations) # Defining the drug combination groups in each cross-validation step
                    ndc_groups = diana_analysis.obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group, me_too_drug_combinations) # Defining the non-drug combination groups in each cross-validation step
                    merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

                    if method == 'random':
                        mean, var, std, list_auc, list_prob = diana_analysis.run_nfold_crossvalidation_dummy(n_fold, merged_groups, classifiers[classifier])
                    else:
                        mean, var, std, list_auc, list_prob = diana_analysis.run_nfold_crossvalidation_scikit_with_prob(n_fold, merged_groups, classifiers[classifier])

                    mean_aucs.append(mean)
                    std_aucs.append(std)
                    all_aucs = all_aucs + list_auc
                    all_probs = all_probs + list_prob

                final_mean = np.mean(all_aucs)
                #final_mean = np.mean(mean_aucs)
                std = np.std(all_aucs)
                mean_std = np.mean(std_aucs)
                std_means = np.std(mean_aucs)
                print('FINAL MEAN: {}'.format(final_mean))
                print('STD: {}\n'.format(std))
                #print('MEAN of STD: {}'.format(mean_std))

            else:
                ndc_repetitions_val = diana_analysis.obtain_n_groups_of_k_length(ndc_data_val, repetitions, num_dc_val)
                num_repetitions=0
                for ndc_data_equal in ndc_repetitions:
                    ndc_data_equal_val = ndc_repetitions_val[num_repetitions]
                    num_repetitions+=1
                    train = pd.concat([dc_data,ndc_data_equal])
                    test = pd.concat([dc_data_val,ndc_data_equal_val])
                    X_train, y_train = train.iloc[:, :-1], train.iloc[:, -1]
                    X_test, y_test = test.iloc[:, :-1], test.iloc[:, -1]

                    if method == 'random':
                        clf = DummyClassifier().fit(X_train, y_train)
                        y_pred = clf.predict(X_test)
                    else:
                        clf = classifiers[classifier].fit(X_train, y_train)
                        y_pred = clf.predict(X_test)
                    auc = metrics.roc_auc_score(y_test, y_pred)
                    all_aucs.append(auc)

                    # Get probabilities of being drug combination
                    prob = clf.predict_proba(X_test) # Get the probability used to classify. This is a list, and there is a probability for each class
                    classes = clf.classes_ # This is the order of the classes. The probabilities are given in this order
                    for index in xrange(len(classes)):
                        if classes[index] == 1:
                            dc_index = index # Obtain in which position is located the probability of being drug combination
                    for p in xrange(len(prob)):
                        dc_prob = prob[p][dc_index] # We use the index to obtain the probability of being drug combination
                        dc_label = y_test[p]
                        dc_name = y_test.index.values[p] # We obtain the name of the drug combination
                        array = [ dc_prob, dc_label, dc_name ] # Create an array with the probability, the label and the name of the pair
                        all_probs.append(array) # Append the array in all_prob

                final_mean = np.mean(all_aucs)
                std = np.std(all_aucs)
                print('FINAL MEAN: {}'.format(final_mean))
                print('STD: {}\n'.format(std))

            # Store the distribution of AUCs in the dictionary
            analysis_results[method] = all_aucs
            method_to_results[method] = (final_mean, std)
            method_to_probs[method] = all_probs


    if options.without_repetition:

        pass

    else:

        #------------------------------#
        #   PLOT DISTRIBUTION OF AUC   #
        #------------------------------#

        methods_without_atc = copy.copy(methods_ordered)
        methods_without_atc.remove('dcatc')
        all_data = [ analysis_results[method] for method in methods_without_atc ]
        data_labels = [ method_to_label[method] for method in methods_without_atc ]

        fig = pylab.figure(dpi=300)
        ax = pylab.axes()
        pos = 1
        all_positions = []

        for x in xrange(len(methods_without_atc)):

            # plot violin plot
            print(data_labels[x])
            print(all_data[x])
            parts = ax.violinplot(all_data[x],
                               positions = [pos],
                               showmeans=False,
                               showmedians=True)
            
            all_positions.append(pos)
            pos+=2

            # Change color of the body
            for pc in parts['bodies']:
                pc.set_facecolor(colors_ordered[x][0])

            # Change color of the segments
            parts['cmedians'].set_color(colors_ordered[x][1])
            parts['cbars'].set_color(colors_ordered[x][1])
            parts['cmins'].set_color(colors_ordered[x][1])
            parts['cmaxes'].set_color(colors_ordered[x][1])

        # adding horizontal grid lines
        ax.yaxis.grid(True)
        ax.set_xticks([y + 1 for y in range(len(all_data))])
        ax.set_ylabel('Distribution of AUC values')
        # add x-tick labels
        plt.setp(ax, xticks=all_positions,
                 xticklabels=data_labels)
        #plt.xticks(rotation=15)
        # Save
        pylab.savefig(plot_name, format=fig_format)
        plt.show()


    #---------------------------------#
    #   PRINT THE RESULTS IN A FILE   #
    #---------------------------------#

    tables_dir = os.path.join(analysis_dir, 'tables')
    create_directory(tables_dir)
    output_file = os.path.join(tables_dir, 'general_performance{}{}.txt'.format(atc_str, pca_str))
    with open(output_file, 'w') as output_fd:
        for method, results in sorted(method_to_results.iteritems(), key=lambda (x, y): y[0], reverse = True):
            output_fd.write('{}\t{}\t{}\n'.format(method, results[0], results[1]))


    #-------------------------------------------------------------------#
    #   TABLE OF COMPARISON OF AUC DISTRIBUTIONS USING MANN WHITNEY U   #
    #-------------------------------------------------------------------#

    mannwhitney_file = os.path.join(tables_dir, 'general_performance_mannwhitney{}{}.txt'.format(atc_str, pca_str))
    with open(mannwhitney_file, 'w') as mannwhitney_fd:

        mann_results = {}

        mannwhitney_fd.write(' ')
        for method in methods_ordered:
            mannwhitney_fd.write('\t{}'.format(method_to_label[method]))
        mannwhitney_fd.write('\n')

        # Perform the comparisons
        for method1 in methods_ordered:
            mann_results.setdefault(method1, {})
            for method2 in methods_ordered:
                if method1 == method2:
                    mann_results[method1][method2] = '-'
                else:
                    method1_dist = analysis_results[method1]
                    method2_dist = analysis_results[method2]
                    stat, pval = scipy.stats.mannwhitneyu(method1_dist, method2_dist)
                    mann_results[method1][method2] = [stat, pval]

        # Write the table of crossings
        for method1 in methods_ordered:
            mannwhitney_fd.write('{}'.format(method_to_label[method1]))
            for method2 in methods_ordered:
                if method1 == method2:
                    mannwhitney_fd.write('\t-')
                else:
                    stat, pval = mann_results[method1][method2]
                    mannwhitney_fd.write('\t{}, {:.2e}'.format(stat,pval))
            mannwhitney_fd.write('\n')


    #-------------------------------------------------------------------------#
    #   PRINT THE MEAN OF PROBABILITIES OF BEING DRUG COMBINATION IN A FILE   #
    #-------------------------------------------------------------------------#

    prob_file = os.path.join(tables_dir, 'general_performance_probabilities{}{}.txt'.format(atc_str, pca_str))
    with open(prob_file, 'w') as prob_fd:
        for method in methods_ordered:
            dc2scoresmean = obtain_drug_combination_scores_mean(method_to_probs[method])
            for dc, mean in sorted(dc2scoresmean.iteritems(), key=lambda (x, y): y, reverse=True):
                drug_id1, drug_id2 = dc.split('---')
                drug1 = diana_id_to_drugbank[drug_id1].upper()
                drug2 = diana_id_to_drugbank[drug_id2].upper()
                atcs_drug1 = ', '.join(sorted(set([ atc[0] for atc in drugbank_to_atcs[drug1] ])))
                atcs_drug2 = ', '.join(sorted(set([ atc[0] for atc in drugbank_to_atcs[drug2] ])))
                prob_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(method, drug1, drug2, atcs_drug1, atcs_drug2, mean))


    #-------------------#
    #   SAVE ALL AUCs   #
    #-------------------#

    auc_file = os.path.join(tables_dir, 'general_performance_aucs{}{}.txt'.format(atc_str, pca_str))
    with open(auc_file, 'w') as auc_fd:
        for method in methods_ordered:
            auc_fd.write('{}\t{}\n'.format(method, ','.join([str(x) for x in analysis_results[method]])))



    # fig = pylab.figure(dpi=300)
    # ax = pylab.axes()
    # #pylab.hold(True)
    # pos = 1
    # col_num = 0

    # xticks = [] # Define the places in which the labels will be
    # xlabels = [] # Define the labels (the names of the methods)
    # #colors = [ ['#9ed0ff, blue'], ['#32f232', 'green'], ['#fbc562', '#d48900'], ['#ff7373', '#b80000'], ['grey', 'black'] ]

    # for method in methods_ordered:

    #     positions = []
    #     positions.append(pos) # Define the positions of the boxplots
    #     pos+=2 # Add separation between boxplots
    #     xlabels.append(method_to_label[method]) # Add the method used at the x axis

    #     # Boxplot group
    #     #bp = boxplot(data, positions = positions, widths = 0.6)
    #     bp = pylab.boxplot(analysis_results[method], positions = positions, widths = 0.6, patch_artist=True)

    #     tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
    #     xticks.append(tick)

    # # Set axes limits and labels
    # pylab.xlim(0,pos-1)
    # pylab.ylim(0,1)
    # ax.set_xticklabels(xlabels)
    # ax.set_xticks(xticks)
    # pylab.xlabel('Type of data')
    # pylab.ylabel('Distribution of AUC values')

    # fig.autofmt_xdate()
    # pylab.savefig(plot_name, format=fig_format)
    # pylab.show()


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


def obtain_drug_combination_scores_mean(all_probs):
    """
    Obtain a dictionary of drug_combination_name => mean of probabilities of being drug combination
    """

    dc2scores = {}
    dc2scoresmean = {}

    for prob_array in all_probs:

        dc_prob, dc_label, dc_name = prob_array

        # Skip the non drug combinations
        if dc_label != 1:
            continue

        dc2scores.setdefault(dc_name, [])
        dc2scores[dc_name].append(dc_prob)

    for dc in dc2scores:
        mean = np.mean(dc2scores[dc])
        dc2scoresmean[dc]=mean

    return dc2scoresmean


if  __name__ == "__main__":
    main()
