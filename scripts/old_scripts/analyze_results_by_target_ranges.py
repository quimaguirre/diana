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
    parser.add_argument('-datc','--different_atc',dest='different_atc',action = 'store_true',
                        help = """ Consider only drug combinations with different ATCs. """)
    parser.add_argument('-pca','--pca',dest='pca',action = 'store_true',
                        help = """" Make a PCA to reduce dimensionality. """)
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



    #-------------------------------------#
    #   EVALUATE PERFORMANCE BY TARGETS   #
    #-------------------------------------#

    img_dir = os.path.join(analysis_dir, 'figures')
    create_directory(img_dir)
    fig_format = 'png'

    tables_dir = os.path.join(analysis_dir, 'tables')
    create_directory(tables_dir)

    # Number of targets
    num_targets = [[1],[2],[3,4,5,6],[7]]

    # Names of the methods
    if consider_se:
        if options.different_atc:
            types_analysis = ['dctargets', 'dcguild', 'dcstructure', 'dcse', 'random']
            types_analysis2 = ['dctargets', 'dcguild', 'dcstructure', 'dcse'] # Without random!!
            #types_analysis_labels = ['dcTargets', 'dcGUILD', 'dcStructure', 'dcSE', 'Random']
            types_analysis_labels = [ 'Target', 'PPI','Structure', 'Side Effects', 'Random']
        else:
            types_analysis = ['dctargets', 'dcguild', 'dcstructure', 'dcatc', 'dcse', 'random']
            types_analysis2 = ['dctargets', 'dcguild', 'dcstructure', 'dcatc', 'dcse'] # Without random!!
            #types_analysis_labels = ['dcTargets', 'dcGUILD', 'dcStructure', 'dcATC', 'dcSE', 'Random']
            types_analysis_labels = [ 'Target', 'PPI','Structure', 'ATC', 'Side Effects', 'Random']
    else:
        types_analysis = ['dctargets', 'dcguild', 'dcstructure', 'random']
        types_analysis2 = ['dctargets', 'dcguild', 'dcstructure'] # Without random!!
        types_analysis_labels = ['dcTargets', 'dcGUILD', 'dcStructure', 'Random']
        types_analysis_labels = [ 'Target', 'PPI','Structure', 'Random']


    # Machine learning parameters
    repetitions = 25 # Number of repetititons
    n_fold = 2     # Number of folds
    min_num_dc_group = 10
    greater_or_smaller = 'greater'
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

    # Plot of distributions of AUC
    plot_auc_distribution = os.path.join(img_dir, 'numtargets_auc_distribution_ranges{}.{}'.format(pca_str, fig_format))

    # Plot of accuracy/sensitivity name
    acc_sens_dctargets = os.path.join(img_dir, 'numtargets_accsens_dctargets_ranges{}.{}'.format(pca_str, fig_format))
    acc_sens_dcguild = os.path.join(img_dir, 'numtargets_accsens_dcguild_ranges{}.{}'.format(pca_str, fig_format))
    acc_sens_dcstructure = os.path.join(img_dir, 'numtargets_accsens_dcstructure_ranges{}.{}'.format(pca_str, fig_format))
    acc_sens_dcatc = os.path.join(img_dir, 'numtargets_accsens_dcatc_ranges{}.{}'.format(pca_str, fig_format))
    acc_sens_dcse = os.path.join(img_dir, 'numtargets_accsens_dcse_ranges{}.{}'.format(pca_str, fig_format))

    # Results table
    results_table = os.path.join(tables_dir, 'numtargets_auc_table_ranges{}.txt'.format(pca_str))

    # Accuracy/Sensitivity results table
    prec_rec_table = os.path.join(tables_dir, 'numtargets_accsens_table_ranges{}.txt'.format(pca_str))

    # File with results of Mann Whitney tests
    mannwhitney_file = os.path.join(tables_dir, 'numtargets_mannwhitney_ranges{}.txt'.format(pca_str))

    # Get the targets file
    drugbank_to_targets_file = os.path.join(toolbox_dir, 'drugbank_to_targets.pcl')
    drugbank_to_targets = cPickle.load(open(drugbank_to_targets_file))

    # Get the DIANA IDs file
    diana_id_to_drugbank_file = os.path.join(toolbox_dir, 'diana_id_to_drugbank.pcl')
    diana_id_to_drugbank = cPickle.load(open(diana_id_to_drugbank_file))


    analysis_results = {} # Defining the dictionary that will store the results

    if consider_se:
        dct_columns, dcg_columns, dcs_columns, dcatc_columns, dcse_columns = diana_analysis.obtain_method_to_columns(threshold_list, ATC_SE=consider_se)
    else:
        dct_columns, dcg_columns, dcs_columns = diana_analysis.obtain_method_to_columns(threshold_list, ATC_SE=consider_se)

    for range_tar in num_targets:

        selected_rows = []

        for index, row in df.iterrows():

            (drug_id1, drug_id2) = index.split('---')
            drug1 = diana_id_to_drugbank[drug_id1].upper()
            drug2 = diana_id_to_drugbank[drug_id2].upper()

            if len(range_tar) == 1:
                # If it is the first of the range
                if range_tar == num_targets[0]:
                    if len(drugbank_to_targets[drug1]) <= range_tar[0] and len(drugbank_to_targets[drug2]) <= range_tar[0]:
                        selected_rows.append(index)
                # If it is the last of the range
                elif range_tar == num_targets[len(num_targets)-1]:
                    if len(drugbank_to_targets[drug1]) >= range_tar[0] and len(drugbank_to_targets[drug2]) >= range_tar[0]:
                        selected_rows.append(index)
                # If it is in the middle of the range
                else:
                    if len(drugbank_to_targets[drug1]) == range_tar[0] and len(drugbank_to_targets[drug2]) == range_tar[0]:
                        selected_rows.append(index)
            else:
                if len(drugbank_to_targets[drug1]) in range_tar and len(drugbank_to_targets[drug2]) in range_tar:
                    selected_rows.append(index)


        df_tar = df.ix[selected_rows]
        dc_data = df_tar[df_tar['combination'] == 1]
        num_dc = len(dc_data.index)
        print('Num drug combinations: {}'.format(num_dc))

        if consider_se:
            list_methods = [ ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['dcatc', dcatc_columns], ['dcse', dcse_columns], ['random', columns] ]
        else:
            list_methods = [ ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['random', columns] ]

        for method, columns_method in list_methods:

            print('Evaluating {} targets with method {}\n'.format(range_tar,method))

            #------------------------------------------------------------------#
            #   SELECT RELEVANT FEATURES / REDUCE DIMENSIONALITY OF THE DATA   #
            #------------------------------------------------------------------#

            if options.pca:

                variance_cut_off = 0.01
                num_components = 0
                df_method = df_tar[columns_method]
                df_raw = df_method.drop('combination', axis=1)
                raw_columns = copy.copy(columns_method)
                raw_columns.remove('combination')
                pca = PCA(n_components=None)
                pca.fit(df_raw)
                values_trans = pca.transform(df_raw)
                explained_variance = pca.explained_variance_ratio_
                for column, var in sorted(zip(raw_columns, explained_variance), key=lambda x: x[1], reverse=True):
                    #print(column, var)
                    if var > variance_cut_off:
                        num_components += 1

                if num_components < len(raw_columns):

                    print('Number of features:\t{}\n'.format(len(raw_columns)))
                    print('Reduction to {} components\n'.format(num_components))

                    pca = PCA(n_components=num_components)
                    pca.fit(df_raw)
                    values_trans = pca.transform(df_raw)
                    indexes = df_method.index.values
                    df_trans = pd.DataFrame.from_records(values_trans, index=indexes)
                    df_comb = df_method[['combination']]
                    df_new = pd.concat([df_trans, df_comb], axis=1)
                    df_method = df_new

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
                    selected_columns = [col for col in selected_columns if col not in dcatc_columns or col == 'combination']

                print('Selected columns: {}\n'.format(', '.join(selected_columns)))
                print('Number of selected features: {}\n'.format(len(selected_columns)-1)) # We take away the combinations column

                # Define the new table with the selected columns
                df_method = df_tar[selected_columns]
                dc_data = df_method[df_method['combination'] == 1]
                ndc_data = df_method[df_method['combination'] == 0]
                num_dc = len(dc_data.index)
                num_ndc = len(ndc_data.index)

            #------------------------------------------------------------------#


            dc_data = df_method[df_method['combination'] == 1]
            ndc_data = df_method[df_method['combination'] == 0]
            num_dc = len(dc_data.index)
            num_ndc = len(ndc_data.index)

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

                if method == 'random':
                    #mean, var, std, list_auc = run_nfold_crossvalidation_random(n_fold, merged_groups, classifiers[classifier])
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

            # Store the distribution of AUCs in the dictionary
            analysis_results.setdefault(range_tar[0], {})
            analysis_results[range_tar[0]].setdefault(method, {})
            analysis_results[range_tar[0]][method]['all_aucs'] = all_aucs
            analysis_results[range_tar[0]][method]['all_probs'] = all_probs
            analysis_results[range_tar[0]][method]['mean'] = final_mean
            analysis_results[range_tar[0]][method]['std'] = std
            analysis_results[range_tar[0]][method]['num_dc'] = num_dc


    #------------------------------------#
    #   PLOT PRECISION VS. SENSITIVITY   #
    #------------------------------------#

    analysis_results = plot_precision_sensitivity(analysis_results, 'dctargets', num_targets, acc_sens_dctargets)
    analysis_results = plot_precision_sensitivity(analysis_results, 'dcguild', num_targets, acc_sens_dcguild)
    analysis_results = plot_precision_sensitivity(analysis_results, 'dcstructure', num_targets, acc_sens_dcstructure)
    if consider_se:
        analysis_results = plot_precision_sensitivity(analysis_results, 'dcatc', num_targets, acc_sens_dcatc)
        analysis_results = plot_precision_sensitivity(analysis_results, 'dcse', num_targets, acc_sens_dcse)


    #----------------------------------------------------#
    #   PLOT DISTRIBUTION OF AUC PER NUMBER OF TARGETS   #
    #----------------------------------------------------#

    plot_auc_distributions(analysis_results, num_targets, types_analysis, types_analysis_labels, plot_auc_distribution, fig_format=fig_format, consider_se=consider_se)


    #--------------------------------------------------------#
    #   TABLE OF DISTRIBUTION OF AUC PER NUMBER OF TARGETS   #
    #--------------------------------------------------------#

    with open(results_table, 'w') as results_table_fd:

        # Header
        results_table_fd.write(' ')
        for method in types_analysis_labels:
            results_table_fd.write('\t{}\t \t '.format(method))
        results_table_fd.write('\n')

        for num in num_targets:
            results_table_fd.write('{}'.format(num))
            for method in types_analysis:
                mean = analysis_results[num[0]][method]['mean']
                std = analysis_results[num[0]][method]['std']
                num_dc = analysis_results[num[0]][method]['num_dc']
                results_table_fd.write('\t{}\t{}\t{}'.format(mean, std, num_dc))
            results_table_fd.write('\n')


    #----------------------------------------#
    #   TABLE OF PRECISION VS. SENSITIVITY   #
    #----------------------------------------#

    with open(prec_rec_table, 'w') as prec_rec_table_fd:

        # Header
        prec_rec_table_fd.write(' ')
        for method in types_analysis2:
            prec_rec_table_fd.write('\t{}\t '.format(method))
        prec_rec_table_fd.write('\n')

        for num in num_targets:
            prec_rec_table_fd.write('{}'.format(num))
            for method in types_analysis2:
                cut_off = analysis_results[num[0]][method]['cut_off']
                value = analysis_results[num[0]][method]['value']
                prec_rec_table_fd.write('\t{}\t{}'.format(cut_off, value))
            prec_rec_table_fd.write('\n')


    #-------------------------------------------------------------------#
    #   TABLE OF COMPARISON OF AUC DISTRIBUTIONS USING MANN WHITNEY U   #
    #-------------------------------------------------------------------#

    with open(mannwhitney_file, 'w') as mannwhitney_fd:

        mann_results = {}

        mannwhitney_fd.write(' \t ')
        for method in types_analysis_labels:
            mannwhitney_fd.write('\t{}'.format(method))
        mannwhitney_fd.write('\n')

        # Perform the comparisons
        for num in num_targets:
            mann_results.setdefault(num[0], {})
            for method1 in types_analysis:
                mann_results[num[0]].setdefault(method1, {})
                for method2 in types_analysis:
                    if method1 == method2:
                        mann_results[num[0]][method1][method2] = '-'
                    else:
                        method1_dist = analysis_results[num[0]][method1]['all_aucs']
                        method2_dist = analysis_results[num[0]][method2]['all_aucs']
                        stat, pval = scipy.stats.mannwhitneyu(method1_dist, method2_dist)
                        mann_results[num[0]][method1][method2] = [stat, pval]

        # Write the table of crossings
        for num in num_targets:
            for method1 in types_analysis:
                mannwhitney_fd.write('{}\t{}'.format(num[0], method1))
                for method2 in types_analysis:
                    if method1 == method2:
                        mannwhitney_fd.write('\t-')
                    else:
                        stat, pval = mann_results[num[0]][method1][method2]
                        mannwhitney_fd.write('\t{}, {:.2e}'.format(stat,pval))
                mannwhitney_fd.write('\n')




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


def plot_precision_sensitivity(analysis_results, type_analysis, array_ordered, name_plot, greater_or_smaller='greater', fig_format='png'):
    """
    Plots the precision vs. sensitivity curve.
    """

    cut_offs = frange(0,1,0.005)
    fig, ax = plt.subplots()
    prec_colors = ['#a6d0ff','#8cc2ff','#4199fd','#0078ff','#005bc1','#004ea5','#003169']
    #prec_colors = ['#a6d0ff','#003169']
    sens_colors = ['#a6ff9e','#51ff41','#15ff00','#12c701','#0d9f00','#0c7e02','#085c00']
    #sens_colors = ['#a6ff9e','#085c00']
    point_colors = ['#ffb1b1','#ff7d7d','#ff3737','#f00000','#d40303','#b20000','#890000']
    #point_colors = ['#ffb1b1','#890000']
    c = 0

    for num in array_ordered:

        precision = []
        sensitivity = []

        all_probs = analysis_results[num[0]][type_analysis]['all_probs']

        for cut_off in cut_offs:

            tp,fp,tn,fn = calculate_tp_fp_tn_fn(all_probs, cut_off)

            try:
                prec_val = float(tp)/(float(tp)+float(fp))
            except:
                prec_val = 1
            sens_val = float(tp)/(float(tp)+float(fn))

            precision.append(prec_val)
            sensitivity.append(sens_val)
        
        ax.plot( cut_offs, precision, '-', color=prec_colors[c])
        ax.plot( cut_offs, sensitivity, '-', color=sens_colors[c])

        # Find the intersection point
        idx = np.argwhere(np.diff(np.sign(np.array(precision) - np.array(sensitivity))) != 0).reshape(-1) + 0

        # Plot the intersection point
        ax.plot(cut_offs[idx[0]], precision[idx[0]], 'o', color=point_colors[c])

        pylab.xlabel('Probability thresholds')
        pylab.ylabel('Precision-Recall')

        analysis_results[num[0]][type_analysis]['cut_off'] = cut_offs[idx[0]]
        analysis_results[num[0]][type_analysis]['value'] = precision[idx[0]]

        c+=1

    labels = []
    for nums in array_ordered:
        if len(nums) == 1:
            if nums == array_ordered[0]:
                labels.append(r'$\leq$' + ' ' +str(nums[0]))
            else:
                labels.append(r'$\geq$' + ' ' +str(nums[0]))
                labels.append(str(nums[0]) + ' ' + r'$\leq$' + ' x ' + r'$\leq$' + ' ' +str(nums[len(nums)-1]))
        else:
            labels.append(str(nums[0]) + ' ' + r'$\leq$' + ' x ' + r'$\leq$' + ' ' +str(nums[len(nums)-1]))
    labels = labels + ['Precision', 'Recall']
    # draw temporary color lines and use them to create a legend
    hB, = pylab.plot([1,1],'o', color='#ffb1b1')
    hG, = pylab.plot([1,1],'o', color='#ff7d7d')
    hY, = pylab.plot([1,1],'o', color='#ff3737')
    hR, = pylab.plot([1,1],'o', color='#f00000')
    hBl, = pylab.plot([1,1],'o', color='#d40303')
    hW, = pylab.plot([1,1],'o', color='#b20000')
    hX, = pylab.plot([1,1],'o', color='#890000')
    hPr, = pylab.plot([1,1],'-', color='#0078ff')
    hRe, = pylab.plot([1,1],'-', color='#0d9f00')
    lgd = ax.legend(handles=(hB, hG, hY, hR, hBl, hW, hX, hPr, hRe), labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hB.set_visible(False)
    hG.set_visible(False)
    hY.set_visible(False)
    hR.set_visible(False)
    hBl.set_visible(False)
    hW.set_visible(False)
    hX.set_visible(False)
    hPr.set_visible(False)
    hRe.set_visible(False)


    pylab.savefig(name_plot, format=fig_format, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    #pylab.show()
    
    return analysis_results


def calculate_tp_fp_tn_fn(all_probs, cut_off):

    tp=0
    fp=0
    tn=0
    fn=0

    for prob, label, _ in all_probs:

        if prob > cut_off:
            if label == 1:
                tp+=1
            else:
                fp+=1
        else:
            if label == 1:
                fn+=1
            else:
                tn+=1

    return tp,fp,tn,fn


def frange(x,y,jump):

    array = []
    x = 0
    y = 1
    jump = 0.005
    while x <= y:
        array.append(x)
        x += jump
    return array


def plot_auc_distributions(analysis_results, num_targets, types_analysis, types_analysis_labels, plot_name, greater_or_smaller='greater', fig_format='png', consider_se=False):

    fig = pylab.figure(dpi=300)
    ax = pylab.axes()
    #pylab.hold(True)
    pos = 2
    xticks = [] # Define the places in which the labels will be

    for num in num_targets:

        positions = []
        for x in xrange(len(types_analysis)):
            positions.append(pos) # Define the positions of the boxplots
            pos+=1
        pos+=1 # Add separation between boxplot groups
        pylab.plt.axvline(x=pos,linewidth=0.3,linestyle='--',color='black',dashes=(1, 1))
        pos+=2 # Add separation between boxplot groups

        data = []
        for method in types_analysis:
            data.append(analysis_results[num[0]][method]['all_aucs']) # Get the groups of plots that we will add

        # Boxplot group
        bp = pylab.boxplot(data, positions = positions, widths = 0.6, patch_artist=True)
        setBoxColors(bp, len(types_analysis), consider_se)

        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # Set axes limits and labels
    pylab.xlim(0,pos-2)
    pylab.ylim(0,1)
    axes_labels = []
    for nums in num_targets:
        if len(nums) == 1:
            if nums == num_targets[0]:
                if nums[0] == 1:
                    axes_labels.append(str(nums[0]))
                else:
                    axes_labels.append(r'$\leq$' + ' ' +str(nums[0]))
            elif nums == num_targets[len(num_targets)-1]:
                axes_labels.append(r'$\geq$' + ' ' +str(nums[0]))
            else:
                axes_labels.append(str(nums[0]))
        else:
            axes_labels.append(str(nums[0]) + ' ' + r'$\leq$' + ' x ' + r'$\leq$' + ' ' +str(nums[len(nums)-1]))
    ax.set_xticklabels(axes_labels)
    ax.set_xticks(xticks)
    ax.yaxis.grid(True)
    pylab.xlabel('Number of targets')
    pylab.ylabel('Distribution of AUC values')

    # draw temporary color lines and use them to create a legend
    hR, = pylab.plot([1,1],'-', color='red')
    hG, = pylab.plot([1,1],'-', color='green')
    hB, = pylab.plot([1,1],'-', color='black')
    hW, = pylab.plot([1,1],'-', color='#aeaeae') # grey
    if consider_se:
        hBl, = pylab.plot([1,1],'-', color='#22a9bd') # blue
        hO, = pylab.plot([1,1],'-', color='#e59600') # orange
        lgd = ax.legend(handles=(hR, hG, hB, hBl, hO, hW), labels=types_analysis_labels,  bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        hBl.set_visible(False)
        hO.set_visible(False)
    else:
        lgd = ax.legend(handles=(hR, hG, hB, hW), labels=types_analysis_labels,  bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hR.set_visible(False)
    hG.set_visible(False)
    hB.set_visible(False)
    hW.set_visible(False)

    pylab.savefig(plot_name, format=fig_format, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pylab.show()

    return


def setBoxColors(bp, bg_color, consider_se):
    """
    Set the colors of the box plots groups
    Code from: http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    """

    # RED (dcTargets)
    pylab.setp(bp['boxes'][0], color='#b80000')
    pylab.setp(bp['caps'][0], color='#b80000')
    pylab.setp(bp['caps'][1], color='#b80000')
    pylab.setp(bp['whiskers'][0], color='#b80000')
    pylab.setp(bp['whiskers'][1], color='#b80000')
    pylab.setp(bp['medians'][0], color='black')
    bp['boxes'][0].set_facecolor('#ff7373') #red

    # GREEN (dcGUILD)
    pylab.setp(bp['boxes'][1], color='green')
    pylab.setp(bp['caps'][2], color='green')
    pylab.setp(bp['caps'][3], color='green')
    pylab.setp(bp['whiskers'][2], color='green')
    pylab.setp(bp['whiskers'][3], color='green')
    pylab.setp(bp['medians'][1], color='black')
    bp['boxes'][1].set_facecolor('#32f232') #green

    # BLACK (dcStructure)
    pylab.setp(bp['boxes'][2], color='black')
    pylab.setp(bp['caps'][4], color='black')
    pylab.setp(bp['caps'][5], color='black')
    pylab.setp(bp['whiskers'][4], color='black')
    pylab.setp(bp['whiskers'][5], color='black')
    pylab.setp(bp['medians'][2], color='#fbc562')
    bp['boxes'][2].set_facecolor('#4f4f4f') #black

    if consider_se:
        # BLUE (dcATC)
        pylab.setp(bp['boxes'][3], color='#0049e5')
        pylab.setp(bp['caps'][6], color='#0049e5')
        pylab.setp(bp['caps'][7], color='#0049e5')
        pylab.setp(bp['whiskers'][6], color='#0049e5')
        pylab.setp(bp['whiskers'][7], color='#0049e5') #dark blue
        pylab.setp(bp['medians'][3], color='black')
        bp['boxes'][3].set_facecolor('#22a9bd') #blue

        # ORANGE (dcATC)
        pylab.setp(bp['boxes'][4], color='#966200')
        pylab.setp(bp['caps'][8], color='#966200')
        pylab.setp(bp['caps'][9], color='#966200')
        pylab.setp(bp['whiskers'][8], color='#966200')
        pylab.setp(bp['whiskers'][9], color='#966200') #brown
        pylab.setp(bp['medians'][4], color='black')
        bp['boxes'][4].set_facecolor('#e59600') #orange

        # GREY (Random)
        pylab.setp(bp['boxes'][5], color='black')
        pylab.setp(bp['caps'][10], color='black')
        pylab.setp(bp['caps'][11], color='black')
        pylab.setp(bp['whiskers'][10], color='black')
        pylab.setp(bp['whiskers'][11], color='black')
        pylab.setp(bp['medians'][5], color='#fbc562')
        bp['boxes'][5].set_facecolor('#aeaeae') #grey
    else:
        # GREY (Random)
        pylab.setp(bp['boxes'][3], color='black')
        pylab.setp(bp['caps'][6], color='black')
        pylab.setp(bp['caps'][7], color='black')
        pylab.setp(bp['whiskers'][6], color='black')
        pylab.setp(bp['whiskers'][7], color='black') #dark blue
        pylab.setp(bp['medians'][3], color='#fbc562')
        bp['boxes'][3].set_facecolor('#aeaeae') #blue

    return




if  __name__ == "__main__":
    main()


