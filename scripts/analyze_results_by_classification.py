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
    parser.add_argument('-cl','--classification',dest='classification',action = 'store',default='dcdb',
                        help = """Define the type of classification that will be used. It can be: dcdb, biological_process, pathway""")
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

    print("\n\t\t-------------------------------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Analysis of results: Analysis by classification\n")
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

    diana_id_to_drugbank_file = os.path.join(toolbox_dir, 'diana_id_to_drugbank.pcl')
    diana_id_to_drugbank = cPickle.load(open(diana_id_to_drugbank_file))

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



    #------------------------------------------------------------------#
    #   EVALUATE PERFORMANCE BY CLASSIFICATION OF DRUGS COMBINATIONS   #
    #------------------------------------------------------------------#

    img_dir = os.path.join(analysis_dir, 'figures')
    create_directory(img_dir)
    fig_format = 'png'

    tables_dir = os.path.join(analysis_dir, 'tables')
    create_directory(tables_dir)

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

    # Define the type of classification
    if options.classification == 'dcdb':
        classifications = ['Different targets in different biological processes', 
                           'Different targets in related biological processes',
                           'Different targets in same biological process',
                           'Same target']
        classifications_labels = ['Class 1', 'Class 2', 'Class 3', 'Class 4']
        type_classification = '_dcdb'

    elif options.classification == 'biological_process':
        classifications = ['different_targets_different_bp', 
                           'different_targets_similar_bp',
                           'similar_targets']
        classifications_labels = ['Class 1', 'Class 2', 'Class 3']
        type_classification = '_bp'
        classification_file = os.path.join(toolbox_dir, 'classification_targets_bp.pcl')
        classification_dict = cPickle.load(open(classification_file))

    elif options.classification == 'pathway':
        classifications = ['different_targets_different_pathways', 
                           'different_targets_similar_pathways',
                           'similar_targets']
        classifications_labels = ['Class 1', 'Class 2', 'Class 3']
        type_classification = '_pathway'
        classification_file = os.path.join(toolbox_dir, 'classification_targets_pathways.pcl')
        classification_dict = cPickle.load(open(classification_file))
    else:
        raise IncorrectClassificationType(options.classification)

    # Machine learning parameters
    repetitions = 25 # Number of repetititons
    n_fold = 2     # Number of folds
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
    plot_auc_distribution = os.path.join(img_dir, 'classification{}_auc_distribution{}{}.{}'.format(type_classification, atc_str, pca_str, fig_format))

    # Plot of accuracy/sensitivity name
    acc_sens_dctargets = os.path.join(img_dir, 'classification{}_accsens_dctargets{}{}.{}'.format(type_classification, atc_str, pca_str, fig_format))
    acc_sens_dcguild = os.path.join(img_dir, 'classification{}_accsens_dcguild{}{}.{}'.format(type_classification, atc_str, pca_str, fig_format))
    acc_sens_dcstructure = os.path.join(img_dir, 'classification{}_accsens_dcstructure{}{}.{}'.format(type_classification, atc_str, pca_str, fig_format))
    acc_sens_dcatc = os.path.join(img_dir, 'classification{}_accsens_dcatc{}{}.{}'.format(type_classification, atc_str, pca_str, fig_format))
    acc_sens_dcse = os.path.join(img_dir, 'classification{}_accsens_dcse{}{}.{}'.format(type_classification, atc_str, pca_str, fig_format))

    # Results table
    results_table = os.path.join(tables_dir, 'classification{}_auc_table{}{}.txt'.format(type_classification, atc_str, pca_str))

    # Accuracy/Sensitivity results table
    prec_rec_table = os.path.join(tables_dir, 'classification{}_accsens_table{}{}.txt'.format(type_classification, atc_str, pca_str))

    # File with results of Mann Whitney tests
    mannwhitney_file = os.path.join(tables_dir, 'classification{}_mannwhitney{}{}.txt'.format(type_classification, atc_str, pca_str))

    # Get the classification files
    drug_int_2_drugs_file = os.path.join(toolbox_dir, 'drug_int_2_drugs.pcl')
    drug_int_2_drugs = cPickle.load(open(drug_int_2_drugs_file))
    drug_int_2_info_file = os.path.join(toolbox_dir, 'drug_int_2_info.pcl')
    drug_int_2_info = cPickle.load(open(drug_int_2_info_file))
    drugbank_to_dcdb_file = os.path.join(toolbox_dir, 'drugbank_to_dcdb.pcl')
    drugbank_to_dcdb = cPickle.load(open(drugbank_to_dcdb_file))
    drugbank_to_atcs_file = os.path.join(toolbox_dir, 'drugbank_to_atcs.pcl')
    drugbank_to_atcs = cPickle.load(open(drugbank_to_atcs_file))



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



    analysis_results = {} # Defining the dictionary that will store the results

    if consider_se:
        dct_columns, dcg_columns, dcs_columns, dcatc_columns, dcse_columns = diana_analysis.obtain_method_to_columns(threshold_list, ATC_SE=consider_se)
    else:
        dct_columns, dcg_columns, dcs_columns = diana_analysis.obtain_method_to_columns(threshold_list, ATC_SE=consider_se)

    for classification in classifications:

        print('\nANALIZING CLASSIFICATION: {}\n'.format(classification))

        class_rows = []
        no_class_rows = []

        for index, row in df.iterrows():

            (drug_id1, drug_id2) = index.split('---')
            drug1 = diana_id_to_drugbank[drug_id1].upper()
            drug2 = diana_id_to_drugbank[drug_id2].upper()
            di1_bool = False
            di2_bool = False
            di_bool = False

            if options.classification == 'dcdb':

                for DI in drug_int_2_drugs:
                    # If it is drug interaction...
                    if drug1 in drugbank_to_dcdb and drug2 in drugbank_to_dcdb:
                        for db_id in drugbank_to_dcdb[drug1]:
                            if db_id in drug_int_2_drugs[DI]:
                                di1_bool = True
                        for db_id in drugbank_to_dcdb[drug2]:
                            if db_id in drug_int_2_drugs[DI]:
                                di2_bool = True
                    if di1_bool and di2_bool:
                        di_bool = True
                        # If it is from the classification of interest we store it in the class group
                        if classification == drug_int_2_info[DI]['classification']:
                            class_rows.append(index)
                            break
                        # If it is NOT from the classification of interest we store it in the no class group
                        else:
                            no_class_rows.append(index)
                            break
                if di_bool == False:
                    # If it not drug interaction, we store it in both groups
                    class_rows.append(index)
                    no_class_rows.append(index)

            elif options.classification == 'pathway' or options.classification == 'biological_process':
                if index in classification_dict:
                    if classification == classification_dict[index]:
                        class_rows.append(index)
                    else:
                        no_class_rows.append(index)
                else:
                    class_rows.append(index)
                    no_class_rows.append(index)


        df_class = df.ix[class_rows] # Create a table with the rows selected (DDIs of the class and non-DDIs)
        df_noclass = df.ix[no_class_rows] # Create a table with the DDIs not of the class and non-DDIs 

        # Get the number of DC for the given classification
        dc_data = df_class[df_class['combination'] == 1]
        num_dc = len(dc_data.index)
        print('Number of DC for {}: {}'.format(classification, num_dc))


        if consider_se:
            if options.different_atc:
                list_methods = [ ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['dcse', dcse_columns], ['random', columns] ]
            else:
                list_methods = [ ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['dcatc', dcatc_columns], ['dcse', dcse_columns], ['random', columns] ]
        else:
            list_methods = [ ['dctargets', dct_columns], ['dcguild', dcg_columns], ['dcstructure', dcs_columns], ['random', columns] ]

        for method, columns_method in list_methods:

            print('Evaluating classification {} with method {}\n'.format(classification, method))

            #------------------------------------------------------------------#
            #   SELECT RELEVANT FEATURES / REDUCE DIMENSIONALITY OF THE DATA   #
            #------------------------------------------------------------------#

            if options.pca:

                variance_cut_off = 0.01
                num_components = 0
                df_method = df_class[columns_method]
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
                df_method = df_class[selected_columns]
                dc_data = df_method[df_method['combination'] == 1]
                ndc_data = df_method[df_method['combination'] == 0]
                num_dc = len(dc_data.index)
                num_ndc = len(ndc_data.index)

            #------------------------------------------------------------------#


            dc_data = df_method[df_method['combination'] == 1]
            ndc_data = df_method[df_method['combination'] == 0]
            num_dc = len(dc_data.index)
            num_ndc = len(ndc_data.index)

            # Stop if the number of drug interactions is smaller than the number of cross-validations!!
            if num_dc < n_fold:
                print('Not possible to do the analysis for classification {}. The number of positive samples is {} and the n-fold is {}\n'.format(classification, num_dc, n_fold))
                analysis_results.setdefault(classification, {})
                analysis_results.setdefault[classification](method, {})
                analysis_results[classification][method] = {'mean':'-','std':'-','num_dc':int(num_dc),'all_aucs':'-'}
                #analysis_results = {'in classification' : {'mean':'-','std':'-','num_dc':int(num_dc)}, 'not in classification' : {'mean':'-','std':'-','num_dc':int(num_dc)}} # Defining the variable
                continue
            # Stop if the number of drug interactions is smaller than the minimum number given!!
            if num_dc < min_num_dc_group:
                print('Not possible to do the analysis for classification {}. The number of positive samples is {} and the minimum number per group is {}\n'.format(classification, num_dc, min_num_dc_group))
                analysis_results.setdefault(classification, {})
                analysis_results[classification].setdefault(method, {})
                analysis_results[classification][method] = {'mean':'-','std':'-','num_dc':int(num_dc),'all_aucs':'-'}
                #analysis_results = {'in classification' : {'mean':'-','std':'-','num_dc':int(num_dc)}, 'not in classification' : {'mean':'-','std':'-','num_dc':int(num_dc)}} # Defining the variable
                continue

            # Obtain the different non-drug combination groups to repeat the analysis
            # We will use the number of drug combinations when classification in order to have same number of samples 
            print('Building {} repetition groups of {} (same) DC and {} (different) non-DC'.format(repetitions,num_dc,num_dc))
            ndc_repetitions = diana_analysis.obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times
            #print(ndc_repetitions)

            mean_aucs = [] # Here we will store the means of AUCs from the cross-validations
            std_aucs = [] # Here we will store the standard deviations of the AUCs from the cross-validations
            all_aucs = [] # Here we will store ALL the AUCs
            all_probs = [] # Here we store all the probabilities and labels

            num_repetitions=0
            for ndc_data_equal in ndc_repetitions:

                num_items_group = int( float(num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation

                num_repetitions+=1
                if num_repetitions == 1:
                    print('Building {} fold groups of {} DC and {} non-DC'.format(n_fold,num_items_group,num_items_group))

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

            analysis_results.setdefault(classification, {})
            analysis_results[classification].setdefault(method, {})
            analysis_results[classification][method]['mean'] = final_mean
            analysis_results[classification][method]['std'] = std
            analysis_results[classification][method]['num_dc'] = int(num_dc)
            analysis_results[classification][method]['all_aucs'] = all_aucs
            analysis_results[classification][method]['all_probs'] = all_probs



    #------------------------------------#
    #   PLOT PRECISION VS. SENSITIVITY   #
    #------------------------------------#

    analysis_results = plot_precision_sensitivity(analysis_results, 'dctargets', classifications, classifications_labels, acc_sens_dctargets)
    analysis_results = plot_precision_sensitivity(analysis_results, 'dcguild', classifications, classifications_labels, acc_sens_dcguild)
    analysis_results = plot_precision_sensitivity(analysis_results, 'dcstructure', classifications, classifications_labels, acc_sens_dcstructure)
    if consider_se:
        if options.different_atc:
            analysis_results = plot_precision_sensitivity(analysis_results, 'dcse', classifications, classifications_labels, acc_sens_dcse)
        else:
            analysis_results = plot_precision_sensitivity(analysis_results, 'dcatc', classifications, classifications_labels, acc_sens_dcatc)
            analysis_results = plot_precision_sensitivity(analysis_results, 'dcse', classifications, classifications_labels, acc_sens_dcse)



    #-------------------------------------------------#
    #   PLOT DISTRIBUTION OF AUC PER CLASSIFICATION   #
    #-------------------------------------------------#

    plot_auc_distributions(analysis_results, classifications, classifications_labels, types_analysis, types_analysis_labels, plot_auc_distribution, fig_format=fig_format, consider_se=consider_se, different_atc=options.different_atc)



    #--------------------------------------------------------#
    #   TABLE OF DISTRIBUTION OF AUC PER NUMBER OF TARGETS   #
    #--------------------------------------------------------#

    with open(results_table, 'w') as results_table_fd:

        # Header
        results_table_fd.write(' ')
        for method in types_analysis_labels:
            results_table_fd.write('\t{}\t \t '.format(method))
        results_table_fd.write('\n')

        for classification in classifications:
            results_table_fd.write('{}'.format(classification))
            for method in types_analysis:
                mean = analysis_results[classification][method]['mean']
                std = analysis_results[classification][method]['std']
                num_dc = analysis_results[classification][method]['num_dc']
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

        for classification in classifications:
            prec_rec_table_fd.write('{}'.format(classification))
            for method in types_analysis2:
                cut_off = analysis_results[classification][method]['cut_off']
                value = analysis_results[classification][method]['value']
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
        for classification in classifications:
            mann_results.setdefault(classification, {})
            for method1 in types_analysis:
                mann_results[classification].setdefault(method1, {})
                for method2 in types_analysis:
                    if method1 == method2:
                        mann_results[classification][method1][method2] = '-'
                    else:
                        method1_dist = analysis_results[classification][method1]['all_aucs']
                        method2_dist = analysis_results[classification][method2]['all_aucs']
                        stat, pval = scipy.stats.mannwhitneyu(method1_dist, method2_dist)
                        mann_results[classification][method1][method2] = [stat, pval]

        # Write the table of crossings
        for classification in classifications:
            for method1 in types_analysis:
                mannwhitney_fd.write('{}\t{}'.format(classification, method1))
                for method2 in types_analysis:
                    if method1 == method2:
                        mannwhitney_fd.write('\t-')
                    else:
                        stat, pval = mann_results[classification][method1][method2]
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


class IncorrectClassificationType(Exception):
    """
    Exception raised when the classification type is not correct.
    """
    def __init__(self, classification):
        self.classification = classification

    def __str__(self):
        return 'The classification {} is not correct.\nPlease, introduce "dcdb", "biological_process", or "pathway".\n'.format(self.classification)


def plot_precision_sensitivity(analysis_results, type_analysis, classifications, classifications_labels, name_plot, fig_format='png'):
    """
    Plots the precision vs. sensitivity curve.
    """

    cut_offs = frange(0,1,0.005)
    fig, ax = plt.subplots()
    prec_colors = ['#a6d0ff','#4199fd','#005bc1','#003169']
    #prec_colors = ['#a6d0ff','#003169']
    sens_colors = ['#a6ff9e','#15ff00','#0d9f00','#085c00']
    #sens_colors = ['#a6ff9e','#085c00']
    point_colors = ['#ffb1b1','#ff3737','#d40303','#890000']
    #point_colors = ['#ffb1b1','#890000']
    c = 0

    for classification in classifications:

        precision = []
        sensitivity = []

        print(type_analysis)
        print(classification)
        all_probs = analysis_results[classification][type_analysis]['all_probs']

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

        analysis_results[classification][type_analysis]['cut_off'] = cut_offs[idx[0]]
        analysis_results[classification][type_analysis]['value'] = precision[idx[0]]

        c+=1

    classifications_labels = classifications_labels + ['Precision', 'Recall']

    # draw temporary color lines and use them to create a legend
    hB, = pylab.plot([1,1],'o', color='#ffb1b1')
    hG, = pylab.plot([1,1],'o', color='#ff3737')
    hY, = pylab.plot([1,1],'o', color='#d40303')
    hR, = pylab.plot([1,1],'o', color='#890000')
    hPr, = pylab.plot([1,1],'-', color='#0078ff')
    hRe, = pylab.plot([1,1],'-', color='#0d9f00')
    lgd = ax.legend(handles=(hB, hG, hY, hR, hPr, hRe), labels=classifications_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hB.set_visible(False)
    hG.set_visible(False)
    hY.set_visible(False)
    hR.set_visible(False)
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


def plot_auc_distributions(analysis_results, classifications, classifications_labels, types_analysis, types_analysis_labels, plot_name, fig_format='png', consider_se=False, different_atc=False):

    fig = pylab.figure(dpi=300)
    ax = pylab.axes()
    #pylab.hold(True)
    pos = 2
    xticks = [] # Define the places in which the labels will be

    for classification in classifications:

        positions = []
        for x in xrange(len(types_analysis)):
            positions.append(pos) # Define the positions of the boxplots
            pos+=1
        pos+=1 # Add separation between boxplot groups
        pylab.plt.axvline(x=pos,linewidth=0.3,linestyle='--',color='black',dashes=(1, 1))
        pos+=2 # Add separation between boxplot groups

        data = []
        for method in types_analysis:
            data.append(analysis_results[classification][method]['all_aucs']) # Get the groups of plots that we will add

        # Boxplot group
        #bp = boxplot(data, positions = positions, widths = 0.6)
        bp = pylab.boxplot(data, positions = positions, widths = 0.6, patch_artist=True)
        setBoxColors(bp, len(types_analysis), consider_se, different_atc)

        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # Set axes limits and labels
    pylab.xlim(0,pos-2)
    pylab.ylim(0,1)
    axes_labels = classifications_labels
    ax.set_xticklabels(axes_labels)
    ax.set_xticks(xticks)
    ax.yaxis.grid(True)
    pylab.ylabel('Distribution of AUC values')
    ax.yaxis.grid(True)

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


def setBoxColors(bp, bg_color, consider_se, different_atc):
    """
    Set the colors of the box plots groups
    Code from: http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    """

    pylab.setp(bp['boxes'][0], color='#b80000')
    pylab.setp(bp['caps'][0], color='#b80000')
    pylab.setp(bp['caps'][1], color='#b80000')
    pylab.setp(bp['whiskers'][0], color='#b80000')
    pylab.setp(bp['whiskers'][1], color='#b80000')
    pylab.setp(bp['medians'][0], color='black')
    bp['boxes'][0].set_facecolor('#ff7373') #red

    pylab.setp(bp['boxes'][1], color='green')
    pylab.setp(bp['caps'][2], color='green')
    pylab.setp(bp['caps'][3], color='green')
    pylab.setp(bp['whiskers'][2], color='green')
    pylab.setp(bp['whiskers'][3], color='green')
    pylab.setp(bp['medians'][1], color='black')
    bp['boxes'][1].set_facecolor('#32f232') #green

    pylab.setp(bp['boxes'][2], color='black')
    pylab.setp(bp['caps'][4], color='black')
    pylab.setp(bp['caps'][5], color='black')
    pylab.setp(bp['whiskers'][4], color='black')
    pylab.setp(bp['whiskers'][5], color='black')
    pylab.setp(bp['medians'][2], color='#fbc562')
    bp['boxes'][2].set_facecolor('#4f4f4f') #black

    if consider_se:
        if not different_atc:
            # BLUE (dcATC)
            pylab.setp(bp['boxes'][3], color='#0049e5')
            pylab.setp(bp['caps'][6], color='#0049e5')
            pylab.setp(bp['caps'][7], color='#0049e5')
            pylab.setp(bp['whiskers'][6], color='#0049e5')
            pylab.setp(bp['whiskers'][7], color='#0049e5') #dark blue
            pylab.setp(bp['medians'][3], color='black')
            bp['boxes'][3].set_facecolor('#22a9bd') #blue

            # ORANGE (dcSE)
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
            # ORANGE (dcSE)
            pylab.setp(bp['boxes'][3], color='#966200')
            pylab.setp(bp['caps'][6], color='#966200')
            pylab.setp(bp['caps'][7], color='#966200')
            pylab.setp(bp['whiskers'][6], color='#966200')
            pylab.setp(bp['whiskers'][7], color='#966200') #brown
            pylab.setp(bp['medians'][3], color='black')
            bp['boxes'][3].set_facecolor('#e59600') #orange

            # GREY (Random)
            pylab.setp(bp['boxes'][4], color='black')
            pylab.setp(bp['caps'][8], color='black')
            pylab.setp(bp['caps'][9], color='black')
            pylab.setp(bp['whiskers'][8], color='black')
            pylab.setp(bp['whiskers'][9], color='black')
            pylab.setp(bp['medians'][4], color='#fbc562')
            bp['boxes'][4].set_facecolor('#aeaeae') #grey
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


