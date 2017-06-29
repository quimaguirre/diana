import cPickle
import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir
from os.path import isfile, isdir, join
import pandas as pd
import re
import sys
import scipy.stats as stats
from scipy import interp

# Scikit Learn modules
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


def main():

    parse_results()

    return


def parse_results():

    ################
    #### INPUTS ####
    ################

    # Folder containing everything about the analysis
    current_dir = '/home/quim/project/diana_results'
    results_dir = current_dir + '/3_targets_analysis'

    # Table with the parsed results of DIANA
    data_frame_file = results_dir + '/all_results_3targets.csv'

    # ATC groups to analyze
    atc_groups = set()
    dump_file = current_dir + "/dcdb2atc_lvl1.pcl"
    dcdb2atc_lvl1 = cPickle.load(open(dump_file))
    for dcdb in dcdb2atc_lvl1:
        atc_groups = atc_groups.union(dcdb2atc_lvl1[dcdb]) # Join all the ATCs in one group

    # Remove always all nan
    remove_all_nan = False

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    type_analysis = 'main'
    print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))

    # Defining the parameters of the cross-validation test
    repetitions = 10 # Number of repetititons
    n_fold = 10     # Number of folds
    classifier = 'SVC'
    classifiers = {
        'KNeighborsClassifier' : KNeighborsClassifier(3),
        'SVC' : SVC(),
        'SVC_0_025' : SVC(kernel="linear", C=0.025),
        'SVC_1' : SVC(gamma=2, C=1),
        'GaussianProcessClassifier' : GaussianProcessClassifier(1.0 * RBF(1.0), warm_start=True),
        'DecisionTreeClassifier' : DecisionTreeClassifier(max_depth=5),
        'RandomForestClassifier' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        'MLPClassifier' : MLPClassifier(alpha=1),
        'AdaBoostClassifier' : AdaBoostClassifier(),
        'GaussianNB' : GaussianNB(),
        'QuadraticDiscriminantAnalysis' : QuadraticDiscriminantAnalysis()}

    # Define the output file
    output_file = results_dir + '/same_ATC_analysis.csv'
    columns_out = ['ATC mean', 'ATC mean of std', 'No-ATC mean', 'No-ATC mean of std', 'Num DC']
    output_df = pd.DataFrame(columns=columns_out)


    ########################
    #### PRE-PROCESSING ####
    ########################

    # Get the columns
    columns = []
    guild_info = []
    link_info = []
    seed_info = []

    for tprof in ('N', 'E'):
        #for top in ('5', '10', '20', '50', '100'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000', '100'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
                guild_info.append(string)
    for tprof in ('F'):
        #for top in ('5', '10', '20', '50'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
                guild_info.append(string)
    for tprof in ('LN', 'LE'):
        for an in ('sp', 'dp'):
            string = tprof+an
            columns.append(string)
            link_info.append(string)
    tprof = 'LF'
    for an in ('sp', 'dp'):
        string = tprof+an
        columns.append(string)
        link_info.append(string)
    for tprof in ('SN', 'SP'):
        for an in ('sp', 'dp', 'ji'):
            string = tprof+an
            columns.append(string)
            if an != 'sp':
                seed_info.append(string)
    columns.append('SFsp')
    seed_info.append('SFsp')
    columns.append('SFdp')
    seed_info.append('SFdp')
    columns.append('Struc')
    columns.append('Comb')

    # MANUAL SELECTION OF THE BEST FEATURES #############################################

    guild_info = ['Nper0.1dp', 'Nper2.5sp', 'E100sp', 'E100dp', 'Fper0.1dp', 'Fper0.1sp']
    link_info = ['LNdp', 'LEsp', 'LFdp']
    seed_info = ['SNdp', 'SPdp', 'SFdp']

    #####################################################################################


    # Load data frame
    df = pd.read_csv(data_frame_file, index_col=0)

    dc_data = df[df['Comb'] == 1]
    num_dc = len(dc_data.index)
    print('Number of DC with at least 3 targets: {}'.format(num_dc))

    # Replace the None values in Struc by nan
    df = df.replace(to_replace={'Struc':{'None':np.nan}})
    # Replace the NA values in Struc by nan
    df = df.replace(to_replace={'Struc':{'NA':np.nan}})

    # Deleting the spearman for seeds. It is useless
    df = df.drop('SNsp', axis=1)
    df = df.drop('SPsp', axis=1)

    #df = df.dropna()

    #dc_data = df[df['Comb'] == 1]
    #num_dc = len(dc_data.index)
    #print('Number of DC after removing nan: {}'.format(num_dc))


    ##########################################
    #### ANALYSIS BY SAME / DIFFERENT ATC ####
    ##########################################

    for atc in atc_groups:

        print('\nANALIZING ATC: {}'.format(atc))

        atc_rows = []
        no_atc_rows = []

        for index, row in df.iterrows():

            (drug1, drug2) = index.split('---')

            if drug1 not in dcdb2atc_lvl1:
                #print('Drug without ATC: {}'.format(drug1))
                continue
            if drug2 not in dcdb2atc_lvl1:
                #print('Drug without ATC: {}'.format(drug2))
                continue
            if atc in dcdb2atc_lvl1[drug1] and atc in dcdb2atc_lvl1[drug2]:
                atc_rows.append(index)
            else:
                no_atc_rows.append(index)

        df_atc = df.ix[atc_rows] # Create a table with the rows selected (the pairs with the same ATC)
        df_noatc = df.ix[no_atc_rows] # Create a table with the pairs with different ATC

        same_atc_num_dc = None # Defining the variable
        results_analysis = {'pair ATC' : {}, 'not pair ATC' : {}} # Defining the variable

        for (df_analysis, analysis) in [(df_atc, 'same ATC'), (df_noatc, 'different ATC')]:

            print('\n{}\n'.format(analysis))

            # Get the number of DC for the given number of targets
            dc_data = df_analysis[df_analysis['Comb'] == 1]
            num_dc = len(dc_data.index)
            print('Number of DC for {} {}: {}'.format(analysis, atc, num_dc))

            # Option of removing all the nan values, no matter the type of analysis
            if remove_all_nan == True:
                df_analysis = df_analysis.dropna()

            df_all = df_analysis[guild_info+seed_info+link_info+['Struc']+['Comb']]
            df_main = df_analysis[guild_info+seed_info+link_info+['Comb']]
            df_guild = df_analysis[guild_info+['Comb']]
            df_seed = df_analysis[seed_info+['Comb']]
            df_link = df_analysis[link_info+['Comb']]
            df_struc = df_analysis[['Struc']+['Comb']]

            # Getting the data frame of the analysis
            if type_analysis == 'all':
                df_ml = df_all
            elif type_analysis == 'main':
                df_ml = df_main
            elif type_analysis == 'guild':
                df_ml = df_guild
            elif type_analysis == 'seeds':
                df_ml = df_seed
            elif type_analysis == 'linkers':
                df_ml = df_link
            elif type_analysis == 'struc':
                df_ml = df_struc

            # Remove nan values
            df_ml = df_ml.dropna()

            # Get the drug combination data and the number of drug combinations
            dc_data = df_ml[df_ml['Comb'] == 1]
            num_dc = len(dc_data.index)
            print('Number of DC after removing nan: {}'.format(num_dc))

            if analysis == 'same ATC':
                same_atc_num_dc = num_dc # Get the number of DC when same ATC

            if analysis == 'different ATC':
                print('Number of DC used: {}' .format(same_atc_num_dc))

            # Stop if the number of drug interactions is smaller than the number of cross-validations!!
            if num_dc < n_fold:
                print('Not possible to do the analysis for ATC {}. The number of positive samples is {} and the n-fold is {}\n'.format(atc, num_dc, n_fold))
                break

            # Get the non-drug combination data
            ndc_data = df_ml[df_ml['Comb'] == 0]

            # Obtain the different non-drug combination groups to repeat the analysis
            # We will use the number of drug combinations when ATCs are the same in order to have same number of samples 
            print('Building {} repetition groups of {} (same) DC and {} (different) non-DC'.format(repetitions,same_atc_num_dc,same_atc_num_dc))
            ndc_repetitions = obtain_n_groups_of_k_length(ndc_data, repetitions, same_atc_num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times
            #print(ndc_repetitions)

            mean_aucs = []
            std_aucs = []

            for ndc_data_equal in ndc_repetitions:

                num_items_group = int( float(same_atc_num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation

                print('Building {} fold groups of {} DC and {} non-DC'.format(n_fold,num_items_group,num_items_group))
                dc_groups = obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group) # Defining the drug combination groups in each cross-validation step
                ndc_groups = obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group) # Defining the non-drug combination groups in each cross-validation step
                #print dc_groups
                #print ndc_groups
                merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

                mean, var, std = run_nfold_crossvalidation_scikit(n_fold, merged_groups, classifiers[classifier])
                mean_aucs.append(mean)
                std_aucs.append(std)
                #run_nfold_crossvalidation_testing_classifiers(n_fold, merged_groups)

            final_mean = np.mean(mean_aucs)
            mean_std = np.mean(std_aucs)
            std_means = np.std(mean_aucs)
            print('FINAL MEAN: {}'.format(final_mean))
            print('MEAN of STD: {}'.format(mean_std))
            print('STD of MEANS: {}'.format(std_means))

    return


def obtain_n_groups_of_k_length(my_df, n, k):
    """"""

    groups = []

    for y in xrange(n):
        selection = pd.DataFrame()
        new_df = my_df.sample(n=k)
        my_df = my_df.loc[~my_df.index.isin(new_df.index)]
        groups.append(new_df)

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
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        mean_tpr += interp(mean_fpr, fpr, tpr)
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

    return mean_auc, var_auc, std_auc


if  __name__ == "__main__":
    main()
