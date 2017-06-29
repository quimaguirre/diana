import cPickle
import math
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
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, ExtraTreesClassifier
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

    # AUC results file
    auc_results_table = results_dir + "/auc_results_table_crossval.csv"
    auc_results_table_tsv = results_dir + '/tables/feature_selection_crossval_v3.tsv'
    auc_results_table_tsv_groups = results_dir + '/tables/feature_selection_crossval_groups_v3.tsv'

    # Remove always all nan
    remove_all_nan = True

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    type_analysis = 'all'
    print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))

    # Defining the parameters of the cross-validation test
    repetitions = 25 # Number of repetititons
    n_fold = 10     # Number of folds
    classifier = 'SVC'
    classifiers = {
        'KNeighbors' : KNeighborsClassifier(3),
        'SVC' : SVC(),
        'SVC linear' : SVC(kernel="linear", C=0.025),
        'SVC rbf' : SVC(gamma=2, C=1),
        'SVC best 1' : SVC(gamma=0.1, C=1000.0, kernel='rbf'),
        'SVC best 2' : SVC(gamma=1.0, C=10.0, kernel='rbf'),
        'GaussianProcess' : GaussianProcessClassifier(1.0 * RBF(1.0), warm_start=True),
        'DecisionTree' : DecisionTreeClassifier(max_depth=5),
        'RandomForest' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        'MLP' : MLPClassifier(alpha=1),
        'AdaBoost' : AdaBoostClassifier(),
        'GaussianNB' : GaussianNB(),
        'QuadraticDiscriminant' : QuadraticDiscriminantAnalysis()
    }


    ########################
    #### PRE-PROCESSING ####
    ########################

    # Get the names of the features
    columns = []
    features = []
    guild_node = []
    guild_edge = []
    guild_func = []
    link_info = []
    seed_info = []
    for tprof in ('N', 'E'):
        #for top in ('5', '10', '20', '50', '100'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000', '100'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
                features.append(string)
                if tprof == 'N':
                    guild_node.append(string)
                elif tprof == 'E':
                    guild_edge.append(string)
    for tprof in ('F'):
        #for top in ('5', '10', '20', '50'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
                features.append(string)
                guild_func.append(string)
    for tprof in ('LN', 'LE'):
        for an in ('sp', 'dp'):
            string = tprof+an
            columns.append(string)
            features.append(string)
            link_info.append(string)
    tprof = 'LF'
    for an in ('sp', 'dp'):
        string = tprof+an
        columns.append(string)
        features.append(string)
        link_info.append(string)
    for tprof in ('SN', 'SP'):
#        for an in ('sp', 'dp', 'ji'):
        for an in ('dp', 'ji'):
            string = tprof+an
            columns.append(string)
            features.append(string)
            if an != 'sp':
                seed_info.append(string)
    columns.append('SFsp')
    features.append('SFsp')
    seed_info.append('SFsp')
    columns.append('SFdp')
    features.append('SFdp')
    seed_info.append('SFdp')
    columns.append('Struc')
    features.append('Struc')
    columns.append('Comb')


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


    ##########################
    #### FEATURE ANALYSIS ####
    ##########################

    #### USING AUC CALCULATION WITH CROSS-VALIDATION ####

    columns_auc = ['AUC']
    results_df = pd.DataFrame(columns=columns_auc)
    feature_results = {}


    all_features = list(columns)
    all_features.remove('Comb')

    for feature in all_features:

        print('\nFEATURE: {}\n'.format(feature))

        # Get the number of DC for the given number of targets
        dc_data = df[df['Comb'] == 1]
        num_dc = len(dc_data.index)
        print('Number of DC for feature {}: {}'.format(feature, num_dc))

        # Option of removing all the nan values, no matter the type of analysis
        if remove_all_nan == True:
            df = df.dropna()

        df_ml = df[[feature]+['Comb']]

        # Remove nan values
        df_ml = df_ml.dropna()

        # Get the drug combination data and the number of drug combinations
        dc_data = df_ml[df_ml['Comb'] == 1]
        num_dc = len(dc_data.index)
        print('Number of DC after removing nan: {}'.format(num_dc))

        # Stop if the number of drug interactions is smaller than the number of cross-validations!!
        if num_dc < n_fold:
            print('Not possible to do the analysis for {} targets. The number of positive samples is {} and the n-fold is {}\n'.format(num, num_dc, n_fold))
            continue

        # Get the non-drug combination data
        ndc_data = df_ml[df_ml['Comb'] == 0]

        # Obtain the different non-drug combination groups to repeat the analysis
        ndc_repetitions = obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times
        #print(ndc_repetitions)

        mean_aucs = []
        std_aucs = []

        for ndc_data_equal in ndc_repetitions:

            num_items_group = int( float(num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation

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
        tdf = pd.DataFrame([[final_mean]], columns=columns_auc, index=[feature])
        results_df = results_df.append(tdf)
        feature_results.setdefault(feature, {})
        feature_results[feature].setdefault('AUC crossval', {})
        feature_results[feature]['AUC crossval'] = final_mean


    #### USING AUC CALCULATION ####

    for feature in all_features:

        #print('\n{}\n'.format(feature))
        df_class = df.copy()
        if feature == 'Struc':
            #df_class = df_class.dropna()
            df_class['Struc'] = pd.to_numeric(df_class['Struc'])

        auc = calculate_roc_curve_scikit(df_class, feature, display=False)
        feature_results.setdefault(feature, {})
        feature_results[feature].setdefault('AUC', {})
        feature_results[feature]['AUC'] = auc
        #tdf = pd.DataFrame([[auc]], columns=columns_auc, index=[feature])
        #results_df = results_df.append(tdf)


    #### USING EXTRA TREES CLASSIFIER ####

    df_class = df.copy()
    df_class = df_class.dropna()
    df_main = df_class[features]
    targets = df_class['Comb']
    model = ExtraTreesClassifier()
    model.fit(df_main, targets)
    scores = model.feature_importances_

    for x in xrange(len(features)):
        #print('{}\t{}'.format(features[x], scores[x]))
        feature_results.setdefault(features[x], {})
        feature_results[features[x]].setdefault('ExtraTrees', {})
        feature_results[features[x]]['ExtraTrees'] = scores[x]


    #### USING Z-SCORE ####

    df_z = df.copy()
    df_z = df_z.dropna()
    df_main = df_z[features]
    targets = df_class['Comb']

    for feature in all_features:
        if feature == 'Struc':
            #df_class = df_class.dropna()
            df_main['Struc'] = pd.to_numeric(df_main['Struc'])

        scores = df_main[feature]
        zscore = calculate_zscore(scores.as_matrix(), targets.as_matrix())

        feature_results.setdefault(feature, {})
        feature_results[feature].setdefault('ZScore', {})
        feature_results[feature]['ZScore'] = zscore


    # Get data frame
    #results_df.to_csv(auc_results_table)


    #### PRINT ####

    fo = open(auc_results_table_tsv, 'w')

    for feature, results in sorted(feature_results.iteritems(), key=lambda (x, y): y['AUC crossval'], reverse = True):

        print(feature, results)
        fo.write('{}\t{}\t{}\t{}\t{}\n'.format( feature, results['AUC crossval'], results['AUC'], results['ExtraTrees'], results['ZScore'] ))

    fo.close()

    fo2 = open(auc_results_table_tsv_groups, 'w')

    for (name, type_data, group) in [ ('GUILD','Node',guild_node), ('GUILD','Edge',guild_edge), ('GUILD','Functional',guild_func), ('Linkers','Linkers', link_info), ('Seeds','Seeds', seed_info), ('Struc.', 'Struc.', ['Struc']) ]:

        #fo2.write('# GROUP {}\n'.format(name))

        for feature, results in sorted(feature_results.iteritems(), key=lambda (x, y): y['AUC crossval'], reverse = True):

            if feature in group:

                if feature[-2:] == 'sp':
                    type_comp = 'sp'
                    comparison = 'Spearman'
                elif feature[-2:] == 'dp':
                    type_comp = 'dp'
                    comparison = 'Dot product'
                elif feature == 'Struc':
                    comparison = '-'

                if feature[1:4] == 'per':
                    type_threshold = 'per'
                    number = feature.split(type_threshold)[1].split(type_comp)[0]
                    threshold = 'top {}%'.format(number)
                elif feature[1:3] == 'it':
                    type_threshold = 'it'
                    number = feature.split(type_threshold)[1].split(type_comp)[0]
                    threshold = 'top {} items'.format(number)
                elif feature[1:4] == '100':
                    threshold = 'complete profile'
                else:
                    threshold = '-'


                fo2.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( feature, name, type_data, threshold, comparison, results['AUC crossval']))

    fo2.close()

    return


def calculate_roc_curve_scikit(data_frame, classifier, display=True):
    """Calculate the ROC curve using scikit learn"""

    df = data_frame[[classifier, 'Comb']]
    df = df[np.isfinite(df[classifier])]
    y = df['Comb']
    scores = df[classifier]
    pd.set_option('display.max_rows', len(scores))
    
    #print scores
    #print np.any(np.isnan(scores))
    #print np.all(np.isfinite(scores))

    fpr, tpr, thresholds = metrics.roc_curve(y, scores)
    auc = metrics.roc_auc_score(y, scores)
    #print('SCIKIT AUC: {}\n'.format(auc))    
    #print('SCIKIT THRESHOLDS: {}'.format(thresholds))

    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curve SciKit Learn ({})'.format(classifier))
    plt.legend(loc="lower right")
    if display == True:
        plt.show()
    plt.clf()

    return auc


def calculate_zscore(scores, targets):
    """Calculate the Z-score test of 2 samples"""

    zscore = float(( np.mean(scores) - np.mean(targets) )) / float(math.sqrt( float(np.var(scores))/float(len(scores)) + float(np.var(targets))/float(len(targets)) ))

    return zscore


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
        y_score = clf.decision_function(X_test)

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
