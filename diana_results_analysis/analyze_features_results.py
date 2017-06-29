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

# Scikit Learn modules
from sklearn import metrics
from sklearn.ensemble import ExtraTreesClassifier


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
    auc_results_table = results_dir + "/auc_results_table.csv"


    ########################
    #### PRE-PROCESSING ####
    ########################

    # Get the names of the features
    columns = []
    features = []
    for tprof in ('N', 'E'):
        #for top in ('5', '10', '20', '50', '100'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000', '100'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
                features.append(string)
    for tprof in ('F'):
        #for top in ('5', '10', '20', '50'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
                features.append(string)
    for tprof in ('LN', 'LE'):
        for an in ('sp', 'dp'):
            string = tprof+an
            columns.append(string)
            features.append(string)
    tprof = 'LF'
    for an in ('sp', 'dp'):
        string = tprof+an
        columns.append(string)
        features.append(string)
    for tprof in ('SN', 'SP'):
        for an in ('sp', 'dp', 'ji'):
            string = tprof+an
            columns.append(string)
            features.append(string)
    columns.append('SFsp')
    features.append('SFsp')
    columns.append('SFdp')
    features.append('SFdp')
    columns.append('Struc')
    features.append('Struc')
    columns.append('Comb')


    # Load data frame
    df = pd.read_csv(data_frame_file, index_col=0)

    dc_data = df[df['Comb'] == 1]
    num_dc = len(dc_data.index)
    print('Number of DC before nan: {}'.format(num_dc))

    # Replace the None values in Struc by nan
    df = df.replace(to_replace={'Struc':{'None':np.nan}})
    # Replace the NA values in Struc by nan
    df = df.replace(to_replace={'Struc':{'NA':np.nan}})

    df = df.dropna()

    dc_data = df[df['Comb'] == 1]
    num_dc = len(dc_data.index)
    print('Number of DC after removing nan: {}'.format(num_dc))


    ##########################
    #### FEATURE ANALYSIS ####
    ##########################

    #### USING AUC CALCULATION ####

    columns_auc = ['AUC']
    results_df = pd.DataFrame(columns=columns_auc)
    feature_results = {}


    all_features = list(columns)
    all_features.remove('Comb')

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


    #print(results_df)

    # Get data frame
    results_df.to_csv(auc_results_table)


    #### USING EXTRA TREES CLASSIFIER ####

    df_class = df.copy()
    df_class = df_class.dropna()
    df_main = df_class[features]
    targets = df_class['Comb']
    model = ExtraTreesClassifier(n_estimators=10,bootstrap=False)
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


    #### PRINT ####

    for feature in sorted(feature_results.iteritems(), key=lambda (x, y): y['AUC'], reverse = True):

        print feature


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


if  __name__ == "__main__":
    main()
