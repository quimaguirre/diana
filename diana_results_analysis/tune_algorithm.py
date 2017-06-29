import cPickle
import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir
from os.path import isfile, isdir, join
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, xlabel, ylabel
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

    # Results file
    results_tuning = results_dir + '/tables/tuning_results.tsv'

    # Number of targets
    num_targets = [3]

    # Remove always all nan
    remove_all_nan = True

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    type_analysis = 'all'
    print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))

    # Defining the parameters of the cross-validation test
    repetition_analysis = 2
    repetitions = 25 # Number of repetititons
    n_fold = 10     # Number of folds


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

    guild_info = ['Nit500sp', 'Nper5dp', 'Eit1000sp', 'Eit50dp', 'Fper0.1sp', 'Fper0.1dp',]
    link_info = ['LNdp', 'LNsp', 'LEsp', 'LEdp', 'LFdp', 'LFsp']
    seed_info = ['SNdp', 'SNji', 'SPdp', 'SPji', 'SFdp', 'SFsp']

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


    #############################
    #### CHOOSE A CLASSIFIER ####
    #############################

    dict_results = {} # Defining the dictionary that will store the results

    for xrep in xrange(repetition_analysis):

        for num in num_targets:

            print('\nNUMBER OF TARGETS: {}\n'.format(num))

            selected_rows = []

    #        if num > 3:

            dump_file = current_dir + "/dcdb2targets.pcl"
            dcdb2targets = cPickle.load(open(dump_file))

            for index, row in df.iterrows():

                (drug1, drug2) = index.split('---')

                if len(dcdb2targets[drug1]) >= num and len(dcdb2targets[drug2]) >= num:
                    selected_rows.append(index)

            df_tar = df.ix[selected_rows]

    #        else:
    #            df_tar = df.copy()


            # Get the number of DC for the given number of targets
            dc_data = df_tar[df_tar['Comb'] == 1]
            num_dc = len(dc_data.index)
            print('Number of DC for {} targets: {}'.format(num, num_dc))

            # Option of removing all the nan values, no matter the type of analysis
            if remove_all_nan == True:
                df_tar = df_tar.dropna()

            df_all = df_tar[guild_info+seed_info+link_info+['Struc']+['Comb']]
            df_main = df_tar[guild_info+seed_info+link_info+['Comb']]
            df_guild = df_tar[guild_info+['Comb']]
            df_seed = df_tar[seed_info+['Comb']]
            df_link = df_tar[link_info+['Comb']]
            df_struc = df_tar[['Struc']+['Comb']]

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

            # Stop if the number of drug interactions is smaller than the number of cross-validations!!
            if num_dc < n_fold:
                print('Not possible to do the analysis for {} targets. The number of positive samples is {} and the n-fold is {}\n'.format(num, num_dc, n_fold))
                continue

            # Get the non-drug combination data
            ndc_data = df_ml[df_ml['Comb'] == 0]

            # Obtain the different non-drug combination groups to repeat the analysis
            ndc_repetitions = obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times
            #ndc_same_length = obtain_n_groups_of_k_length(ndc_data, 1, num_dc)[0] # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times


            for ndc_same_length in ndc_repetitions:

                merged_groups = pd.concat([dc_data, ndc_same_length])

                X_train, y_train = merged_groups.iloc[:, :-1], merged_groups.iloc[:, -1]

                pipe_svc = Pipeline([('slc', StandardScaler()),
                                     ('clf', SVC(random_state=1))])

                param_range = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]
                param_grid = [{'clf__C': param_range,
                               'clf__kernel': ['linear']},
                              {'clf__C': param_range,
                               'clf__gamma': param_range,
                               'clf__kernel': ['rbf']}]

                gs = GridSearchCV(estimator=pipe_svc,
                                  param_grid=param_grid,
                                  scoring='accuracy',
                                  cv=10,
                                  n_jobs=-1)


                gs = gs.fit(X_train, y_train)
                result = str(gs.best_params_)
                print(result)
                print(gs.best_score_)
                dict_results.setdefault(result, 0)
                dict_results[result] += 1

            print(dict_results)

    print('\nFINAL RESULT\n')

    print(dict_results)

    fo = open(results_tuning, 'w')

    for param_comb in sorted(dict_results, reverse = True):

        print('{}\t{}\n'.format(param_comb, dict_results[param_comb]))
        fo.write('{}\t{}\n'.format(param_comb, dict_results[param_comb]))

    fo.close()

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

    return mean_auc, var_auc, std_auc, all_auc


def setBoxColors(bp, bg_color, edge_color):
    """
    Set the colors of the box plots groups
    Code from: http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    """

    setp(bp['boxes'][0], color=edge_color)
    setp(bp['caps'][0], color=edge_color)
    setp(bp['caps'][1], color=edge_color)
    setp(bp['whiskers'][0], color=edge_color)
    setp(bp['whiskers'][1], color=edge_color)
    setp(bp['medians'][0], color='black')
    bp['boxes'][0].set_facecolor(bg_color)

    return


if  __name__ == "__main__":
    main()
