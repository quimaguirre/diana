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
    #auc_results_table = results_dir + "/classifier_results_table.csv"

    # Plot name
    fig_format = 'eps'
    #plot_name = results_dir + '/figures/classifiers_selection_before_tuning.'+fig_format
    plot_name = results_dir + '/figures/classifiers_selection_after_tuning.'+fig_format

    # Number of targets
    num_targets = [3]

    # Remove always all nan
    remove_all_nan = True

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    type_analysis = 'all'
    print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))

    # Defining the parameters of the cross-validation test
    repetitions = 25 # Number of repetitions
    n_fold = 10     # Number of folds
    # classifiers = {
    #     'KNei' : KNeighborsClassifier(3),
    #     'SVC' : SVC(),
    #     'SVC0.025' : SVC(kernel="linear", C=0.025),
    #     'SVC1' : SVC(gamma=2, C=1),
    #     'SVCB1' : SVC(gamma=0.1, C=1000.0, kernel='rbf'),
    #     'SVCB2' : SVC(gamma=1.0, C=10.0, kernel='rbf'),
    #     'GauPr' : GaussianProcessClassifier(1.0 * RBF(1.0), warm_start=True),
    #     'DTree' : DecisionTreeClassifier(max_depth=5),
    #     'RFor' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
    #     'MLP' : MLPClassifier(alpha=1),
    #     'AdaB' : AdaBoostClassifier(),
    #     'GauNB' : GaussianNB(),
    #     'QDiscr' : QuadraticDiscriminantAnalysis()
    # }

    classifiers = {
        'KNeighbors' : KNeighborsClassifier(3),
        'SVC' : SVC(),
        'SVC linear' : SVC(kernel="linear", C=0.025),
        'SVC rbf' : SVC(gamma=2, C=1),
        'SVC best 1' : SVC(gamma=1.0, C=1.0, kernel='rbf'),
        'SVC best 2' : SVC(gamma=0.1, C=10.0, kernel='rbf'),
        'SVC best 3' : SVC(gamma=0.01, C=1000.0, kernel='rbf'),
        'SVC best 4' : SVC(gamma=0.1, C=1.0, kernel='rbf'),
        'DecisionTree' : DecisionTreeClassifier(max_depth=5),
        'RandomForest' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        'MLP' : MLPClassifier(alpha=1),
        'AdaBoost' : AdaBoostClassifier(),
        'GaussianNB' : GaussianNB(),
        'QuadraticDiscr.' : QuadraticDiscriminantAnalysis()
    }
#    classifiers_order = ['AdaB', 'DTree', 'GauNB', 'GauPr', 'KNei', 'MLP', 'QDiscr', 'RFor', 'SVC', 'SVC0.025', 'SVC1', 'SVCB1', 'SVCB2']
#    classifiers_order = ['AdaB', 'DTree', 'GauNB', 'GauPr', 'KNei', 'MLP', 'QDiscr', 'RFor', 'SVC', 'SVC0.025', 'SVC1']
#    classifiers_order = ['AdaBoost', 'DecisionTree', 'GaussianNB', 'KNeighbors', 'MLP', 'QuadraticDiscr.', 'RandomForest', 'SVC', 'SVC linear', 'SVC rbf']
    classifiers_order = ['AdaBoost', 'DecisionTree', 'GaussianNB', 'KNeighbors', 'MLP', 'QuadraticDiscr.', 'RandomForest', 'SVC', 'SVC linear', 'SVC rbf', 'SVC best 1', 'SVC best 2', 'SVC best 3', 'SVC best 4']
#    classifiers_order = ['SVC', 'SVC0.025', 'SVC1', 'SVCB1', 'SVCB2']

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

    #guild_info = ['Nit500sp', 'Nit1000dp', 'Eit1000sp', 'Eit50dp', 'Fit500dp', 'Fper0.1sp']
    #link_info = ['LNdp', 'LNsp', 'LEsp', 'LEdp', 'LFdp', 'LFsp']
    #seed_info = ['SNdp', 'SNji', 'SPdp', 'SPji', 'SFdp', 'SFsp']

    guild_info = ['Nper0.5sp', 'Nper0.1dp', 'Eper10sp', 'Eper2.5dp', 'Fper0.1sp', 'Fper0.1dp',]
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

    analysis_results = {} # Defining the dictionary that will store the results

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
        #print(ndc_repetitions)

        for classifier in classifiers:

            print('\nCLASSIFIER: {}\n'.format(classifier))

            mean_aucs = [] # Here we will store the means of AUCs from the cross-validations
            std_aucs = [] # Here we will store the standard deviations of the AUCs from the cross-validations
            all_aucs = [] # Here we will store ALL the AUCs

            for ndc_data_equal in ndc_repetitions:

                num_items_group = int( float(num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation

                dc_groups = obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group) # Defining the drug combination groups in each cross-validation step
                ndc_groups = obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group) # Defining the non-drug combination groups in each cross-validation step
                #print dc_groups
                #print ndc_groups
                merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

                mean, var, std, list_auc = run_nfold_crossvalidation_scikit(n_fold, merged_groups, classifiers[classifier])
                mean_aucs.append(mean)
                std_aucs.append(std)
                all_aucs = all_aucs + list_auc
                #run_nfold_crossvalidation_testing_classifiers(n_fold, merged_groups)

            final_mean = np.mean(mean_aucs)
            mean_std = np.mean(std_aucs)
            std_means = np.std(mean_aucs)
            print('FINAL MEAN: {}'.format(final_mean))
            print('MEAN of STD: {}'.format(mean_std))
            print('STD of MEANS: {}'.format(std_means))

            # Store the distribution of AUCs in the dictionary
            analysis_results[classifier] = all_aucs


    ####### PLOT GRAPHIC OF DISTRIBUTION OF AUCS VS. NUM. TARGETS AND METHOD
    #print(analysis_results)

    fig = figure(dpi=300)
    ax = axes()
    #hold(True)
    pos = 1
    col_num = 0

    xticks = [] # Define the places in which the labels will be
    xlabels = [] # Define the labels (the names of the classifiers)
    #colors = [ ['#9ed0ff, blue'], ['#32f232', 'green'], ['#fbc562', '#d48900'], ['#ff7373', '#b80000'], ['grey', 'black'] ]

    for classifier in classifiers_order:

        positions = []
        positions.append(pos) # Define the positions of the boxplots
        pos+=2 # Add separation between boxplots
        xlabels.append(classifier) # Add the classifier used at the x axis

        # Boxplot group
        #bp = boxplot(data, positions = positions, widths = 0.6)
        bp = boxplot(analysis_results[classifier], positions = positions, widths = 0.6, patch_artist=True)
        #setBoxColors(bp, colors[col_num][0], colors[col_num][1])

        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # Set axes limits and labels
    xlim(0,pos-1)
    ylim(0,1)
    ax.set_xticklabels(xlabels)
    ax.set_xticks(xticks)
    xlabel('Classifiers')
    ylabel('Distribution of AUC values')

    fig.autofmt_xdate()
    savefig(plot_name, format=fig_format)
    show()


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
