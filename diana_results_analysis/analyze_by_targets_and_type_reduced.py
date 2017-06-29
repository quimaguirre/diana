import cPickle
import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir
from os.path import isfile, isdir, join
import pandas as pd
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, xlabel, ylabel
import re
import sys
import scipy.stats as stats
from scipy import interp

# Scikit Learn modules
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

    # Get greater or smaller number of targets
    greater_or_smaller = 'greater'

    # Plot of results name
    plot_name = results_dir + '/figures/targets_and_type_{}.eps'.format(greater_or_smaller)

    # File with results of Mann Whitney tests:
    mannwhitney_file = results_dir + '/tables/targets_and_type_{}_mannwhitney.txt'.format(greater_or_smaller)
    mannwhitney_struc_file = results_dir + '/tables/targets_and_type_{}_mannwhitneystruc.txt'.format(greater_or_smaller)

    # Results table:
    results_table = results_dir + '/tables/targets_and_type_{}_table.txt'.format(greater_or_smaller)


    # Number of targets
    num_targets = [3,4,5,6,7,8,9]
    #num_targets = [3,4]

    # Remove always all nan
    remove_all_nan = True

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    #type_analysis = 'seeds'
    #print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))
    types_analysis = ['node', 'edge', 'function', 'struc', 'random']
    types_analysis_labels = ['Nodes', 'Edges', 'Functions', 'Structure', 'Random']

    # Defining the parameters of the cross-validation test
    repetitions = 25 # Number of repetititons
    n_fold = 2     # Number of folds
    min_num_dc_group = 10
    classifier = 'SVC best 1'
    classifiers = {
        'KNeighbors' : KNeighborsClassifier(3),
        'SVC' : SVC(),
        'SVC linear' : SVC(kernel="linear", C=0.025),
        'SVC rbf' : SVC(gamma=2, C=1),
        'SVC best 1' : SVC(gamma=1.0, C=1.0, kernel='rbf', probability=True),
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

    node_info = ['Nit500sp', 'Nper5dp', 'LNdp', 'LNsp', 'SNdp', 'SNji']
    edge_info = ['Eit1000sp', 'Eit50dp', 'LEsp', 'LEdp']
    func_info = ['Fper0.1sp', 'Fper0.1dp', 'LFdp', 'LFsp', 'SFdp', 'SFsp']

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
    #### ANALYSIS BY TARGETS ####
    #############################

    analysis_results = {} # Defining the dictionary that will store the results

    for num in num_targets:

        print('\nNUMBER OF TARGETS: {}\n'.format(num))
        if greater_or_smaller == 'greater':
            print('Analyzing drugs with number of targets greater or equal to {}\n'.format(num))
        elif greater_or_smaller == 'smaller':
            print('Analyzing drugs with number of targets smaller or equal to {}\n'.format(num))

        selected_rows = []

#        if num > 3:

        dump_file = current_dir + "/dcdb2targets.pcl"
        dcdb2targets = cPickle.load(open(dump_file))

        for index, row in df.iterrows():

            (drug1, drug2) = index.split('---')

            if greater_or_smaller == 'greater':
                if len(dcdb2targets[drug1]) >= num and len(dcdb2targets[drug2]) >= num:
                    selected_rows.append(index)
            elif greater_or_smaller == 'smaller':
                if len(dcdb2targets[drug1]) <= num and len(dcdb2targets[drug2]) <= num:
                    selected_rows.append(index)
            else:
                print('\nERROR: Please, for the parameter greater_or_smaller, select \'greater\' or \'smaller\'\n')
                sys.exit(10)

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

        ######## ANALYSIS BY TYPE OF PROGRAM
        for type_analysis in types_analysis:

            print('\nType of analysis: {}\n'.format(type_analysis))

            df_all = df_tar[node_info+edge_info+func_info+['Struc']+['Comb']]
            df_node = df_tar[node_info+['Comb']]
            df_edge = df_tar[edge_info+['Comb']]
            df_func = df_tar[func_info+['Comb']]
            df_struc = df_tar[['Struc']+['Comb']]

            # Getting the data frame of the analysis
            if type_analysis == 'all':
                df_ml = df_all
            elif type_analysis == 'node':
                df_ml = df_node
            elif type_analysis == 'edge':
                df_ml = df_edge
            elif type_analysis == 'function':
                df_ml = df_func
            elif type_analysis == 'struc':
                df_ml = df_struc
            elif type_analysis == 'random':
                #df_ml = df_tar[['Comb']]
                df_ml = df_all

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

            # Stop if the number of drug interactions is smaller than the minimum number given!!
            if num_dc < min_num_dc_group:
                print('Not possible to do the analysis for num of targets {} and type of analysis {}. The number of positive samples is {} and the minimum number per group is {}\n'.format(num, type_analysis, num_dc, min_num_dc_group))
                analysis_results.setdefault(num, {})
                analysis_results[num][type_analysis] = '-'
                continue

            # Get the non-drug combination data
            ndc_data = df_ml[df_ml['Comb'] == 0]

            # Obtain the different non-drug combination groups to repeat the analysis
            print('Building {} repetition groups of {} (same) DC and {} (different) non-DC'.format(repetitions,num_dc,num_dc))
            ndc_repetitions = obtain_n_groups_of_k_length(ndc_data, repetitions, num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times
            #print(ndc_repetitions)

            mean_aucs = [] # Here we will store the means of AUCs from the cross-validations
            std_aucs = [] # Here we will store the standard deviations of the AUCs from the cross-validations
            all_aucs = [] # Here we will store ALL the AUCs

            num_repetitions=0
            for ndc_data_equal in ndc_repetitions:

                num_repetitions+=1

                num_items_group = int( float(num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation
                if num_repetitions == 1:
                    print('Building {} fold groups of {} DC and {} non-DC x {} repetitions'.format(n_fold,num_items_group,num_items_group, repetitions))

                dc_groups = obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group) # Defining the drug combination groups in each cross-validation step
                ndc_groups = obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group) # Defining the non-drug combination groups in each cross-validation step
                #print dc_groups
                #print ndc_groups
                merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

                if type_analysis == 'random':
                    #mean, var, std, list_auc = run_nfold_crossvalidation_random(n_fold, merged_groups, classifiers[classifier])
                    mean, var, std, list_auc, list_prob = run_nfold_crossvalidation_dummy(n_fold, merged_groups, classifiers[classifier])
                else:
                    mean, var, std, list_auc = run_nfold_crossvalidation_scikit(n_fold, merged_groups, classifiers[classifier])

                mean_aucs.append(mean)
                std_aucs.append(std)
                all_aucs = all_aucs + list_auc
                #run_nfold_crossvalidation_testing_classifiers(n_fold, merged_groups)

            final_mean = np.mean(mean_aucs)
            mean_std = np.mean(std_aucs)
            std_means = np.std(mean_aucs)
            std = np.std(all_aucs)
            print('FINAL MEAN: {}'.format(final_mean))
            print('MEAN of STD: {}'.format(mean_std))
            print('STD: {}'.format(std))

            # Store the distribution of AUCs in the dictionary
            analysis_results.setdefault(num, {})
            analysis_results[num].setdefault(type_analysis, {})
            analysis_results[num][type_analysis]['all_aucs'] = all_aucs
            analysis_results[num][type_analysis]['mean'] = final_mean
            analysis_results[num][type_analysis]['std'] = std
            analysis_results[num][type_analysis]['num_dc'] = num_dc


    ####### PLOT GRAPHIC OF DISTRIBUTION OF AUCS VS. NUM. TARGETS AND METHOD
    #print(analysis_results)

    fig = figure(dpi=300)
    ax = axes()
    #hold(True)
    pos = 2

    xticks = [] # Define the places in which the labels will be
    for num in num_targets:

        positions = []
        for x in xrange(len(types_analysis)):
            positions.append(pos) # Define the positions of the boxplots
            pos+=1
        pos+=1 # Add separation between boxplot groups
        plt.axvline(x=pos,linewidth=0.3,linestyle='--',color='black',dashes=(1, 1))
        pos+=2 # Add separation between boxplot groups

        data = []
        for type_analysis in types_analysis:
            data.append(analysis_results[num][type_analysis]['all_aucs']) # Get the groups of plots that we will add

        # Boxplot group
        #bp = boxplot(data, positions = positions, widths = 0.6)
        bp = boxplot(data, positions = positions, widths = 0.6, patch_artist=True)
        setBoxColors(bp, len(types_analysis))

        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # Set axes limits and labels
    xlim(0,pos-2)
    ylim(0,1)
    if greater_or_smaller == 'greater':
        axes_labels = [r'$\geq$' + ' ' + str(x) for x in num_targets]
    elif greater_or_smaller == 'smaller':
        axes_labels = [r'$\leq$' + ' ' + str(x) for x in num_targets]
    ax.set_xticklabels(axes_labels)
    ax.set_xticks(xticks)
    xlabel('Number of targets')
    ylabel('Distribution of AUC values')

    # draw temporary color lines and use them to create a legend
    hG, = plot([1,1],'-', color='green')
    hR, = plot([1,1],'-', color='red')
    hB, = plot([1,1],'-', color='blue')
    hBl, = plot([1,1],'-', color='black')
    hW, = plot([1,1],'-', color='#aeaeae')
    lgd = ax.legend(handles=(hG, hR, hB, hBl, hW), labels=types_analysis_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hG.set_visible(False)
    hR.set_visible(False)
    hB.set_visible(False)
    hBl.set_visible(False)
    hW.set_visible(False)

    savefig(plot_name, format='eps', bbox_extra_artists=(lgd,), bbox_inches='tight')
    show()



    ################ RESULTS TABLE
    fr = open(results_table, 'w')

    # Header
    fr.write(' ')
    for method in types_analysis_labels:
        fr.write('\t{}\t \t '.format(method))
    fr.write('\n')

    for num in num_targets:
        fr.write('{}'.format(num))
        for method in types_analysis:
            mean = analysis_results[num][method]['mean']
            std = analysis_results[num][method]['std']
            num_dc = analysis_results[num][method]['num_dc']
            fr.write('\t{}\t{}\t{}'.format(mean, std, num_dc))
        fr.write('\n')
    fr.close()


    ################ MANN WHITNEY U
    fo = open(mannwhitney_file, 'w')

    mann_results = {}

    fo.write(' \t ')
    for method in types_analysis:
        fo.write('\t{}'.format(method))
    fo.write('\n')

    # Perform the comparisons
    for num in num_targets:
        mann_results.setdefault(num, {})
        for method1 in types_analysis:
            mann_results[num].setdefault(method1, {})
            for method2 in types_analysis:
                if method1 == method2:
                    mann_results[num][method1][method2] = '-'
                else:
                    method1_dist = analysis_results[num][method1]['all_aucs']
                    method2_dist = analysis_results[num][method2]['all_aucs']
                    stat, pval = stats.mannwhitneyu(method1_dist, method2_dist)
                    mann_results[num][method1][method2] = [stat, pval]

    # Write the table of crossings
    for num in num_targets:
        for method1 in types_analysis:
            fo.write('{}\t{}'.format(num, method1))
            for method2 in types_analysis:
                if method1 == method2:
                    fo.write('\t-')
                else:
                    stat, pval = mann_results[num][method1][method2]
                    fo.write('\t{}, {:.2e}'.format(stat,pval))
            fo.write('\n')

    fo.close()


    ################ MANN WHITNEY U BETWEEN STRUCTURAL VALUES
    fs = open(mannwhitney_struc_file, 'w')
    fs.write(' ')

    mann_struc_results = {}

    for num1 in num_targets:
        mann_struc_results.setdefault(num1, {})
        fs.write('\t{}'.format(num1))
        for num2 in num_targets:
            if num1 == num2:
                mann_struc_results[num1][num2] = '-'
            else:
                num1_dist = analysis_results[num1]['struc']['all_aucs']
                num2_dist = analysis_results[num2]['struc']['all_aucs']
                stat, pval = stats.mannwhitneyu(num1_dist, num2_dist)
                mann_struc_results[num1][num2] = [stat, pval]
    fs.write('\n')

    # Write the table of crossings
    for num1 in num_targets:
        fs.write('{}'.format(num1))
        for num2 in num_targets:
            if num1 == num2:
                fs.write('\t-')
            else:
                stat, pval = mann_struc_results[num1][num2]
                fs.write('\t{}, {:.2e}'.format(stat,pval))
        fs.write('\n')

    fs.close()



    return


def obtain_n_groups_of_k_length(my_df, n, k):
    """Obtain randomly and without repetition n groups of k length"""

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

    return mean_auc, var_auc, std_auc, all_auc


def run_nfold_crossvalidation_random(n, groups, classifier):
    """
    Run a cross-validation assigning randomly the prediction
    n = number of folds
    groups = list with the balanced groups in each fold of the cross-validation
    classifier = classifier used in the machine learning approach
    """

    for group in groups:
        group['Random'] = shuffle(group)

    all_auc = []
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    stdsc = StandardScaler()
    #pca = PCA()

    for x in xrange(n):

        test = groups[x]
        train_groups = [item for index,item in enumerate(groups) if index != x]
        train = pd.concat(train_groups)

        y_test = test['Comb']
        y_pred = test['Random']

        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        auc = metrics.roc_auc_score(y_test, y_pred)
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


def run_nfold_crossvalidation_dummy(n, groups, classifier):
    """
    n = number of folds
    groups = list with the balanced groups in each fold of the cross-validation
    classifier = classifier used in the machine learning approach
    """

    all_auc = []
    all_prob = []
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

        clf = DummyClassifier().fit(X_train, y_train)
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

        prob = clf.predict_proba(X_test) # Get the probability used to classify. This is a list, and there is a probability for each class
        
        classes = clf.classes_ # This is the order of the classes. The probabilities are given in this order

        for index in xrange(len(classes)):
            if classes[index] == 1:
                dc_index = index # Obtain in which position is located the probability of being drug combination

        for p in xrange(len(prob)):
            dc_prob = prob[p][dc_index] # We use the index to obtain the probability of being drug combination
            dc_label = y_test[p]
            array = [ dc_prob, dc_label ] # Create an array with the probability and the label
            all_prob.append(array) # Append the array in all_prob

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

    return mean_auc, var_auc, std_auc, all_auc, all_prob


def shuffle(dfinit, n=1, axis=0):     
    dfcopy = dfinit.copy()
    for _ in range(n):
        dfcopy.apply(np.random.shuffle, axis=axis)
    return dfcopy


def setBoxColors(bp, length):
    """
    Set the colors of the box plots groups
    Code from: http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    """
    if length != 5:
        print('\nERROR: This function is thought to work on groups of 5 types of analysis. Redo the function!!\n')
        sys.exit(10)

    setp(bp['boxes'][0], color='green')
    setp(bp['caps'][0], color='green')
    setp(bp['caps'][1], color='green')
    setp(bp['whiskers'][0], color='green')
    setp(bp['whiskers'][1], color='green')
    setp(bp['medians'][0], color='black')
    bp['boxes'][0].set_facecolor('#32f232') #green

    setp(bp['boxes'][1], color='#b80000')
    setp(bp['caps'][2], color='#b80000')
    setp(bp['caps'][3], color='#b80000')
    setp(bp['whiskers'][2], color='#b80000')
    setp(bp['whiskers'][3], color='#b80000')
    setp(bp['medians'][1], color='black')
    bp['boxes'][1].set_facecolor('#ff7373') #red

    setp(bp['boxes'][2], color='blue')
    setp(bp['caps'][4], color='blue')
    setp(bp['caps'][5], color='blue')
    setp(bp['whiskers'][4], color='blue')
    setp(bp['whiskers'][5], color='blue')
    setp(bp['medians'][2], color='black')
    bp['boxes'][2].set_facecolor('#9ed0ff') #blue

    setp(bp['boxes'][3], color='black')
    setp(bp['caps'][6], color='black')
    setp(bp['caps'][7], color='black')
    setp(bp['whiskers'][6], color='black')
    setp(bp['whiskers'][7], color='black')
    setp(bp['medians'][3], color='#fbc562')
    bp['boxes'][3].set_facecolor('#4f4f4f') #light black

    setp(bp['boxes'][4], color='black')
    setp(bp['caps'][8], color='black')
    setp(bp['caps'][9], color='black')
    setp(bp['whiskers'][8], color='black')
    setp(bp['whiskers'][9], color='black')
    setp(bp['medians'][4], color='#fbc562')
    bp['boxes'][4].set_facecolor('#aeaeae') # grey

    return


if  __name__ == "__main__":
    main()
