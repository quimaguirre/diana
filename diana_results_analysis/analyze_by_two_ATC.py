import cPickle
import matplotlib.pyplot as plt
import mysql.connector
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

    # Plot of results name
    plot_name = results_dir + '/figures/atc_pairs.eps'

    # ATC groups to analyze
    atc_groups = set()
    dump_file = current_dir + "/dcdb2atc_lvl1.pcl"
    dcdb2atc_lvl1 = cPickle.load(open(dump_file))
    for dcdb in dcdb2atc_lvl1:
        atc_groups = atc_groups.union(dcdb2atc_lvl1[dcdb]) # Join all the ATCs in one group

    atc_list = sorted(list(atc_groups))

    # Define the pairs of ATC groups to analyze
    pairs_atc = []
    for atc1 in atc_groups:
        for atc2 in atc_groups:
            #if atc1 != atc2:
                comb1 = '{}---{}'.format(atc1, atc2)
                comb2 = '{}---{}'.format(atc2, atc1)
                if comb1 not in pairs_atc and comb2 not in pairs_atc:
                    pairs_atc.append(comb1)

    # Remove always all nan
    remove_all_nan = False

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    type_analysis = 'guild'
    print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))

    # Best predictions list
    drug_comb_scores_file = results_dir + '/tables/drug_comb_scores_list_atcs.txt'
    method_best_predictions = 'guild'
    cnx = mysql.connector.connect(user='quim', password='',
                                  host='localhost',
                                  database='BIANA_JAN_2017')
    dcdb2name = obtain_dcdb2name(cnx)
    dcdb2atc = obtain_dcdb2atc(cnx)

    # Defining the parameters of the cross-validation test
    repetitions = 10 # Number of repetitions
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

    # Define the output file
    output_file = results_dir + '/tables/ATC_pairs.csv'
    columns_out = ['Pair mean', 'No-pair mean', 'Pair mean of std', 'No-pair mean of std', 'Num DC']
    output_df = pd.DataFrame(columns=columns_out)

    # Results table
    results_table = results_dir + '/tables/ATC_results_table.txt'


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

    guild_info = ['Nit500sp', 'Nper5dp', 'Eit1000sp', 'Eit50dp', 'Fper0.1sp', 'Fper0.1dp']
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


    ##########################################
    #### ANALYSIS BY SAME / DIFFERENT ATC ####
    ##########################################

    results_analysis = {} # Defining the variable that will store the results

    for atc_pair in pairs_atc:

        (atc1, atc2) = atc_pair.split('---')

        print('\nANALIZING ATC: {} and ATC: {}'.format(atc1, atc2))

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
            if atc1 in dcdb2atc_lvl1[drug1] and atc2 in dcdb2atc_lvl1[drug2]:
                atc_rows.append(index)
            elif atc2 in dcdb2atc_lvl1[drug1] and atc1 in dcdb2atc_lvl1[drug2]:
                atc_rows.append(index)
            else:
                no_atc_rows.append(index)

        df_atc = df.ix[atc_rows] # Create a table with the rows selected (the pairs with the same ATC)
        df_noatc = df.ix[no_atc_rows] # Create a table with the pairs with different ATC

        same_atc_num_dc = None # Defining the variable

        results_analysis.setdefault(atc_pair, {})
        #results_analysis[atc_pair] = {'pair ATC' : {}, 'not pair ATC' : {}} # Defining the variable
        results_analysis[atc_pair] = {'pair ATC' : {}, 'random' : {}} # Defining the variable

        #for (df_analysis, analysis) in [(df_atc, 'pair ATC'), (df_noatc, 'not pair ATC')]:
        for (df_analysis, analysis) in [(df_atc, 'pair ATC'), (df_atc, 'random')]:

            print('\n{}\n'.format(analysis))

            # Get the number of DC for the given number of targets
            dc_data = df_analysis[df_analysis['Comb'] == 1]
            num_dc = len(dc_data.index)
            print('Number of DC for {} {}: {}'.format(analysis, atc_pair, num_dc))

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

            if analysis == 'random':
                df_ml = df_guild

            # Remove nan values
            df_ml = df_ml.dropna()

            # Get the drug combination data and the number of drug combinations
            dc_data = df_ml[df_ml['Comb'] == 1]
            num_dc = len(dc_data.index)
            print('Number of DC after removing nan: {}'.format(num_dc))

            if analysis == 'pair ATC':
                same_atc_num_dc = num_dc # Get the number of DC when same ATC

            #if analysis == 'not pair ATC':
            if analysis == 'random':
                print('Number of DC used: {}'.format(same_atc_num_dc))

            # Stop if the number of drug interactions is smaller than the number of cross-validations!!
            if num_dc < n_fold:
                print('Not possible to do the analysis for ATC pair {}. The number of positive samples is {} and the n-fold is {}\n'.format(atc_pair, num_dc, n_fold))
                #results_analysis[atc_pair] = {'pair ATC' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}, 'not pair ATC' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}} # Defining the variable
                results_analysis[atc_pair] = {'pair ATC' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}, 'random' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}} # Defining the variable
                break
            # Stop if the number of drug interactions is smaller than the minimum number given!!
            if num_dc < min_num_dc_group:
                print('Not possible to do the analysis for ATC pair {}. The number of positive samples is {} and the minimum number per group is {}\n'.format(atc_pair, num_dc, min_num_dc_group))
                #results_analysis[atc_pair] = {'pair ATC' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}, 'not pair ATC' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}} # Defining the variable
                results_analysis[atc_pair] = {'pair ATC' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}, 'random' : {'mean':'-','std':'-','num_dc':int(same_atc_num_dc)}} # Defining the variable
                break

            # Get the non-drug combination data
            ndc_data = df_ml[df_ml['Comb'] == 0]

            # Obtain the different non-drug combination groups to repeat the analysis
            # We will use the number of drug combinations when ATCs are the same in order to have same number of samples 
            print('Building {} repetition groups of {} (same) DC and {} (different) non-DC'.format(repetitions,same_atc_num_dc,same_atc_num_dc))
            ndc_repetitions = obtain_n_groups_of_k_length(ndc_data, repetitions, same_atc_num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times
            #print(ndc_repetitions)

            mean_aucs = [] # Here we will store the means of AUCs from the cross-validations
            std_aucs = [] # Here we will store the standard deviations of the AUCs from the cross-validations
            all_aucs = [] # Here we will store ALL the AUCs
            all_probs = [] # Here we store all the probabilities and labels

            for ndc_data_equal in ndc_repetitions:

                num_items_group = int( float(same_atc_num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation

                print('Building {} fold groups of {} DC and {} non-DC'.format(n_fold,num_items_group,num_items_group))
                dc_groups = obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group) # Defining the drug combination groups in each cross-validation step
                ndc_groups = obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group) # Defining the non-drug combination groups in each cross-validation step
                #print dc_groups
                #print ndc_groups
                merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

                if analysis == 'random':
                    #mean, var, std, list_auc = run_nfold_crossvalidation_random(n_fold, merged_groups, classifiers[classifier])
                    mean, var, std, list_auc, list_prob = run_nfold_crossvalidation_dummy(n_fold, merged_groups, classifiers[classifier])
                else:
                    mean, var, std, list_auc, list_prob = run_nfold_crossvalidation_scikit(n_fold, merged_groups, classifiers[classifier])
                mean_aucs.append(mean)
                std_aucs.append(std)
                all_aucs = all_aucs + list_auc
                all_probs = all_probs + list_prob
                #run_nfold_crossvalidation_testing_classifiers(n_fold, merged_groups)

            final_mean = np.mean(mean_aucs)
            mean_std = np.mean(std_aucs)
            std_means = np.std(mean_aucs)
            std = np.std(all_aucs)
            results_analysis[atc_pair][analysis]['mean'] = final_mean
            results_analysis[atc_pair][analysis]['std'] = std
            results_analysis[atc_pair][analysis]['num_dc'] = int(same_atc_num_dc)
            results_analysis[atc_pair][analysis]['all_aucs'] = all_aucs
            results_analysis[atc_pair][analysis]['all_probs'] = all_probs
            print('MEAN AUCs: {}'.format(final_mean))
            print('STD AUCs: {}'.format(std))

        #row = [ results_analysis[atc_pair]['pair ATC']['mean'], results_analysis[atc_pair]['not pair ATC']['mean'], results_analysis[atc_pair]['pair ATC']['std'], results_analysis[atc_pair]['not pair ATC']['std'], results_analysis[atc_pair]['pair ATC']['num_dc'] ]
        row = [ results_analysis[atc_pair]['pair ATC']['mean'], results_analysis[atc_pair]['random']['mean'], results_analysis[atc_pair]['pair ATC']['std'], results_analysis[atc_pair]['random']['std'], results_analysis[atc_pair]['pair ATC']['num_dc'] ]
        odf = pd.DataFrame([row], columns=columns_out, index=[atc_pair])
        output_df = output_df.append(odf)

    print(output_df)
    output_df.to_csv(output_file)


    ####### RANKING OF DRUG COMBINATIONS BY THE MEAN OF THE SCORES OF THE PREDICTIONS
    scores_fd = open(drug_comb_scores_file, 'w')

    for atc_pair in pairs_atc:

        if 'all_probs' not in results_analysis[atc_pair]['pair ATC']:
            continue
        scores_fd.write('##### {} #####\n'.format(atc_pair))
        scores_best_predictions = results_analysis[atc_pair]['pair ATC']['all_probs']
        dc2scoresmean = obtain_drug_combination_scores_mean(scores_best_predictions)
        for dc, mean in sorted(dc2scoresmean.iteritems(), key=lambda (x, y): y, reverse=True):
            dcdb1, dcdb2 = dc.split('---')
            name1 = dcdb2name[dcdb1]
            try:
                atcs1 = ','.join(list(obtain_set_of_first_letter_ATC(dcdb2atc[dcdb1])))
            except:
                atcs1 = '-'
            name2 = dcdb2name[dcdb2]
            try:
                atcs2 = ','.join(list(obtain_set_of_first_letter_ATC(dcdb2atc[dcdb2])))
            except:
                atcs2 = '-'
            scores_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(name1, atcs1, name2, atcs2, mean))

    scores_fd.close()


    ####### PLOT GRAPHIC OF DISTRIBUTION OF AUCS VS. CLASSIFICATION AND METHOD

    figsize = (len(atc_pair)*2, 8)
    fig = figure(figsize=figsize, dpi=300)
    ax = axes()
    #hold(True)
    pos = 2
    axes_labels = [] # The labels that we will add

    xticks = [] # Define the places in which the labels will be
    for atc_pair in sorted(pairs_atc):

        if results_analysis[atc_pair]['pair ATC']['mean'] == '-': # If the analysis has not been done, skip
            continue

        (atc1, atc2) = atc_pair.split('---')
        label = '{} vs. {}'.format(atc1, atc2)
        axes_labels.append(label) # If we have result, we add the label in the list


        positions = []
        #for x in xrange(len(['pair ATC', 'not pair ATC'])):
        for x in xrange(len(['pair ATC', 'random'])):
            positions.append(pos) # Define the positions of the boxplots
            pos+=1
        pos+=1 # Add separation between boxplot groups
        plt.axvline(x=pos,linewidth=0.3,linestyle='--',color='black',dashes=(1, 1))
        pos+=2 # Add separation between boxplot groups

        data = []
        #for analysis in ['pair ATC', 'not pair ATC']:
        for analysis in ['pair ATC', 'random']:
            data.append(results_analysis[atc_pair][analysis]['all_aucs']) # Get the data that we will add at the plot

        # Boxplot group
        #bp = boxplot(data, positions = positions, widths = 0.6)
        bp = boxplot(data, positions = positions, widths = 0.6, patch_artist=True)
        #setBoxColors(bp, len(['pair ATC', 'not pair ATC']))
        setBoxColors(bp, len(['pair ATC', 'random']))

        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # Set axes limits and labels
    xlim(0,pos-2)
    ylim(0,1)
    ax.set_xticklabels(axes_labels)
    ax.set_xticks(xticks)
    xlabel('ATC pairs')
    ylabel('Distribution of AUC values')

    # draw temporary color lines and use them to create a legend
    hG, = plot([1,1],'-', color='green')
    hR, = plot([1,1],'-', color='#aeaeae')
    #legend([hG, hR],['pair ATC', 'not pair ATC'])
    lgd= ax.legend(handles=[hG, hR], labels=['Guild', 'Random'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hG.set_visible(False)
    hR.set_visible(False)

    fig.autofmt_xdate()
    savefig(plot_name, format='eps', bbox_extra_artists=(lgd,), bbox_inches='tight')
    show()


    ####### RESULTS TABLE

    fo = open(results_table, 'w')

    # Header
    fo.write(' ')
    for atc in atc_list:
        fo.write('\t{}\t \t '.format(atc))
    fo.write('\n')

    for atc1 in atc_list:
        fo.write('{}'.format(atc1))
        for atc2 in atc_list:
            atc_pair1 = '{}---{}'.format(atc1, atc2)
            atc_pair2 = '{}---{}'.format(atc2, atc1)
            if atc_pair1 in results_analysis:
                atc_pair = atc_pair1
            else:
                atc_pair = atc_pair2
            mean = results_analysis[atc_pair]['pair ATC']['mean']
            num_dc = results_analysis[atc_pair]['pair ATC']['num_dc']

            if mean != '-':
                auc_dist = results_analysis[atc_pair]['pair ATC']['all_aucs']
                random_dist = results_analysis[atc_pair]['random']['all_aucs']
                stat, pval = stats.mannwhitneyu(auc_dist, random_dist)
            else:
                stat = '-'
                pval = '-'

            fo.write('\t{}\t{}\t{}'.format(mean, num_dc, pval))
        fo.write('\n')

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

        clf = classifier.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        y_score = clf.decision_function(X_test)

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


def obtain_drug_combination_scores_mean(all_probs):
    """
    Obtain a dictionary of drug_combination_name => scores_mean
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


def obtain_dcdb2name(cnx):
    """
    Obtain dictionary DCDB_drugID : drug_name
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, N.value FROM externalEntityDCDB_drugID D, externalEntityName N WHERE D.externalEntityID = N.externalEntityID
             ''')

    cursor.execute(query)

    dcdb2name = {}
    for dcdb, name in cursor:
        dcdb2name[dcdb] = name

    cursor.close()

    return dcdb2name

def obtain_dcdb2atc(cnx):
    """
    Obtain dictionary DCDB_drugID : set(ATCs)
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, A.value FROM externalEntityDCDB_drugID D, externalEntityATC A WHERE D.externalEntityID = A.externalEntityID
             ''')

    cursor.execute(query)

    dcdb2atc = {}
    for dcdb, atc in cursor:
        dcdb2atc.setdefault(dcdb, set())
        dcdb2atc[dcdb].add(atc)

    cursor.close()

    return dcdb2atc


def obtain_set_of_first_letter_ATC(ATC_set):
    """ set(['C10AB05','J01XE01']) --> set(['C','J']) """
    one_letter = set()
    for atc in ATC_set:
        one_letter.add(atc[0])
    return one_letter


def setBoxColors(bp, length):
    """
    Set the colors of the box plots groups
    Code from: http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    """
    if length != 2:
        print('\nERROR: This function is thought to work on groups of 2 types of analysis, not {}. Redo the function!!\n'.format(length))
        sys.exit(10)

    setp(bp['boxes'][0], color='green')
    setp(bp['caps'][0], color='green')
    setp(bp['caps'][1], color='green')
    setp(bp['whiskers'][0], color='green')
    setp(bp['whiskers'][1], color='green')
    setp(bp['medians'][0], color='black')
    bp['boxes'][0].set_facecolor('#32f232') # green

    setp(bp['boxes'][1], color='black')
    setp(bp['caps'][2], color='black')
    setp(bp['caps'][3], color='black')
    setp(bp['whiskers'][2], color='black')
    setp(bp['whiskers'][3], color='black')
    setp(bp['medians'][1], color='#fbc562')
    bp['boxes'][1].set_facecolor('#aeaeae') #light grey

    return


if  __name__ == "__main__":
    main()
