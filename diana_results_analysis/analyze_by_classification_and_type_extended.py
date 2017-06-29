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
    plot_name = results_dir + '/figures/classification_and_type_extended.eps'

    # File with results of Mann Whitney tests:
    mannwhitney_file = results_dir + '/tables/classification_and_type_mannwhitney_extended.txt'

    # classifications to analyze
    classifications = ['Different targets in different biological processes', 
                       'Different targets in related biological processes',
                       'Different targets in same biological process',
                       'Same target']
    classifications_short = ['Diff tar diff BP', 
                       'Diff tar rel BP',
                       'Diff tar same BP',
                       'Same tar']

    dump_file = current_dir + "/drug_int_2_drugs.pcl"
    drug_int_2_drugs = cPickle.load(open(dump_file))
    dump_file = current_dir + "/drug_int_2_info.pcl"
    drug_int_2_info = cPickle.load(open(dump_file))

    # Remove always all nan
    remove_all_nan = True

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    types_analysis = ['all', 'guild-node', 'guild-edge', 'guild-function', 'link-node', 'link-edge', 'link-function', 'seed-node', 'seed-pfam', 'seed-function', 'struc', 'random']
    #print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))
    legend_labels = ['all', 'guild', 'linkers', 'seeds', 'struc', 'random']

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
        'SVC best 1' : SVC(gamma=1.0, C=1.0, kernel='rbf'),
        'SVC best 2' : SVC(gamma=0.01, C=100.0, kernel='rbf'),
        'SVC best 3' : SVC(gamma=0.01, C=1000.0, kernel='rbf'),
        'SVC best 4' : SVC(gamma=0.1, C=100.0, kernel='rbf'),
        'DecisionTree' : DecisionTreeClassifier(max_depth=5),
        'RandomForest' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        'MLP' : MLPClassifier(alpha=1),
        'AdaBoost' : AdaBoostClassifier(),
        'GaussianNB' : GaussianNB(),
        'QuadraticDiscr.' : QuadraticDiscriminantAnalysis()
    }

    # Define the output file
    #output_file = results_dir + '/classification_analysis.csv'
    #columns_out = ['Class mean', 'No-class mean', 'Class mean of std', 'No-class mean of std', 'Num DC']
    #output_df = pd.DataFrame(columns=columns_out)


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

    guild_node_info = ['Nit500sp', 'Nit1000dp', 'SNdp', 'SNji']
    guild_edge_info = ['Eit1000sp', 'Eit50dp', 'LEsp', 'LEdp']
    guild_func_info = ['Fit500dp', 'Fper0.1sp', 'LFdp', 'LFsp', 'SFdp', 'SFsp']
    link_node_info = ['LNdp', 'LNsp']
    link_edge_info = ['LEsp', 'LEdp']
    link_func_info = ['LFdp', 'LFsp']
    seed_node_info = ['SNdp', 'SNji']
    seed_pfam_info = ['SPdp', 'SPji']
    seed_func_info = ['SFdp', 'SFsp']

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



    ####################################
    #### ANALYSIS BY CLASSIFICATION ####
    ####################################

    results_analysis = {} # Defining the variable that will store the results

    for classification in classifications:

        print('\nANALIZING CLASSIFICATION: {}\n'.format(classification))

        class_rows = []
        no_class_rows = []

        for index, row in df.iterrows():

            (drug1, drug2) = index.split('---')
            di_bool = False

            for DI in drug_int_2_drugs:
                # If it is drug interaction...
                if drug1 in drug_int_2_drugs[DI] and drug2 in drug_int_2_drugs[DI]:
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


        df_class = df.ix[class_rows] # Create a table with the rows selected (DDIs of the class and non-DDIs)
        df_noclass = df.ix[no_class_rows] # Create a table with the DDIs not of the class and non-DDIs 

        same_class_num_dc = None # Defining the variable

        for (df_analysis, analysis) in [(df_class, 'in classification'), (df_noclass, 'not in classification')]:

            print('\n{}\n'.format(analysis))

            # Get the number of DC for the given classification
            dc_data = df_analysis[df_analysis['Comb'] == 1]
            num_dc = len(dc_data.index)
            print('Number of DC for {} {}: {}'.format(analysis, classification, num_dc))

            # Option of removing all the nan values, no matter the type of analysis
            if remove_all_nan == True:
                df_analysis = df_analysis.dropna()


            ######## ANALYSIS BY TYPE OF PROGRAM
            for type_analysis in types_analysis:

                print('\nType of analysis: {}\n'.format(type_analysis))

                # Getting the data frame of the analysis
                if type_analysis == 'all':
                    df_all = df_analysis[guild_node_info+guild_edge_info+guild_func_info+link_node_info+link_edge_info+link_func_info+seed_node_info+seed_pfam_info+seed_func_info+['Struc']+['Comb']]
                    df_ml = df_all
                elif type_analysis == 'guild-node':
                    df_node = df_analysis[guild_node_info+['Comb']]
                    df_ml = df_node
                elif type_analysis == 'guild-edge':
                    df_edge = df_analysis[guild_edge_info+['Comb']]
                    df_ml = df_edge
                elif type_analysis == 'guild-function':
                    df_func = df_analysis[guild_func_info+['Comb']]
                    df_ml = df_func
                elif type_analysis == 'link-node':
                    df_node = df_analysis[link_node_info+['Comb']]
                    df_ml = df_node
                elif type_analysis == 'link-edge':
                    df_edge = df_analysis[link_edge_info+['Comb']]
                    df_ml = df_edge
                elif type_analysis == 'link-function':
                    df_func = df_analysis[link_func_info+['Comb']]
                    df_ml = df_func
                elif type_analysis == 'seed-node':
                    df_node = df_analysis[seed_node_info+['Comb']]
                    df_ml = df_node
                elif type_analysis == 'seed-edge':
                    df_edge = df_analysis[seed_edge_info+['Comb']]
                    df_ml = df_edge
                elif type_analysis == 'seed-function':
                    df_func = df_analysis[seed_func_info+['Comb']]
                    df_ml = df_func
                elif type_analysis == 'struc':
                    df_struc = df_analysis[['Struc']+['Comb']]
                    df_ml = df_struc
                elif type_analysis == 'random':
                    df_ml = df_analysis[['Comb']]

                # Remove nan values
                df_ml = df_ml.dropna()

                # Get the drug combination data and the number of drug combinations
                dc_data = df_ml[df_ml['Comb'] == 1]
                num_dc = len(dc_data.index)
                print('Number of DC after removing nan: {}'.format(num_dc))

                if analysis == 'in classification':
                    same_class_num_dc = num_dc # Get the number of DC when classification

                if analysis == 'not in classification':
                    print('Number of DC used: {}' .format(same_class_num_dc))

                # Stop if the number of drug interactions is smaller than the number of cross-validations!!
                if num_dc < n_fold:
                    print('Not possible to do the analysis for classification {}. The number of positive samples is {} and the n-fold is {}\n'.format(classification, num_dc, n_fold))
                    results_analysis.setdefault(classification, {})
                    results_analysis.setdefault[classification](type_analysis, {})
                    results_analysis.setdefault[classification][type_analysis](analysis, {})
                    results_analysis[classification][type_analysis][analysis] = {'mean':'-','std':'-','num_dc':int(same_class_num_dc),'all_aucs':'-'}
                    #results_analysis = {'in classification' : {'mean':'-','std':'-','num_dc':int(same_class_num_dc)}, 'not in classification' : {'mean':'-','std':'-','num_dc':int(same_class_num_dc)}} # Defining the variable
                    continue
                # Stop if the number of drug interactions is smaller than the minimum number given!!
                if num_dc < min_num_dc_group:
                    print('Not possible to do the analysis for classification {}. The number of positive samples is {} and the minimum number per group is {}\n'.format(classification, num_dc, min_num_dc_group))
                    results_analysis.setdefault(classification, {})
                    results_analysis[classification].setdefault(type_analysis, {})
                    results_analysis[classification][type_analysis].setdefault(analysis, {})
                    results_analysis[classification][type_analysis][analysis] = {'mean':'-','std':'-','num_dc':int(same_class_num_dc),'all_aucs':'-'}
                    #results_analysis = {'in classification' : {'mean':'-','std':'-','num_dc':int(same_class_num_dc)}, 'not in classification' : {'mean':'-','std':'-','num_dc':int(same_class_num_dc)}} # Defining the variable
                    continue

                # Get the non-drug combination data
                ndc_data = df_ml[df_ml['Comb'] == 0]

                # Obtain the different non-drug combination groups to repeat the analysis
                # We will use the number of drug combinations when classification in order to have same number of samples 
                print('Building {} repetition groups of {} (same) DC and {} (different) non-DC'.format(repetitions,same_class_num_dc,same_class_num_dc))
                ndc_repetitions = obtain_n_groups_of_k_length(ndc_data, repetitions, same_class_num_dc) # Obtain n number of groups containing different non-drug combinations to repeat the analysis n times
                #print(ndc_repetitions)

                mean_aucs = [] # Here we will store the means of AUCs from the cross-validations
                std_aucs = [] # Here we will store the standard deviations of the AUCs from the cross-validations
                all_aucs = [] # Here we will store ALL the AUCs

                for ndc_data_equal in ndc_repetitions:

                    num_items_group = int( float(same_class_num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation

                    print('Building {} fold groups of {} DC and {} non-DC'.format(n_fold,num_items_group,num_items_group))
                    dc_groups = obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group) # Defining the drug combination groups in each cross-validation step
                    ndc_groups = obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group) # Defining the non-drug combination groups in each cross-validation step
                    #print dc_groups
                    #print ndc_groups
                    merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

                    if type_analysis == 'random':
                        mean, var, std, list_auc = run_nfold_crossvalidation_random(n_fold, merged_groups, classifiers[classifier])
                    else:
                        mean, var, std, list_auc = run_nfold_crossvalidation_scikit(n_fold, merged_groups, classifiers[classifier])
                    mean_aucs.append(mean)
                    std_aucs.append(std)
                    all_aucs = all_aucs + list_auc
                    #run_nfold_crossvalidation_testing_classifiers(n_fold, merged_groups)

                final_mean = np.mean(mean_aucs)
                mean_std = np.mean(std_aucs)
                std_means = np.std(mean_aucs)
                results_analysis.setdefault(classification, {})
                results_analysis[classification].setdefault(type_analysis, {})
                results_analysis[classification][type_analysis].setdefault(analysis, {})
                results_analysis[classification][type_analysis][analysis]['mean'] = final_mean
                results_analysis[classification][type_analysis][analysis]['std'] = mean_std
                results_analysis[classification][type_analysis][analysis]['num_dc'] = int(same_class_num_dc)
                results_analysis[classification][type_analysis][analysis]['all_aucs'] = all_aucs
                print('FINAL MEAN: {}'.format(final_mean))
                print('MEAN of STD: {}'.format(mean_std))
                print('STD of MEANS: {}'.format(std_means))

        #row = [ results_analysis['in classification']['mean'], results_analysis['not in classification']['mean'], results_analysis['in classification']['std'], results_analysis['not in classification']['std'], results_analysis['in classification']['num_dc'] ]
        #odf = pd.DataFrame([row], columns=columns_out, index=[classification])
        #output_df = output_df.append(odf)

    #print(output_df)
    #output_df.to_csv(output_file)

    ####### PLOT GRAPHIC OF DISTRIBUTION OF AUCS VS. CLASSIFICATION AND METHOD
    print(results_analysis)

    fig = figure(dpi=300)
    ax = axes()
    #hold(True)
    pos = 1

    xticks = [] # Define the places in which the labels will be
    for classification in classifications:

        positions = []
        for x in xrange(len(types_analysis)):
            positions.append(pos) # Define the positions of the boxplots
            pos+=1
        pos+=1 # Add separation between boxplot groups

        data = []
        for type_analysis in types_analysis:
            data.append(results_analysis[classification][type_analysis]['in classification']['all_aucs']) # Get the groups of plots that we will add

        # Boxplot group
        #bp = boxplot(data, positions = positions, widths = 0.6)
        bp = boxplot(data, positions = positions, widths = 0.6, patch_artist=True)
        setBoxColors(bp, len(types_analysis))

        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # Set axes limits and labels
    xlim(0,pos-1)
    ylim(0,1)
    axes_labels = classifications_short
    ax.set_xticklabels(axes_labels)
    ax.set_xticks(xticks)
    xlabel('Classification')
    ylabel('Distribution of AUC values')

    # draw temporary color lines and use them to create a legend
    hB, = plot([1,1],'-', color='blue')
    hG, = plot([1,1],'-', color='green')
    hY, = plot([1,1],'-', color='orange')
    hR, = plot([1,1],'-', color='red')
    hBl, = plot([1,1],'-', color='black')
    hW, = plot([1,1],'-', color='#e5e5e5')
    legend([hB, hG, hY, hR, hBl, hW],legend_labels)
    hB.set_visible(False)
    hG.set_visible(False)
    hY.set_visible(False)
    hR.set_visible(False)
    hBl.set_visible(False)
    hW.set_visible(False)

    #fig.autofmt_xdate()
    savefig(plot_name, format='eps')
    show()

    sys.exit(0)
    ################ MANN WHITNEY U
    fo = open(mannwhitney_file, 'w')

    all_diff = results_analysis['Different targets in different biological processes']['all']['in classification']['all_aucs']
    guild_diff = results_analysis['Different targets in different biological processes']['guild']['in classification']['all_aucs']
    link_diff = results_analysis['Different targets in different biological processes']['linkers']['in classification']['all_aucs']
    seed_diff = results_analysis['Different targets in different biological processes']['seeds']['in classification']['all_aucs']

    all_seed = stats.mannwhitneyu(all_diff, seed_diff)
    guild_seed = stats.mannwhitneyu(guild_diff, seed_diff)
    link_seed = stats.mannwhitneyu(link_diff, seed_diff)

    print('\nRESULT OF ALL VS SEEDS COMPARISON FOR DIFFERENT PATHWAYS: {}\n'.format(all_seed))
    print('RESULT OF GUILD VS SEEDS COMPARISON FOR DIFFERENT PATHWAYS: {}\n'.format(guild_seed))
    print('RESULT OF LINKERS VS SEEDS COMPARISON FOR DIFFERENT PATHWAYS: {}\n'.format(link_seed))
    fo.write('# Different targets in different biological processes\n')
    fo.write('All vs. Seeds:\t{}\n'.format(all_seed))
    fo.write('GUILD vs. Seeds:\t{}\n'.format(guild_seed))
    fo.write('Linkers vs. Seeds:\t{}\n'.format(link_seed))

    all_rel = results_analysis['Different targets in related biological processes']['all']['in classification']['all_aucs']
    guild_rel = results_analysis['Different targets in related biological processes']['guild']['in classification']['all_aucs']
    link_rel = results_analysis['Different targets in related biological processes']['linkers']['in classification']['all_aucs']
    seed_rel = results_analysis['Different targets in related biological processes']['seeds']['in classification']['all_aucs']

    all_seed_rel = stats.mannwhitneyu(all_rel, seed_rel)
    guild_seed_rel = stats.mannwhitneyu(guild_rel, seed_rel)
    link_seed_rel = stats.mannwhitneyu(link_rel, seed_rel)

    print('RESULT OF ALL VS SEEDS COMPARISON FOR DIFFERENT PATHWAYS: {}\n'.format(all_seed_rel))
    print('RESULT OF GUILD VS SEEDS COMPARISON FOR DIFFERENT PATHWAYS: {}\n'.format(guild_seed_rel))
    print('RESULT OF LINKERS VS SEEDS COMPARISON FOR DIFFERENT PATHWAYS: {}\n'.format(link_seed_rel))
    fo.write('\n# Different targets in related biological processes\n')
    fo.write('All vs. Seeds:\t{}\n'.format(all_seed_rel))
    fo.write('GUILD vs. Seeds:\t{}\n'.format(guild_seed_rel))
    fo.write('Linkers vs. Seeds:\t{}\n'.format(link_seed_rel))

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
    if length != 12:
        print('\nERROR: This function is thought to work on groups of 12 types of analysis. Redo the function!!\n')
        sys.exit(10)

    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['medians'][0], color='black')
    bp['boxes'][0].set_facecolor('#9ed0ff')

    setp(bp['boxes'][1], color='green')
    setp(bp['caps'][2], color='green')
    setp(bp['caps'][3], color='green')
    setp(bp['whiskers'][2], color='green')
    setp(bp['whiskers'][3], color='green')
    setp(bp['medians'][1], color='black')
    bp['boxes'][1].set_facecolor('#32f232')

    setp(bp['boxes'][2], color='green')
    setp(bp['caps'][4], color='green')
    setp(bp['caps'][5], color='green')
    setp(bp['whiskers'][4], color='green')
    setp(bp['whiskers'][5], color='green')
    setp(bp['medians'][2], color='black')
    bp['boxes'][2].set_facecolor('#32f232')

    setp(bp['boxes'][3], color='green')
    setp(bp['caps'][6], color='green')
    setp(bp['caps'][7], color='green')
    setp(bp['whiskers'][6], color='green')
    setp(bp['whiskers'][7], color='green')
    setp(bp['medians'][3], color='black')
    bp['boxes'][3].set_facecolor('#32f232')

    setp(bp['boxes'][4], color='#d48900')
    setp(bp['caps'][8], color='#d48900')
    setp(bp['caps'][9], color='#d48900')
    setp(bp['whiskers'][8], color='#d48900')
    setp(bp['whiskers'][9], color='#d48900')
    setp(bp['medians'][4], color='black')
    bp['boxes'][4].set_facecolor('#fbc562') #orange

    setp(bp['boxes'][5], color='#d48900')
    setp(bp['caps'][10], color='#d48900')
    setp(bp['caps'][11], color='#d48900')
    setp(bp['whiskers'][10], color='#d48900')
    setp(bp['whiskers'][11], color='#d48900')
    setp(bp['medians'][5], color='black')
    bp['boxes'][5].set_facecolor('#fbc562') #orange

    setp(bp['boxes'][6], color='#d48900')
    setp(bp['caps'][12], color='#d48900')
    setp(bp['caps'][13], color='#d48900')
    setp(bp['whiskers'][12], color='#d48900')
    setp(bp['whiskers'][13], color='#d48900')
    setp(bp['medians'][6], color='black')
    bp['boxes'][6].set_facecolor('#fbc562') #orange

    setp(bp['boxes'][7], color='#b80000')
    setp(bp['caps'][14], color='#b80000')
    setp(bp['caps'][15], color='#b80000')
    setp(bp['whiskers'][14], color='#b80000')
    setp(bp['whiskers'][15], color='#b80000')
    setp(bp['medians'][7], color='black')
    bp['boxes'][7].set_facecolor('#ff7373') #red

    setp(bp['boxes'][8], color='#b80000')
    setp(bp['caps'][16], color='#b80000')
    setp(bp['caps'][17], color='#b80000')
    setp(bp['whiskers'][16], color='#b80000')
    setp(bp['whiskers'][17], color='#b80000')
    setp(bp['medians'][8], color='black')
    bp['boxes'][8].set_facecolor('#ff7373') #red

    setp(bp['boxes'][9], color='#b80000')
    setp(bp['caps'][18], color='#b80000')
    setp(bp['caps'][19], color='#b80000')
    setp(bp['whiskers'][18], color='#b80000')
    setp(bp['whiskers'][19], color='#b80000')
    setp(bp['medians'][9], color='black')
    bp['boxes'][9].set_facecolor('#ff7373') #red

    setp(bp['boxes'][10], color='black')
    setp(bp['caps'][20], color='black')
    setp(bp['caps'][21], color='black')
    setp(bp['whiskers'][20], color='black')
    setp(bp['whiskers'][21], color='black')
    setp(bp['medians'][10], color='#fbc562')
    bp['boxes'][10].set_facecolor('grey')

    setp(bp['boxes'][11], color='black')
    setp(bp['caps'][22], color='black')
    setp(bp['caps'][23], color='black')
    setp(bp['whiskers'][22], color='black')
    setp(bp['whiskers'][23], color='black')
    setp(bp['medians'][11], color='grey')
    bp['boxes'][11].set_facecolor('#e5e5e5')


    return


if  __name__ == "__main__":
    main()
