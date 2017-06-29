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
    plot_name = results_dir + '/figures/structural_similarity_in_drug_combinations.eps'

    me_too_drugs_file = current_dir + '/me_too_drug_pairs.pcl'
    me_too_drug_comb_file = current_dir + '/me_too_drug_comb_pairs.pcl'
    me_too_drugs_table = results_dir + '/tables/me_too_drug_pairs.txt'
    me_too_drug_combs_table = results_dir + '/tables/me_too_drug_comb_pairs.txt'

    cnx = mysql.connector.connect(user='quim', password='',
                                  host='localhost',
                                  database='BIANA_JAN_2017')

    # # Plot of accuracy/sensitivity name
    # acc_sens_all = results_dir + '/figures/classification_and_method_accsens_all.eps'
    # acc_sens_guild = results_dir + '/figures/classification_and_method_accsens_guild.eps'
    # acc_sens_linkers = results_dir + '/figures/classification_and_method_accsens_linkers.eps'
    # acc_sens_seeds = results_dir + '/figures/classification_and_method_accsens_seeds.eps'
    # acc_sens_struc = results_dir + '/figures/classification_and_method_accsens_struc.eps'
    # acc_sens_random = results_dir + '/figures/classification_and_method_accsens_random.eps'
    # #acc_sens_name_method = results_dir + '/figures/classification_and_method_accsens_method.eps'

    # # File with results of Mann Whitney tests:
    # mannwhitney_file = results_dir + '/tables/classification_and_method_mannwhitney.txt'
    # mannwhitney_struc_file = results_dir + '/tables/classification_and_method_mannwhitneystruc.txt'

    # # Results table:
    # results_table = results_dir + '/tables/classification_and_method_table.txt'

    # # Results table:
    # prec_rec_table = results_dir + '/tables/classification_and_method_precrec_table.txt'

    # # classifications to analyze
    # classifications = ['Different targets in different biological processes', 
    #                    'Different targets in related biological processes',
    #                    'Different targets in same biological process',
    #                    'Same target']
    # classifications_short = ['Diff tar diff BP', 
    #                    'Diff tar rel BP',
    #                    'Diff tar same BP',
    #                    'Same tar']

    dump_file = current_dir + "/drug_int_2_drugs.pcl"
    drug_int_2_drugs = cPickle.load(open(dump_file))
    dump_file = current_dir + "/drug_int_2_info.pcl"
    drug_int_2_info = cPickle.load(open(dump_file))

    # Remove always all nan
    remove_all_nan = True

    # Type of machine learning analysis
    # (all, main, guild, seeds, linkers, struc)
    types_analysis = ['all', 'guild', 'linkers', 'seeds', 'struc', 'random']
    #print('\nTYPE OF ANALYSIS: {}\n'.format(type_analysis))

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


    ###########################################
    #### ANALYSIS OF STRUCTURAL SIMILARITY ####
    ###########################################

    df_struc = df[['Struc']+['Comb']]
    df_struc = df_struc.dropna()
    df_struc_all = df_struc[['Struc']]
    df_struc_all=df_struc_all.astype(float)

    dc_data = df_struc[df_struc['Comb'] == 1]
    df_struc_dc = dc_data[['Struc']]
    df_struc_dc=df_struc_dc.astype(float)
    
    ndc_data = df_struc[df_struc['Comb'] == 0]
    df_struc_ndc = ndc_data[['Struc']]
    df_struc_ndc=df_struc_ndc.astype(float)

    # Plot of me-too drugs

    df_struc_dc.plot.hist(alpha=0.5,rwidth = 0.8,bins=10)
    xlabel('Structural similarity score')
    ylabel('Frequency of drug combinations')
    ax = axes()
    ax.legend_.remove()
    savefig(plot_name, format='eps', dpi=300)
    show()
    plt.clf()

    num_metoo_dc = 0
    num_nonmetoo_dc = 0
    me_too_drug_pairs = set()
    me_too_drug_comb_pairs = set()

    done_pairings = []

    me_too_drugs_dict = {}
    me_too_drug_combs_dict = {}


    for dc_index1, dc_row1 in df_struc_dc.iterrows():

        score = dc_row1['Struc']
        if score >= 0.7:
            me_too_drug_pairs.add(dc_index1)
            me_too_drugs_dict[dc_index1] = score

        for dc_index2, dc_row2 in df_struc_dc.iterrows():

            if dc_index1 == dc_index2:
                continue
            combpair1 = '___'.join([dc_index1, dc_index2])
            combpair2 = '___'.join([dc_index2, dc_index1])
            if combpair1 in done_pairings or combpair2 in done_pairings:
                continue
            done_pairings.append(combpair1)
            done_pairings.append(combpair2)

            (drug11, drug12) = dc_index1.split('---')
            (drug21, drug22) = dc_index2.split('---')

            pairing11_1 = '---'.join([drug11, drug21])
            pairing11_2 = '---'.join([drug21, drug11])
            pairing12_1 = '---'.join([drug12, drug22])
            pairing12_2 = '---'.join([drug22, drug12])

            pairing21_1 = '---'.join([drug11, drug22])
            pairing21_2 = '---'.join([drug22, drug11])
            pairing22_1 = '---'.join([drug12, drug21])
            pairing22_2 = '---'.join([drug21, drug12])

            group1 = []
            no_pairing = False
            for possib1, possib2 in [ (pairing11_1,pairing11_2), (pairing12_1, pairing12_2) ]:

                if possib1 in df_struc_all.index:
                    pairing = df_struc_all.loc[[possib1]]
                    group1.append(pairing)
                elif possib2 in df_struc_all.index:
                    pairing = df_struc_all.loc[[possib2]]
                    group1.append(pairing)
                else:
                    #print('No pairing found!')
                    num_nonmetoo_dc+=1
                    no_pairing = True
            group2 = []
            for possib1, possib2 in [ (pairing21_1,pairing21_2), (pairing22_1, pairing22_2) ]:

                if possib1 in df_struc_all.index:
                    pairing = df_struc_all.loc[[possib1]]
                    group2.append(pairing)
                elif possib2 in df_struc_all.index:
                    pairing = df_struc_all.loc[[possib2]]
                    group2.append(pairing)
                else:
                    #print('No pairing found!')
                    if no_pairing == False:
                        num_nonmetoo_dc+=1
                    no_pairing = True

            if no_pairing:
                continue

            score11 = group1[0].iloc[0]['Struc']
            score12 = group1[1].iloc[0]['Struc']

            score21 = group2[0].iloc[0]['Struc']
            score22 = group2[1].iloc[0]['Struc']

            #print(score11, score12, score21, score22)

            if (score11 < 0.7 and score12 < 0.7) or (score21 < 0.7 and score22 < 0.7):
                num_nonmetoo_dc+=1
            else:
                num_metoo_dc+=1
                me_too_drug_comb_pairs.add(combpair1)

                if (score11 >= 0.7 and score12 >= 0.7):
                    me_too_drug_combs_dict.setdefault(combpair1, {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_1', {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_2', {})
                    me_too_drug_combs_dict[combpair1]['me_too_1'][group1[0].index[0]] = score11
                    me_too_drug_combs_dict[combpair1]['me_too_2'][group1[1].index[0]] = score12
                elif (score21 < 0.7 and score22 < 0.7):
                    me_too_drug_combs_dict.setdefault(combpair2, {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_1', {})
                    me_too_drug_combs_dict[combpair1].setdefault('me_too_2', {})
                    me_too_drug_combs_dict[combpair1]['me_too_1'][group2[0].index[0]] = score21
                    me_too_drug_combs_dict[combpair1]['me_too_2'][group2[1].index[0]] = score22

    print('Number of me-too drug combinations: {}'.format(num_metoo_dc))
    print('Number of non me-too drug combinations: {}'.format(num_nonmetoo_dc))

    dcdb2name = obtain_dcdb2name(cnx)

    me_too_drugs_fd = open(me_too_drugs_table, 'w')
    me_too_drug_comb_fd = open(me_too_drug_combs_table, 'w')

    # Write the results of me-too drug pairs
    for drug_pair, score in sorted(me_too_drugs_dict.iteritems(), key=lambda (x, y): y, reverse=True):
        (drug1, drug2) = drug_pair.split('---')
        name1 = dcdb2name[drug1]
        name2 = dcdb2name[drug2]
        me_too_drugs_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(drug1, name1, drug2, name2, score))

    # Write the results of me-too drug combination pairs
    for drug_comb_pair in me_too_drug_combs_dict:
        (dc1, dc2) = drug_comb_pair.split('___')
        (drug1, drug2) = dc1.split('---')
        name1 = dcdb2name[drug1]
        name2 = dcdb2name[drug2]
        (drug3, drug4) = dc2.split('---')
        name3 = dcdb2name[drug3]
        name4 = dcdb2name[drug4]

        me_too_1 = me_too_drug_combs_dict[drug_comb_pair]['me_too_1'].keys()[0]
        score1 = me_too_drug_combs_dict[drug_comb_pair]['me_too_1'][me_too_1]
        (mtd1, mtd2) = me_too_1.split('---') # mtd = me too drug
        mtn1 = dcdb2name[mtd1] # mtn = me too drug name
        mtn2 = dcdb2name[mtd2]

        me_too_2 = me_too_drug_combs_dict[drug_comb_pair]['me_too_2'].keys()[0]
        score2 = me_too_drug_combs_dict[drug_comb_pair]['me_too_2'][me_too_2]
        (mtd3, mtd4) = me_too_2.split('---')
        mtn3 = dcdb2name[mtd3]
        mtn4 = dcdb2name[mtd4]

        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name1,name2,name3,name4,mtn1,mtn2,score1,mtn3,mtn4,score2))
        me_too_drug_comb_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name1,name2,name3,name4,mtn1,mtn2,score1,mtn3,mtn4,score2))

    me_too_drugs_fd.close()
    me_too_drug_comb_fd.close()

    cPickle.dump(me_too_drug_comb_pairs, open(me_too_drug_comb_file, 'w')) 
    cPickle.dump(me_too_drug_pairs, open(me_too_drugs_file, 'w')) 

    sys.exit(0)


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

        for (df_analysis, analysis) in [(df_class, 'in classification')]:

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
                elif type_analysis == 'random':
                    #df_ml = df_analysis[['Comb']]
                    df_ml = df_guild

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
                all_probs = [] # Here we store all the probabilities and labels

                for ndc_data_equal in ndc_repetitions:

                    num_items_group = int( float(same_class_num_dc) / float(n_fold) ) # Calculate the number of items in each group of the cross-validation

                    print('Building {} fold groups of {} DC and {} non-DC'.format(n_fold,num_items_group,num_items_group))
                    dc_groups = obtain_n_groups_of_k_length(dc_data, n_fold, num_items_group) # Defining the drug combination groups in each cross-validation step
                    ndc_groups = obtain_n_groups_of_k_length(ndc_data_equal, n_fold, num_items_group) # Defining the non-drug combination groups in each cross-validation step
                    #print dc_groups
                    #print ndc_groups
                    merged_groups = [pd.concat([x,y]) for x,y in zip(dc_groups, ndc_groups)]

                    if type_analysis == 'random':
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
                results_analysis.setdefault(classification, {})
                results_analysis[classification].setdefault(type_analysis, {})
                results_analysis[classification][type_analysis].setdefault(analysis, {})
                results_analysis[classification][type_analysis][analysis]['mean'] = final_mean
                results_analysis[classification][type_analysis][analysis]['std'] = std
                results_analysis[classification][type_analysis][analysis]['num_dc'] = int(same_class_num_dc)
                results_analysis[classification][type_analysis][analysis]['all_aucs'] = all_aucs
                results_analysis[classification][type_analysis][analysis]['all_probs'] = all_probs
                print('FINAL MEAN: {}'.format(final_mean))
                print('STD: {}'.format(std))


        #row = [ results_analysis['in classification']['mean'], results_analysis['not in classification']['mean'], results_analysis['in classification']['std'], results_analysis['not in classification']['std'], results_analysis['in classification']['num_dc'] ]
        #odf = pd.DataFrame([row], columns=columns_out, index=[classification])
        #output_df = output_df.append(odf)

    #print(output_df)
    #output_df.to_csv(output_file)

    results_analysis = plot_precision_sensitivity(results_analysis, 'all', classifications, acc_sens_all)
    results_analysis = plot_precision_sensitivity(results_analysis, 'guild', classifications, acc_sens_guild)
    results_analysis = plot_precision_sensitivity(results_analysis, 'linkers', classifications, acc_sens_linkers)
    results_analysis = plot_precision_sensitivity(results_analysis, 'seeds', classifications, acc_sens_seeds)
    results_analysis = plot_precision_sensitivity(results_analysis, 'struc', classifications, acc_sens_struc)
    #results_analysis = plot_precision_sensitivity(results_analysis, 'random', classifications, acc_sens_random)
    #_ = plot_precision_sensitivity_by_method(results_analysis, 'Different targets in different biological processes', types_analysis, acc_sens_name_method)

    ####### PLOT GRAPHIC OF DISTRIBUTION OF AUCS VS. CLASSIFICATION AND METHOD
    #print(results_analysis)

    fig = figure(dpi=300)
    ax = axes()
    #hold(True)
    pos = 2

    xticks = [] # Define the places in which the labels will be
    for classification in classifications:

        positions = []
        for x in xrange(len(types_analysis)):
            positions.append(pos) # Define the positions of the boxplots
            pos+=1
        pos+=1 # Add separation between boxplot groups
        plt.axvline(x=pos,linewidth=0.3,linestyle='--',color='black',dashes=(1, 1))
        pos+=2 # Add separation between boxplot groups

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
    xlim(0,pos-2)
    ylim(0,1)
    axes_labels = classifications_short
    ax.set_xticklabels(axes_labels)
    ax.set_xticks(xticks)
    #xlabel('Classification')
    ylabel('Distribution of AUC values')

    # draw temporary color lines and use them to create a legend
    hB, = plot([1,1],'-', color='blue')
    hG, = plot([1,1],'-', color='green')
    hY, = plot([1,1],'-', color='orange')
    hR, = plot([1,1],'-', color='red')
    hBl, = plot([1,1],'-', color='black')
    hW, = plot([1,1],'-', color='#aeaeae')
    lgd = ax.legend(handles=(hB, hG, hY, hR, hBl, hW), labels=types_analysis, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hB.set_visible(False)
    hG.set_visible(False)
    hY.set_visible(False)
    hR.set_visible(False)
    hBl.set_visible(False)
    hW.set_visible(False)

    savefig(plot_name, format='eps', bbox_extra_artists=(lgd,), bbox_inches='tight')
    show()


    ################ RESULTS TABLE
    fr = open(results_table, 'w')

    # Header
    fr.write(' ')
    for method in types_analysis:
        fr.write('\t{}\t \t '.format(method))
    fr.write('\n')

    for classification in classifications:
        fr.write('{}'.format(classification))
        for method in types_analysis:
            mean = results_analysis[classification][method]['in classification']['mean']
            std = results_analysis[classification][method]['in classification']['std']
            num_dc = results_analysis[classification][method]['in classification']['num_dc']
            fr.write('\t{}\t{}\t{}'.format(mean, std, num_dc))
        fr.write('\n')
    fr.close()

    ################ PRECISION/RECALL TABLE
    fp = open(prec_rec_table, 'w')

    # Header
    fp.write(' ')
    types_analysis2 = ['all', 'guild', 'linkers', 'seeds', 'struc'] # Without random!
    for method in types_analysis2:
        fp.write('\t{}\t '.format(method))
    fp.write('\n')

    for classification in classifications:
        fp.write('{}'.format(classification))
        for method in types_analysis2:
            cut_off = results_analysis[classification][method]['in classification']['cut_off']
            value = results_analysis[classification][method]['in classification']['value']
            fp.write('\t{}\t{}'.format(cut_off, value))
        fp.write('\n')
    fp.close()

    ################ MANN WHITNEY U
    fo = open(mannwhitney_file, 'w')

    mann_results = {}

    # Header
    fo.write(' \t ')
    for method in types_analysis:
        fo.write('\t{}'.format(method))
    fo.write('\n')

    # Perform the comparisons
    for classification in classifications:
        mann_results.setdefault(classification, {})
        for method1 in types_analysis:
            mann_results[classification].setdefault(method1, {})
            for method2 in types_analysis:
                if method1 == method2:
                    mann_results[classification][method1][method2] = '-'
                else:
                    method1_dist = results_analysis[classification][method1]['in classification']['all_aucs']
                    method2_dist = results_analysis[classification][method2]['in classification']['all_aucs']
                    stat, pval = stats.mannwhitneyu(method1_dist, method2_dist)
                    mann_results[classification][method1][method2] = [stat, pval]

    # Write the table of crossings
    for classification in classifications:
        for method1 in types_analysis:
            fo.write('{}\t{}'.format(classification, method1))
            for method2 in types_analysis:
                if method1 == method2:
                    fo.write('\t-')
                else:
                    stat, pval = mann_results[classification][method1][method2]
                    fo.write('\t{}, {:.2e}'.format(stat,pval))
            fo.write('\n')

    fo.close()


    ################ MANN WHITNEY U BETWEEN STRUCTURAL VALUES
    fs = open(mannwhitney_struc_file, 'w')
    fs.write(' ')

    mann_struc_results = {}

    for classification1 in classifications:
        mann_struc_results.setdefault(classification1, {})
        fs.write('\t{}'.format(classification1))
        for classification2 in classifications:
            if classification1 == classification2:
                mann_struc_results[classification1][classification2] = '-'
            else:
                class1_dist = results_analysis[classification1]['struc']['in classification']['all_aucs']
                class2_dist = results_analysis[classification2]['struc']['in classification']['all_aucs']
                stat, pval = stats.mannwhitneyu(class1_dist, class2_dist)
                mann_struc_results[classification1][classification2] = [stat, pval]
    fs.write('\n')

    # Write the table of crossings
    for classification1 in classifications:
        fs.write('{}'.format(classification1))
        for classification2 in classifications:
            if classification1 == classification2:
                fs.write('\t-')
            else:
                stat, pval = mann_struc_results[classification1][classification2]
                fs.write('\t{}, {:.2e}'.format(stat,pval))
        fs.write('\n')

    fs.close()


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


def plot_precision_sensitivity(results_analysis, type_analysis, array_ordered, name_plot):

    cut_offs = frange(0,1,0.005)
    fig, ax = plt.subplots()
    prec_colors = ['#a6d0ff','#4199fd','#005bc1','#003169']
    #prec_colors = ['#a6d0ff','#003169']
    sens_colors = ['#a6ff9e','#15ff00','#0d9f00','#085c00']
    #sens_colors = ['#a6ff9e','#085c00']
    point_colors = ['#ffb1b1','#ff3737','#d40303','#890000']
    #point_colors = ['#ffb1b1','#890000']
    c = 0

    for classification in array_ordered:

        precision = []
        sensitivity = []

        all_probs = results_analysis[classification][type_analysis]['in classification']['all_probs']
        for cut_off in cut_offs:

            tp,fp,tn,fn = calculate_tp_fp_tn_fn(all_probs, cut_off)

            try:
                prec_val = float(tp)/(float(tp)+float(fp))
            except:
                prec_val = 1
            sens_val = float(tp)/(float(tp)+float(fn))
            #print('PRECISION: {}\tSENSITIVITY: {}\n'.format(prec_val, sens_val))

            precision.append(prec_val)
            sensitivity.append(sens_val)
        
        ax.plot( cut_offs, precision, '-', color=prec_colors[c])
        ax.plot( cut_offs, sensitivity, '-', color=sens_colors[c])

        # Find the intersection point
        idx = np.argwhere(np.diff(np.sign(np.array(precision) - np.array(sensitivity))) != 0).reshape(-1) + 0

        # Plot the intersection point
        ax.plot(cut_offs[idx[0]], precision[idx[0]], 'o', color=point_colors[c])

        results_analysis[classification][type_analysis]['in classification']['cut_off'] = cut_offs[idx[0]]
        results_analysis[classification][type_analysis]['in classification']['value'] = precision[idx[0]]

        print('Classification: {}\tValue: {}\tCut-off: {}'.format(classification, precision[idx[0]], cut_offs[idx[0]]))

        # Plot vertical line from point
        #ax.plot( [cut_offs[idx[0]], cut_offs[idx[0]]], [precision[idx[0]], 0], 'r--')
        #plt.annotate(' {}'.format(cut_offs[idx[0]]), xy=(cut_offs[idx[0]], 0))
        # Plot horizontal line from point
        #ax.plot( [cut_offs[idx[0]], 0], [precision[idx[0]], precision[idx[0]]], 'r--')
        #plt.annotate('{:.3f}'.format(precision[idx[0]]), xy=(0, precision[idx[0]]+0.02))

        c+=1

    xlabel('Probability thresholds')
    ylabel('Precision-Recall')
    
    # draw temporary color lines and use them to create a legend
    hB, = plot([1,1],'o', color='#ffb1b1')
    hG, = plot([1,1],'o', color='#ff3737')
    hY, = plot([1,1],'o', color='#d40303')
    hR, = plot([1,1],'o', color='#890000')
    hPr, = plot([1,1],'-', color='#0078ff')
    hRe, = plot([1,1],'-', color='#0d9f00')
    classifications_short = ['Diff tar diff BP', 
                       'Diff tar rel BP',
                       'Diff tar same BP',
                       'Same tar']
    classifications_short = classifications_short + ['Precision', 'Recall']
    lgd = ax.legend(handles=(hB, hG, hY, hR, hPr, hRe), labels=classifications_short, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hB.set_visible(False)
    hG.set_visible(False)
    hY.set_visible(False)
    hR.set_visible(False)
    hPr.set_visible(False)
    hRe.set_visible(False)

    savefig(name_plot, format='eps', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    show()
    
    return results_analysis

def plot_precision_sensitivity_by_method(results_analysis, classification, array_ordered, name_plot):

    cut_offs = frange(0,1,0.005)
    fig, ax = plt.subplots()
    prec_colors = ['blue','green','#d48900','#b80000','black','#a5a2a2']
    #prec_colors = ['#a6d0ff','#003169']
    sens_colors = ['#9ed0ff','#32f232','#fbc562','#ff7373','#5c5c5c','#c7c7c7']
    #sens_colors = ['#a6ff9e','#085c00']
    point_colors = ['#3a6ffd','#3bbf1e','#ffb037','#ff0000','black','#a5a2a2']
    #point_colors = ['#ffb1b1','#890000']
    c = 0

    for type_analysis in array_ordered:

        precision = []
        sensitivity = []

        all_probs = results_analysis[classification][type_analysis]['in classification']['all_probs']
        for cut_off in cut_offs:

            tp,fp,tn,fn = calculate_tp_fp_tn_fn(all_probs, cut_off)

            try:
                prec_val = float(tp)/(float(tp)+float(fp))
            except:
                prec_val = 1
            sens_val = float(tp)/(float(tp)+float(fn))
            #print('PRECISION: {}\tSENSITIVITY: {}\n'.format(prec_val, sens_val))

            precision.append(prec_val)
            sensitivity.append(sens_val)
        
        ax.plot( cut_offs, precision, '-', color=prec_colors[c])
        ax.plot( cut_offs, sensitivity, '-', color=sens_colors[c])

        # Find the intersection point
        idx = np.argwhere(np.diff(np.sign(np.array(precision) - np.array(sensitivity))) != 0).reshape(-1) + 0

        # Plot the intersection point
        ax.plot(cut_offs[idx[0]], precision[idx[0]], 'o', color=point_colors[c])

        xlabel('Probability cut-offs')

        results_analysis[classification][type_analysis]['in classification']['cut_off'] = cut_offs[idx[0]]
        results_analysis[classification][type_analysis]['in classification']['value'] = precision[idx[0]]

        print('Method: {}\tValue: {}\tCut-off: {}'.format(type_analysis, precision[idx[0]], cut_offs[idx[0]]))

        # Plot vertical line from point
        #ax.plot( [cut_offs[idx[0]], cut_offs[idx[0]]], [precision[idx[0]], 0], 'r--')
        #plt.annotate(' {}'.format(cut_offs[idx[0]]), xy=(cut_offs[idx[0]], 0))
        # Plot horizontal line from point
        #ax.plot( [cut_offs[idx[0]], 0], [precision[idx[0]], precision[idx[0]]], 'r--')
        #plt.annotate('{:.3f}'.format(precision[idx[0]]), xy=(0, precision[idx[0]]+0.02))

        c+=1

    # draw temporary color lines and use them to create a legend
    hB, = plot([1,1],'-', color='blue')
    hG, = plot([1,1],'-', color='green')
    hY, = plot([1,1],'-', color='#d48900')
    hR, = plot([1,1],'-', color='#b80000')
    hBl, = plot([1,1],'-', color='black')
    hW, = plot([1,1],'-', color='#a5a2a2')
    lgd = ax.legend(handles=(hB, hG, hY, hR, hBl, hW), labels=array_ordered, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    hB.set_visible(False)
    hG.set_visible(False)
    hY.set_visible(False)
    hR.set_visible(False)
    hBl.set_visible(False)
    hW.set_visible(False)

    savefig(name_plot, format='eps', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    show()
    
    return results_analysis

def calculate_tp_fp_tn_fn(all_probs, cut_off):

    tp=0
    fp=0
    tn=0
    fn=0

    for prob, label in all_probs:

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
    all_prob = []

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
        prob = clf.predict_proba(X_test)
        
        classes = clf.classes_

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

    #print(all_prob)

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


def setBoxColors(bp, length):
    """
    Set the colors of the box plots groups
    Code from: http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    """
    if length != 6:
        print('\nERROR: This function is thought to work on groups of 6 types of analysis. Redo the function!!\n')
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

    setp(bp['boxes'][2], color='#d48900')
    setp(bp['caps'][4], color='#d48900')
    setp(bp['caps'][5], color='#d48900')
    setp(bp['whiskers'][4], color='#d48900')
    setp(bp['whiskers'][5], color='#d48900')
    setp(bp['medians'][2], color='black')
    bp['boxes'][2].set_facecolor('#fbc562') #orange

    setp(bp['boxes'][3], color='#b80000')
    setp(bp['caps'][6], color='#b80000')
    setp(bp['caps'][7], color='#b80000')
    setp(bp['whiskers'][6], color='#b80000')
    setp(bp['whiskers'][7], color='#b80000')
    setp(bp['medians'][3], color='black')
    bp['boxes'][3].set_facecolor('#ff7373') #red

    setp(bp['boxes'][4], color='black')
    setp(bp['caps'][8], color='black')
    setp(bp['caps'][9], color='black')
    setp(bp['whiskers'][8], color='black')
    setp(bp['whiskers'][9], color='black')
    setp(bp['medians'][4], color='#fbc562')
    bp['boxes'][4].set_facecolor('#4f4f4f')

    setp(bp['boxes'][5], color='black')
    setp(bp['caps'][10], color='black')
    setp(bp['caps'][11], color='black')
    setp(bp['whiskers'][10], color='black')
    setp(bp['whiskers'][11], color='black')
    setp(bp['medians'][5], color='#fbc562')
    bp['boxes'][5].set_facecolor('#aeaeae')

    return

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

if  __name__ == "__main__":
    main()
