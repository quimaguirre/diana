import argparse
import copy
import cPickle
import matplotlib.pyplot as plt
import ntpath
import numpy as np
import pandas as pd
import pylab
import scipy
import time
import sys, os, re

from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler
from sklearn.dummy import DummyClassifier



def obtain_header_fields(first_line, separator='\t'):
    """ 
    Obtain a dictionary: "field_name" => "position" 
    """
    fields_dict = {}
    header_fields = first_line.strip().split(separator)
    for x in xrange(0, len(header_fields)):
        fields_dict[header_fields[x].lower()] = x
    return fields_dict


def obtain_columns(threshold_list, ATC_SE=False):
    """ 
    Obtain the name of the column from the values of the method, data_type and threshold.
    """

    columns = []

    # Get dcTargets columns
    for data_type in ['target', 'pfam', 'function']:
        for scoring_function in ['dot_product', 'spearman', 'jaccard']:
            col = 'dct'+'_'+data_type+'_'+scoring_function
            columns.append(col)

    # Get dcGUILD columns
    for top_threshold in threshold_list:
        for data_type in ['node', 'edge', 'function']:
            for scoring_function in ['dot_product', 'spearman', 'jaccard']:
                col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
                columns.append(col)

    # Get dcStructure columns
    columns.append('dcstructure')

    if ATC_SE:
        # Get dcATC columns
        for scoring_function in ['dot_product', 'spearman', 'jaccard']:
            col = 'dcatc'+'_'+scoring_function
            columns.append(col)

    if ATC_SE:
        # Get dcSE columns
        for scoring_function in ['dot_product', 'spearman', 'jaccard']:
            col = 'dcse'+'_'+scoring_function
            columns.append(col)

    columns.append('combination')

    return columns


def obtain_columns_best_features(guild_thresholds, rank_scoring, list_scoring, ATC_SE):
    """ 
    Introduce LISTS containing the values that will be used.
    """

    columns = []

    # Get dcTargets columns
    for data_type in ['target', 'pfam', 'function']:
        if data_type == 'function':
            for scoring_function in rank_scoring:
                    col = 'dct'+'_'+data_type+'_'+scoring_function
                    columns.append(col)
        else:
            for scoring_function in list_scoring:
                col = 'dct'+'_'+data_type+'_'+scoring_function
                columns.append(col)

    # Get dcGUILD columns
    for top_threshold in guild_thresholds:
        for data_type in ['node', 'edge', 'function']:
            for scoring_function in rank_scoring:
                col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
                columns.append(col)

    # Get dcStructure columns
    columns.append('dcstructure')

    if ATC_SE:
        # Get dcATC columns
        for scoring_function in list_scoring:
            col = 'dcatc'+'_'+scoring_function
            columns.append(col)

    if ATC_SE:
        # Get dcSE columns
        for scoring_function in list_scoring:
            col = 'dcse'+'_'+scoring_function
            columns.append(col)

    columns.append('combination')

    return columns


def obtain_method_to_columns_best_features(guild_thresholds, rank_scoring, list_scoring, ATC_SE):
    """ 
    Introduce LISTS containing the values that will be used.
    """

    # Get dcTargets columns
    dct_columns = []
    for data_type in ['target', 'pfam', 'function']:
        if data_type == 'function':
            for scoring_function in rank_scoring:
                    col = 'dct'+'_'+data_type+'_'+scoring_function
                    dct_columns.append(col)
        else:
            for scoring_function in list_scoring:
                col = 'dct'+'_'+data_type+'_'+scoring_function
                dct_columns.append(col)

    # Get dcGUILD columns
    dcg_columns = []
    for top_threshold in guild_thresholds:
        for data_type in ['node', 'edge', 'function']:
            for scoring_function in rank_scoring:
                col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
                dcg_columns.append(col)

    # Get dcStructure columns
    dcs_columns = []
    dcs_columns.append('dcstructure')

    dct_columns.append('combination')
    dcg_columns.append('combination')
    dcs_columns.append('combination')

    if ATC_SE:
        # Get dcATC columns
        dcatc_columns = []
        for scoring_function in list_scoring:
            col = 'dcatc'+'_'+scoring_function
            dcatc_columns.append(col)

        # Get dcSE columns
        dcse_columns = []
        for scoring_function in list_scoring:
            col = 'dcse'+'_'+scoring_function
            dcse_columns.append(col)

        dcatc_columns.append('combination')
        dcse_columns.append('combination')

        return dct_columns, dcg_columns, dcs_columns, dcatc_columns, dcse_columns
    else:
        return dct_columns, dcg_columns, dcs_columns


def obtain_columns_best_features_for_specific_method(method, guild_thresholds, rank_scoring, list_scoring):
    """ 
    Introduce LISTS containing the values that will be used.
    """
    if method == 'dctargets':
        # Get dcTargets columns
        dct_columns = []
        for data_type in ['target', 'pfam', 'function']:
            if data_type == 'function':
                for scoring_function in rank_scoring:
                        col = 'dct'+'_'+data_type+'_'+scoring_function
                        dct_columns.append(col)
            else:
                for scoring_function in list_scoring:
                    col = 'dct'+'_'+data_type+'_'+scoring_function
                    dct_columns.append(col)
        dct_columns.append('combination')
        return dct_columns

    elif method == 'dcguild':
        # Get dcGUILD columns
        dcg_columns = []
        for top_threshold in guild_thresholds:
            for data_type in ['node', 'edge', 'function']:
                for scoring_function in rank_scoring:
                    col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
                    dcg_columns.append(col)
        dcg_columns.append('combination')
        return dcg_columns

    elif method == 'dcstructure':
        # Get dcStructure columns
        dcs_columns = []
        dcs_columns.append('dcstructure')
        dcs_columns.append('combination')
        return dcs_columns

    elif method == 'dcatc':
        # Get dcATC columns
        dcatc_columns = []
        for scoring_function in list_scoring:
            col = 'dcatc'+'_'+scoring_function
            dcatc_columns.append(col)
        dcatc_columns.append('combination')
        return dcatc_columns

    elif method == 'dcse':
        # Get dcSE columns
        dcse_columns = []
        for scoring_function in list_scoring:
            col = 'dcse'+'_'+scoring_function
            dcse_columns.append(col)
        dcse_columns.append('combination')
        return dcse_columns


def obtain_method_to_columns(threshold_list, ATC_SE=False):
    """ 
    Obtain the name of the column from the values of the method, data_type and threshold.
    """

    dct_columns = []

    for data_type in ['target', 'pfam', 'function']:
        for scoring_function in ['dot_product', 'spearman', 'jaccard']:
            col = 'dct'+'_'+data_type+'_'+scoring_function
            dct_columns.append(col)
    
    dcg_columns = []

    for top_threshold in threshold_list:
        for data_type in ['node', 'edge', 'function']:
            for scoring_function in ['dot_product', 'spearman', 'jaccard']:
                col = 'dcg'+'_'+data_type+'_'+str(top_threshold)+'_'+scoring_function
                dcg_columns.append(col)

    dcs_columns = []
    dcs_columns.append('dcstructure')

    if ATC_SE:

        # Get dcATC columns
        dcatc_columns = []
        for scoring_function in ['dot_product', 'spearman', 'jaccard']:
            col = 'dcatc'+'_'+scoring_function
            dcatc_columns.append(col)

        # Get dcSE columns
        dcse_columns = []
        for scoring_function in ['dot_product', 'spearman', 'jaccard']:
            col = 'dcse'+'_'+scoring_function
            dcse_columns.append(col)


    dct_columns.append('combination')
    dcg_columns.append('combination')
    dcs_columns.append('combination')

    if ATC_SE:
        dcatc_columns.append('combination')
        dcse_columns.append('combination')
        return dct_columns, dcg_columns, dcs_columns, dcatc_columns, dcse_columns
    else:
        return dct_columns, dcg_columns, dcs_columns


def get_results_from_table(results_table, columns, combination_field):
    """ 
    Obtain the results from the results table of the comparison.
    """

    column_to_results = {}
    with open(results_table, 'r') as results_table_fd:

        # method    data_type   threshold   dot_product spearman    jaccard
        first_line = results_table_fd.readline()
        fields_dict = obtain_header_fields(first_line, separator='\t')

        for line in results_table_fd:

            fields = line.strip().split('\t')
            method = fields[ fields_dict['method'] ]
            data_type = fields[ fields_dict['data_type'] ]
            threshold = fields[ fields_dict['threshold'] ]
            dot_product = fields[ fields_dict['dot_product'] ]
            spearman = fields[ fields_dict['spearman'] ]
            jaccard = fields[ fields_dict['jaccard'] ]

            # Get dcTargets
            if method == 'dctargets':
                for scoring_function, result in [['dot_product', dot_product], ['spearman', spearman], ['jaccard', jaccard]]:
                    col = 'dct'+'_'+data_type+'_'+scoring_function
                    column_to_results[col] = result
            # Get dcGUILD
            elif method == 'dcguild':
                for scoring_function, result in [['dot_product', dot_product], ['spearman', spearman], ['jaccard', jaccard]]:
                    col = 'dcg'+'_'+data_type+'_'+str(threshold)+'_'+scoring_function
                    column_to_results[col] = result
            # Get dcStructure
            elif method == 'dcstructure':
                column_to_results['dcstructure'] = dot_product
            # Get dcATC
            elif method == 'dcatc':
                for scoring_function, result in [['dot_product', dot_product], ['spearman', spearman], ['jaccard', jaccard]]:
                    col = 'dcatc'+'_'+scoring_function
                    column_to_results[col] = result
            # Get dcSE
            elif method == 'dcse':
                for scoring_function, result in [['dot_product', dot_product], ['spearman', spearman], ['jaccard', jaccard]]:
                    col = 'dcse'+'_'+scoring_function
                    column_to_results[col] = result

    results = []
    for column in columns:
        if column in column_to_results:
            results.append(column_to_results[column])
        elif column == 'combination':
            results.append(combination_field)
        else:
            print('The column {} is not among the result columns!'.format(column))
            print('Predefined columns: {}'.format(sorted(columns)))
            print('Result columns: {}\n'.format(sorted(column_to_results.keys())))
            sys.exit(10)

    return results


def obtain_me_too_drugs_and_combinations(df, columns, me_too_drugs_table, me_too_drug_combs_table):
    """ 
    Obtain me-too drugs and me-to drug combinations in the dataset.
    """

    df_me_too = pd.DataFrame(columns=columns)

    me_too_drug_pairs = set()
    me_too_drug_comb_pairs = set()

    me_too_drugs_dict = {}
    me_too_drug_combs_dict = {}

    num_metoo_dc = 0
    num_nonmetoo_dc = 0

    done_pairings = []

    for index, row in df.iterrows():

        score = row['dcstructure']

        if score >= 0.7:
            me_too_drug_pairs.add(index)
            me_too_drugs_dict[index] = score

            df2 = pd.DataFrame([row], columns=columns, index=[index])
            df_me_too = df_me_too.append(df2)


    for index1, row1 in df_me_too.iterrows():

        score = row1['dcstructure']


        for index2, row2 in df_me_too.iterrows():

            if index1 == index2:
                continue

            combpair1 = '___'.join([index1, index2])
            combpair2 = '___'.join([index2, index1])

            if combpair1 in done_pairings or combpair2 in done_pairings:
                continue

            done_pairings.append(combpair1)
            done_pairings.append(combpair2)

            (drug11, drug12) = index1.split('---')
            (drug21, drug22) = index2.split('---')

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

                if possib1 in df_me_too.index:
                    pairing = df_me_too.loc[[possib1]]
                    group1.append(pairing)
                elif possib2 in df_me_too.index:
                    pairing = df_me_too.loc[[possib2]]
                    group1.append(pairing)
                else:
                    #print('No pairing found!')
                    num_nonmetoo_dc+=1
                    no_pairing = True

            group2 = []
            for possib1, possib2 in [ (pairing21_1,pairing21_2), (pairing22_1, pairing22_2) ]:

                if possib1 in df_me_too.index:
                    pairing = df_me_too.loc[[possib1]]
                    group2.append(pairing)
                elif possib2 in df_me_too.index:
                    pairing = df_me_too.loc[[possib2]]
                    group2.append(pairing)
                else:
                    #print('No pairing found!')
                    if no_pairing == False:
                        num_nonmetoo_dc+=1
                    no_pairing = True

            if no_pairing:
                continue

            score11 = group1[0].iloc[0]['dcstructure']
            score12 = group1[1].iloc[0]['dcstructure']

            score21 = group2[0].iloc[0]['dcstructure']
            score22 = group2[1].iloc[0]['dcstructure']


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

    print('Number of me-too drug combinations:\t{}\n'.format(num_metoo_dc))
    print('Number of non me-too drug combinations:\t{}\n'.format(num_nonmetoo_dc))

    me_too_drugs_fd = open(me_too_drugs_table, 'w')
    me_too_drug_comb_fd = open(me_too_drug_combs_table, 'w')

    # Write the results of me-too drug pairs
    for drug_pair, score in sorted(me_too_drugs_dict.iteritems(), key=lambda (x, y): y, reverse=True):
        (drug1, drug2) = drug_pair.split('---')
        name1 = '-'
        name2 = '-'
        #name1 = dcdb2name[drug1]
        #name2 = dcdb2name[drug2]
        me_too_drugs_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(drug1, name1, drug2, name2, score))

    # Write the results of me-too drug combination pairs
    for drug_comb_pair in me_too_drug_combs_dict:
        (dc1, dc2) = drug_comb_pair.split('___')
        (drug1, drug2) = dc1.split('---')
        name1 = drug1
        name2 = drug2
        #name1 = dcdb2name[drug1]
        #name2 = dcdb2name[drug2]
        (drug3, drug4) = dc2.split('---')
        name3 = drug3
        name4 = drug4
        #name3 = dcdb2name[drug3]
        #name4 = dcdb2name[drug4]

        me_too_1 = me_too_drug_combs_dict[drug_comb_pair]['me_too_1'].keys()[0]
        score1 = me_too_drug_combs_dict[drug_comb_pair]['me_too_1'][me_too_1]
        (mtd1, mtd2) = me_too_1.split('---') # mtd = me too drug
        mtn1 = mtd1
        mtn2 = mtd2
        #mtn1 = dcdb2name[mtd1] # mtn = me too drug name
        #mtn2 = dcdb2name[mtd2]

        me_too_2 = me_too_drug_combs_dict[drug_comb_pair]['me_too_2'].keys()[0]
        score2 = me_too_drug_combs_dict[drug_comb_pair]['me_too_2'][me_too_2]
        (mtd3, mtd4) = me_too_2.split('---')
        mtn3 = mtd3
        mtn4 = mtd4
        #mtn4 = dcdb2name[mtd3]
        #mtn4 = dcdb2name[mtd4]

        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name1,name2,name3,name4,mtn1,mtn2,score1,mtn3,mtn4,score2))
        me_too_drug_comb_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name1,name2,name3,name4,mtn1,mtn2,score1,mtn3,mtn4,score2))

    me_too_drugs_fd.close()
    me_too_drug_comb_fd.close()

    return me_too_drug_pairs, me_too_drug_comb_pairs


def obtain_n_groups_of_k_length(my_df, n, k, me_too_drug_combinations=None):
    """
    Obtain n number of groups of length k.
    If me_too_drug_combinations contains a list of the me-too drug combination pairs,
    it will check that there are no me-too drug combinations in different groups.
    If me_too_drug_combinations = None, it will not do anything.
    """

    repeat = True

    while repeat:
        groups = []
        repeat = False
        curr_df = my_df.copy() # Copy the dataframe, so that if we repeat we have the initial one
        for y in xrange(n):
            new_df = curr_df.sample(n=k) # Get a random sample of length k from the main dataframe
            curr_df = curr_df.loc[~curr_df.index.isin(new_df.index)]  # Remove the sample that we have taken from the main dataframe
            groups.append(new_df) # Append the sample in the list gropus

        # Check if two me-too drug combinations are part of two different groups
        # If this happens, we will repeat the process (because they could be used in training / testing at the same time)
        if me_too_drug_combinations:
            for pair in me_too_drug_combinations:
                drug_comb1, drug_comb2 = pair
                comb1_group = None
                comb2_group = None
                for x in xrange(len(groups)):
                    indexes = groups[x].index.values
                    if drug_comb1 in groups[x].index.values:
                        comb1_group = x
                    if drug_comb2 in groups[x].index.values:
                        comb2_group = x
                if comb1_group and comb2_group and comb1_group != comb2_group:
                    repeat = True
                    break

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
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
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


def run_nfold_crossvalidation_scikit_with_prob(n, groups, classifier):
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
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        auc = metrics.roc_auc_score(y_test, y_pred)
        all_auc.append(auc)

    mean_tpr /= n
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    mean_auc2 = np.mean(all_auc)

    var_auc = np.var(all_auc)
    std_auc = np.std(all_auc)

    return mean_auc, var_auc, std_auc, all_auc, all_prob


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
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
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

