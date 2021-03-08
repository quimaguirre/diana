import argparse
import configparser
import networkx as nx
import pandas as pd
import pymysql
import time
import sys, os, re


def main():
    """
    module load Python/3.6.2
    python /home/quim/PHD/Projects/DIANA/diana/scripts/filter_drugs_without_metab.py
    """

    outputs_dir = '/home/quim/PHD/Projects/DIANA/diana/outputs'
    outputs_metab_dir = '/home/quim/PHD/Projects/DIANA/diana/outputs_targets_with_metab'
    drugbank_benchmark_file = os.path.join(outputs_dir, 'DIANA_drugbank_benchmark.txt')
    metab_drugbank_benchmark_file = os.path.join(outputs_metab_dir, 'DIANA_drugbank_benchmark.txt')
    drug_to_targets_file = os.path.join(outputs_dir, 'drugbank_geneid_drug_to_targets.txt')
    metab_drug_to_targets_file = os.path.join(outputs_metab_dir, 'drugbank_geneid_drug_to_targets.txt')
    old_results_dir = '/sbi/users/quim/DIANA_data/diana/workspace/profiles_with_metab'
    new_results_dir = '/sbi/users/quim/DIANA_data/diana/workspace/profiles'
    drugbank_drug_changes_file = os.path.join(outputs_dir, 'DIANA_drugbank_benchmark_changes.txt')

    # Gather the drug identifiers
    drugs = set()
    with open(drugbank_benchmark_file, 'r') as input_fd:
        for line in input_fd:
            drug = line.strip()
            drugs.add(drug)

    # Gather the drug identifiers in metab benchmark
    metab_drugs = set()
    with open(metab_drugbank_benchmark_file, 'r') as input_fd:
        for line in input_fd:
            drug = line.strip()
            metab_drugs.add(drug)

    drugs_to_remove = metab_drugs - drugs


    # Gather drugs with different targets
    drug_to_targets_df = pd.read_csv(drug_to_targets_file, sep='\t', index_col=None)
    metab_drug_to_targets_df = pd.read_csv(metab_drug_to_targets_file, sep='\t', index_col=None)
    drugs_with_different_targets = set()
    drugs_with_same_targets = set()
    for drug in drugs:
        targets = set()
        metab_targets = set()
        drugbankids = set()
        drugbankids.add(drug)
        for group_targets in drug_to_targets_df.loc[drug_to_targets_df['#drugbankid'].isin(drugbankids), 'geneids'].tolist():
            targets = targets | set(group_targets.split('; '))

        for group_targets in metab_drug_to_targets_df.loc[metab_drug_to_targets_df['#drugbankid'].isin(drugbankids), 'geneids'].tolist():
            metab_targets = metab_targets | set(group_targets.split('; '))

        if metab_targets != targets:
            drugs_with_different_targets.add(drug)
        else:
            drugs_with_same_targets.add(drug)
            #Copy the results in the new results folder
            old_drug_dir = os.path.join(old_results_dir, drug)
            command = 'cp -r {} {}'.format(old_drug_dir, new_results_dir)
            print(command)
            os.system(command)

    print('Drugs including metab targets: {}'.format(len(metab_drugs)))
    print('Drugs: {}'.format(len(drugs)))
    print('Drugs to remove: {}'.format(len(drugs_to_remove)))
    print('Drugs with different targets: {}'.format(len(drugs_with_different_targets)))
    print('Drugs with same targets: {}'.format(len(drugs_with_same_targets)))

    with open(drugbank_drug_changes_file, 'w') as out_fd:
        for drug in drugs_with_different_targets:
            out_fd.write('{}\n'.format(drug))


if  __name__ == "__main__":
    main()
