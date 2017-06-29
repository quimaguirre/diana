import pandas as pd
import cPickle
import os
from os import listdir
from os.path import isfile, isdir, join
import re
import sys

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

    # File where we will store the results of the analysis
    #data_frame_file = results_dir + '/all_results_5targets.csv'
    data_frame_file = results_dir + '/all_results_3targets.csv'

    # Crossings file
    analysis_file = results_dir + "/diana_crossings_cluster_3targets.txt"



    ###########################
    #### PARSE THE RESULTS ####
    ###########################

    main_directory = results_dir + "/results"

    dump_file = current_dir + "/drug_int_2_drugs.pcl"
    drug_int_2_drugs = cPickle.load(open(dump_file))
    dump_file = current_dir + "/pair2comb.pcl"
    pair2comb = cPickle.load(open(dump_file))

    ddi = sum(1 for x in pair2comb.values() if x == 1)
    non_ddi = sum(1 for x in pair2comb.values() if x == 0)

    print('NUMBER OF DRUG COMBINATIONS:\t\t{}'.format(ddi))
    print('NUMBER OF NON-DRUG COMBINATIONS:\t{}'.format(non_ddi))


    #columns = ['N5sp', 'N5dp', 'N10sp', 'N10dp', 'N20sp', 'N20dp', 'N50sp', 'N50dp', 'N100sp', 'N100dp', 'E5sp', 'E5dp', 'E10sp', 'E10dp', 'E20sp', 'E20dp', 'E50sp', 'E50dp', 'E100sp', 'E100dp', 'F5sp', 'F5dp', 'F10sp', 'F10dp', 'F20sp', 'F20dp', 'F50sp', 'F50dp', 'LNsp', 'LNdp', 'LEsp', 'LEdp', 'LFsp', 'LFdp', 'SNsp', 'SNdp', 'SNji', 'SPsp', 'SPdp', 'SPji', 'SFsp', 'SFdp', 'Struc', 'Comb']
    #print(len(columns))
    #columns = ['per0.1','per0.5','per1','per2.5','per5','per10','per20', 'it50', 'it100', 'it250', 'it500', 'it1000']

    # Define the names of the columns in the future data frame
    columns = []
    for tprof in ('N', 'E'):
        #for top in ('5', '10', '20', '50', '100'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000', '100'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
    for tprof in ('F'):
        #for top in ('5', '10', '20', '50'):
        for top in ('per0.1', 'per0.5', 'per1', 'per2.5', 'per5', 'per10', 'per20', 'it50', 'it100', 'it250', 'it500', 'it1000'):
            for an in ('sp', 'dp'):
                string = tprof+top+an
                columns.append(string)
    for tprof in ('LN', 'LE'):
        for an in ('sp', 'dp'):
            string = tprof+an
            columns.append(string)
    tprof = 'LF'
    for an in ('sp', 'dp'):
        string = tprof+an
        columns.append(string)
    for tprof in ('SN', 'SP'):
        for an in ('sp', 'dp', 'ji'):
            string = tprof+an
            columns.append(string)
    columns.append('SFsp')
    columns.append('SFdp')
    columns.append('Struc')
    columns.append('Comb')



    # Create a data frame to store the results
    df = pd.DataFrame(columns=columns)


    # Obtain all the results subfolders of the results main folder
    results_dir_list = [join(main_directory, f) for f in listdir(main_directory) if isdir(join(main_directory, f))]

    for directory in results_dir_list:

        # Obtain all the result files inside each result subfolder
        results_file_list = [join(directory, f) for f in listdir(directory) if isfile(join(directory, f))]

        table_regex = re.compile('/table_results')
        results_regex = re.compile('results_.+/results_.+.txt$')
        extended_regex = re.compile('/extended_analysis')

        for results_file in results_file_list:

            # The file matches with the TABLE RESULTS pattern
            if table_regex.search(results_file) != None:

                file_name = os.path.basename(results_file)
                (drug1, drug2) = results_file.split('.txt')[0].split('table_results_')[1].split('_')
                pair = '{}---{}'.format(drug1, drug2)

                fg = open(results_file,"r")

                for line in fg:

                    # Skip the comments
                    if line[0] == "#":
                        continue

                    # Obtain the fields
                    fields = line.strip().split('\t')
                    # Add the Comb field (if it is drug combination or not)
                    fields.append(pair2comb[pair])

                    # If the number of fields does not correspond to the expected, the script is stopped
                    #if len(fields) != 44:
                    if len(fields) != 92:
                    #if len(fields) != 148:
                        #print('The length is not 44')
                        print('The length does not correspond to the number of columns')
                        sys.exit(10)

                    df2 = pd.DataFrame([fields], columns=columns, index=[pair])
                    # Add the information to the main data frame
                    df = df.append(df2)

    #print(df)

    # Get data frame
    df.to_csv(data_frame_file)

    return


if  __name__ == "__main__":
    main()