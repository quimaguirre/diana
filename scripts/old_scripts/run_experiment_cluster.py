import os, sys, re
import ConfigParser 
import optparse
import shutil
import subprocess
import difflib
import collections
#import numpy as np

# Alberto Meseguer file; 18/11/2016
# Modified by Quim Aguirre; 13/03/2017

# This file is the master coordinator of the DIANA project. It is used to run multiple DIANA commands in parallel in the cluster

#-------------#
# Functions   #
#-------------#

#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("pddi.py [--dummy=DUMMY_DIR] -i INPUT_FILE [-o OUTPUT_DIR]  [-v]")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input crossings file", metavar="INPUT_FILE")
    parser.add_option("-s", action="store", type="string", dest="sif_file", help="Input SIF file")
    parser.add_option("-t", action="store", type="string", dest="type_of_analysis", help="Type of analysis: 'profile_creation' or 'comparison'")
    parser.add_option("--dummy_dir", default="dummy/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = ./)", metavar="DUMMY_DIR")
    parser.add_option('-ws','--worspace',dest='workspace',action = 'store',default=os.path.join(os.path.dirname(__file__), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory will be created""")


    (options, args) = parser.parse_args()

    if options.input_file is None or options.sif_file is None or options.type_of_analysis is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#-------------#
# Main        #
#-------------#

# Add "." to sys.path #
src_path =  os.path.abspath(os.path.dirname(__file__))
sys.path.append(src_path)

# Read configuration file #     
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, "config_marvin.ini")
config.read(config_file)

import hashlib

# Imports my functions #
import functions


# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Arguments & Options #
options = parse_options()
# Directory arguments
input_file = os.path.abspath(options.input_file)
dummy_dir = os.path.abspath(options.dummy_dir)


# Create directories if necessary

logs_dir = src_path + "/logs"
if not os.path.exists(logs_dir):
    os.mkdir(logs_dir)


f = open(input_file, "r")

# Depending on the type of analysis, we will submit different commands

if options.type_of_analysis == 'profile_creation':

    analysis = '-prof'
    all_drugs = set()

    for line in f:

        (drug1, drug2) = line.strip().split('---')
        all_drugs.add(drug1)
        all_drugs.add(drug2)

    f.close()

    for drug in all_drugs:

        # Check if the p-value file is already created. If so, skip
        pvalue_file = data_dir + "/" + drug + "/guild_results_using_sif/output_scores.sif.netcombo.pval"
        if os.path.exists(pvalue_file):
            continue        

        guild_path = '/gpfs42/robbyfs/homes/users/qaguirre/guild/scoreN'
        command = 'python {}/diana_cluster/scripts/generate_profiles.py -d {} -pt geneid -sif {} -gu {}'.format( src_path, drug, options.sif_file, guild_path )
        print(command)
        # python /home/quim/project/diana_cluster/scripts/generate_profiles.py -d 'DCC0303' -pt 'geneid' -sif /home/quim/project/diana_cluster/workspace/sif/human_eAFF_geneid_2017.sif -gu /home/quim/project/diana_cluster/diana/toolbox/scoreN

        # To run the command at the local machine
        #os.system(command)

        #To run in the cluster submitting files to queues
        functions.submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file="command_queues_marvin.txt", dummy_dir=dummy_dir)


elif options.type_of_analysis == 'comparison':

    analysis = '-comp'

    for line in f:

        (drug1, drug2) = line.strip().split('---')

        # Check if the results are already done
        comp_results_dir = res_dir + "/results_" + drug1 + "_" + drug2
        table_file = comp_results_dir + '/table_results_' + drug1 + '_' + drug2 + '.txt'
        if os.path.exists(table_file):
            continue

        command = 'python {}/diana_cluster/scripts/compare_profiles.py -d1 {} -d2 {} -pt geneid'.format( src_path, drug1, drug2 )
        print(command)
        # python /home/quim/project/diana_cluster/scripts/compare_profiles.py -d1 'DCC0303' -d2 'DCC1743' -pt 'geneid'
        
        # To run the command at the local machine
        #os.system(command)

        #To run in the cluster submitting files to queues
        functions.submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file="command_queues_marvin.txt", dummy_dir=dummy_dir)

    f.close()

else:
    print('The type of analysis has been wrongly defined. Introduce \'profile_creation\' or \'comparison\'')
    sys.exit(10)



