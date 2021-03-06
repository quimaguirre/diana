#!/usr/bin/env python
import imp
import sys
import os
from context import diana
import diana.toolbox.guild_utilities as guild_utilities

def main():
    # Set the seed and network files
    if len(sys.argv)>6:
     data_dir       =sys.argv[1]+"/"
     seed_file      =sys.argv[2]
     network_file   =sys.argv[3]
     scoring_folder =sys.argv[4]+"/"
     random_networks_folder =sys.argv[5]+"/"
     executable_path = sys.argv[6]
    elif len(sys.argv)==6:
     data_dir       =sys.argv[1]+"/"
     seed_file      =sys.argv[2]
     network_file   =sys.argv[3]
     scoring_folder =sys.argv[4]+"/"
     random_networks_folder =sys.argv[5]+"/"
    elif len(sys.argv)==5:
     data_dir       =sys.argv[1]+"/"
     seed_file      =sys.argv[2]
     network_file   =sys.argv[3]
     scoring_folder =sys.argv[4]+"/"
     random_networks_folder = scoring_folder
    elif len(sys.argv)==4:
     data_dir       =sys.argv[1]+"/"
     seed_file = data_dir + "seeds.txt"
     network_file = data_dir + "network_guild.sif"
     scoring_folder = "./guild/"
     random_networks_folder = scoring_folder
    elif len(sys.argv)==3:
     data_dir       =sys.argv[1]+"/"
     seed_file      =sys.argv[2]
     network_file = data_dir + "network_guild.sif"
     scoring_folder = "./guild/"
     random_networks_folder = scoring_folder
    elif len(sys.argv)==2:
     data_dir       =sys.argv[1]+"/"
     seed_file      =sys.argv[2]
     network_file   =sys.argv[3]
     scoring_folder = "./guild/"
     random_networks_folder = scoring_folder
    else:
     data_dir = "./"
     seed_file = data_dir + "seeds.txt"
     network_file = data_dir + "network_guild.sif"
     scoring_folder = "./guild/"
     random_networks_folder = scoring_folder

    script_path = os.path.abspath(os.path.dirname(__file__))
    #executable_path = os.path.join(script_path, 'scoreN')

    # Create input files for scoring
    guild_utilities.prepare_scoring(network_file, seed_file, scoring_folder, random_networks_folder, non_seed_score=0.01,seed_score=1.0,
                                    edge_score=1.0, n_sample=100, delim="\t")

    # Run GUILD and create output files, the case for Netcombo
    guild_utilities.run_scoring(scoring_folder, random_networks_folder, executable_path, scoring_type="netcombo", calculate_pvalue=False)
    #run_scoring(scoring_folder, executable_path, scoring_type="netzcore", parameters={"n_iteration":5,
    #"n_sample":100, "sampling_prefix":scoring_folder+"sampled_graph."}, qname=None)
    return

if __name__ == "__main__":
    main()

