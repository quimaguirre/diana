import argparse
import matplotlib.pyplot as plt
import numpy as np
import pylab
import sys, os, re



def main():

    options = parse_user_arguments()
    plot_general_performance(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of the input drug",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-a1','--auc_file1',dest='auc_file1',action = 'store',
                        help = """Define the file containing the AUCs for all data types for gold standard 1 (all pairs).""")
    parser.add_argument('-a2','--auc_file2',dest='auc_file2',action = 'store',
                        help = """Define the file containing the AUCs for all data types for gold standard 2 (different ATC pairs).""")
    parser.add_argument('-o','--output_plot',dest='output_plot',action = 'store',
                        help = """Define the output plot.""")

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def plot_general_performance(options):
    """
    Analyzes the results of the comparisons
    """
    methods_ordered = ['Combination', 'dctargets', 'dcguild', 'dcstructure', 'dcatc', 'dcse', 'random']
    method_to_label = {'Combination':'All', 'dctargets':'Target', 'dcguild':'PPI', 'dcstructure':'Structure', 'dcatc':'ATC', 'dcse':'Side Eff.', 'random':'Random'}
    colors_ordered1 = [['yellow', '#baba26'], ['#ff7373', '#ff7373'], ['#32f232', '#32f232'], ['#4f4f4f', '#4f4f4f'], ['#22a9bd', '#22a9bd'], ['#e59600', '#e59600'], ['#aeaeae', '#aeaeae']] # yellow, red, green, black, orange, grey            
    colors_ordered2 = [['#baba26', '#a1a13c'], ['#a11a1a', '#a11a1a'], ['#128312', '#128312'], ['#2d2b2b', 'black'], ['blue', 'blue'], ['#a96f00', '#a96f00'], ['#838181', '#838181']] # yellow, red, green, black, orange, grey            
    method_to_aucs1 = read_auc_file(options.auc_file1)
    method_to_aucs2 = read_auc_file(options.auc_file2)


    #------------------------------#
    #   PLOT DISTRIBUTION OF AUC   #
    #------------------------------#

    all_data1 = [ method_to_aucs1[method] for method in methods_ordered ]
    all_data2 = [ method_to_aucs2[method] for method in methods_ordered ]
    data_labels = [ method_to_label[method] for method in methods_ordered ]

    fig = pylab.figure(dpi=300)
    ax = pylab.axes()
    pos = 1
    xticks = []

    for x in xrange(len(methods_ordered)):


        ### Plot general data ###

        positions = []
        print(data_labels[x])
        print(all_data1[x])
        parts = ax.violinplot(all_data1[x],
                           positions = [pos],
                           showmeans=False,
                           showmedians=True)
        
        positions.append(pos)

        # Change color of the body
        for pc in parts['bodies']:
            pc.set_facecolor(colors_ordered1[x][0])

        # Change color of the segments
        parts['cmedians'].set_color(colors_ordered1[x][1])
        parts['cbars'].set_color(colors_ordered1[x][1])
        parts['cmins'].set_color(colors_ordered1[x][1])
        parts['cmaxes'].set_color(colors_ordered1[x][1])

        pos+=1

        ### Plot different ATC data ###

        parts = ax.violinplot(all_data2[x],
                           positions = [pos],
                           showmeans=False,
                           showmedians=True)
        
        positions.append(pos)

        # Change color of the body
        for pc in parts['bodies']:
            pc.set_facecolor(colors_ordered2[x][0])

        # Change color of the segments
        parts['cmedians'].set_color(colors_ordered2[x][1])
        parts['cbars'].set_color(colors_ordered2[x][1])
        parts['cmins'].set_color(colors_ordered2[x][1])
        parts['cmaxes'].set_color(colors_ordered2[x][1])

        pos+=2
        tick = np.mean(positions) # The label will be at the mean of the positions (in the middle)
        xticks.append(tick)

    # adding horizontal grid lines
    ax.yaxis.grid(True)
    ax.set_xticks([y + 1 for y in range(len(all_data1))])
    ax.set_ylabel('Distribution of AUC values')
    # add x-tick labels
    plt.setp(ax, xticks=xticks,
             xticklabels=data_labels)
    #plt.xticks(rotation=15)
    # Save
    pylab.savefig(options.output_plot, format='png')
    plt.show()


    return

#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)

def read_auc_file(auc_file):
    """
    Checks if a file exists AND is a file
    """
    method_to_aucs = {}
    with open(auc_file, 'r') as auc_fd:
        for line in auc_fd:
            method, aucs = line.strip().split('\t')
            aucs = aucs.split(',')
            aucs = [float(x) for x in aucs]
            method_to_aucs[method] = aucs
    return method_to_aucs


if  __name__ == "__main__":
    main()
