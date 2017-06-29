import sys
import argparse
import os
import re

def main():

    options = parse_user_arguments()
    weight(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Weight the network edges by tissue or subcellular-specific localization of the interactions",
        epilog      = "@oliva's lab 2016")
    parser.add_argument('-tissue','--seeds_input_file',dest='seed',action = 'store',default='input_seed',
                        help = 'Seeds Input file (default is input_seed)')
    options=parser.parse_args()

    return options

def weight(options):







if  __name__ == "__main__":
    main()
