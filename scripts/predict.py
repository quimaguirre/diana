import cPickle
import pandas as pd
import time
import sys, os, re




def main():

    options = parse_user_arguments()
    predict(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Predict if a pair of drugs is a drug combination",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-d1','--drug_name1',dest='drug_name1',action = 'store',
                        help = """ Name of the drug number 1. If you do not provide targets for this drug or the number of targets is not large enough,
                        the program will use this name to search for targets in BIANA database. If targets are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between 
                        double quotes. """)
    parser.add_argument('-d2','--drug_name2',dest='drug_name2',action = 'store',
                        help = """ Name of the drug number 2. If you do not provide targets for this drug or the number of targets is not large enough,
                        the program will use this name to search for targets in BIANA database. If targets are provided, this field will be only used
                        for naming purposes and will be completely optional.
                        If the name of the drug has more than one word or special characters (parentheses, single quotes), introduce the name between 
                        double quotes. """)
    parser.add_argument('-t1','--targets1',dest='targets1',action = 'store',
                        help = 'Input file with the targets of the drug 1. Each target must be separated by a newline character.')
    parser.add_argument('-t2','--targets2',dest='targets2',action = 'store',
                        help = 'Input file with the targets of the drug 2. Each target must be separated by a newline character.')
    parser.add_argument('-pt','--proteins_type_id',dest='proteins_type_id',action = 'store', default='geneid',
                        help = 'Input the type of ID of the targets introduced / proteins of the network. It must be the same! (default is geneid).')
    parser.add_argument('-th','--threshold_list',dest='threshold_list',action = 'store',
                        help = """List of percentages that will be used as cut-offs to define the profiles of the drugs. It has to be a file containing:
                        - Different numbers that will be the threshold values separated by newline characters. 
                        For example, a file called "top_threshold.list" containing:
                        0.1
                        0.5
                        1
                        5
                        10
                        """)
    parser.add_argument('-ws','--worspace',dest='workspace',action = 'store',default=os.path.join(os.path.join(os.path.dirname(__file__), '..'), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory will be created""")
    parser.add_argument('-db','--database',dest='database',action = 'store',default='BIANA_JUN_2017',
                        help = """Define the database to use for the search of targets: 
                        (default is BIANA_JUN_2017)""")
    parser.add_argument('-dbu','--db_user',dest='db_user',action = 'store',default='quim',
                        help = """Define the MySQL user to access to the database: 
                        (default is quim)""")
    parser.add_argument('-dbp','--db_pass',dest='db_pass',action = 'store',default='',
                        help = """Define the MySQL password to access to the database: 
                        (default is '')""")
    parser.add_argument('-dbh','--db_host',dest='db_host',action = 'store',default='localhost',
                        help = """Define the MySQL host to access to the database: 
                        (default is localhost)""")
    parser.add_argument('-up','--unification',dest='unification_protocol',action = 'store',default='geneid_seqtax_v1',
                        help = """Define the unification protocol used in BIANA database (default is BIANA_JUN_2017)""")

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def predict(options):
    """
    Predicts if a pair of drugs is a drug combination
    """

    # Start marker for time measure
    start = time.time()

    print("\n\t\t------------------------------------------------------------------------------------------------------\n")
    print("\t\tStarting Drug Interactions ANAlysis (DIANA), a program created by @OLIVA'S LAB. Third part: Prediction\n")
    print("\t\t------------------------------------------------------------------------------------------------------\n")

    # Get the script path
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    # Check that the other_data directory exists
    other_data_dir = os.path.join(main_path, 'workspace/other_data')
    check_directory(other_data_dir)

    # Table with the results of DIANA
    data_frame_file = os.path.join(other_data_dir ,'all_results_3targets.csv')

    # Load data frame
    df = pd.read_csv(data_frame_file, index_col=0)

    dc_data = df[df['Comb'] == 1]
    num_dc = len(dc_data.index)
    print('Number of drug combinations with at least 3 targets: {}\n'.format(num_dc))

    # Replace the None values in Struc by nan
    df = df.replace(to_replace={'Struc':{'None':np.nan}})
    # Replace the NA values in Struc by nan
    df = df.replace(to_replace={'Struc':{'NA':np.nan}})

    # Deleting the spearman for seeds. It is useless
    df = df.drop('SNsp', axis=1)
    df = df.drop('SPsp', axis=1)



    # End marker for time
    end = time.time()
    print('\n  DIANA INFO:\tTIME OF EXECUTION: {:.3f} seconds or {:.3f} minutes.\n'.format(end - start, (end - start) / 60))

    return



#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


if  __name__ == "__main__":
    main()

