import sys, os, re
import collections
import ConfigParser 
import cPickle
import difflib
import hashlib
import optparse
import shutil
import subprocess

# Personal modules
import functions


def main():

    options = parse_options()
    run_comparison_profiles_cluster(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("pddi.py [--dummy=DUMMY_DIR] -i INPUT_FILE [-o OUTPUT_DIR]  [-v]")

    # Directory arguments
    parser.add_option("-c", action="store", type="string", dest="crossings_file", help="File with all the crossings that we want to compare", metavar="INPUT_FILE")
    parser.add_option("-t", action="store", type="string", dest="drug2targets", help="Pickle file containing the dictionary drug2targets", metavar="INPUT_FILE")
    parser.add_option('-w','--worspace',dest='workspace',action = 'store',default=os.path.join(os.path.dirname(__file__), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory will be created""", metavar="INPUT_DIR")
    parser.add_option("--dummy_dir", default="dummy/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = ./)", metavar="DUMMY_DIR")


    (options, args) = parser.parse_args()

    if options.crossings_file is None or options.drug2targets is None or options.workspace is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options




#################
#################
# MAIN FUNCTION #
#################
#################

def run_comparison_profiles_cluster(options):
    """
    Compare profiles in the cluster.
    """

    # Add "." to sys.path #
    src_path =  os.path.abspath(os.path.dirname(__file__))
    sys.path.append(src_path)

    # Read configuration file #     
    config = ConfigParser.ConfigParser()
    config_file = os.path.join(src_path, "config_marvin.ini")
    config.read(config_file)

    # Define which python to be used #
    python = os.path.join(config.get("Paths", "python_path"), "python")


    # Arguments & Options #
    options = parse_options()
    crossings_file = options.crossings_file
    workspace_dir = os.path.abspath(options.workspace)
    create_directory(workspace_dir)
    comparisons_dir = os.path.join(workspace_dir, 'comparisons')
    dummy_dir = os.path.abspath(options.dummy_dir)
    create_directory(dummy_dir)


    # Create log directory
    logs_dir = os.path.join(src_path, 'logs')
    create_directory(logs_dir)


    # Get the names of the drugs that we want to analyze
    crossings = set()
    with open(crossings_file, 'r') as crossings_file_fd:
        for line in crossings_file_fd:
            crossing = line.strip()
            crossings.add(crossing)

    # Open the drug2targets file
    drug2targets = cPickle.load(open(options.drug2targets))


    # Generate the profiles for each drug
    for crossing in crossings:

        drug1, drug2 = crossing.split('---')

        drug_name1 = drug1.lower()
        targets1 = list(drug2targets[drug1.upper()])
        drug_id1 = generate_drug_id(drug_name1, targets1)

        drug_name2 = drug2.lower()
        targets2 = list(drug2targets[drug2.upper()])
        drug_id2 = generate_drug_id(drug_name2, targets2)

        # Check if the results table file is already created. If so, skip
        comp_dir = os.path.join(comparisons_dir, '{}---{}'.format(drug_id1, drug_id2))
        results_table = os.path.join(comp_dir, 'results_table.tsv')
        if fileExist(results_table):
            continue


        # Define the command
        command = 'python {}/diana_cluster/scripts/compare_profiles.py -d1 {} -d2 {} -pt geneid -ws {}'.format( src_path, drug_name1, drug_name2, workspace_dir )
        print(command)

        # To run the command at the local machine
        #os.system(command)

        # To run in the cluster submitting files to queues
        functions.submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file="command_queues_marvin.txt", dummy_dir=dummy_dir)

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

def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return

def generate_drug_id(drug_name, targets):
    """
    Generates an id for the drug containing the drug name followed
    by an ID made from the processing of all the targets of the drug
    """
    m = hashlib.md5()
    targets = [str(x) for x in targets] # Transform the targets to strings
    targets_str = ''.join(sorted(targets)) # Join them in one string
    m.update(targets_str) # Introduce the string in the hashlib instance
    targets_id = m.hexdigest()[:12] # Obtain a unique ID from the targets string. Only get the first 12 characters
    drug_str = ''.join(drug_name.split('\s')) # Obtain the drug name, and remove any space in the name
    unique_id = drug_str+'_'+targets_id # Add the drug name with the targets ID creating a unique ID for the drug
    return unique_id


if  __name__ == "__main__":
    main()

