import sys, os, re
import collections
import cPickle
import ConfigParser 
import difflib
import hashlib
import ntpath
import optparse
import shutil
import subprocess

# Personal modules
import functions


def main():

    options = parse_options()
    run_all_profiles_in_cluster(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("pddi.py [--dummy=DUMMY_DIR] -i INPUT_FILE [-o OUTPUT_DIR]  [-v]")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_file", help="File with all the input drugs from which we want the profiles", metavar="INPUT_FILE")
    parser.add_option("-s", action="store", type="string", dest="sif_file", help="Input SIF file")
    parser.add_option("-t", action="store", type="string", dest="drug2targets", help="Pickle file containing the dictionary drug2targets", metavar="INPUT_FILE")
    parser.add_option('-w','--workspace',dest='workspace',action = 'store',default=os.path.join(os.path.dirname(__file__), 'workspace'),
                        help = """Define the workspace directory where the data directory and the results directory will be created""")
    parser.add_option("--dummy_dir", default="dummy/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = ./)", metavar="DUMMY_DIR")


    (options, args) = parser.parse_args()

    if options.input_file is None or options.sif_file is None or options.drug2targets is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options




#################
#################
# MAIN FUNCTION #
#################
#################

def run_all_profiles_in_cluster(options):
    """
    Generates all the profiles in the cluster
    """
    # Limit of jobs
    limit = None
    l = 1

    # Add "." to sys.path #
    src_path =  os.path.abspath(os.path.dirname(__file__))
    sys.path.append(src_path)

    # Read configuration file #     
    config = ConfigParser.ConfigParser()
    config_file = os.path.join(src_path, "config_hydra.ini")
    config.read(config_file)

    # Define which python to be used #
    python = os.path.join(config.get("Paths", "python_path"), "python")


    # Arguments & Options #
    options = parse_options()
    input_file = options.input_file
    sif_file = options.sif_file
    workspace_dir = os.path.abspath(options.workspace)
    create_directory(workspace_dir)
    profiles_dir = os.path.join(workspace_dir, 'profiles')

    # Create the dummy directory (to store the commands)
    dummy_dir = os.path.abspath(options.dummy_dir)
    create_directory(dummy_dir)

    # Create log directory (to store the log files of the jobs)
    logs_dir = os.path.join(src_path, 'logs')
    create_directory(logs_dir)

    # Get the names of the drugs that we want to analyze
    all_drugs = set()
    with open(input_file, 'r') as input_file_fd:
        for line in input_file_fd:
            drug = line.strip()
            all_drugs.add(drug)

    # Open the drug2targets file
    drug2targets = cPickle.load(open(options.drug2targets))

    # Generate the profiles for each drug
    for drug in all_drugs:

        # Check if the p-value file is already created. If so, skip
        drug_name = drug.lower()
        targets = list(drug2targets[drug.upper()])
        network_filename = ntpath.basename(options.sif_file)
        drug_id = generate_drug_id(drug_name, targets, network_filename)
        drug_dir = os.path.join(profiles_dir, drug_id)
        dcguild_dir = os.path.join(drug_dir, 'dcguild_profiles')
        pvalue_file = os.path.join(dcguild_dir, 'output_scores.sif.netcombo.pval')
        dcstructure_dir = os.path.join(drug_dir, 'dcstructure_profiles')
        structure_file = os.path.join(dcstructure_dir, 'structure_profile.txt')
        if fileExist(structure_file):
            pass
        else:
            continue

        if limit: # Break the loop if a limit of jobs is introduced
            if l > limit:
                break

        # Define the command
        command = 'python {}/diana/scripts/generate_profiles.py -d {} -pt geneid -sif {} -ws {}'.format( src_path, drug_name, options.sif_file, workspace_dir )
        print(command)

        #To run in the cluster submitting files to queues
        functions.submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file="command_queues_hydra.txt", dummy_dir=dummy_dir)

        l += 1

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

def generate_drug_id(drug_name, targets, network_name):
    """
    Generates an ID for the drug using the drug name, the sorted targets and
    the network name.
    """
    id_str = ''
    drug_str = ''.join(drug_name.split('\s')) # Obtain the drug name, and remove any space in the name
    id_str += drug_str.lower() # Add the drug name in the ID string
    targets = [str(x) for x in targets] # Transform the targets to strings
    targets_str = ''.join(sorted(targets)) # Join the targets in one string
    id_str += targets_str.lower() # Add the targets in the ID string
    id_str += network_name.lower() # Add the network name in the ID string
    m = hashlib.md5()
    m.update(id_str) # Introduce the string in the hashlib instance
    unique_id = m.hexdigest()[:12] # Obtain a unique ID from the string. Only get the first 12 characters
    return unique_id

def old_generate_drug_id(drug_name, targets):
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

