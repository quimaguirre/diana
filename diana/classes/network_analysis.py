import pickle
import networkx as nx
import numpy as np
import scipy.stats
import os, sys

import diana.classes.network_translation as NT
import diana.classes.tissue_specificity as TS
import diana.classes.top_scoring as TOP

"""
    NetworkAnalysis
    2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
    Universitat Pompeu Fabra
"""

#################
#### CLASSES ####
#################

class Network(object):
    """ 
    Class defining a network object 
    """

    def __init__(self, network_file, type_id, network_format):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    type_id
        @pdef:     Type of IDs in the network
        @ptype:    {String}

        @param:    network_format
        @pdef:     Format of the network
        @ptype:    {String} {'sif' or 'multi-fields'}

        @raises: {IncorrectNetworkFormat} if the network_format is not in
        self.formats.
        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        self.network_file = network_file
        self.type_id = type_id
        self.network_format = network_format
        self.formats = ['sif', 'multi-fields', 'raw']
        self.main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        self.pickles_path = os.path.join(self.main_path, 'data')

        self._tissue_specific = False

        self.network = self.parse_network(tissue_specific=self._tissue_specific)

    ###########
    # METHODS #
    ###########

    def get_edges(self):
        """
        Obtain a list of all the edges, where each edge is a tuple with the two
        nodes of the interaction
        """
        return self.network.edges()

    def get_nodes(self):
        """
        Obtain a list of all the nodes in the network
        """
        return self.network.nodes()

    def get_targets_in_network(self, targets):
        """
        Get the targets that are inside the network given by the user
        """
        targets_in_network = set()
        str_tar = [str(x) for x in targets]
        nodes = self.get_nodes()
        for target in str_tar:
            if target in nodes:
                targets_in_network.add(target)
        return list(targets_in_network)

    def get_methodid2interactions(self, verbose=False):
        """
        Obtains a dictionary containing all the method IDs that appear in the 
        network and the number of times that they appear.
        If verbose:
        Prints the methods and number of interactions in order (higher to lower)
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        psimi2method_file = os.path.join(self.pickles_path,'psimi2method.pcl')
        self.psimi2method = pickle.load(open(psimi2method_file, 'rb'))

        methodid2interactions = {}

        for u,v,d in self.network.edges_iter(data=True):

            for method_id in d['method_ids']:
                method_id = int(method_id)
                methodid2interactions.setdefault(method_id, 0)
                methodid2interactions[method_id] += 1

        if verbose:
            for psimi, interactions in sorted(methodid2interactions.items(), key=lambda x: x[1], reverse=True):
                if psimi in self.psimi2method:
                    name = self.psimi2method[psimi]
                else:
                    name = '-'
                print('Method ID: {}\tName: {}\tNumber of interactions: {}'.format(psimi, name, interactions))

        return methodid2interactions

    def get_method2interactions(self, verbose=False):
        """
        Obtains a dictionary containing all the method names that appear in the 
        network and the number of times that they appear.
        If verbose:
        Prints the methods and number of interactions in order (higher to lower)
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        method2interactions = {}

        for u,v,d in self.network.edges_iter(data=True):

            for method in d['method_names']:
                method2interactions.setdefault(method, 0)
                method2interactions[method] += 1

        if verbose:
            for method, interactions in sorted(method2interactions.iteritems(), key=lambda x: x[1], reverse=True):
                print('Method: {}\tNumber of interactions: {}'.format(method, interactions))

        return method2interactions

    def get_numpmids2interactions(self, verbose=False):
        """
        Obtains a dictionary containing the number of interactions containing a 
        given number of pubmed IDs
        If verbose:
        Prints the number of pubmeds and number of interactions in order of
        number of pubmeds (higher to lower)
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        numpmids2interactions = {}

        for u,v,d in self.network.edges_iter(data=True):

            num_pmids = len(d['pmids'])
            numpmids2interactions.setdefault(num_pmids, 0)
            numpmids2interactions[num_pmids] += 1

        if verbose:
            for num_pmids, interactions in sorted(numpmids2interactions.iteritems(), key=lambda x: x[1]):
                print('Number of pubmeds: {}\tNumber of interactions: {}'.format(num_pmids, interactions))

        return numpmids2interactions

    def get_database2interactions(self, verbose=False):
        """
        Obtains a dictionary containing all the databases that appear in the 
        network and the number of times that they appear.
        If verbose:
        Prints the databases and number of interactions in order (higher to 
        lower)
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        database2interactions = {}

        for u,v,d in self.network.edges_iter(data=True):

            for source in d['sources']:
                database2interactions.setdefault(source, 0)
                database2interactions[source] += 1

        if verbose:
            for source, interactions in sorted(database2interactions.iteritems(), key=lambda x: x[1], reverse=True):
                print('Database: {}\tNumber of interactions: {}'.format(source, interactions))

        return database2interactions

    def parse_network(self, tissue_specific=False, verbose=False):
        """
        Parse the network file using the Python module NetworkX.
        It is possible to parse the network in two formats:
            - 'sif' : <node1>\t<score>\t<node2>
            - 'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>
            - 'multi-fields' + tissue_specific=True : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>\t<tissue_db>\t<tissue_additional>
        """

        if verbose:
            print('Parsing network...\n')

        G=nx.Graph()

        network_fd = open(self.network_file, 'r')

        for line in network_fd:

            if self.network_format == 'sif':
                fields = line.strip().split('\t')
                if len(fields) == 3:
                    (node1, score, node2) = fields
                    G.add_edge(node1,node2,score=score)
                else:
                    print('The SIF input network does not have 3 fields: {}'.format(line.strip()))
                    print('Please, provide a SIF network that has 3 fields in each line')
                    sys.exit(10)

            elif self.network_format == 'multi-fields':
                if tissue_specific:
                    (node1, node2, sources, method_ids, method_names, pubmeds, databases, additional) = line.strip().split('\t')
                    sources = sources.split(';')
                    method_ids = method_ids.split(';')
                    method_names = method_names.split(';')
                    pubmeds = pubmeds.split(';')
                    databases = databases.split(';')
                    additional = additional.split(';')
                    G.add_edge(node1, node2, sources=sources, method_ids=method_ids, method_names=method_names, pmids=pubmeds, tissue_db=databases, tissue_additional=additional)
                else:
                    (node1, node2, sources, method_ids, method_names, pubmeds) = line.strip().split('\t')
                    sources = sources.split(';')
                    method_ids = method_ids.split(';')
                    method_names = method_names.split(';')
                    pubmeds = pubmeds.split(';')
                    G.add_edge(node1, node2, sources=sources, method_ids=method_ids, method_names=method_names, pmids=pubmeds)
            elif self.network_format == 'raw':
                (node1, node2) = line.strip().split('\t')
                G.add_edge(node1,node2)
            else:
                raise IncorrectNetworkFormat(self.network_format, self.formats)

        network_fd.close()

        if verbose:
            print('Parsing network... finished!\n')

        return G


    def score_network(self, node_to_score, output_file):
        """
        Uses a BIANA translation file to return a new Network object with the 
        desired type of IDs and format

        @param:    node_to_score
        @pdef:     Dictionary of all the nodes to a tuple containing (score, pvalue)
        @ptype:    String

        @param:    output_file
        @pdef:     Output file where the scored network will be written
        @ptype:    String
        """

        # For now, it is only possible to score the network in sif format
        # It has to be implemented for the other formats
        if self.network_format != 'sif':
            raise IncorrectNetworkFormat(self.network_format)

        TOP.score_network(self.network, node_to_score, output_file)

        return EdgeProfile(network_file=output_file, type_id=self.type_id, network_format=self.network_format, top=100)


    def filter_network_by_tissue(self, filtered_network_file, filtered_nodes_file, tissue_object, permission, verbose=False):
        """
        Filter the network using only interactions where the proteins are 
        present in a Tissue object

        @param:    filtered_network_file
        @pdef:     File where the tissue-specific network will be written
        @ptype:    String

        @param:    filtered_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @param:    tissue_object
        @pdef:     Tissue class object used to filter the network
        @ptype:    Tissue class object

        @param:    permission
        @pdef:     Level of permission to create the network
        @ptype:    Integer {0,1,2}
        """

        if verbose:
            print('Filtering network by tissue...\n')

        TS.filter_network_tissue_specific(self.network_file, tissue_object, permission, filtered_network_file, filtered_nodes_file)

        if verbose:
            print('Filtering network by tissue... finished!\n')

        return TissueSpecificNetwork(filtered_network_file, filtered_nodes_file, self.type_id, self.network_format, tissue_object, permission)


    def filter_network_by_method(self, methods_excluded=None, method_ids_excluded=None, methods_included=None, method_ids_included=None, output_network_file=None, output_nodes_file=None, verbose=False):
        """
        Filter the network: 
        - By excluding interactions only reported by methods in 'methods_excluded'
          and 'method_ids_excluded' list
        - By including interactions at least reported by one ofthe methods in 
          'methods_included' list and 'method_ids_included' list

        @param:    methods_excluded
        @pdef:     Method names list which will exclude interactions if they are
                   constituted by methods in this list
        @ptype:    List

        @param:    method_ids_excluded
        @pdef:     Same as methods_excluded for psi-mi IDs list
        @ptype:    List

        @param:    methods_included
        @pdef:     Method names list which will only include interactions that
                   at least contain one of the methods in this list
        @ptype:    List

        @param:    method_ids_included
        @pdef:     Same as methods_included for psi-mi IDs list
        @ptype:    List

        @param:    output_network_file
        @pdef:     File where the network will be written
        @ptype:    String

        @param:    output_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @param:    verbose
        @pdef:     If True, informs about when the parsing starts and ends
        @ptype:    Bool (Default: False)

        @return:   Network object with the filtered network
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        if verbose:
            print('Filtering network by method...\n')

        fnet=nx.Graph()

        for u,v,d in self.network.edges_iter(data=True):
            skip_inc = False
            skip_exc = False

            # Check if at least one of the imprescindible methods is included
            if methods_included != None:
                skip_inc = True
                for method in d['method_names']:
                    if method in methods_included:
                        skip_inc = False
            if method_ids_included != None:
                skip_inc = True
                for method in d['method_ids']:
                    if method in methods_included:
                        skip_inc = False

            # Check if the interaction has at least another method apart from
            # the ones in methods_excluded/method_ids_excluded
            if methods_excluded != None:
                skip_exc = True
                for method in d['method_names']:
                    if method not in methods_excluded:
                        skip_exc = False

            if method_ids_excluded != None:
                method_ids_excluded = [str(x) for x in method_ids_excluded]
                skip_exc = True
                for method in d['method_ids']:
                    if method not in method_ids_excluded:
                        skip_exc = False

            if skip_inc == False and skip_exc == False:
                fnet.add_edge(u,v,d)

        write_network_file_from_networkx_graph(fnet, output_network_file, output_nodes_file, self.network_format, tissue_specific=self._tissue_specific)

        if verbose:
            print('Filtering network by method... finished!\n')

        return Network(output_network_file, output_nodes_file, self.type_id, self.network_format)

    def filter_network_by_number_pubmeds(self, min_num_pubmeds, output_network_file, output_nodes_file, verbose=False):
        """
        Filter the network by minimum number of pubmed ids

        @param:    min_num_pubmeds
        @pdef:     Minimum number of pubmeds which an interaction has to have
        @ptype:    Integer

        @param:    output_network_file
        @pdef:     File where the network will be written
        @ptype:    String

        @param:    output_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @param:    verbose
        @pdef:     If True, informs about when the parsing starts and ends
        @ptype:    Bool (Default: False)

        @return:   Network object with the filtered network
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        if verbose:
            print('Filtering network by number of pubmeds...\n')

        fnet=nx.Graph()

        for u,v,d in self.network.edges_iter(data=True):
            skip_inc = False
            skip_exc = False

            # Check if it has at least the minimum number of pubmeds required
            number_pubmeds = len(d['pmids'])
            if number_pubmeds >= min_num_pubmeds:
                fnet.add_edge(u,v,d)

        write_network_file_from_networkx_graph(fnet, output_network_file, output_nodes_file, self.network_format, tissue_specific=self._tissue_specific)

        if verbose:
            print('Filtering network by number of pubmeds... finished!\n')

        return Network(output_network_file, output_nodes_file, self.type_id, self.network_format)

    def filter_network_by_database(self, databases_included, output_network_file, output_nodes_file, verbose=False):
        """
        Filter the network by interactions included only in certain databases.
        If the interaction is from one of the databases, it is included

        @param:    databases_included
        @pdef:     Databases from which the interaction must come from
        @ptype:    List

        @param:    output_network_file
        @pdef:     File where the network will be written
        @ptype:    String

        @param:    output_nodes_file
        @pdef:     File where the nodes file will be written
        @ptype:    String

        @param:    verbose
        @pdef:     If True, informs about when the parsing starts and ends
        @ptype:    Bool (Default: False)

        @return:   Network object with the filtered network
        """
        if self.network_format != 'multi-fields':
            print('It is only possible to use this method with a multi-fields network\n')
            sys.exit(10)

        if verbose:
            print('Filtering network by databases...\n')

        fnet=nx.Graph()

        for u,v,d in self.network.edges_iter(data=True):
            skip_inc = False
            skip_exc = False

            # Check if at least one of the databases is included in the provided
            # list
            for database in d['sources']:
                if database in databases_included:
                    fnet.add_edge(u,v,d)
                    break

        write_network_file_from_networkx_graph(fnet, output_network_file, output_nodes_file, self.network_format, tissue_specific=self._tissue_specific)

        if verbose:
            print('Filtering network by databases... finished!\n')

        return Network(output_network_file, output_nodes_file, self.type_id, self.network_format)

    def get_number_of_housekeeping_genes(self, housekeeping_genes):
        """
        Returns the number of housekeeping genes in the network

        @param:    housekeeping_genes
        @pdef:     Set containing the dataset of housekeeping genes
        @ptype:    Set
        """
        housekeeping_genes = set([ str(x) for x in housekeeping_genes ])
        return len( set(self.network.nodes()) & housekeeping_genes )


class TissueSpecificNetwork(Network):
    """ 
    Child class of Network defining a tissue-specific network object 
    """

    def __init__(self, network_file, type_id, network_format, tissue_object, permission):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    type_id
        @pdef:     Type of IDs in the network
        @ptype:    {String}

        @param:    network_format
        @pdef:     Format of the network
        @ptype:    {String} {'sif' or 'multi-fields'}

        @param:    tissue_object
        @pdef:     Instance from Tissue Class associated to the Network
        @ptype:    {String}

        @param:    permission
        @pdef:     Level of permissivity of the network filtering
        @ptype:    {String}

        @raises: {IncorrectNetworkFormat} if the network_format is not in
        self.formats.
        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        #super(TissueSpecificNetwork, self).__init__(network_file, type_id, network_format)

        self.network_file = network_file
        self.type_id = type_id
        self.network_format = network_format
        self.formats = ['sif', 'multi-fields']
        self.pickles_path = '/home/quim/project/tissue_specificity/scripts/pickles'

        self.tissue_object = tissue_object
        self.permission = permission

        self._tissue_specific = True 

        self.network = self.parse_network(tissue_specific=self._tissue_specific)

        self.hpa_edges = self.get_hpa_edges()
        self.jensen_edges = self.get_jensen_edges()

    ###########
    # METHODS #
    ###########

    def get_hpa_edges(self):
        """
        Obtain the edges in the tissue-specific network according to Human
        Protein Atlas 
        """
        hpa=nx.Graph()
        for u,v,d in self.network.edges_iter(data=True):
            if 'hpa' in d['tissue_db']:
                hpa.add_edge(u,v,d)
        return hpa.edges()

    def get_jensen_edges(self):
        """
        Obtain the edges in the tissue-specific network according to Tissues
        (Jensen Lab)
        """
        jensen=nx.Graph()
        for u,v,d in self.network.edges_iter(data=True):
            if 'jensen' in d['tissue_db']:
                jensen.add_edge(u,v,d)
        return jensen.edges()

    def get_union(self):
        """
        Obtain the common tissue-specific interactions between Human Protein
        Atlas and Tissues (Jensen Lab)
        """
        union=nx.Graph()
        for u,v,d in self.network.edges_iter(data=True):
            if 'jensen' in d['tissue_db'] or 'hpa' in d['tissue_db']:
                union.add_edge(u,v,d)
        return union.edges()

    def get_intersection(self):
        """
        Obtain the tissue-specific interactions in Human Protein Atlas, Tissues
        (Jensen Lab) or both
        """
        intersection=nx.Graph()
        for u,v,d in self.network.edges_iter(data=True):
            if 'jensen' in d['tissue_db'] and 'hpa' in d['tissue_db']:
                intersection.add_edge(u,v,d)
        return intersection.edges()


class GUILDProfile(object):
    """ 
    Class defining a GUILD profile object 
    """

    def __init__(self, scores_file, type_id='geneid', top=100, top_type='percentage'):
        """ 
        @param:    scores_file
        @pdef:     File resulting from running GUILD, which contains the network scored
        @ptype:    {String}

        @param:    type_id
        @pdef:     Type of IDs of the nodes in the pvalue_file (i.e. geneid, biana...)
        @ptype:    {String}

        @param:    top
        @pdef:     Percentage of the top-scoring nodes to parse (100, 10...)
        @ptype:    {String}

        @param:    top_type
        @pdef:     Type of top threshold (it can be either 'percentage' or 'number_of_nodes')
        @ptype:    {String}

        @raises: {IncorrectThresholdType} if a wrong type of threshold is introduced
        """

        self.scores_file = scores_file
        self.type_id = type_id
        self.top = float(top)
        self.possible_threshold_types = ['percentage', 'number_of_nodes']
        if top_type not in self.possible_threshold_types:
            raise TOP.IncorrectThresholdType(wrong_top_type=top_type, top_types=self.possible_threshold_types)
        self.top_type = top_type
        self.node_to_score = self.parse_scores_file()


    ###########
    # METHODS #
    ###########

    def parse_scores_file(self):
        """
        Obtains the score of every node in the GUILD scores file
        """
        node_to_score = {}
        with open(self.scores_file, 'r') as f:
            for line in f:
                if line[0] == '#' or line.strip().split('\t')[0] == 'Id':
                    continue
                node, score = line.strip().split('\t')
                node_to_score[node] = float(score)
        return node_to_score


    def create_node_profile(self, threshold, threshold_type, output_file=None):
        """
        Select the most relevant nodes of the pvalue_profile using a threshold
        to select the top scoring nodes. Returns a new GUILDProfile containing the
        selected nodes.
        """
        TOP.node_top_scoring(self.node_to_score, threshold, threshold_type, output_file)

        return GUILDProfile(output_file, self.type_id, threshold, threshold_type)


    def create_functional_profile(self, type_correction, output_file, associations_file):
        """
        Perform a functional enrichment analysis of the nodes of the top scoring profile,
        creating a functional profile.
        It requires all the nodes of the network as well, which will be used as background genes.
        """
        TOP.functional_top_scoring(top_geneids=self.node_to_score.keys(), type_correction=type_correction, output_file=output_file, associations_file=associations_file)

        return FunctionalProfile(output_file, self.top, self.scores_file)




class EdgeProfile(Network):
    """ 
    Class defining a GUILD profile object 
    """

    def __init__(self, network_file, type_id, network_format, top, top_type='percentage'):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    type_id
        @pdef:     Type of IDs in the network
        @ptype:    {String}

        @param:    network_format
        @pdef:     Format of the network
        @ptype:    {String} {'sif' or 'multi-fields'}

        @param:    top
        @pdef:     Percentage of the nodes with respect to the initial GUILD file (100, 10...)
        @ptype:    {String}

        @param:    top_type
        @pdef:     Type of top threshold (it can be either 'percentage' or 'number_of_nodes')
        @ptype:    {String}

        @raises: {IncorrectNetworkFormat} if the network_format is not in
        self.formats.
        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        @raises: {IncorrectThresholdType} if a wrong type of threshold is introduced
        """

        self.network_file = network_file
        self.type_id = type_id
        self.network_format = network_format
        self.formats = ['sif', 'multi-fields', 'raw']

        self.top = float(top)
        self.possible_threshold_types = ['percentage', 'number_of_nodes']
        if top_type not in self.possible_threshold_types:
            raise TOP.IncorrectThresholdType(wrong_top_type=top_type, top_types=self.possible_threshold_types)
        self.top_type = top_type

        self._tissue_specific = False

        self.network = self.parse_network(tissue_specific=self._tissue_specific)
        self.edge_to_score = self.parse_edge_profile()


    ###########
    # METHODS #
    ###########

    def parse_edge_profile(self):
        """
           Obtains from the edges profile a dict containing:
           - Key = frozenset([node1, node2])
           - Value = (Score, "useless")
        """
        edge_to_score = {}
        with open(self.network_file,'r') as f:
            for line in f:
                node1, score, node2 = line.strip().split('\t')
                edge = frozenset([node1, node2])
                edge_to_score[edge] = float(score)
        return edge_to_score


    def create_edge_profile(self, node_to_score, threshold, threshold_type, output_file):
        """
        Selects the most relevant nodes of the network selecting the top scoring ones.
        Generates a subnetwork with them.
        Scores the edges calculating the mean of the nodes.
        Returns a new EdgeProfile with the most relevant edges.
        """
        TOP.edge_top_scoring(self.network_file, node_to_score, threshold, threshold_type, output_file)

        return EdgeProfile(network_file=output_file, type_id=self.type_id, network_format=self.network_format, top=threshold, top_type=threshold_type)




class FunctionalProfile(object):
    """ 
    Class defining a GUILD profile object 
    """

    def __init__(self, functional_file, top, node_file):
        """ 
        @param:    functional_file
        @pdef:     File resulting from the functional enrichment analysis, which contains the enriched functions
        @ptype:    {String}

        @param:    top
        @pdef:     Percentage of the nodes with respect to the initial GUILD file (100, 10...)
        @ptype:    {String}

        @param:    node_file
        @pdef:     Node profile file from which the functional enrichment has been done
        @ptype:    {String}

        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        self.functional_file = functional_file
        self.top = top
        self.node_file = node_file

        self.term_id_to_values = self.parse_functional_profile()


    ###########
    # METHODS #
    ###########

    def parse_functional_profile(self):
        """
        Parses a functional profile obtained from an functional enrichment analysis of a list of nodes.
        Returns a dictionary with:
        - Key = GO/Reactome term id
        - Value = Set --> (Log of odds ratio, Adjusted p-value)
        """
        term_id_to_values = {}
        with open(self.functional_file, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                term_id, term_name, num_genes, total_genes, pval, adj_pval = line.strip().split('\t')
                term_id_to_values[term_id] = (float(num_genes)/float(total_genes), float(adj_pval))
        return term_id_to_values




class Tissue(object):
    """ Class defining a tissue object """

    def __init__(self, tissue_terms_hpa, tissue_terms_jensen, jensen_conf=3, hpa_level='medium', hpa_rel='approved', verbose=False):
        """ 
        @param:    tissue_terms_hpa
        @pdef:     Tissue terms of interest in Human Protein Atlas
        @ptype:    {List or String if only one term}

        @param:    tissue_terms_jensen
        @pdef:     Brenda Tissue Ontology names of interest in Tissues (Jensen Lab)
        @ptype:    {List or String if only one term}

        @param:    jensen_conf
        @pdef:     Level of confidence cut-off in Tissues (Jensen Lab) proteins
        @pdefault: 3
        @ptype:    {Integer or Float} {from 0 to 4}

        @param:    hpa_level
        @pdef:     Expression level cut-off in Human Protein Atlas proteins
        @pdefault: 'medium'
        @ptype:    {String} {'not detected','low','medium','high'}

        @param:    hpa_rel
        @pdef:     Reliability cut-off in Human Protein Atlas proteins
        @pdefault: 'approved'
        @ptype:    {String} {'uncertain','approved','supported'}

        @param:    verbose
        @pdef:     If True, informs about when the parsing starts and ends
        @pdefault: False
        @ptype:    Bool
        """

        self.tissue_terms_hpa = self.check_tissue_terms(tissue_terms_hpa)
        self.tissue_terms_jensen = self.check_tissue_terms(tissue_terms_jensen)
        self.jensen_conf = jensen_conf
        self.hpa_level = hpa_level
        self.hpa_rel = hpa_rel
        self.pickles_path = '/home/quim/project/tissue_specificity/scripts/pickles'

        # The scales of Human Protein Atlas levels and reliability:
        self.hpa_scale_level = {
        'not detected' : 0,
        'low' : 1,
        'medium' : 2,
        'high' : 3
        }
        self.hpa_scale_rel = {
        'uncertain' : 0,
        'approved' : 1,
        'supported' : 2
        }

        # Mapping files from tissue terms to biana user entities of tissues

        if verbose:
            print('Generating tissue...\n')

        BTOname_file = os.path.join(self.pickles_path,'BTOname2uE.pcl')
        HPA_tissue_file = os.path.join(self.pickles_path,'tissue2uEs.pcl')
        self.BTOname2uE = pickle.load(open(BTOname_file, 'rb'))
        self.tissue2uEs = pickle.load(open(HPA_tissue_file, 'rb'))

        self.user_entities_hpa = self.get_user_entities_hpa()
        self.user_entities_jensen = self.get_user_entities_jensen()
        self.all_tissue_user_entities = self.get_all_tissue_user_entities()

        # Mapping files from user entities of proteins to user entities of
        # tissues
        prot2tissues_file = os.path.join(self.pickles_path,'UEprot2UETissues.pcl')
        prot2HPA_file = os.path.join(self.pickles_path,'UEprot2UEHPA.pcl')
        self.UEprot2UETissues = pickle.load(open(prot2tissues_file, 'rb'))
        self.UEprot2UEHPA = pickle.load(open(prot2HPA_file, 'rb'))

        # Get all the tissue-specific proteins
        self.proteins_hpa = self.get_tissue_specific_proteins_hpa()
        self.proteins_jensen = self.get_tissue_specific_proteins_jensen()
        self.all_tissue_proteins = self.proteins_hpa | self.proteins_jensen

        if verbose:
            print('Generating tissue... finished!\n')

    ###########
    # METHODS #
    ###########

    def check_tissue_terms(self, tissue_terms):
        """ If a string is introduced, it is added in a list object """
        if isinstance(tissue_terms, list):
            return tissue_terms
        elif isinstance(tissue_terms, str):
            termlist = []
            termlist.append(tissue_terms)
            return termlist
        else:
            print('Introduce a list in the tissue terms parameters!\n')
            raise ValueError

    def get_user_entities_hpa(self): 
        """ Returns user entities from Human Protein Atlas """
        self.user_entities_hpa = set([])
        for term in self.tissue_terms_hpa:
            if term in self.tissue2uEs:
                for uE in self.tissue2uEs[term]:
                    self.user_entities_hpa.add(uE)
        return self.user_entities_hpa

    def get_user_entities_jensen(self): 
        """ Returns user entities from Tissues (Jensen lab) """
        self.user_entities_jensen = set([])
        for term in self.tissue_terms_jensen:
            if term in self.BTOname2uE:
                self.user_entities_jensen.add(self.BTOname2uE[term])
        return self.user_entities_jensen

    def get_all_tissue_user_entities(self): 
        """ Returns user entities from Tissues (Jensen lab) """
        return self.user_entities_hpa | self.user_entities_jensen

    def get_tissue_specific_proteins_hpa(self): 
        """ 
        Returns all the user entities from Human Protein Atlas expressed in the
        tissues of interest above a certain level of expression and reliability
        """
        proteins_hpa = set()
        prot_hpa_to_values = {}
        for uE_protein in self.UEprot2UEHPA:
            for uE_tissue in self.UEprot2UEHPA[uE_protein]:
                if uE_tissue in self.user_entities_hpa:
                    level = self.UEprot2UEHPA[uE_protein][uE_tissue]['level']
                    reliability = self.UEprot2UEHPA[uE_protein][uE_tissue]['reliability']

                    # If our reliability is bigger or equal than the cut-off,
                    # we continue
                    if self.hpa_scale_rel[reliability] >= self.hpa_scale_rel[self.hpa_rel]:

                        # If our level of expression is bigger or equal than
                        # the cut-off, we continue
                        if self.hpa_scale_level[level] >= self.hpa_scale_level[self.hpa_level]:
                            proteins_hpa.add(uE_protein)
                            prot_hpa_to_values.setdefault(uE_protein, {})
                            prot_hpa_to_values[uE_protein].setdefault(uE_tissue, {})
                            prot_hpa_to_values[uE_protein][uE_tissue]['level'] = level
                            prot_hpa_to_values[uE_protein][uE_tissue]['reliability'] = reliability

        return proteins_hpa

    def get_tissue_specific_proteins_jensen(self): 
        """ 
        Returns all the user entities from Tissues (Jensen lab) expressed in the
        tissues of interest above a certain level of confidence
        """
        proteins_jensen = set()
        for uE_protein in self.UEprot2UETissues:
            for uE_tissue in self.UEprot2UETissues[uE_protein]:
                if uE_tissue in self.user_entities_jensen:
                    conf = self.UEprot2UETissues[uE_protein][uE_tissue]['confidence']
                    src = self.UEprot2UETissues[uE_protein][uE_tissue]['source']
                    evidence = self.UEprot2UETissues[uE_protein][uE_tissue]['evidence']

                    # If our confidence is bigger or equal than the cut-off, we
                    # continue
                    if float(conf) >= float(self.jensen_conf):
                        proteins_jensen.add(uE_protein)

        return proteins_jensen


class IncorrectNetworkFormat(NameError):
    """
    Subclass exception of the NameError which raises when a network format is
    not provided correctly
    """
    def __init__(self, network_format, formats):
        self.network_format = network_format
        self.formats = formats

    def __str__(self):
        return 'The network format {} is not valid. It must be one of the following formats: {}'.format(self.network_format, self.formats)

class IncorrectTypeID(Exception):
    """
    Exception that raises when a translation is done without having biana codes
    as initial type of IDs
    """
    def __init__(self, type_id):
        self.type_id = type_id

    def __str__(self):
        return 'The initial type of IDs of the network is not biana, it is {}. It is only possible to translate from BIANA codes'.format(self.type_id)


###################
#### FUNCTIONS ####
###################

def calculate_contingency_table(contingency_table):
    """Gets a complete network and filters by tissue interactions in the same tissue"""
    chi2, pval, dof, expected = scipy.stats.chi2_contingency(contingency_table)
    return chi2, pval, dof, expected

def write_network_file_from_networkx_graph(input_network, output_network_file, output_nodes_file, output_network_format, tissue_specific=False):
    """
    Writes a network file and a nodes file from a networkx graph object.
    (Currently only available for multi-fields networks)
    """

    output_network_fd = open(output_network_file, 'w')

    if output_network_format == 'multi-fields':
        for u,v,d in input_network.edges_iter(data=True):
            sources = ';'.join(d['sources'])
            method_ids = ';'.join(d['method_ids'])
            method_names = ';'.join(d['method_names'])
            pmids = ';'.join(d['pmids'])
            if tissue_specific:
                tissue_db = ';'.join(d['tissue_db'])
                tissue_additional = ';'.join(d['tissue_additional'])
                output_network_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format( u,v,sources,method_ids,method_names,pmids,tissue_db,tissue_additional ))
            else:
                output_network_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format( u,v,sources,method_ids,method_names,pmids ))
    elif output_network_format == 'raw':
        for u,v in input_network.edges_iter():
            output_network_fd.write('{}\t{}\n'.format( u,v ))
    else:
        print('Not available format\n')
        sys.exit(10)

    output_network_fd.close()
    output_nodes_fd = open(output_nodes_file, 'w')

    for node in input_network.nodes():
        output_nodes_fd.write('{}\n'.format(node))

    output_nodes_fd.close()

    return

def get_nodes_intersection_of_two_networks(network1, network2):
    """Gets the intersecting nodes of two networks"""
    network1_nodes = set(network1.get_nodes())
    network2_nodes = set(network2.get_nodes())
    return network1_nodes & network2_nodes

# def get_edges_intersection_of_two_networks(network1, network2):
#     """Gets the intersecting edges of two networks"""
#     intersection = set()
#     network2_edges = network2.get_edges()
#     for u,v in network1.network.edges_iter():
#         comb1 = (u,v)
#         comb2 = (v,u)
#         if comb1 in network2_edges or comb2 in network2_edges:
#             if comb1 in network2_edges:
#                 intersection.add(comb1)
#             else:
#                 intersection.add(comb2)
#     return intersection
def get_edges_intersection_of_two_networks(network1, network2):
    """Gets the intersecting edges of two networks"""
    network1_edges = network1.get_edges()
    network2_edges = network2.get_edges()

    network1_set = set(frozenset(i) for i in network1_edges)
    network2_set = set(frozenset(i) for i in network2_edges)

    return network1_set & network2_set