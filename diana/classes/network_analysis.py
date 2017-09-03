import cPickle
import networkx as nx
import numpy as np
import scipy.stats
import os, sys

import network_translation as NT
import tissue_specificity as TS
import top_scoring as TOP

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

    def __init__(self, network_file, node_file, type_id, network_format):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    node_file
        @pdef:     Path to the file containing the nodes of the network
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
        self.node_file = node_file
        self.type_id = type_id
        self.network_format = network_format
        self.formats = ['sif', 'multi-fields', 'raw']
        self.pickles_path = '/home/quim/project/tissue_specificity/scripts/pickles'

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
        self.psimi2method = cPickle.load(open(psimi2method_file))

        methodid2interactions = {}

        for u,v,d in self.network.edges_iter(data=True):

            for method_id in d['method_ids']:
                method_id = int(method_id)
                methodid2interactions.setdefault(method_id, 0)
                methodid2interactions[method_id] += 1

        if verbose:
            for psimi, interactions in sorted(methodid2interactions.iteritems(), key=lambda (x, y): y, reverse=True):
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
            for method, interactions in sorted(method2interactions.iteritems(), key=lambda (x, y): y, reverse=True):
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
            for num_pmids, interactions in sorted(numpmids2interactions.iteritems(), key=lambda (x, y): x):
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
            for source, interactions in sorted(database2interactions.iteritems(), key=lambda (x, y): y, reverse=True):
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
                (node1, score, node2) = line.strip().split('\t')
                G.add_edge(node1,node2,score=score)

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


    def translate_network(self, translation_file, translation_id, translation_format, translated_network, translated_nodes):
        """
        Uses a BIANA translation file to return a new Network object with the 
        desired type of IDs and format

        @param:    translation_file
        @pdef:     File containing the translations from biana codes to the translation id
        @ptype:    String

        @param:    translation_id
        @pdef:     Type of ID in which the network will be translated
        @ptype:    String

        @param:    translation_format
        @pdef:     Format of the final network
        @ptype:    String

        @param:    translated_network
        @pdef:     File where the translated network will be written
        @ptype:    String

        @param:    translated_nodes
        @pdef:     File where the nodes file will be written
        @ptype:    String
        """

        if self.type_id != 'biana':
            raise IncorrectTypeID(self.type_id)

        NT.translate(self.network_file, self.node_file, translation_file, self.network_format, translation_format, translated_network, translated_nodes)

        return Network(translated_network, translated_nodes, translation_id, translation_format)


    def score_network(self, node_to_values, output_file):
        """
        Uses a BIANA translation file to return a new Network object with the 
        desired type of IDs and format

        @param:    node_to_values
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

        TOP.score_network(self.network_file, node_to_values, output_file)

        return EdgeProfile(output_file, None, self.type_id, self.network_format, 100)


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

    def __init__(self, network_file, node_file, type_id, network_format, tissue_object, permission):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    node_file
        @pdef:     Path to the file containing the nodes of the network
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

        #super(TissueSpecificNetwork, self).__init__(network_file, node_file, type_id, network_format)

        self.network_file = network_file
        self.node_file = node_file
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

    def translate_network(self, translation_file, translation_id, translation_format, translated_network, translated_nodes, verbose=False):
        """
        Uses a BIANA translation file to return a new Network object with the 
        desired type of IDs and format
        """
        if verbose:
            print('Translating network to {}...\n'.format(translation_id))

        if self.type_id != 'biana':
            raise IncorrectTypeID(self.type_id)

        NT.translate(self.network_file, self.node_file, translation_file, self.network_format, translation_format, translated_network, translated_nodes)

        if verbose:
            print('Translating network to {}... finished!\n'.format(translation_id))

        return TissueSpecificNetwork(translated_network, translated_nodes, translation_id, translation_format, self.tissue_object, self.permission)

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

    def __init__(self, pvalue_file, type_id, top):
        """ 
        @param:    pvalue_file
        @pdef:     File resulting from running GUILD, which contains the network scored
        @ptype:    {String}

        @param:    type_id
        @pdef:     Type of IDs of the nodes in the pvalue_file (i.e. geneid, biana...)
        @ptype:    {String}

        @param:    top
        @pdef:     Percentage of the nodes with respect to the initial GUILD file (100, 10...)
        @ptype:    {String}

        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        self.pvalue_file = pvalue_file
        self.type_id = type_id
        self.top = float(top)

        self.node_to_values = self.parse_pvalue_file()


    ###########
    # METHODS #
    ###########

    def parse_pvalue_file(self):
        """
        Obtains the score and p-value for every node in the p-value file
        """
        node_to_values = {}
        f = open(self.pvalue_file, 'r')
        for line in f:
            if line[0] == '#':
                continue
            node, score, pval = line.strip().split(' ')
            node_to_values[node] = (float(score), float(pval))
        f.close()
        return node_to_values


    def translate_pvalue_file(self, translation_file, translation_id, output_file, verbose=False):
        """
        Uses a BIANA translation file to return a new GUILDProfile object with the 
        desired type of IDs and format
        """
        if verbose:
            print('Translating network to {}...\n'.format(translation_id))

        if self.type_id != 'biana':
            raise IncorrectTypeID(self.type_id)

        NT.translate(None, self.pvalue_file, translation_file, 'guild', 'guild', None, output_file)

        if verbose:
            print('Translating network to {}... finished!\n'.format(translation_id))

        return GUILDProfile(output_file, translation_id, self.top)


    def create_node_profile(self, threshold, output_file):
        """
        Select the most relevant nodes of the pvalue_profile using a threshold
        to select the top scoring nodes. Returns a new GUILDProfile containing the
        selected nodes.
        """
        TOP.node_top_scoring(self.node_to_values, threshold, output_file)

        return GUILDProfile(output_file, self.type_id, threshold)


    def create_functional_profile(self, obodag, geneid2gos, all_node_to_vals, output_file):
        """
        Perform a functional enrichment analysis of the nodes of the top scoring profile,
        creating a functional profile.
        It requires all the nodes of the network as well, which will be used as background genes.
        """
        TOP.functional_top_scoring(obodag, geneid2gos, all_node_to_vals.keys(), self.node_to_values.keys(), output_file)

        return FunctionalProfile(output_file, self.top, self.pvalue_file)




class EdgeProfile(Network):
    """ 
    Class defining a GUILD profile object 
    """

    def __init__(self, network_file, node_file, type_id, network_format, top):
        """ 
        @param:    network_file
        @pdef:     Path to the file containing the edges of the network
        @ptype:    {String}

        @param:    node_file
        @pdef:     Path to the file containing the nodes of the network
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

        @raises: {IncorrectNetworkFormat} if the network_format is not in
        self.formats.
        @raises: {IncorrectTypeID} if the method translate_network is used with
        a network of type_id different from 'biana'
        """

        self.network_file = network_file
        self.node_file = node_file
        self.type_id = type_id
        self.network_format = network_format
        self.formats = ['sif', 'multi-fields', 'raw']

        self.top = float(top)

        self._tissue_specific = False

        self.network = self.parse_network(tissue_specific=self._tissue_specific)


    ###########
    # METHODS #
    ###########

    def translate_network(self, translation_file, translation_id, output_file, verbose=False):
        """
        Uses a BIANA translation file to return a new Network object with the 
        desired type of IDs and format
        """
        if verbose:
            print('Translating network to {}...\n'.format(translation_id))

        if self.type_id != 'biana':
            raise IncorrectTypeID(self.type_id)

        NT.translate(self.network_file, None, translation_file, self.network_format, self.network_format, output_file, None)

        if verbose:
            print('Translating network to {}... finished!\n'.format(translation_id))

        return EdgeProfile(output_file, None, translation_id, self.network_format, self.top)


    def create_edge_profile(self, node_to_values, threshold, output_file):
        """
        Selects the most relevant nodes of the network selecting the top scoring ones.
        Generates a subnetwork with them.
        Scores the edges calculating the mean of the nodes.
        Returns a new EdgeProfile with the most relevant edges.
        """
        TOP.edge_top_scoring(self.network_file, node_to_values, threshold, output_file)

        return EdgeProfile(output_file, None, self.type_id, self.network_format, threshold)




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

        self.go_id_to_values = self.parse_functional_profile()


    ###########
    # METHODS #
    ###########

    def parse_functional_profile(self):
        """
        Obtains the score and p-value for every node in the p-value file
        """
        go_id_to_values = {}
        f = open(self.functional_file, 'r')
        for line in f:
            if line[0] == '#':
                continue
            num_genes, total_genes, log_odds_ratio, pval, adj_pval, go_id, go_name, go_type = line.strip().split('\t')
            go_id_to_values[go_id] = (log_odds_ratio, adj_pval, go_name, go_type)
        f.close()
        return go_id_to_values




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
        self.BTOname2uE = cPickle.load(open(BTOname_file))
        self.tissue2uEs = cPickle.load(open(HPA_tissue_file))

        self.user_entities_hpa = self.get_user_entities_hpa()
        self.user_entities_jensen = self.get_user_entities_jensen()
        self.all_tissue_user_entities = self.get_all_tissue_user_entities()

        # Mapping files from user entities of proteins to user entities of
        # tissues
        prot2tissues_file = os.path.join(self.pickles_path,'UEprot2UETissues.pcl')
        prot2HPA_file = os.path.join(self.pickles_path,'UEprot2UEHPA.pcl')
        self.UEprot2UETissues = cPickle.load(open(prot2tissues_file))
        self.UEprot2UEHPA = cPickle.load(open(prot2HPA_file))

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