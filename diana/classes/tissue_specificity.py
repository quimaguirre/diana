import argparse
import cPickle
import mysql.connector
import os, sys

"""
    NetworkAnalysis
    2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
    Universitat Pompeu Fabra
"""

#################
#### CLASSES ####
#################

class HouseKeepingGenes(object):
    """ Class defining a HouseKeepingGenes object """

    def __init__(self, type_id):
        """ 
        @param:    type_id
        @pdef:     Type of IDs in the network
        @ptype:    {String}
        """

        self.pickles_path = '/home/quim/project/tissue_specificity/scripts/pickles'
        self.type_id = type_id

        self.hpa_hk_genes = self.get_hpa_hk_genes()
        self.eisenberg_hk_genes = self.get_eisenberg_hk_genes()
        self.all_hk_genes = self.get_all_hk_genes()


    def get_hpa_hk_genes(self):
        """
        Obtains a set containing the Human Protein Atlas house keeping genes
        """

        if self.type_id == 'biana':
            hpa_uE_dump = os.path.join(self.pickles_path,'hpa_hk_uEs.pcl')
            hpa_hk_genes = cPickle.load(open(hpa_uE_dump))
        elif self.type_id == 'geneid':
            hpa_geneid_dump = os.path.join(self.pickles_path,'hpa_hk_geneIDs.pcl')
            hpa_hk_genes = cPickle.load(open(hpa_geneid_dump))
        else:
            print('The input_type must be \'geneid\' or \'biana\'\n')
            sys.exit(10)

        return hpa_hk_genes

    def get_eisenberg_hk_genes(self):
        """
        Obtains a set containing the Eisenberg data of house keeping genes
        """

        if self.type_id == 'biana':
            elieis_uE_dump = os.path.join(self.pickles_path,'eisenberg_hk_uEs.pcl')
            eisenberg_hk_genes = cPickle.load(open(elieis_uE_dump))
        elif self.type_id == 'geneid':
            elieis_geneid_dump = os.path.join(self.pickles_path,'eisenberg_hk_geneIDs.pcl')
            eisenberg_hk_genes = cPickle.load(open(elieis_geneid_dump))
        else:
            print('The input_type must be \'geneid\' or \'biana\'\n')
            sys.exit(10)

        return eisenberg_hk_genes

    def get_all_hk_genes(self):
        """
        Obtains a set containing all the house keeping genes in Human Protein
        Atlas and Eisenberg datasets
        """
        return self.hpa_hk_genes | self.eisenberg_hk_genes


###################
#### FUNCTIONS ####
###################


def filter_network_tissue_specific(input_network, tissue, permission, output_network, output_nodes):
    """Gets a complete network and filters by tissue interactions in the same tissue"""

    permission = int(permission)

    # Check that the permission parameter is well introduced
    if permission != 0 and permission != 1 and permission != 2:
        print('In the parameter --permission, introduce 0, 1 or 2\n')
        sys.exit(10)

    # Open files
    network_fd = open(input_network, 'r')
    output_fd = open(output_network, 'w')
    output_nodes_fd = open(output_nodes, 'w')

    ue2hpainfo = {}
    ue2jenseninfo = {}
    tissue_nodes = set()

    # Process the network
    for line in network_fd:

        fields = line.strip().split('\t')
        uE1 = int(fields[0])
        uE2 = int(fields[1])
        other = fields[2:len(fields)]
        other = '\t'.join(other)
        #print('\nChecking interaction {} and {}...\n'.format(uE1, uE2))

        for uE in [uE1, uE2]:

            if uE not in ue2jenseninfo:

                ue2jenseninfo.setdefault(uE, {'info':False, 'specific':False})

                # Report if there is info about the protein in Jensen or not
                if uE not in tissue.UEprot2UETissues:
                    ue2jenseninfo[uE]['info'] = False
                else:
                    ue2jenseninfo[uE]['info'] = True

                    # Report if there is specificity in Jensen or not
                    if uE in tissue.proteins_jensen:
                        ue2jenseninfo[uE]['specific'] = True
                    else:
                        ue2jenseninfo[uE]['specific'] = False

            if uE not in ue2hpainfo:

                ue2hpainfo.setdefault(uE, {'info':False, 'specific':False})

                # Report if there is info about the protein in HPA or not
                if uE not in tissue.UEprot2UEHPA:
                    ue2hpainfo[uE]['info'] = False
                else:
                    ue2hpainfo[uE]['info'] = True

                    # Report if there is specificity in HPA or not
                    if uE in tissue.proteins_hpa:
                        ue2hpainfo[uE]['specific'] = True
                    else:
                        ue2hpainfo[uE]['specific'] = False


        ############ Order the information obtained from Tissues and HPA ############

        result_hpa = { 'info' : False, 'specific' : False, 'db' : 'hpa' }
        result_jensen = { 'info' : False, 'specific' : False, 'db' : 'jensen' }

        for database, result in [ [ue2jenseninfo, result_jensen], [ue2hpainfo, result_hpa] ]:

            # If we have info in both proteins, it is True
            if database[uE1]['info'] == True and database[uE1]['info'] == True:
                result['info'] = True

            # If we only have info about one of the proteins...
            if (database[uE1]['info'] == True and database[uE2]['info'] == False) or (database[uE1]['info'] == False and database[uE2]['info'] == True):
                for uE in [uE1, uE2]:
                    if database[uE]['info'] == True: # Identify the protein in which we have info
                        if database[uE]['specific'] == True: # If the protein is tissue-specific, we say that the interaction is partially tissue-specific
                            result['info'] = 'partial'

            # If one of the proteins is specific and the other not, or both are not, we consider the interaction as not specific!
            if (database[uE1]['specific'] == True and database[uE2]['specific'] == False) or (database[uE1]['specific'] == False and database[uE2]['specific'] == True) or (database[uE1]['specific'] == False and database[uE2]['specific'] == False):
                result['specific'] = False
            # If both are specific, then it is specific!
            elif database[uE1]['specific'] == True and database[uE2]['specific'] == True:
                result['specific'] = True


        #print('\nHPA info: {}\tHPA specificity: {}'.format(result_hpa['info'], result_hpa['specific']))
        #print('JENSEN info: {}\tJENSEN specificity: {}'.format(result_jensen['info'], result_jensen['specific']))


        ############ Decide if they are tissue specific or not... ############

        # If both specificity results are True, we consider tissue-specificity

        if result_hpa['specific'] == True and result_jensen['specific'] == True:

            #print('\n... tissue specific!\n')

            database = ';'.join([result_hpa['db'], result_jensen['db']])
            additional = '-'
            output_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(uE1, uE2, other, database, additional))
            tissue_nodes.add(uE1)
            tissue_nodes.add(uE2)


        # If there is True info in only one of the databases and it is tissue-specific, we consider tissue-specificity

        if (result_hpa['info'] == True and result_jensen['info'] == False) or (result_hpa['info'] == False and result_jensen['info'] == True):

            for result in [result_jensen, result_hpa]:

                if result['info'] == True and result['specific'] == True:

                    #print('\n... tissue specific!\n')

                    database = result['db']
                    additional = '-'
                    output_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(uE1, uE2, other, database, additional))
                    tissue_nodes.add(uE1)
                    tissue_nodes.add(uE2)


        # If there is contradiction (one of them is True and the other is False), we add the one which is True

        if (result_hpa['specific'] == True and result_jensen['info'] == True and result_jensen['specific'] == False) or (result_hpa['info'] == True and result_hpa['specific'] == False and result_jensen['specific'] == True):

            #print('\n... contradiction!\n')

            for result in [result_jensen, result_hpa]:

                if result['info'] == True and result['specific'] == True:

                    database = result['db']
                    additional = 'contradiction'
                    output_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(uE1, uE2, other, database, additional))
                    tissue_nodes.add(uE1)
                    tissue_nodes.add(uE2)


        # If there is no info, we check the permissivity

        if result_hpa['info'] == False and result_jensen['info'] == False:

            #print('\n... no info!\n')

            # If the level of permissivity is 2, we include it
            if permission == 2:
                #print('\n... as permission is level {}, we include it!\n'.format(permission))
                database = '-'
                additional = 'no info'
                output_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(uE1, uE2, other, database, additional))
                tissue_nodes.add(uE1)
                tissue_nodes.add(uE2)
            else:
                continue


        # If there is only partial info in one of the databases or both, we check the level of permission

        if (result_hpa['info'] == 'partial' and result_jensen['info'] == False) or (result_hpa['info'] == False and result_jensen['info'] == 'partial') or (result_hpa['info'] == 'partial' and result_jensen['info'] == 'partial'):

            # If the permission is medium or high, we include it!
            if permission == 1 or permission == 2:

                database = []
                for result in [result_jensen, result_hpa]:

                    if result['info'] == 'partial':

                        database.append(result['db'])

                #print('\n... only partial info in one of the databases... tissue specific!\n')
                database = ';'.join(database)
                additional = 'partial info'
                output_fd.write('{}\t{}\t{}\t{}\t{}\n'.format(uE1, uE2, other, database, additional))
                tissue_nodes.add(uE1)
                tissue_nodes.add(uE2)

        #print('\n')


    # Print the nodes in the nodes file
    for node in tissue_nodes:
        output_nodes_fd.write('{}\n'.format(node))
        # if nodes[node]['tissue_specific'] == True:
        #     if 'level' not in nodes[node]:
        #         nodes[node]['level'] = '-'
        #         nodes[node]['reliability'] = '-'
        #     if 'confidence' not in nodes[node]:
        #         nodes[node]['confidence'] = '-'
        #     output_nodes_fd.write('{}\t{}\t{}\t{}\n'.format(node, nodes[node]['level'], nodes[node]['reliability'], nodes[node]['confidence']))


    network_fd.close()
    output_fd.close()
    output_nodes_fd.close()

    return


