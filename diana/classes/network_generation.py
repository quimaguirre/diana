import sys, os, re
import biana
try: from biana import *
except: sys.exit(10)


"""
    NetworkAnalysis
    2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
    Universitat Pompeu Fabra
"""

def generate_network(targets_list, targets_type_id, radius, taxid, translation_file, translation_type_id, node_file, edge_file,
                     restricted_to_TAP = False, restricted_to_Y2H = False, restricted_to_user = None,
                     except_TAP = False, except_Y2H = False, except_user = None,
                     database = 'BIANA_JUN_2017', unification_protocol = 'geneid_seqtax_v1',
                     output_format = 'sif', verbose = False):
    """
    Generates a complete interactome network if the radius = 0.
    Generates a network of expansion if the radius > 0. This means that the center of the network will be the targets,
    and the network will be expanded as many levels as indicated in the parameter radius.

    Parameters:
        @targets_list:          List of targets from which the network of expansion starts
        @targets_type_id:       Type of ID of the targets introduced (i.e. geneid, uniprotaccession, uniprotentry)
        @radius:                Number of levels of expansion of the network. If 0, the complete interactome will be generated
        @taxid:                 Restrict the proteins of the network to a concrete species, indicated by its Taxonomy ID
        @translation_file:      The function will generate a network in BIANA codes, but also a translation file containing the
                                translations from BIANA codes to the desired type of ID
        @translation_type_id:   Indicate the type of ID in which the function will translate the BIANA codes (i.e. geneid,
                                uniprotaccession, uniprotentry...)
        @node_file:             Name of the output file with all the nodes of the network
        @edge_file:             Name of the output file with all the edges (interactions) of the network
        @restricted_to_TAP:     Restrict the interactions of the network to the ones that have been reported by Tandem Affinity
                                Purification (TAP) methods
        @restricted_to_Y2H:     Restrict the interactions of the network to the ones that have been reported by Yeast to Hybrid
                                (Y2H) methods
        @restricted_to_user:    File containing methods IDs separated by newline characters. It will restrict the interactions of
                                the network to the ones reported by the methods indicated in the file
        @except_TAP:            Restrict the interactions of the network to the ones that have been reported by at least one
                                method which is not a Tandem Affinity Purification method
        @except_Y2H:            Restrict the interactions of the network to the ones that have been reported by at least one
                                method which is not a Yeast to Hybrid method
        @except_user:           File containing methods IDs separated by newline characters. Restrict the interactions of the 
                                network to the ones that have been reported by at least one method which is not in the file
        @database:              Name of the BIANA database used.
                                default: 'BIANA_JUN_2017'
        @unification_protocol:  Name of the unification protocol used in the BIANA database selected.
                                default: 'geneid_seqtax_v1'
        @output_format:         Output format of the network. It can be:
                                    'sif': <node1>\tscore\t<node2>
                                    'netscore': <node1>\t<node2>\t<score>
                                    'multi-fields' : <node1>\t<node2>\t<sources>\t<method_ids>\t<method_names>\t<pmids>
                                default: 'sif'
        @verbose:               default: False
    """

    # Parameters that I have decided to fix
    restricted_to_seeds = False
    minimum_number_of_methods = 1
    minimum_number_of_db = 1

    # Restricted by the user:
    if not fileExist(restricted_to_user):
        print "No restriction on methods selected by the user"
        user_selection=False
    else:
        use_methods=[]
        input_method=open(restricted_to_user)
        for line in input_method:
            fields = line.strip().split("\t")
            use_methods.append(fields[0])
        input_method.close()
        user_selection=True
        #print "Input to use only Methods:",repr(use_methods)
    if not fileExist(except_user):
        print "No rejection of methods selected by the user"
        user_rejection=False
    else:
        no_methods=[]
        input_method=open(except_user)
        for line in input_method:
            fields = line.strip().split("\t")
            no_methods.append(fields[0])
        input_method.close()
        user_rejection=True
        #print "Input of rejected Methods:",repr(no_methods)

    session = create_new_session(sessionID="biana_session",dbname=database,dbhost="localhost",
                                 dbuser="quim",dbpassword=None,
                                 unification_protocol=unification_protocol)

    # Create network network of expansion if the radius is larger than 0
    if radius > 0:

        if len(targets_list) < 3:
            print "There are not enough targets to perform the analysis"
            sys.exit(10)
        else:
            level=radius

            proteome = session.create_new_user_entity_set(  identifier_description_list =targets_list,
                                                          attribute_restriction_list=[("taxid",taxid)],
                                                          id_type=targets_type_id,new_user_entity_set_id="proteome",
                                                          negative_attribute_restriction_list=[] )

    # Create complete network or interactome
    else:
        level=0
        proteome = session.create_new_user_entity_set( identifier_description_list = [("taxid",taxid)],
                                                      attribute_restriction_list=[], id_type="embedded",
                                                      new_user_entity_set_id="proteome",
                                                      negative_attribute_restriction_list=[] )


    # Select interactions

    if restricted_to_TAP:
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               relation_attribute_restriction_list = [("Method_id",400)],
                                #relation_attribute_restriction_list  = [("psimi_name","affinity technology")],
                               include_relations_last_level = (not restricted_to_seeds) , use_self_relations = False)
    elif restricted_to_Y2H:
        #print "restricting to y2h"
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               relation_attribute_restriction_list = [("Method_id",18)],
                               #relation_attribute_restriction_list = [("psimi_name","y2h2")],
                               include_relations_last_level = (not restricted_to_seeds) , use_self_relations = False)
    else:
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               include_relations_last_level = (not restricted_to_seeds) , use_self_relations = False)



    # Summary of interactions

    out_network = open(edge_file,'w')
    all_interactions = proteome.getRelations()
    print "Num interactions:", len(all_interactions)

    nodes=set()

    # Get all the user entity ids from the user entity set 'proteome'
    all_uEs = proteome.get_user_entity_ids()
    # Obtain a dictionary user entity ID => type
    uEId_to_type = session.dbAccess.get_user_entity_type(unification_protocol, all_uEs)

    skip_interactions=0
    for (uE_id1, uE_id2) in all_interactions:

    #self.dbAccess.get_external_entities_dict( externalEntityIdsList = [external_entity_relation_id] )

        # Get TYPE of user entity
        uE1_type = uEId_to_type[uE_id1]
        uE2_type = uEId_to_type[uE_id2]
        # If type is not protein, we skip the interaction
        if uE1_type != 'protein' or uE2_type != 'protein':
            if verbose:
                print('Skipping interaction because the type of one of the user entities is not protein!')
                print('Node 1: {}\tType: {}'.format(uE_id1, uE1_type))
                print('Node 2: {}\tType: {}'.format(uE_id2, uE2_type))
            skip_interactions=skip_interactions+1
            continue

        eErIDs_list = proteome.get_external_entity_relation_ids(uE_id1, uE_id2)

        method_names = set()
        method_ids = set()
        source_databases = set()
        use_method_ids=set()
        pubmed_ids = set()


        relationObj_dict = session.dbAccess.get_external_entities_dict(
                                externalEntityIdsList = eErIDs_list, attribute_list = [],
                                relation_attribute_list = ["method_id","psimi_name","pubmed"], participant_attribute_list = [] )

        num_methods=0
        for current_eErID in eErIDs_list:
            relationObj = relationObj_dict[current_eErID]
            if verbose:
                print "Interaction: (",uE_id1,",",uE_id2,")"
                print relationObj

            #if relationObj.get_attribute(attribute_identifier="psimi_name") is not None:
            #    print "\t".join([ x.value for x in relationObj.get_attribute(attribute_identifier="psimi_name") ])
            #if relationObj.get_attribute(attribute_identifier="method_id") is not None:
            #print "\t".join([ x.value for x in relationObj.get_attribute(attribute_identifier="method_id") ])
            #print relationObj.get_attributes_dict()
            #print [ x.value for x in relationObj.get_attributes_dict()["psimi_name"] ]
            #print [ x.value for x in relationObj.get_attributes_dict()["method_id"] ]

            if "psimi_name" in relationObj.get_attributes_dict():
                method_names.update([ str(x.value) for x in relationObj.get_attributes_dict()["psimi_name"] ])
            if "method_id" in relationObj.get_attributes_dict():
                method_ids.update([ x.value for x in relationObj.get_attributes_dict()["method_id"]])
            if "pubmed" in relationObj.get_attributes_dict():
                pubmed_ids.update([ x.value for x in relationObj.get_attributes_dict()["pubmed"]])
            source_databases.add(str(session.dbAccess.get_external_database(
                                    database_id = relationObj.get_source_database()) ))
            if except_TAP:
                affinity = get_affinity_method_ids()
                for m in method_ids:
                    if m not in affinity:
                        use_method_ids.add(m)
                        #print "Add", m
            elif except_Y2H:
                complementation = get_complementation_method_ids()
                #print "check Y2H"
                for m in method_ids:
                    if m not in complementation:
                        use_method_ids.add(m)
                        #print "Add", m
            elif user_rejection:
                for m in method_ids:
                    if m not in no_methods:
                        use_method_ids.add(m)
            elif user_selection:
                for m in method_ids:
                    #print "Check",repr(use_methods)
                    if m in set(use_methods):
                        use_method_ids.add(m)
                    elif verbose:
                        print "Not among selected methods ",m
            else:
                use_method_ids.update(method_ids)

        if len(source_databases) > 0:
            info_sources=";".join([str(x) for x in source_databases])
        else:
            if verbose:
                print('Skipping interaction it has no source database!')
                print('Node 1: {}\tNode 2: {}'.format(uE_id1, uE_id2))
            skip_interactions=skip_interactions+1
            continue
        if len(method_names) > 0:
            info_methods=";".join([str(x) for x in method_names])
        else:
            info_methods='-'
        if len(use_method_ids) > 0:
            info_methods_ids=";".join([str(x) for x in use_method_ids])
        else:
            if verbose:
                print('Skipping interaction it has no method!')
                print('Node 1: {}\tNode 2: {}'.format(uE_id1, uE_id2))
            skip_interactions=skip_interactions+1
            continue
        if len(pubmed_ids) > 0:
            info_pubmed_ids=";".join([str(x) for x in pubmed_ids])
        else:
            info_pubmed_ids='-'
        num_databases=len(source_databases)
        num_methods=len(use_method_ids)
        num_pubmeds = len(pubmed_ids)

        if verbose:
            print "Methods",num_methods,info_methods,"\tSelected:",info_methods_ids
            print "Databases",num_databases,info_sources
            print "Pubmeds",num_pubmeds,info_pubmed_ids

        if num_methods >= minimum_number_of_methods:
            use=True
        else:
            use=False

        if use and num_databases >= minimum_number_of_db:
            use=True
        else:
            use=False

        if not use:
            skip_interactions=skip_interactions+1
        #print method_names, method_ids, source_databases


        # Output edge file

        if use:
            nodes.add(uE_id1)
            nodes.add(uE_id2)

            if verbose:
                if output_format == 'multi-fields' :
                    out_network.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".
                                  format(uE_id1,uE_id2,info_sources,info_methods_ids,info_methods,info_pubmed_ids))
                elif output_format == 'netscore':
                    out_network.write('\t{}\t{}\t{:.2f}\n'.format(uE_id1,uE_id2,1))
                else:
                    out_network.write("{}\t{:.2f}\t{}\n".format(uE_id1,1.,uE_id2))
            else:
                if output_format == 'multi-fields' :
                    out_network.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".
                                  format(uE_id1,uE_id2,info_sources,info_methods_ids,info_methods,info_pubmed_ids))
                elif output_format == 'netscore':
                    out_network.write('\t{}\t{}\t{:.2f}\n'.format(uE_id1,uE_id2,1.))
                else:
                    out_network.write("{}\t{:.2f}\t{}\n".format(uE_id1,1.,uE_id2))

    print "Num neglected interactions:", skip_interactions
    out_network.close()



    # If we wanted the complete interactome, the translation will be done differently
    if radius == 0:

        # Output node file

        out_proteins = open(node_file,'w')
        for protein in nodes:
            if output_format == 'multi-fields':
                out_proteins.write("{0}\t{1:10.2f}\t{2:10.2f}\t{3:10.2f}\n".format(protein,1.,1.,0.1))
            elif output_format == 'netscore':
                out_proteins.write("{0}\t{1:10.2f}\t{2:10.2f}\t{3:10.2f}\n".format(protein,1.,1.,0.1))
            else:
                out_proteins.write("{0}\t{1:10.2f}\n".format(protein,0.1))
        out_proteins.close()


        ################################# TRANSLATION ####################################
        out_translation = open(translation_file,'w')
        for protein in nodes:
            uE = session.get_user_entity(protein)
            translate=set()
            translate_uni=set()
            if translation_type_id == "proteinsequence":
                maxlen=0;
                for current_id in uE.get_attribute(attribute_identifier=translation_type_id):
                    if maxlen < len(current_id.value.get_sequence().upper()):
                        maxlen=len(current_id.value.get_sequence().upper())
                translation=",".join([str(current_id.value.get_sequence().upper()) for current_id in uE.get_attribute(attribute_identifier=translation_type_id) if len(str(current_id.value.get_sequence().upper())) == maxlen ] )
                #print "Translation",protein,translation
                #print("{0}\t'{1}'\n".format(protein,translation))
            else:
                ##### TRANSLATION TO 'translation_type_id'
                for current_id in uE.get_attribute(attribute_identifier=translation_type_id):
                    translate.add(current_id.value.upper())
                translation="','".join(["{0}".format(x) for x in translate])
            out_translation.write("{0}\t'{1}'\n".format(protein,translation))
        out_translation.close()
        ####################################################################################


    # If we wanted a network of expansion, the translation will be done differently
    elif radius > 0:

        seeds=set()
        for seed in targets_list:
            seeds.add(seed.lower())


        # Output node file

        out_proteins = open(node_file,'w')
        translate={}
        for protein in nodes:
            score=0.1
            uE = session.get_user_entity(protein)
            for current_id in uE.get_attribute(attribute_identifier=targets_type_id):
                if current_id.value.lower() in seeds:
                    translate.setdefault(current_id.value.lower(),[])
                    translate[current_id.value.lower()].append(protein)
                    score=1.0
            if output_format == 'multi-fields':
                out_proteins.write("{0}\t{1:10.2f}\t{2:10.2f}\t{3:10.2f}\n".format(protein,1.,1.,score))
            elif output_format == 'netscore':
                out_proteins.write("{0}\t{1:10.2f}\t{2:10.2f}\t{3:10.2f}\n".format(protein,1.,1.,score))
            else:
                out_proteins.write("{0}\t{1:10.2f}\n".format(protein,score))
        out_proteins.close()


        # Get the IDS of single nodes that were not previously found in the network

        single=set()
        for uE_id in proteome.get_unconnected_nodes():
            single.add(uE_id)
        for protein in single:
            uE = session.get_user_entity(protein)
            for current_id in uE.get_attribute(attribute_identifier=targets_type_id):
                if current_id.value.lower() in seeds:
                    translate.setdefault(current_id.value.lower(),[])
                    translate[current_id.value.lower()].append(protein)


        # Get all IDS of SEEDS, defined as "proteome", and check missing codes to be
        # added for translation

        allseed=set()
        for uE_id in proteome.get_user_entity_ids():
            allseed.add(uE_id)
        for protein in allseed:
            if protein not in single and protein not in nodes:
                uE = session.get_user_entity(protein)
                for current_id in uE.get_attribute(attribute_identifier=targets_type_id):
                    if current_id.value.lower() in seeds:
                        translate.setdefault(current_id.value.lower(),[])
                        translate[current_id.value.lower()].append(protein)


        ################################# TRANSLATION ####################################
        out_translation = open("translation_seeds_to_BIANA_codes.txt",'w')
        for s in seeds:
            if s == '': continue
            if s in translate:
                codes=set(translate[s])
                translation="','".join([str(x) for x in codes])
                #out_translation.write("%s\t'%s'\n" % (s.upper(),translation))
                out_translation.write("{0}\t'{1}'\n".format(s.upper(),translation))
            else:
                out_translation.write("{0}\t'Unknown'\n".format(s.upper()))
        out_translation.close()
        ####################################################################################


        # Output translation file

        out_translation = open(translation_file,'w')
        for protein in nodes:
            uE = session.get_user_entity(protein)
            translate=set()
            if translation_type_id == "proteinsequence":
                maxlen=0;
                for current_id in uE.get_attribute(attribute_identifier=translation_type_id):
                    if maxlen < len(current_id.value.get_sequence().upper()):
                        maxlen=len(current_id.value.get_sequence().upper())
                translation=",".join([str(current_id.value.get_sequence().upper()) for current_id in uE.get_attribute(attribute_identifier=translation_type_id) if len(str(current_id.value.get_sequence().upper())) == maxlen ] )
            else:
                for current_id in uE.get_attribute(attribute_identifier=translation_type_id):
                    translate.add(current_id.value.upper())
                translation="','".join(["{0}".format(x) for x in translate])
            out_translation.write("{0}\t'{1}'\n".format(protein,translation))
        out_translation.close()

    return

def fileExist (file):               #Checks if a file exists AND is a file
    if file is not None: return os.path.exists(file) and os.path.isfile(file)
    else: return False

def get_affinity_method_ids():
    affinity_dict={
    '0':    'molecular interaction',
    '4':    'affinity chromatography technology',
    '6':    'anti bait coimmunoprecipitation',
    '7':    'anti tag coimmunoprecipitation',
    '8':    'array technology',
    '9':    'bacterial display',
    '19':   'coimmunoprecipitation',
    '28':   'cosedimentation in solution',
    '29':   'cosedimentation through density gradient',
    '34':   'display technology',
    '47':   'far western blotting',
    '48':   'filamentous phage display',
    '49':   'filter binding',
    '50':   'flag tag coimmunoprecipitation',
    '60':   'ha tag coimmunoprecipitation',
    '62':   'his tag coimmunoprecipitation',
    '66':   'lambda phage display',
    '71':   'molecular sieving',
    '73':   'mrna display',
    '75':   'myc tag coimmunoprecipitation',
    '81':   'peptide array',
    '84':   'phage display',
    '89':   'protein array',
    '92':   'protein in situ array',
    '95':   'proteinchip(r) on a surface-enhanced laser desorption/ionization',
    '96':   'pull down',
    '98':   'ribosome display',
    '108':  't7 phage display',
    '109':  'tap tag coimmunoprecipitation',
    '115':  'yeast display',
    '225':  'chromatin immunoprecipitation array',
    '400':  'affinity technology',
    '402':  'chromatin immunoprecipitation assay',
    '405':  'competition binding',
    '411':  'enzyme linked immunosorbent assay',
    '412':  'electrophoretic mobility supershift assay',
    '413':  'electrophoretic mobility shift assay',
    '440':  'saturation binding',
    '492':  'in vitro',
    '493':  'in vivo',
    '657':  'systematic evolution of ligands by exponential enrichment',
    '676':  'tandem affinity purification',
    '678':  'antibody array',
    '695':  'sandwich immunoassay',
    '729':  'luminescence based mammalian interactome mapping',
    '858':  'immunodepleted coimmunoprecipitation',
    '892':  'solid phase assay',
    '899':  'p3 filamentous phage display',
    '900':  'p8 filamentous phage display',
    '921':  'surface plasmon resonance array',
    '946':  'miniaturized immunoprecipitation',
    '947':  'bead aggregation assay',
    '963':  'interactome parallel affinity capture',
    '1017': 'rna immunoprecipitation',
    '1028': 'modified chromatin immunoprecipitation',
    '1029': 'proteomics of isolated chromatin segments',
    '1031': 'protein folding/unfolding',
    '1087': 'monoclonal antibody blockade',
    }
    affinity=set(affinity_dict.keys())
    return affinity

def get_complementation_method_ids():
    complementation_dict={
    '0':    'molecular interaction',
    '10':   'beta galactosidase complementation',
    '11':   'beta lactamase complementation',
    '14':   'adenylate cyclase complementation',
    '18':   'two hybrid',
    '80':   'partial DNA sequence identification by hybridization',
    '90':   'protein complementation assay',
    '97':   'reverse ras recruitment system',
    '111':  'dihydrofolate reductase reconstruction',
    '112':  'ubiquitin reconstruction',
    '228':  'cytoplasmic complementation assay',
    '230':  'membrane bound complementation assay',
    '231':  'mammalian protein protein interaction trap',
    '232':  'transcriptional complementation assay',
    '369':  'lex-a dimerization assay',
    '370':  'tox-r dimerization assay',
    '397':  'two hybrid array',
    '398':  'two hybrid pooling approach',
    '399':  'two hybrid fragment pooling approach',
    '432':  'one hybrid',
    '437':  'protein three hybrid',
    '438':  'rna three hybrid',
    '492':  'in vitro',
    '493':  'in vivo',
    '588':  'three hybrid',
    '655':  'lambda repressor two hybrid',
    '726':  'reverse two hybrid',
    '727':  'lexa b52 complementation',
    '728':  'gal4 vp16 complementation',
    '809':  'bimolecular fluorescence complementation',
    '895':  'protein kinase A complementation',
    '916':  'lexa vp16 complementation',
    '1037': 'Split renilla luciferase complementation',
    '1111': 'two hybrid bait or prey pooling approach',
    '1112': 'two hybrid prey pooling approach',
    '1113': 'two hybrid bait and prey pooling approach',
    '1203': 'split luciferase complementation',
    '1204': 'split firefly luciferase complementation',
    '1320': 'membrane yeast two hybrid',
    }
    complementation=set(complementation_dict.keys())
    return complementation

def check_restriction(restriction):
    """
    Checks if the restriction has been correctly introduced.
    Returns the values for:
    restricted_to_TAP, restricted_to_Y2H, restricted_to_user, except_TAP, except_Y2H, except_user

    """
    if not restriction:
        return [False, False, None, False, False, None]
    else:
        res = restriction.lower()
        if res == 'aff':
            return [True, False, None, False, False, None]
        elif res == 'y2h':
            return [False, True, None, False, False, None]
        elif res == 'eaff':
            return [False, False, None, True, False, None]
        elif res == 'ey2h':
            return [False, False, None, False, True, None]
        else:
            raise IncorrectRestrictionType(res)



class IncorrectRestrictionType(Exception):
    """
    Exception that raises when a the restriction introduced to generate the
    network is incorrect.
    """
    def __init__(self, restriction):
        self.restriction = restriction
        self.restriction_types = ['aff', 'y2h', 'eaff', 'ey2h']
    def __str__(self):
        return 'The restriction of the network introduced ({}) is not admitted.\nThe types of restriction admitted in DIANA are: {}\n'.format(self.restriction, ', '.join(self.restriction_types))


