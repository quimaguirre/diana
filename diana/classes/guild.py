

class GUILDifier(object):
    """
    Class defining a GUILDifier object
    """

    def __init__(self, drug_name):
        """
        """


    ###########
    # METHODS #
    ###########

    def prepare(self, user_entity_ids):
        """
        Create files needed to run GUILD.
        """

        # Create initial node scores file
        with open(self.node_scores_file, 'w') as node_scores_fd:
            for user_entity_id in user_entity_ids_all:
                if user_entity_id in user_entity_ids:
                    score = GUILDifier.DEFAULT_SEED_SCORE
                    seed_to_score[user_entity_id] = GUILDifier.DEFAULT_SEED_SCORE
                else:
                    score = GUILDifier.DEFAULT_NON_SEED_SCORE
                node_scores_fd.write("%s %f\n" % (user_entity_id, score))

        # Create initial edge scores as node scores file for NetShort
        self.create_edge_scores_as_node_scores_file(edges = edges, node_to_score = seed_to_score, edge_to_score = edge_to_score, edge_scores_file = self.edge_scores_as_node_scores_file, ignored_nodes = None, default_score = GUILDifier.DEFAULT_NON_SEED_SCORE)

        # Create sampled graphs for NetZcore
        if not os.path.exists(self.sampled_graph_prefix + "1"): 
            print "Creating sampled networks"
            self.sample_network_preserving_topology(self.edge_scores_file, GUILDifier.N_SAMPLE_GRAPH, self.sampled_graph_prefix)

        # Write seed user entity information
        fetcher = BIANAInfoFetcher(self.settings)
        user_entity_to_values = fetcher.get_user_entities_info(self.node_info_file, user_entity_ids)
        f = open(self.seed_file, 'w')
        f.write("BIANA ID\tUniProt ID\tGene Symbol\tDescription\tGene ID\tEquivalent Entries\n") #% "\t".join(BIANAQuery.NODE_ATTRIBUTES[1:]))
        for user_entity_id in user_entity_ids:
            entry_ids, gene_symbols, descriptions, gene_ids, equivalent_entries = user_entity_to_values[user_entity_id]
            if entry_ids == "-":
                continue
            entry_ids = entry_ids.strip("; ")
            inner_values = [user_entity_id]
            inner_values.extend([entry_ids, gene_symbols, descriptions, gene_ids, equivalent_entries])
            f.write("%s\n" % "\t".join(inner_values))
        f.close()
        return
