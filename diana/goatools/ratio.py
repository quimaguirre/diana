#!/usr/bin/env python
# -*- coding: UTF-8 -*-

__copyright__ = "Copyright (C) 2010-2017, H Tang et al., All rights reserved."
__author__ = "various"

from collections import defaultdict, Counter


def count_terms(geneset, assoc, obo_dag):
    """count the number of terms in the study group
    """
    term_cnt = Counter()
    for gene in (g for g in geneset if g in assoc):
        for x in assoc[gene]:
            if x in obo_dag:
                term_cnt[obo_dag[x].id] += 1

    return term_cnt

def get_terms(desc, geneset, assoc, obo_dag, log):
    """Get the terms in the study group
    """
    _chk_gene2go(assoc)
    term2itemids = defaultdict(set)
    genes = [g for g in geneset if g in assoc]
    for gene in genes:
        for x in assoc[gene]:
            if x in obo_dag:
                term2itemids[obo_dag[x].id].add(gene)
    log.write("{N:>6,} out of {M:>6,} {DESC} items found in association\n".format(
        DESC=desc, N=len(genes), M=len(geneset)))
    return term2itemids

def is_ratio_different(min_ratio, study_go, study_n, pop_go, pop_n):
    """
    check if the ratio go /n is different between the study group and
    the population
    """
    if min_ratio is None:
        return True
    s = float(study_go) / study_n
    p = float(pop_go) / pop_n
    if s > p:
        return s / p > min_ratio
    return p / s > min_ratio

def _chk_gene2go(assoc):
    """Check that associations is gene2go, not go2gene."""
    if not assoc:
        raise RuntimeError("NO ITEMS FOUND IN ASSOCIATIONS")
    for key, val in assoc.items():
        if isinstance(key, str) and key[:3] == "GO:":
            raise Exception("ASSOCIATIONS EXPECTED TO BE gene2go, NOT go2gene: {EX}".format(
                EX=assoc.items()[:2]))
        return

# Copyright (C) 2010-2017, H Tang et al., All rights reserved.
