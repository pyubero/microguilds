# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 14:19:05 2022

@author: Pablo
"""
import numpy as np
from functional_clustering_utils import verboseprint
import functional_clustering_utils as fcutils

VERBOSE = True
FILENAME_ENRICHMENT = 'data_enrichment_potF.npz'
FILENAME_CLADE_DATA = 'data_clades_potF.npz'
FILENAME_TREE = 'tree_potF.newick'
# ...
FILENAME_OUT = 'significant_nodes.tsv'  # exporta por nodo MRCA su bcode
FILENAME_NPZ = 'significant_nodes.npz'  # in numpy format

# Accumulators
Z_THRESHOLD = 3
SIGNIFICANT = []


# Load tree adjacency matrix
ADJACENCY = fcutils.get_adjacency_matrix(FILENAME_TREE)

# Load clade data
data = fcutils.get_clade_data(FILENAME_TREE, treetype="newick")
clade_ids, clade_lfs, clade_dpt, leaf_names = data
# ...
CLADE_NLEAFS = np.array([len(leafs) for leafs in clade_lfs])

# Load enrichment results
data = np.load(FILENAME_ENRICHMENT, allow_pickle=True)
SPECIES_NAMES = data['species_names']
FEATURES_NAMES = data['features_names']
FEATURES = data['features']
UNIVOCITY = data['univocity']
ZSCORES = data['zscores']
verboseprint(f"Loaded enrichment data with {data['mcmax']} randomizations.",
             VERBOSE)

# ... filter out nans, substitute with 0s
X = ZSCORES.copy()
X[np.isnan(X)] = 0

# Clear output file
with open(FILENAME_OUT, 'w+', encoding="utf-8") as file:
    pass

# For every feature
for feature_idx in range(FEATURES.shape[1]):

    # Local accumulator
    mrca_nodes = []

    # Find nodes that are significant
    nodes_sign = np.argwhere(np.abs(X[:, feature_idx]) > Z_THRESHOLD)[:, 0]

    for node_idx in nodes_sign:

        # Search parent
        parent_idx = np.argwhere(ADJACENCY[:, node_idx])[:, 0]
        assert len(parent_idx) == 1

        # Add this node only when its parent is not significant...
        if np.abs(X[parent_idx, feature_idx]) < Z_THRESHOLD:
            mrca_nodes.append(node_idx)

    verboseprint(f"---- {FEATURES_NAMES[feature_idx]} ----", VERBOSE)
    verboseprint(f"Significant nodes: {len(nodes_sign)}", VERBOSE)
    verboseprint("Significant nodes without significant parent:", VERBOSE)
    _ = [verboseprint(f"  {_:4d} has {CLADE_NLEAFS[_]} leafs", VERBOSE)
         for _ in mrca_nodes]
    verboseprint("")

    # Show nodes that are contained in other significant nodes
    final_mrca_nodes = mrca_nodes.copy()
    for jj in mrca_nodes:
        jj_leafs = set(clade_lfs[jj])

        for kk in mrca_nodes:
            if jj == kk:
                continue
            kk_leafs = set(clade_lfs[kk])

            if kk_leafs.issubset(jj_leafs):
                if kk in final_mrca_nodes:
                    final_mrca_nodes.remove(kk)
                verboseprint(f'>> node {jj:4d} contains node {kk:4d}', VERBOSE)

    SIGNIFICANT.append(mrca_nodes)

    # with open(FILENAME_OUT, 'a+') as file:
    #    file.write(f'[{feature_idx}] {FEATURES_NAMES[feature_idx]}: \n')

    #    for idx in mrca_nodes:
    #        message = f"Node {idx:4d} contains {CLADE_NLEAFS[idx]:4d} leafs
    # #with Z={ZSCORES[idx, feature_idx]:+2.3f}"
    #        verboseprint(message, VERBOSE)
    #        file.write(message+"\n")

    #       jj_leafs = set(clade_lfs[idx])

    #        for kk in mrca_nodes:
    #            if idx == kk:
    #                continue
    #            kk_leafs = set(clade_lfs[kk])

    #            if kk_leafs.issubset(jj_leafs):
    #                if kk in final_mrca_nodes:
    #                    final_mrca_nodes.remove(kk)
    #                file.write(f'  contains node {kk:4d}\n')

    #    file.write('\n')

# Export numpy file
np.savez(FILENAME_NPZ, significant_nodes=SIGNIFICANT)
verboseprint(f"{FILENAME_NPZ} saved.", VERBOSE)

# Export table
idc_nodes = np.unique([_ for a in SIGNIFICANT for _ in a])
with open(FILENAME_OUT, 'w+', encoding="utf-8") as file:
    file.write('node\tnleafs\t' + '\t'.join(list(FEATURES_NAMES)) + '\n')
    for idx in idc_nodes:
        message = f"{idx:5d}\t" + "\t  ".join([f"{_:.2f}" for _ in X[idx, :]])
        file.write(message + "\n")
