# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 14:19:05 2022

@author: Pablo
"""
import numpy as np
from functional_clustering_utils import verboseprint
import functional_clustering_utils as fcutils

GENE = "potF"
FILENAME_ENRICHMENT = f'data_enrichment_{GENE}.npz'
FILENAME_CLADE_DATA = f'data_clades_{GENE}.npz'
FILENAME_TREE = f'tree_{GENE}.newick'
# ...
FILENAME_OUT = f'significant_nodes_{GENE}.tsv'  # exporta por nodo su bcode
FILENAME_NPZ = f'significant_nodes_{GENE}.npz'  # in numpy format
VERBOSE = False
# Accumulators
Z_THRESHOLD = 2.9
SIGNIFICANT = []


# Load tree adjacency matrix
ADJACENCY = fcutils.get_adjacency_matrix(FILENAME_TREE)

# Load clade data
data = fcutils.get_clade_data(FILENAME_TREE,
                              treetype="newick",
                              filename_out=FILENAME_CLADE_DATA)
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

    verboseprint(f"\n---- {FEATURES_NAMES[feature_idx]} ----", VERBOSE)
    verboseprint(f"Significant nodes: {len(nodes_sign)}", VERBOSE)
    verboseprint("Significant nodes without significant parent:", VERBOSE)
    _ = [verboseprint(f"  {_:4d} has {CLADE_NLEAFS[_]} leafs", VERBOSE)
         for _ in mrca_nodes]

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

    # Export final_mrca_nodes if you wish to remove significant nodes whose
    # direct parent is not significant, but perhaps some other "secondary" or
    # more distant parent is again significant. From 70 to 50 entries in potF.
    SIGNIFICANT.append(mrca_nodes)
    # SIGNIFICANT.append(final_mrca_nodes)


# Export table
idc_nodes = np.unique([_ for a in SIGNIFICANT for _ in a])
with open(FILENAME_OUT, 'w+', encoding="utf-8") as file:
    file.write('node\t' + '\t'.join(list(FEATURES_NAMES)) + '\n')
    for idx in idc_nodes:
        message = f"{idx:5d}\t" + "\t  ".join([f"{_:.2f}" for _ in X[idx, :]])
        file.write(message + "\n")

# Export numpy file
np.savez(FILENAME_NPZ, significant_nodes=np.array(SIGNIFICANT, dtype=object))
verboseprint(f"\n{FILENAME_NPZ} saved.", VERBOSE)

print(len(idc_nodes))
