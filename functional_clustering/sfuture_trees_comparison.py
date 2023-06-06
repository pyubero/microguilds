# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:07:26 2022

@author: Pablo
"""
import warnings

from Bio import Phylo
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import functional_clustering_utils as fcutils


def standardize_name(name):
    '''Adapts and lower-cases name iff it is a string.'''
    if isinstance(name, str):
        name = name.replace("[", "")
        name = name.replace("]", "")
        return name.lower()
    else:
        return "IGNORE"


def leaf2id(treeleaf):
    '''Converts a leaf name to an id.'''
    species = standardize_name(fcutils.get_species(treeleaf.name))
    genus = standardize_name(fcutils.get_genus(treeleaf.name))
    return f"{genus}{species}"


def tree_mapping(tree1, tree2):
    '''Maps the leafs from one tree to the leafs in the other.
    It returns a matrix mapping with shape (numleafs1, numleafs2) such that,
    - leaf 0 in tree 1 can be found in leafs np.nonzero(mapping[0,:]) of tree 2
    - leaf 7 in tree 2 can be found in leafs np.nonzero(mapping[:,7]) of tree 1

    Important: leaf identifiers with "IGNORE" are effectively unmapped.
     '''

    # Compute map from one tree to the other
    leafs1 = tree1.get_terminals()
    leafs2 = tree2.get_terminals()

    mapping = np.zeros((len(leafs1), len(leafs2)))

    # Build lists of identifiers, we currently use genus+species
    tree1_ids = np.array([leaf2id(leaf) for leaf in leafs1])
    tree2_ids = np.array([leaf2id(leaf) for leaf in leafs2])

    # Find mapping
    for row, leaf_id in enumerate(tree1_ids):
        if "IGNORE" in leaf_id:
            continue
        cols = np.argwhere(tree2_ids == leaf_id)[:, 0]
        mapping[row, cols] = 1
    return mapping, (tree1_ids, tree2_ids)


GENEX = "16S"
GENEY = "gyrB"
EXPORT_DATA = True
# ...
FILENAME_TREE_X = f'./data/tree_{GENEX}.newick'
FILENAME_TREE_Y = f'./data/contree_{GENEY}_labelled.newick'
FILENAME_OUT = f'./data/data_tree_comparison_{GENEX}_{GENEY}.npz'


# Load trees
treeX = Phylo.read(FILENAME_TREE_X, 'newick')
leafsX = treeX.get_terminals()
radiusX = np.max([treeX.distance(treeX.root, leaf) for leaf in leafsX])
# ...
treeY = Phylo.read(FILENAME_TREE_Y, 'newick')
leafsY = treeY.get_terminals()
radiusY = np.max([treeY.distance(treeY.root, leaf) for leaf in leafsY])

# Prepare output data accumulators
distance_x, distance_y = [], []
nof_leafs = []

# Precompute tree mappings
MAP, (ids_x, ids_y) = tree_mapping(treeX, treeY)


for clade in tqdm(treeY.get_nonterminals()):

    leafs_in_x = []
    leafs_in_y = []

    for leaf in clade.get_terminals():
        # Build leaf id:
        leafid = leaf2id(leaf)
        map_col = np.nonzero(leafid == ids_y)[0][0]

        # Find y_species and y_genus in the other tree,
        # ... and accumulate
        map_rows = np.nonzero(MAP[:, map_col])[0]

        # Accumulate
        leafs_in_y.append(map_col)
        _ = [leafs_in_x.append(row) for row in map_rows]

    # Raise warning
    if len(leafs_in_x) != len(leafs_in_y):
        warnings.warn(
            '''\nThere might be some leafs in tree Y present more
            than once in tree X; or some leafs in tree Y that
            tree X is lacking.'''
        )

    # Compute distances between clades
    clade_terminals_x = [leafsX[ii] for ii in np.unique(leafs_in_x)]
    clade_terminals_y = [leafsY[ii] for ii in np.unique(leafs_in_y)]

    mrca_x = treeX.common_ancestor(clade_terminals_x)
    mrca_y = treeY.common_ancestor(clade_terminals_y)

    rel_x = np.mean(
        [treeX.distance(mrca_x, ctx) for ctx in clade_terminals_x])/radiusX
    rel_y = np.mean(
        [treeY.distance(mrca_y, cty) for cty in clade_terminals_y])/radiusY

    distance_x.append(rel_x)
    distance_y.append(rel_y)
    nof_leafs.append(len(leafs_in_y))


# Remember that
if EXPORT_DATA:
    print('')
    print(f'Number of points: {len(distance_x)}')
    print(f'Data saved to {FILENAME_OUT}.')
    np.savez(FILENAME_OUT,
             nof_leafs=nof_leafs,
             distance_x=distance_x,
             distance_y=distance_y
             )

plt.plot(distance_x, distance_y, '.')
plt.xlabel(f"Relatedness of {FILENAME_TREE_X}")
plt.ylabel(f"Relatedness of {FILENAME_TREE_Y}")
plt.show()
