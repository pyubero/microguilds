# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:17:35 2022

@author: Pablo

This file loads a phylogenetic tree and exports:
    - leaf names
    - (internal) node identifiers
    - the leafs associated to each (internal) node
    - (internal) node "depth", "height" or simply, "distance"

---> convert this file into  i) analyzer and ii) loader functions (to/from some temporal folder?)
"""
from Bio import Phylo
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
import os


def create_temporal_folder():
    '''Creates a temp folder in the current working directory if it doesn't already axist'''
    cwd = os.get_cwd()
    temporal_folder = os.path.join(cwd, "temp")
    if not os.path.is_dir(temporal_folder):
        os.mkdir(temporal_folder)
    return temporal_folder


def get_clade_data(filename_in, treetype="newick", filename_out=None):
    # Import tree
    t = Phylo.read(filename_in, treetype)
    depths = t.depths()

    # Obtain (terminal) leaf data
    leaf_names= np.array([leaf.name for leaf in t.get_terminals()])
    nleafs= len(leaf_names)

    # Create accumulators...
    clade_ids = []  # clade identifier
    clade_lfs = []  # clade leaf identifiers
    clade_dpt = []  # clade's MRCA distance to the root

    # For every clade, extract relevant data
    for clade in tqdm(t.get_nonterminals()):
        _leaf_names = [leaf.name for leaf in clade.get_terminals()]
        _leaf_idc = [np.argwhere(name == leaf_names)[:, 0][0] for name in _leaf_names]
        # _leaf_depths = [ depths[leaf] for leaf in clade.get_terminals()] #<---- what is this doing here??

        clade_ids.append(clade.name)
        clade_lfs.append(_leaf_idc )
        clade_dpt.append(depths[clade])

    if filename_out is not None:
        np.savez(filename_out,
                leaf_names=np.array(leaf_names),
                clade_ids=np.array(clade_ids),
                clade_lfs=np.array(clade_lfs),
                clade_dpt=np.array(clade_dpt))
        print(f"Data exported to {filename_out}.")

    return clade_ids, clade_lfs, clade_dpt

filename = "data_clades_potF.npz"
def load_clade_data(filename):
    # Load all clade data
    data = np.load(filename, allow_pickle=True)
    leaf_names = data['leaf_names']
    clade_ids  = data['clade_ids']
    clade_lfs  = data['clade_lfs']
    clade_dpt  = data['clade_dpt']

    print(f"Data loaded from {filename}.")
    print(f"Found data for > {len(leaf_names)} leafs,")
    print(f"               > {len(clade_ids)} internal nodes.")
    return leaf_names, clade_ids, clade_lfs, clade_dpt



