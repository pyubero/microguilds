# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:17:35 2022

@author: Pablo Yubero
"""
import os
from Bio import Phylo
from tqdm import tqdm
import numpy as np


def verboseprint(msg, verbose=True, condition=True):
    '''Prinst a message if VERBOSE is True.'''
    if verbose and condition:
        print(msg)


def verbosebar(iterable, verbose=True):
    '''Generates a progress bar with tqdm if VERBOSE is True.'''
    if verbose:
        return tqdm(iterable)
    else:
        return iterable


def create_temporal_folder():
    '''Creates a temp folder in the current working directory if it doesn't
    exist already.'''
    cwd = os.getcwd()
    temporal_folder = os.path.join(cwd, "temp")
    if not os.path.isdir(temporal_folder):
        os.mkdir(temporal_folder)
    return temporal_folder


def get_clade_data(filename, treetype="newick", filename_out=None):
    '''Reads a phylogenetic tree and extracts relevant data from clades,
    which is how we call internal nodes, as opposed to leafs (external)'''

    # If FILENAME_OUT already exists, then simply load it.
    if isinstance(filename_out, str) and os.path.exists(filename_out):
        return load_clade_data(filename_out)

    # Import tree
    tree = Phylo.read(filename, treetype)
    depths = tree.depths()

    # Obtain (terminal) leaf data
    leaf_names = np.array([leaf.name for leaf in tree.get_terminals()])

    # Create accumulators...
    clade_ids = []  # clade identifiers (names)
    clade_lfs = []  # clade leaf identifiers (indices)
    clade_dpt = []  # clade's MRCA distance to the root

    # For every clade, extract relevant data
    for clade in tqdm(tree.get_nonterminals()):
        _leaf_names = [leaf.name for leaf in clade.get_terminals()]
        _leaf_idc = [np.argwhere(name == leaf_names)[:, 0][0]
                     for name in _leaf_names]

        clade_ids.append(clade.name)
        clade_lfs.append(_leaf_idc)
        clade_dpt.append(depths[clade])

    if filename_out is not None:
        np.savez(filename_out,
                 leaf_names=np.array(leaf_names, dtype=object),
                 clade_ids=np.array(clade_ids),
                 clade_lfs=np.array(clade_lfs, dtype=object),
                 clade_dpt=np.array(clade_dpt))
        print(f"\nData exported to {filename_out}.")

    print(f"\nData loaded from {filename}.")
    print(f"Found data for > {len(leaf_names)} leafs,")
    print(f"               > {len(clade_ids)} internal nodes.")

    return clade_ids, clade_lfs, clade_dpt, leaf_names


def load_clade_data(filename):
    '''Load clade data from npz file.'''
    data = np.load(filename, allow_pickle=True)
    leaf_names = data['leaf_names']
    clade_ids = data['clade_ids']
    clade_lfs = data['clade_lfs']
    clade_dpt = data['clade_dpt']

    print(f"\nData loaded from {filename}.")
    print(f"Found data for > {len(leaf_names)} leafs,")
    print(f"               > {len(clade_ids)} internal nodes.")
    return clade_ids, clade_lfs, clade_dpt, leaf_names


def get_species(name):
    '''Extract species name from a string formatted like:
    d_Bacteria_
    p_Proteobacteria_
    c_Gammaproteobacteria_
    o_Pseudomonadales_
    f_Pseudomonadaceae_
    g_Pseudomonas_F_
    s_Pseudomonas_F_alcaligenes <--- would get this tag
    '''
    if "_s_" not in name:
        return None
    sp_name = name.split('_s_')[-1]
    if sp_name == '':
        # print('Species name not found.')
        return None
    else:
        return sp_name.split('_')[-1]


def get_genus(name):
    '''Extract species name from a string formatted like:
    d_Bacteria_
    p_Proteobacteria_
    c_Gammaproteobacteria_
    o_Pseudomonadales_
    f_Pseudomonadaceae_
    g_Pseudomonas_F_ <---- would get this tag
    s_Pseudomonas_F_alcaligenes
    '''
    sp_name = name.split('_g_')[-1].split('_')[0]
    if sp_name == '':
        return None
    else:
        return sp_name


def get_adjacency_matrix(tree_filename):
    # Load tree
    tree = Phylo.read(tree_filename, 'newick')
    nodes = tree.get_nonterminals()

    # Accumulator
    connectivity = np.zeros((len(nodes), len(nodes)))

    # Compute connectivity matrix
    for node in nodes:
        parent_idx = int(node.name.split('_')[1])
        children = node.clades
        for child in children:
            if not child.is_terminal():
                child_idx = int(child.name.split('_')[1])
                connectivity[parent_idx, child_idx] = 1

    return connectivity
