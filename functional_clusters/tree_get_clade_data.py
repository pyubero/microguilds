# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:17:35 2022

@author: Pablo
"""
from Bio import Phylo
from tqdm import tqdm 
import numpy as np
from matplotlib import pyplot as plt


# Load tree
ifile = 'tree_potF_labelled.newick' 
#...
t = Phylo.read(ifile, 'newick')
depths = t.depths()

# Obtain (terminal) leaf data
leafs = t.get_terminals()
leaf_names= np.array( [ leaf.name for leaf in leafs ] )
nleafs= len(leafs)

# Obtain internal node data
all_internal = t.get_nonterminals()


# I would like to output a list with:
clade_ids = []
clade_dpt = []
clade_lfs = []


for clade in all_internal:
    _id = clade.name
    _depth = depths[clade]
    _leaf_names = [ leaf.name for leaf in clade.get_terminals() ]
    _leaf_idc = [np.argwhere(name==leaf_names)[:,0][0] for name in _leaf_names]
    
    clade_ids.append( _id )
    clade_dpt.append( _depth )
    clade_lfs.append(_leaf_idc )
    
    
clade_ids = np.array(clade_ids)
clade_dpt = np.array(clade_dpt)
clade_lfs = np.array(clade_lfs)

np.savez('all_clade_data.npz', leaf_names=leaf_names, clade_ids=clade_ids,
          clade_dpt = clade_dpt,
          clade_lfs = clade_lfs)


# Load all clade data
data = np.load('data_all_clades.npz', allow_pickle=True)
leaf_names = data['leaf_names']
clade_ids  = data['clade_ids']
clade_dpt  = data['clade_dpt'] 
clade_lfs  = data['clade_lfs']






