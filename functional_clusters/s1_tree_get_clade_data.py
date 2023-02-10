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
# FILENAME_TREE = 'tree_rplB.newick' 
FILENAME_TREE = 'tree_potF_labelled.newick' 
# FILENAME_TREE = 'tree_16s_ref.newick' 
FILENAME_OUT  = 'data_clades_potF.npz'
#...
t = Phylo.read(FILENAME_TREE, 'newick')
depths = t.depths()

# Obtain (terminal) leaf data
leafs = t.get_terminals()
leaf_names= np.array( [ leaf.name for leaf in leafs ] )
nleafs= len(leafs)

# Obtain internal node data
all_internal = t.get_nonterminals()


# I would like to output a list with:
clade_ids = []
clade_lfs = []
clade_dpt = [] #... internal node's distance to the root


for clade in all_internal:
    _id = clade.name
    _depth = depths[clade]
    _leafs = clade.get_terminals()
    _leaf_names = [ leaf.name for leaf in _leafs ]
    _leaf_idc = [np.argwhere(name==leaf_names)[:,0][0] for name in _leaf_names]
    
    _leaf_depths = [ depths[leaf] for leaf in _leafs ]
    # _depth2 = np.mean(_leaf_depths)   - depths[clade]
    # _depth3 = np.median(_leaf_depths) - depths[clade]
    # _depth4 = np.max(_leaf_depths)    - depths[clade]
    
    
    clade_ids.append( _id )
    clade_lfs.append(_leaf_idc )
    clade_dpt.append(  _depth )
    # clade_dpt2.append( _depth2 )
    # clade_dpt3.append( _depth3 )
    # clade_dpt4.append( _depth4 )
    
    
clade_ids = np.array(clade_ids)
clade_lfs = np.array(clade_lfs)
clade_dpt = np.array(clade_dpt)
# clade_dpt2= np.array(clade_dpt2)
# clade_dpt3= np.array(clade_dpt3)
# clade_dpt4= np.array(clade_dpt4)

np.savez(FILENAME_OUT, leaf_names=leaf_names, clade_ids=clade_ids,
          clade_lfs = clade_lfs,
          clade_dpt = clade_dpt)
          # clade_dpt2 = clade_dpt2,
          # clade_dpt3 = clade_dpt3,
          # clade_dpt4 = clade_dpt4,
          # description = description )


# Load all clade data
data = np.load(FILENAME_OUT, allow_pickle=True)
leaf_names = data['leaf_names']
clade_ids  = data['clade_ids']
clade_dpt  = data['clade_dpt'] 
clade_lfs  = data['clade_lfs']

print('Number of terminal leafs: ', len(leaf_names) )
print('Number of internal nodes: ', len(clade_dpt) )




