# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 14:19:05 2022

@author: Pablo
"""

import numpy as np
from matplotlib import pyplot as plt
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch
from sklearn.metrics import calinski_harabasz_score
from Bio import Phylo


FILENAME_ZSCORES = 'data_clade_enrichment_potF.npz'
FILENAME_CLADE_DATA = 'data_clades_potF.npz'
FILENAME_TREE = 'tree_potF_labelled.newick' 
FILENAME_OUT ='para_juan_nuevo.txt'                 # exporta por features el MRCA
FILENAME_OUT3= 'para_juan_per_node_nuevo.txt'       # exporta por nodo MRCA su bcode
FILENAME_OUT2= 'data_nodes_sign_filtered.npz' # exporta los nodos MRCA


# Load tree
tree = Phylo.read( FILENAME_TREE, 'newick')
clades = tree.get_nonterminals()

# Compute the connectivity matrix
C = np.zeros(( len(clades), len(clades) ) )
for node in clades:
    parent_idx = int( node.name.split('_')[1] )
    for child in node.clades:
        if not child.is_terminal():
            child_idx = int( child.name.split('_')[1] )
            C[parent_idx, child_idx] = 1
            



# Load all clade data
data = np.load(FILENAME_CLADE_DATA, allow_pickle=True)
leaf_names = data['leaf_names']
clade_ids  = data['clade_ids']
clade_dpt  = data['clade_dpt'] 
clade_lfs  = data['clade_lfs']
clade_nleafs= np.array( [len(leafs) for leafs in clade_lfs])


# Load enrichment results
data = np.load(FILENAME_ZSCORES, allow_pickle=True)
F = data['F']
S = data['S']
ZSCORES= data['ZSCORES']
names = data['names']
features= data['features']
MCMAX = data['MCMAX']

# ... substitute nans with 0s
X = ZSCORES.copy()
X[ np.isnan(X)] = 0

# ... the Z-threshold 
threshold = 3
all_nodes = []
all_barcodes=[]

with open(FILENAME_OUT,'w+') as file: pass


# For every feature
for feature_idx in range(15):
    
    #... find nodes that are significant
    nodes_sign = np.argwhere( np.abs(X[:,feature_idx])>threshold)[:,0]
    
    
    MRCA_NODES = []
    for node_idx in nodes_sign:
        
        # Search parent
        parent_idx = np.argwhere(C[:, node_idx])[:,0][0]
    
        # Add this node only when its parent is not significant...
        if np.abs(X[parent_idx,feature_idx])<threshold:
            MRCA_NODES.append(node_idx)
            
            
    print( len(nodes_sign), len(MRCA_NODES))        
    print(MRCA_NODES)
    print( [clade_nleafs[_] for _ in MRCA_NODES])
    
    
    
    
    # Across all MRCA nodes, discard those nodes that are contained in other
    # nodes that are also significant
    FINAL_MRCA_NODES = MRCA_NODES.copy()
    for jj in MRCA_NODES:
        jj_leafs = set(clade_lfs[jj])
    
        for kk in MRCA_NODES:
            if jj==kk:
                continue
            kk_leafs = set(clade_lfs[kk])
            
            if kk_leafs.issubset( jj_leafs):
                if kk in FINAL_MRCA_NODES:
                    FINAL_MRCA_NODES.remove(kk)
                print('OJO el nodo %d contiene al nodo %d' % (jj, kk) )
    
    print(MRCA_NODES)
    print(FINAL_MRCA_NODES)
    all_nodes.append(FINAL_MRCA_NODES)
    
    with open(FILENAME_OUT,'a+') as file:
        file.write('Feature index: %d\n' % feature_idx )
        
        for idx in FINAL_MRCA_NODES:
            print( "Node %3d contains %4d leafs with Z = %2.2f" % (idx, clade_nleafs[idx], ZSCORES[idx, feature_idx]) )
            file.write("Node %3d contains %4d leafs with Z = %2.2f\n" % (idx, clade_nleafs[idx], ZSCORES[idx, feature_idx]) )
        file.write('\n')



np.savez( FILENAME_OUT2, sign_nodes=all_nodes)

    
idc_nodes = np.unique([_ for a in all_nodes for _ in a ])
with open(FILENAME_OUT3, 'w+') as file:
    file.write('Node;nleafs;'+';'.join(features) + '\n')
    for idx in idc_nodes:
        # string = ' '.join([ features[_] for _ in range(len(features)) if np.abs(X[idx,_])> threshold])
        # file.write('Node '+str(idx)+': '+string+'\n')
        file.write("%d;"%idx +';'.join(["%1.2f" % _ for _ in X[idx,:]])+'\n')  
        
        
        
        


