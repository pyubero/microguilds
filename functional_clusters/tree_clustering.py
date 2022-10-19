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


FILENAME_ZSCORES = 'data_clade_enrichment.npz'
FILENAME_CLADE_DATA = 'data_potF_all_clades.npz'
FILENAME_TREE = 'tree_potF_labelled.newick' 


# Load tree
tree = Phylo.read( FILENAME_TREE, 'newick')
clades = tree.get_nonterminals()

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



data = np.load(FILENAME_ZSCORES, allow_pickle=True)
F = data['F']
S = data['S']
ZSCORES= data['ZSCORES']
# names = data['names']
# features= data['features']

X = ZSCORES.copy()
X[ np.isnan(X)] = 0

threshold = 2

for feature_idx in range(15):
    
    nodes_sign = np.argwhere( np.abs(X[:,feature_idx])>threshold)[:,0]
    
    
    MRCA_NODES = []
    for node_idx in nodes_sign:
        
        # Buscar padre
        parent_idx = np.argwhere(C[:, node_idx])[:,0][0]
    
        # Si el padre no es significativo:
        if np.abs(X[parent_idx,feature_idx])<threshold:
            MRCA_NODES.append(node_idx)
            
    print( len(nodes_sign), len(MRCA_NODES))        
    print(MRCA_NODES)
    print( [clade_nleafs[_] for _ in MRCA_NODES])
    
    
    
    
    
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
    
    with open('para_juan.txt','a+') as file:
        file.write('Feature index: %d\n' % feature_idx )
        
        for idx in FINAL_MRCA_NODES:
            print( "Node %3d contains %4d leafs with Z = %2.2f" % (idx, clade_nleafs[idx], ZSCORES[idx, feature_idx]) )
            file.write("Node %3d contains %4d leafs with Z = %2.2f\n" % (idx, clade_nleafs[idx], ZSCORES[idx, feature_idx]) )
        file.write('\n')





    

def poisson_nevents(lam=1.0, ntrials=1):
    def poisson_pdf(lam=1.0, k=1.0):
        return lam**k*np.exp(-lam)/np.math.factorial(k)
    
    
    ref = np.cumsum( [poisson_pdf(lam,_) for _ in range(5*lam)] )
    nevents = []
    for _ in range(ntrials):
        rr = np.random.rand()       
        nevents.append( np.argwhere(rr<ref)[:,0][0]    )
    return nevents
    
    
    
    
    
    
