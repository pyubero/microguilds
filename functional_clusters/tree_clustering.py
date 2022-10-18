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
X = np.clip(X, -3,3)


order = np.argsort( clade_dpt )


t = 1

clades[1]





# # Find all nodes closer than threshold to the root
# idx = np.argwhere(clade_dpt<threshold)[:,0]
# print(idx)

# # Find clades that compose a partition, that is
# partition_idc =  [ idx[-1],]
# partition_leafs= set( [_ for _ in clade_lfs[idx[-1]]] )

# for jj in np.arange( len(idx)-2, -1,-1):
    
#     if not set( clade_lfs[ idx[jj]] ).issubset( partition_leafs ):
#         partition_idc[-1] = jj
#         partition_leafs |= set ( clade_lfs[jj] )
# print(partition_idc)


# clade_lfs[8] = [8,9]








# # first detect features with all nan values
# idx1_allnan = [_ for _ in range( X.shape[1]) if not (_ in np.argwhere(np.all(np.isnan(X), axis=0))[:,0])]
# X = X[:,idx1_allnan]

# # then detect samples with all nan values
# idx0_allnan = [_ for _ in range( X.shape[0]) if not (_ in np.argwhere(np.any(np.isnan(X), axis=1))[:,0])]
# X = X[idx0_allnan,:]

# CH = []
# # for nclusters in range(30):
# for nclusters in (1,):
        
#     model = AgglomerativeClustering(n_clusters=nclusters+2, affinity='cosine', linkage='complete')
#     model.fit(X)
#     clabels = model.labels_
    
    
#     LABELS = []
#     C_BCODE= []
    
#     for icluster in np.unique(clabels):
#         idx = np.argwhere( clabels == icluster)[:,0]
        
#         LABELS.append('C%d, %d' % (icluster, len(idx)) )
#         C_BCODE.append( np.mean(X[idx,:], axis=0) )
    
#     C_BCODE = np.array(C_BCODE)
    
#     plt.imshow(C_BCODE, aspect='auto', vmin=-2, vmax=2)
#     plt.yticks( ticks=range(clabels.max()+1), labels= LABELS)
#     plt.colorbar()

# # dendrogram = sch.dendrogram(sch.linkage(X, method='ward'))



    
#     ch_index = calinski_harabasz_score(X, clabels)
#     CH.append(ch_index)
#     # print(ch_index)



# # plt.plot( 2+np.array(range(len(CH))), CH,'.-')
# # plt.xlabel('Number of clusters')
# # plt.ylabel('Calinski-Harabsz score')
# # plt.grid()
