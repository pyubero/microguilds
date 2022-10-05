# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:36:26 2022

@author: Pablo
"""

import numpy as np
import pandas as pd
import pickle as pkl
from tqdm import tqdm
from matplotlib import pyplot as plt




def entropy(f):
    f = np.array(f)
    f = f[ np.nonzero(f) ]
    p = f/np.sum(f) # convert frequencies to probabilities
    return -np.sum( p*np.log(p) )



# name1 = M_names[0]
# name2 = M_names[1]

# sp_name = name1.split('_g_')[-1]


def get_species(name):
    sp_name = name.split('_s_')[-1]
    if sp_name == '':
        # print('Species name not found.')
        return None
    else:
        return sp_name.split('_')[-1]
    
def get_genus(name):
    sp_name = name.split('_g_')[-1].split('_')[0]
    if sp_name=='':
        return None
    else: 
        return sp_name

def bin_count(array):
    unique_values = np.unique(array)
    return np.array( [ np.sum(array==value) for value in unique_values] )




# Load tree data (functional tree, not phylogenetic)
data = np.load('tree_distance.npz')
M = data['M']
M_names=data['names']

# Load feature matrix
df = pd.read_csv("potfclusters_final.csv", sep=",", encoding='latin1')
ft = df.to_numpy()
ft_names = np.array( [ line[0] for line in ft ] )
col_names= np.array( [ col for col in df.columns] )[1:16]    
ft = np.array([ line[1:16] for line in ft ])




# ... create true feature matrix
F = np.zeros((M.shape[0], ft.shape[1]))*np.nan

for idx_ft, species in enumerate(ft_names):
    _name = species.split(' ')
    for idx_M, full_name in enumerate(M_names):
        if (_name[0] in full_name) and (_name[1] in full_name):
            F[idx_M,:] = ft[idx_ft,:]

print('Data obtained for %d tree entries.' % (1158-np.sum(np.all( np.isnan(F), axis=1))))            




# First obtain clusters
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import cophenet, linkage
from sklearn.metrics  import calinski_harabasz_score
import warnings
warnings.filterwarnings('ignore')


S = [] 
MCR=[]
PVAL_L=[]
MCMAX = 9999

genus_list  = np.array( [ get_genus(name) for name in M_names] )
species_list= np.array( [ get_species(name) for name in M_names] )


# NCLUSTERS = 4
cluster_array = np.arange(1,202,2)

#...
for i_nclusters in tqdm(range(len(cluster_array))):
    NCLUSTERS = cluster_array[i_nclusters]

    
    
    # Automatic clustering
    # ... assign cluster labels to each tree entry
    clustering_api = KMeans(n_clusters= NCLUSTERS, # Number of clusters to form
                            n_init = 20,  # Number of times the k-means algorithm will be run with different seeds.
                            random_state=1337)
    
    kmeans = clustering_api.fit(M)
    labels = kmeans.labels_
        
    
    
    # Compute the taxonomic diversity in each cluster and append to S
    # ... tax. div. being equal to the inverse of Shannon's entropy on genus
    # ... note however that we could also use species' diversity.
    
    _s=[]
    for cluster_idx in range(labels.max()+1):
        idx = np.argwhere( labels==cluster_idx )[:,0]
        cluster_taxonomy = np.array( [ genus for genus in genus_list[idx] if genus is not None] )
        h = bin_count( cluster_taxonomy )
        _s.append( 1/entropy( h ) )
    S.append(_s)
    
    
    
    _obs    = np.zeros((NCLUSTERS, ft.shape[1] ) )
    _pval_l = np.zeros((NCLUSTERS, ft.shape[1] ) )
    _mc     = np.zeros((NCLUSTERS, ft.shape[1], MCMAX) )
    
    
    for imc in range(MCMAX):
        
        # Do not randomize in the first montecarlo, thus, true values are in position 0.
        if imc==0:
            rr_labels = labels.copy()
        else:    
            rr_labels = np.random.permutation(labels)
        
        # For every cluster...
        for i_cluster in range( NCLUSTERS ):
            _rr  = np.argwhere( rr_labels== i_cluster)[:,0]   # find randomized labels 
            _mc[i_cluster,:,imc] = np.nanmean( F[_rr,:], axis=0)        # save randomized values
        
    pval = np.sum( _mc[:,:,:1] <= _mc , axis=2)/MCMAX
    
    
    PVAL_L.append(pval)
    MCR.append( _mc )
# >>>>>>> AQUI >>>>>>>


np.savez('alotofdata.npz', cluster_array=cluster_array, S=S, MCR=MCR, PVAL_L=PVAL_L )





# Print taxon univocity
s = [np.array(_s) for _s in S]
s = [ _s[_s>0] for _s in s ]


x_axis = 1 + 2*np.arange(len(s))


plt.figure(figsize=(6,4), dpi=300)
plt.plot(x_axis, [np.mean(_s) for _s in s],'.r',zorder=99999 )
for jj in range(len(s)):
    plt.plot( (1+2*jj)*np.ones((len(s[jj]))), s[jj],'.', color=np.ones((3,))*0.5, alpha=0.4)
# plt.plot(x_axis, [np.median(_s) for _s in s],'g' )

plt.yscale('log')

plt.xscale('log')
plt.xlabel('Number of clusters')
plt.ylabel('Taxonomic univocity')
plt.legend(('Mean', 'Raw'))



# for cluster_idx in range( labels.max()+1 ):
#     idx = np.argwhere( labels==cluster_idx )[:,0]       # Find elements in cluster cluster_idx
#     group_size = len(idx)                               # Number of elements in cluster cluster_idx
#     obs_mean = np.nanmean(F[idx,:], axis=0)             # Observed mean features across cluster cluster_idx
#     for imc in tqdm(range(MCMAX)):
#         rr = np.random.permutation( M.shape[0] )[:group_size] # Randomize elements in cluster_idx
#         obs_mc[imc,:] = np.nanmean(F[rr,:], axis=0)           # Observed mean features across randomized cluster
        
#     left_pvals = (np.sum( obs_mean <= obs_mc, axis=0)+1)/(MCMAX+1)
#     right_pvals= (np.sum( obs_mean >= obs_mc, axis=0)+1)/(MCMAX+1)
    
#     plt.figure( figsize=(12,8) )
#     for jj in range(15):
#         plt.subplot(3,5,jj+1)
#         plt.hist( obs_mc[:,jj] , bins=11, density=True)
#         xlims, ylims = plt.xlim(), plt.ylim()
#         plt.vlines( obs_mean[jj], 0, 1e3,'r')
#         if (right_pvals[jj]<0.05) or (right_pvals[jj]>0.95):
#             col = 'r'
#         else:
#             col='k'
#         _pval = np.min( [right_pvals[jj], left_pvals[jj] ])
#         plt.text( xlims[0], ylims[1]*0.9, 'pval=%1.4f' % _pval, color=col, fontsize=12)
#         plt.xlabel(col_names[jj])
#         plt.ylim(ylims)
#     plt.suptitle('Cluster %d; Size %d' % (cluster_idx, group_size) )
#     plt.tight_layout()


