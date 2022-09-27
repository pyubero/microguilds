# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:36:26 2022

@author: Pablo
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import pickle as pkl


# Load tree data (functional tree, not phylogenetic)
data = np.load('tree_distance.npz')
M = data['M']
M_names=data['names']


# Load feature matrix
df = pd.read_csv("test_potf.csv", sep=",", encoding='latin1')
ft = df.to_numpy()
ft_names = np.array( [ line[0] for line in ft ] )
col_names= np.array( [ col for col in df.columns] )[1:16]    
ft = np.array([ line[1:16] for line in df ])

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
clustering_api = KMeans(n_clusters=8, # Number of clusters to form
                        n_init = 10,  # Number of times the k-means algorithm will be run with different seeds.
                        random_state=1337)

kmeans = clustering_api.fit(M)
labels = kmeans.labels_


mc_max = 9999
obs_mc = np.zeros((mc_max, ft.shape[1]))



for cluster_idx in range( labels.max()+1 ):
    idx = np.argwhere( labels==cluster_idx )[:,0]       # Find elements in cluster cluster_idx
    group_size = len(idx)                               # Number of elements in cluster cluster_idx
    obs_mean = np.nanmean(F[idx,:], axis=0)             # Observed mean features across cluster cluster_idx
    for imc in tqdm(range(mc_max)):
        rr = np.random.permutation( M.shape[0] )[:group_size] # Randomize elements in cluster_idx
        obs_mc[imc,:] = np.nanmean(F[rr,:], axis=0)           # Observed mean features across randomized cluster
        
    left_pvals = (np.sum( obs_mean <= obs_mc, axis=0)+1)/(mc_max+1)
    right_pvals= (np.sum( obs_mean >= obs_mc, axis=0)+1)/(mc_max+1)
    
    plt.figure( figsize=(12,8) )
    for jj in range(15):
        plt.subplot(3,5,jj+1)
        plt.hist( obs_mc[:,jj] , density=True)
        xlims, ylims = plt.xlim(), plt.ylim()
        plt.vlines( obs_mean[jj], 0, 1e3,'r')
        if (right_pvals[jj]<0.05) or (right_pvals[jj]>0.95):
            col = 'r'
        else:
            col='k'
        plt.text( xlims[0], ylims[1]*0.9, 'pval=%1.4f' % right_pvals[jj], color=col)
        plt.xlabel(col_names[jj])
        plt.ylim(ylims)
    plt.suptitle('Cluster %d' % cluster_idx)
    plt.tight_layout()


