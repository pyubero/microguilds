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



name1 = M_names[0]
name2 = M_names[1]

sp_name = name1.split('_g_')[-1]


def get_species(name):
    sp_name = name.split('_s_')[-1]
    if sp_name == '':
        print('Species name not found.')
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
    return [ np.sum(array==value) for value in unique_values]




# Load tree data (functional tree, not phylogenetic)
data = np.load('tree_distance.npz')
M = data['M']
M_names=data['names']

# Load feature matrix
df = pd.read_csv("test2_potf.csv", sep=",", encoding='latin1')
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

clustering_api = KMeans(n_clusters= 3, # Number of clusters to form
                        n_init = 20,  # Number of times the k-means algorithm will be run with different seeds.
                        random_state=1337)


kmeans = clustering_api.fit(M)
labels = kmeans.labels_
    



#
genus_list = np.array( [ get_genus(name) for name in M_names] )
species_list= np.array( [ get_species(name) for name in M_names] )
S=[]

for cluster_idx in range(labels.max()+1):
    idx = np.argwhere( labels==cluster_idx )[:,0]
    
    
    cluster_genus = np.array( [ genus for genus in genus_list[idx] if genus is not None] )
    h = bin_count(cluster_genus)
    S.append( 1/entropy( h ) )










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
    plt.suptitle('Cluster %d; Size %d' % (cluster_idx, group_size) )
    plt.tight_layout()


