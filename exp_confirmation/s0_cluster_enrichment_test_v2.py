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



# Load all clade data
data = np.load('all_clade_data.npz', allow_pickle=True)
leaf_names = data['leaf_names']
clade_ids  = data['clade_ids']
clade_dpt  = data['clade_dpt'] 
clade_lfs  = data['clade_lfs']


# Load feature matrix
df = pd.read_csv("test_potf.csv", sep=",", encoding='latin1')
ft = df.to_numpy()
ft_sp_names = np.array( [ line[0] for line in ft ] )
col_names= np.array( [ col for col in df.columns] )[1:16]    
ft = np.array([ line[1:16] for line in ft ])
#...
nleafs    = len(leaf_names)
nclusters = len(clade_ids)
nfeatures = ft.shape[1]

# ... create true feature matrix
F = np.zeros((nleafs, nfeatures))*np.nan

for idx_ft_sp, species in enumerate(ft_sp_names):
    _name = species.split(' ')
    for idx_leaf, full_name in enumerate(leaf_names):
        if (_name[0] in full_name) and (_name[1] in full_name):
            F[idx_leaf,:] = ft[idx_ft_sp,:]

print('Data obtained for %d tree entries.' % (1158-np.sum(np.all( np.isnan(F), axis=1))))            





MCMAX = 999 # 999 takes 220s
PVALS_L = np.zeros( (nclusters, nfeatures))*np.nan
PVALS_R = np.zeros( (nclusters, nfeatures))*np.nan

# First obtain clusters
for ii, cluster_leafs in tqdm(enumerate(clade_lfs)):
    if len(cluster_leafs)<5:
        continue
    
    # Compute true observed value
    obs_mean = np.nanmean(F[ cluster_leafs,:], axis=0)             # Observed mean features across cluster cluster_idx
    
    # Randomize values
    obs_mc   = np.zeros((MCMAX, nfeatures))
    for imc in range(MCMAX):
        rr = np.random.permutation(nleafs)[:len(cluster_leafs)]
        obs_mc[imc,:] = np.nanmean(F[ rr ,:], axis=0)

    PVALS_L[ii,:] = (np.sum( obs_mean <= obs_mc, axis=0)+1)/(MCMAX+1)
    PVALS_R[ii,:] = (np.sum( obs_mean >= obs_mc, axis=0)+1)/(MCMAX+1)