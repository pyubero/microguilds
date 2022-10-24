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


FILENAME_CLADE_DATA = 'data_potF_all_clades.npz'
FILENAME_ENV_DATA   = 'environmental_data.csv'
FILENAME_OUT = 'data_clade_enrichment.npz'
MCMAX = 99 # 999 takes 220s


# Load all clade data
data = np.load(FILENAME_CLADE_DATA, allow_pickle=True)
leaf_names = data['leaf_names']
clade_ids  = data['clade_ids']
clade_dpt  = data['clade_dpt'] 
clade_lfs  = data['clade_lfs']


# Load environmental feature matrix
df = pd.read_csv( FILENAME_ENV_DATA , sep=",", encoding='latin1')
ft = df.to_numpy()
ft_sp_names = np.array( [ line[0] for line in ft ] )
col_names= np.array( [ col for col in df.columns] )[1:16]    
ft = np.array([ line[1:16] for line in ft ])
#...
nleafs    = len(leaf_names)
nclusters = len(clade_ids)
nfeatures = ft.shape[1]
#...
clade_nleafs= np.array( [len(leafs) for leafs in clade_lfs])
genus_list  = np.array( [ get_genus(name) for name in leaf_names] )
species_list= np.array( [ get_species(name) for name in leaf_names] )


#%%
# ... create true feature matrix
F = np.zeros((nleafs, nfeatures))*np.nan

for idx_ft_sp, species in enumerate(ft_sp_names):
    _name = species.split(' ')
    for idx_leaf, full_name in enumerate(leaf_names):
        if (_name[0] in full_name) and (_name[1] in full_name):
            F[idx_leaf,:] = ft[idx_ft_sp,:]

print('Data obtained for %d tree entries.' % (1158-np.sum(np.all( np.isnan(F), axis=1))))            

import warnings
warnings.filterwarnings('ignore')




ZSCORES = np.zeros( (nclusters, nfeatures))*np.nan
PVALS   = np.zeros( (nclusters, nfeatures))*np.nan


# First obtain clusters
for ii, cluster_leafs in tqdm(enumerate(clade_lfs), total=len(clade_lfs)):


    # Compute true observed value
    obs_mean = np.nanmean(F[ cluster_leafs,:], axis=0)             # Observed mean features across cluster cluster_idx
    
    # Randomize values
    obs_mc   = np.zeros((MCMAX, nfeatures))
    for imc in range(MCMAX):
        rr = np.random.permutation(nleafs)[:len(cluster_leafs)]
        obs_mc[imc,:] = np.nanmean(F[ rr ,:], axis=0)

    mc_mean = np.nanmean(obs_mc, axis=0)
    mc_std  = np.nanstd(obs_mc, axis=0)
    mc_std[ np.abs(mc_std)<=1e-8]=np.nan
    ZSCORES[ii,:] = (mc_mean-obs_mean)/mc_std
    PVALS[ii,:]
    
    
# Compute univocity
S = np.zeros((nclusters,))*np.nan
for ii, cluster_leafs in enumerate(clade_lfs):
    cluster_taxonomy = np.array( [ genus_list[idx] for idx in cluster_leafs if genus_list[idx] is not None ] )
    h = bin_count( cluster_taxonomy )
    S[ii] = entropy( h )

    
    
    
# np.savez(FILENAME_OUT,
#           names = ft_sp_names,
#           features= col_names,
#           F = F,
#           S = S,
#           MCMAX = MCMAX,
#           ZSCORES=ZSCORES)  
    
#%%

# Entropy/Univocity correlates with the number of leafs within clade
# plt.figure( figsize=(6,4), dpi=300)
# plt.plot(S, clade_nleafs,'.')
# plt.xlabel('Diversity (S)')
# plt.ylabel('# of terminal leafs')    



# Massive plot of z-scores
# plt.figure( figsize=(15,8), dpi=300)    
# for jj in range(15):
#     plt.subplot(3,5,jj+1)
#     plt.scatter( clade_nleafs, ZSCORES[:,jj], s=8); 
#     xlim = plt.xlim()
    
#     plt.hlines(2, xlim[0], xlim[1] ,'r', lw=0.7)
#     plt.hlines(-2, xlim[0], xlim[1] ,'r', lw=0.7)
#     plt.grid()
    
#     if jj>9:
#         plt.xlabel('# terminal leafs')
        
#     if np.mod(jj,5)==0:
#         plt.ylabel('abs(Z-score)')
        
#     plt.title('%s' % col_names[jj])
#     # plt.xlim(-0.05, 5.5 )
#     # plt.ylim(0,5)

# plt.suptitle('Clade enrichment')
# plt.tight_layout()


data = np.load(FILENAME_OUT, allow_pickle=True)
F = data['F']
S = data['S']
ZSCORES= data['ZSCORES']

# Find unique patterns
THOLD_Z = 3
barcode = 1*(ZSCORES>THOLD_Z) - 1*(ZSCORES<-THOLD_Z)


unq_bc = np.unique(barcode, axis=0)
for jj, ref_bc in enumerate(unq_bc):
    # ref_bc = unq_bc[23]
    idc = np.argwhere( np.all(barcode==ref_bc, axis=1))[:,0]
    



# for jj, ref_bc in enumerate(unq_bc):
# ref_bc = unq_bc[23]
# idc = np.argwhere( np.all(barcode==ref_bc, axis=1))[:,0]


# start_idx = idc[0]
# for test_idx in idc[1:]:
#   if set( clade_lfs[start_idx]) > set(clade_lfs[test_idx]):
#       # print('yay', test_idx)
#       start_idx = test_idx
# print('Ref_bc ', jj, ',node idx ',start_idx, ', # of leafs:', len(clade_lfs[start_idx]),
#       '\t',ref_bc)      

  
    
  
    
  
    
  