# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:36:26 2022

@author: Pablo
"""
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from matplotlib import pyplot as plt
from functional_clustering_utils import verboseprint
import functional_clustering_utils as fcutils

# def accession_to_name(query, fastafile):
#     with open(fastafile,'r') as file:
#         for line in file.readlines():
#             if np.all( [ _ in line for _ in query ]):
#                 return line
#     return None


# # Create dictionary
# acc2names = {}
# with open('seq_recA_std.fasta','r') as file:
#     for line in file.readlines():
#         if line[0]=='>':
#             line = line.replace('>','').split(' ')
#             accession = line[0]
#             name = ' '.join( line[1:3])
#             acc2names.update( {accession : name} )

def entropy(f):
    '''Computes the entropy of a categorical histogram.'''
    f = np.array(f)
    idc = np.argwhere(f > 0)[:, 0]

    # Normalize ocurrences to probabilities
    p = f[idc]/np.sum(f[idc])

    # Return entropy
    _entropy = np.zeros(f.shape)
    _entropy[idc] = -p*np.log(p)
    return np.sum(_entropy)


def bin_count(array):
    '''Coun bins as to build a categorical histogram.'''
    unique_values = np.unique(array)
    return np.array([np.sum(array == value) for value in unique_values])


FILENAME_TREE = 'tree_potF.newick'
FILENAME_ENV_DATA = 'environmental_data.csv'
FILENAME_OUT = 'data_enrichment_potF.npz'
MCMAX = 1999  # 999 takes 220s; 99 takes 30s
VERBOSE = True
DISPLAY_PLOTS = True
IGNORE_WARNINGS = True

if IGNORE_WARNINGS:
    warnings.filterwarnings('ignore')

# Load environmental feature matrix
df = pd.read_csv(FILENAME_ENV_DATA, sep=",", encoding='latin1')
ft = df.to_numpy()
# ...
SPECIES_NAMES = np.array([line[0] for line in ft])
FEATURES_NAMES = np.array([col for col in df.columns])[1:16]
FEATURES_SPECIES = np.array([line[1:16] for line in ft])
NFEATURES = FEATURES_SPECIES.shape[1]

# Load clade data
data = fcutils.get_clade_data(FILENAME_TREE, treetype="newick")
clade_ids, clade_lfs, clade_dpt, leaf_names = data
# ...
NLEAFS = len(leaf_names)
NCLADES = len(clade_ids)
CLADE_NLEAFS = np.array([len(leafs) for leafs in clade_lfs])

# Declare accumulators
# PVALS = np.zeros((nclades, NFEATURES))*np.nan
FEATURES = np.zeros((NLEAFS, NFEATURES))*np.nan
ZSCORES = np.zeros((NCLADES, NFEATURES))*np.nan
UNIVOCITY = np.zeros((NCLADES,))*np.nan

# Uncomment for potF
GENUS = np.array([fcutils.get_genus(name) for name in leaf_names])
SPECIES = np.array([fcutils.get_species(name) for name in leaf_names])

# Uncomment for recA
# Create dictionary
# acc2names = {}
# with open('seq_rplB_std.fasta','r') as file:
#     for line in file.readlines():
#         if line[0]=='>':
#             line = line.replace('>','').split(' ')
#             accession = line[0]
#             name = ' '.join( line[1:3])
#             acc2names.update( {accession : name} )
# GENUS   = np.array( [acc2names[name].split(' ')[0] for name in leaf_names] )
# SPECIES = np.array( [acc2names[name].split(' ')[1] for name in leaf_names] )


# Fill feature matrix of the organisms found in the tree
for idx_sp, species in tqdm(enumerate(SPECIES_NAMES)):
    _speacies_name = species.split(' ')
    for idx_leaf, leaf in enumerate(leaf_names):
        # If "species" is matched to the tree "leaf"
        if (_speacies_name[0] in leaf) and (_speacies_name[1] in leaf):
            FEATURES[idx_leaf, :] = FEATURES_SPECIES[idx_sp, :]

        # Uncomment this when analysing recA data
        # ... or comment this when analysing potF
        # full_name = acc2names[full_name]

_missing_data = np.sum(np.all(np.isnan(FEATURES), axis=1))
_not_missing_data = FEATURES.shape[0] - _missing_data
verboseprint('')
verboseprint(f'Data obtained for {_not_missing_data} tree entries.', VERBOSE)
verboseprint(f'Data missing  for {_missing_data} tree entries.', VERBOSE)

# Uncomment to print missing species <- UNAVAILABLE
# idx_missing = np.argwhere( np.all( np.isnan(F), axis=1))[:,0]
# _=[print('  -%s' % acc2names[ leaf_names[_] ]) for _ in idx_missing]

for idx_clade, _leafs in tqdm(enumerate(clade_lfs), total=NCLADES):
    # Compute true observed values for each feature
    obs_mean = np.nanmean(FEATURES[_leafs, :], axis=0)

    # Randomize values
    values_mc = np.zeros((MCMAX, NFEATURES))
    for imc in range(MCMAX):
        rr = np.random.permutation(NLEAFS)[:len(_leafs)]
        values_mc[imc, :] = np.nanmean(FEATURES[rr, :], axis=0)

    mc_mean = np.nanmean(values_mc, axis=0)
    mc_std = np.nanstd(values_mc, axis=0)
    mc_std[np.abs(mc_std) <= 1e-8] = np.nan
    ZSCORES[idx_clade, :] = (obs_mean-mc_mean)/mc_std

    # Compute univocity
    cluster_taxonomy = np.array(
        [GENUS[idx] for idx in _leafs if GENUS[idx] is not None])
    UNIVOCITY[idx_clade] = entropy(bin_count(cluster_taxonomy))

# Export data
np.savez(FILENAME_OUT,
         species_names=SPECIES_NAMES,
         features_names=FEATURES_NAMES,
         features=FEATURES,
         univocity=UNIVOCITY,
         mcmax=MCMAX,
         zscores=ZSCORES)

if DISPLAY_PLOTS:
    # Load data
    data = np.load(FILENAME_OUT, allow_pickle=True)
    SPECIES_NAMES = data['species_names']
    FEATURES_NAMES = data['features_names']
    FEATURES = data['features']
    UNIVOCITY = data['univocity']
    ZSCORES = data['zscores']

    # Entropy/Univocity correlates with the number of leafs within clade
    plt.figure(figsize=(6, 4))
    plt.plot(UNIVOCITY, CLADE_NLEAFS, '.')
    plt.xlabel('Diversity (S)')
    plt.ylabel('# of terminal leafs')
    plt.show()

    # Massive plot of z-scores
    plt.figure(figsize=(12, 6))
    for jj in range(15):
        plt.subplot(3, 5, jj+1)
        plt.scatter(CLADE_NLEAFS, ZSCORES[:, jj], s=8)
        xlim = plt.xlim()

        plt.hlines(2, xlim[0], xlim[1], 'r', lw=0.7)
        plt.hlines(-2, xlim[0], xlim[1], 'r', lw=0.7)
        plt.grid()

        if jj > 9:
            plt.xlabel('# terminal leafs')

        if np.mod(jj, 5) == 0:
            plt.ylabel('abs(Z-score)')

        plt.title(f'{FEATURES_NAMES[jj]}')
        # plt.xlim(-0.05, 5.5 )
        # plt.ylim(0,5)

    plt.suptitle('Clade enrichment')
    plt.tight_layout()
    plt.show()
