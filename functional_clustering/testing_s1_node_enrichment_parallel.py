# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:36:26 2022

@author: Pablo
"""
# %%
import warnings
import numpy as np
import pandas as pd
from time import sleep
from tqdm import tqdm
from matplotlib import pyplot as plt
from multiprocessing import Pool
from functional_clustering_utils import verboseprint
import functional_clustering_utils as fcutils


def entropy(frequencies):
    '''Computes the entropy of a categorical histogram.'''
    frequencies = np.array(frequencies)
    idc = np.argwhere(frequencies > 0)[:, 0]

    # Normalize ocurrences to probabilities
    p = frequencies[idc]/np.sum(frequencies[idc])

    # Return entropy
    _entropy = np.zeros(frequencies.shape)
    _entropy[idc] = -p*np.log(p)
    return np.sum(_entropy)


def bin_count(array):
    '''Coun bins as to build a categorical histogram.'''
    unique_values = np.unique(array)
    return np.array([np.sum(array == value) for value in unique_values])


# General parameters
GENE = "potF"
FILENAME_ENV_DATA = './data/environmental_data.csv'
MCMAX = 99  # 999 takes 90s; 99 takes 10s
VERBOSE = True
DISPLAY_PLOTS = True
IGNORE_WARNINGS = True

FILENAME_TREE = f'./data/new_{GENE}.newick'
FILENAME_CLADE_DATA = f'./data/data_clades_{GENE}.npz'
FILENAME_OUT = f'./data/data_enrichment_{GENE}.npz'

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
data = fcutils.get_clade_data(FILENAME_TREE,
                              treetype="newick",
                              filename_out=FILENAME_CLADE_DATA)
clade_ids, clade_lfs, clade_dpt, leaf_names = data
# ...
GENUS = np.array([fcutils.get_genus(name) for name in leaf_names])
SPECIES = np.array([fcutils.get_species(name) for name in leaf_names])
# ...
NLEAFS = len(leaf_names)
NCLADES = len(clade_ids)
CLADE_NLEAFS = np.array([len(leafs) for leafs in clade_lfs])

# Declare accumulators
FEATURES = np.zeros((NLEAFS, NFEATURES)) * np.nan
ZSCORES = np.zeros((NCLADES, NFEATURES)) * np.nan
UNIVOCITY = np.zeros((NCLADES,)) * np.nan

# Fill feature matrix of the organisms found in the tree
for idx_sp, species in tqdm(enumerate(SPECIES_NAMES)):
    _speacies_name = species.split(' ')
    for idx_leaf, leaf in enumerate(leaf_names):
        # If "species" is matched to the tree "leaf"
        if (_speacies_name[0] in leaf) and (_speacies_name[1] in leaf):
            FEATURES[idx_leaf, :] = FEATURES_SPECIES[idx_sp, :]

_missing_data = np.sum(np.all(np.isnan(FEATURES), axis=1))
_not_missing_data = FEATURES.shape[0] - _missing_data
idx_missing = np.argwhere(np.all(np.isnan(FEATURES), axis=1))[:, 0]

verboseprint(f'\nData obtained for {_not_missing_data} tree entries.', VERBOSE)
verboseprint(f'Data missing  for {_missing_data} tree entries.', VERBOSE)
_ = [verboseprint(f"  -[{_:4d}] {GENUS[_]} {SPECIES[_]}", VERBOSE)
     for _ in idx_missing]


def parallel_function(ii_clade):
    '''Null model / Randomization function for a given clade idx'''
    _leafs = np.array(clade_lfs[ii_clade])
    _obsmean = np.nanmean(FEATURES[_leafs, :], axis=0)

    # Randomize values
    _mcvalues = np.zeros((MCMAX, NFEATURES))
    for imc in range(MCMAX):
        rr = np.random.permutation(NLEAFS)[:len(_leafs)]
        _mcvalues[imc, :] = np.nanmean(FEATURES[rr, :], axis=0)

    _mcmean = np.nanmean(_mcvalues, axis=0)
    _mcstd = np.nanstd(_mcvalues, axis=0)

    return ii_clade, _obsmean, _mcmean, _mcstd

# from time import sleep
# def parallel_function2(dummy):
#     sleep(1)
#     return True


parallel_function(50)

# %%
# with Pool(processes=4) as p:
#     results = p.imap_unordered(parallel_function2, [10,10,10,10])

#     for res in results:
#         pass

# %%
with Pool() as p:
    results = p.imap_unordered(parallel_function, range(NCLADES), chunksize=20)

    for ii, obs_mean, mc_mean, mc_std in tqdm(results, total=NCLADES):
        # pass
        mc_std[np.abs(mc_std) <= 1e-8] = np.nan
        ZSCORES[ii, :] = (obs_mean - mc_mean) / mc_std

        # Compute univocity
        # _leafs = np.array(clade_lfs[ii])
        # cluster_taxonomy = np.array(
        #     [GENUS[idx] for idx in _leafs if GENUS[idx] is not None])
        # UNIVOCITY[ii] = entropy(bin_count(cluster_taxonomy))


# %% Export and plot data

verboseprint(f"Data exported to {FILENAME_OUT}.", VERBOSE)
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

#%%
# for idx_clade, _leafs in tqdm(enumerate(clade_lfs), total=NCLADES):
    # Compute true observed values for each feature
#     obs_mean = np.nanmean(FEATURES[_leafs, :], axis=0)

    # Randomize values
#     values_mc = np.zeros((MCMAX, NFEATURES))
#     for imc in range(MCMAX):
#         rr = np.random.permutation(NLEAFS)[:len(_leafs)]
#         values_mc[imc, :] = np.nanmean(FEATURES[rr, :], axis=0)

#     mc_mean = np.nanmean(values_mc, axis=0)
#     mc_std = np.nanstd(values_mc, axis=0)
#     mc_std[np.abs(mc_std) <= 1e-8] = np.nan
#     ZSCORES[idx_clade, :] = (obs_mean-mc_mean)/mc_std

    # Compute univocity
#     cluster_taxonomy = np.array(
#         [GENUS[idx] for idx in _leafs if GENUS[idx] is not None])
#     UNIVOCITY[idx_clade] = entropy(bin_count(cluster_taxonomy))
