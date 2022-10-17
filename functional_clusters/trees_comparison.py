# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:07:26 2022

@author: Pablo
"""

from Bio import Phylo
from tqdm import tqdm 
import numpy as np
from matplotlib import pyplot as plt


FILENAME_GENE_SEQS = 'pot16s.fasta'
FILENAME_CLADE_DATA= 'data_potF_all_clades.npz'
FILENAME_TREE_16s  = 'tree_16s_labelled.newick'
FILENAME_TREE_GENE = 'tree_potF_labelled.newick' 
FILENAME_OUT = 'data_tree_comparison.npz'

def search_accession(name):
    if (name[0] == None) or (name[1] == None):
        return None
    
    with open(FILENAME_GENE_SEQS,'r') as file:
        for line in file.readlines():    
            if (name[0] in line) and (name[1] in line):
                return line.split(' ')[0][1:]

def get_species(name):
    sp_name = name.split('_s_')[-1]
    if sp_name == '':
        return None
    else:
        return sp_name.split('_')[-1]
    
def get_genus(name):
    sp_name = name.split('_g_')[-1].split('_')[0]
    if sp_name=='':
        return None
    else: 
        return sp_name
    
#%%
# Load all clade data
data = np.load(FILENAME_CLADE_DATA, allow_pickle=True)
leaf_names = data['leaf_names']
clade_ids  = data['clade_ids']
clade_dpt  = data['clade_dpt'] 
clade_lfs  = data['clade_lfs']




# Load trees
t_16s = Phylo.read(FILENAME_TREE_16s , 'newick')
depths_16s = t_16s.depths()
leafs_16s  = t_16s.get_terminals()
leaf_16s_names= np.array( [ leaf.name for leaf in leafs_16s ] )
depths_16s = t_16s.depths()


#...
t_gene = Phylo.read( FILENAME_TREE_GENE, 'newick')
depths_gene = t_gene.depths()
leafs_gene  = t_gene.get_terminals()
leaf_gene_names= np.array( [ leaf.name for leaf in leafs_gene ] )
depths_gene = t_gene.depths()

# Output data
depth_lca_16s = []
depth_lca_gene = []
leafs_idx_lca  = []
n_bichos = []



# For every internal node in the gene tree...
for clade_idx in tqdm( range( len(clade_ids) )):
    
    # Encontrar leafs del clado clade_idx
    subleafs = clade_lfs[clade_idx]
    leafs_lca_16s = []
    leafs_lca_gene= []
    
    #... for every terminal leaf in the the clade_idx
    for leaf_idx  in subleafs:
        name_leaf= leaf_names[leaf_idx]
        name_bicho = [ get_genus(name_leaf), get_species(name_leaf) ]
        accession  = search_accession(name_bicho)
        
        if accession is not None:
            idx_in_16s = np.argwhere( leaf_16s_names == accession)[:,0]
            leafs_lca_16s.append( idx_in_16s[0]  )
            leafs_lca_gene.append(leaf_idx )            
        
    leafs_lca_16s  = np.unique(leafs_lca_16s)
    leafs_lca_gene = np.unique(leafs_lca_gene)
    
    if len(leafs_lca_16s) > 0:    
        lca_16s = t_16s.common_ancestor( [ leafs_16s[idx] for idx in leafs_lca_16s ] )
        lca_gene= t_gene.common_ancestor([ leafs_gene[idx] for idx in leafs_lca_gene ] )
        d1 = depths_16s[lca_16s] 
        d2 = depths_gene[lca_gene] 
    else:
        d1, d2 = np.nan, np.nan
        
    
    depth_lca_16s.append( d1 )
    depth_lca_gene.append( d2)
    leafs_idx_lca.append( leafs_lca_gene )
    n_bichos.append( len( subleafs) )
    # print('16S:  #leafs: %d \t depth: %1.4f' % ( len(leafs_lca_16s), d1 ) )
    # print('Gene: #leafs: %d \t depth: %1.4f' % ( len(leafs_lca_gene), d2 ) ) 
    # print('-----')
    
n_bichos = np.array(n_bichos)
depth_lca_16s = np.array( depth_lca_16s)
depth_lca_gene= np.array( depth_lca_gene)
leafs_idx_lca = np.array( leafs_idx_lca )

np.savez(FILENAME_OUT, n_bichos = n_bichos, 
                       depth_lca_16s = depth_lca_16s,
                       depth_lca_gene = depth_lca_gene,
                       leafs_idx_lca = leafs_idx_lca)


#%% PLOTS

# Load tree data
data = np.load( FILENAME_OUT, allow_pickle=True)
n_bichos = data['n_bichos']
depth_lca_16s = data['depth_lca_16s']
depth_lca_gene = data['depth_lca_gene']
leafs_idx_lca = data['leafs_idx_lca']

# Load enrichment data
data = np.load('data_clade_enrichment.npz')
F = data['F']
S = data['S']
ZSCORES = data['ZSCORES']
idx_significant = np.argwhere( np.any(np.abs(ZSCORES)>3, axis=1))[:,0]


depth_lca_16s += np.random.randn( *depth_lca_16s.shape )/50000
depth_lca_gene+= np.random.randn( *depth_lca_gene.shape )/50000

plt.figure( figsize=(6,4), dpi=300)      
plt.plot(  depth_lca_16s, depth_lca_gene ,'.', color=np.ones((3,))*0.5, ms=5, zorder=0 )
plt.scatter( depth_lca_16s[idx_significant], depth_lca_gene[idx_significant],c = np.log10(1+np.array(n_bichos[idx_significant])), s=6)    


plt.colorbar( label=' log10 Number of leafs in potF')
plt.xlabel('Taxonomic relatedness')
plt.ylabel('Relatedness')
