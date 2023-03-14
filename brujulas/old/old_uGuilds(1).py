# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:48:47 2022

@author: Pablo Yubero

This script reads a master table and outputs a file with K values for the 
specified GENE_NAME per CONTEXT and LEVEL (taxonomic). Values of K are 
computed as a function of the abundances (A) and counts of sequences (C) 
as per Kfun(A,C).
    
The mastertable should have AT LEAST the following columns:
    gene_fun : gene function by SQM or by who?? <<<
    
    cluster_id : cluster identifier according to who??? <<<
    
    Species_GTDB : Although this is the defualt, it can change depending on the
        taxonomic level to group K values desired. (Related to LEVEL_NAME input)
    
    
    
* INPUTS *

FILENAME  : str
    Input of master table (as per Juan Rivas)

GENE_NAME : str 
    Gene of interest to compute the K values

LEVEL_NAME: str 
    Taxonomic level to generate K values, eg, Species_GTDB. It NEEDS to be the
    name of a column of master table.
    
SEP : str
    default : "\t" 
    
"""

import pandas as pd
import numpy as np

def verboseprint(msg):
    if VERBOSE:
        print(msg)
      
def verbosebar(iterable):
    if VERBOSE:
        return tqdm(iterable)
    else:
        return iterable

def Kfun(A,D,U):
    # Abundance, Diversity, Univocity
    return A*U

def compute_adu(df, inputs, loc):
    # Find where the interesting abundances are in the df/master table
    idc = ( df["taxonomic_classification_level"] == inputs[0][loc[0]] ) & \
          ( df["Context"] ==  inputs[1][loc[1]]) & \
          ( df["cluster_id"] ==  inputs[2][loc[2]] ) 
          
    # Extract value
    subtable = df[ idc ]
    
    # Obtain abundances and relevant parameters
    abundances = subtable["value"]
    diversity = len(abundances)
    abundance = sum(abundances)
    univocity = 1.0
    return abundance, diversity, univocity

    
# How to nicely load our tensor?
def from_df_to_ktensor(df, column="k_value"):
    # you should call before:
    # df = pd.read_csv(f"kValuesPerTaxon_{GENE_NAME}.tsv", sep="\t")
    ntaxons = len(df["Taxon"].unique())
    ncontexts = len(df["Context"].unique())
    nclusters = len(df["Cluster"].unique())
    K = df[column].to_numpy().reshape(nclusters, ntaxons, ncontexts)

    K = np.moveaxis(K2, 0,2)
    return K




FILENAME    = 'mastertable_w_ctxt.tsv'
GENE_NAME   = 'amoa'
LEVEL_NAME  = 'Species_GTDB'
VERBOSE = True
# ...
out_filename = f'kMatrixPerTaxon_{GENE_NAME}_{LEVEL_NAME}_v2.csv'



if VERBOSE:
    from tqdm import tqdm


# Import mastertable
master_table = pd.read_csv(FILENAME, sep="\t")
verboseprint( f"Loaded mastertable with {len(master_table)} rows." )

#... and standardize taxonomic column name
master_table = master_table.rename( columns = {LEVEL_NAME : "taxonomic_classification_level"} )


# Filter by gene name
gene_table = master_table[ master_table ["gene_fun"] == GENE_NAME ]
verboseprint( f"Subtable for gene *{GENE_NAME}* has {len(gene_table)} rows." )

# Find all clusters in gene subtable
clusters   = gene_table['cluster_id'].unique()
n_clusters = len(clusters)
verboseprint( f"Found {n_clusters} clusters in gene subtable:" )
_=[ verboseprint( "\t {}".format(clustername)) for clustername in clusters]


# Find all taxons in gene subtable
taxons = gene_table["taxonomic_classification_level"].unique()
n_taxons= len(taxons)
verboseprint( f"There are {n_taxons} taxonomic levels according to column *{LEVEL_NAME}*.")


# Find all contexts in gene subtable
contexts = gene_table["Context"].unique()
n_contexts = len(contexts)
verboseprint( f"There are {n_contexts} contexts in gene subtable.")
_=[ verboseprint( "\t {}".format(ctxt)) for ctxt in contexts ]


# Process data
data = [ ]
idc = np.array(np.meshgrid( range(n_taxons), range(n_contexts), range(n_clusters))).T.reshape(-1,3)
included = gene_table["taxonomic_classification_level"] == np.random.rand() # to keep track of included and missing samples

for j_tx, j_ct, j_cl in verbosebar(idc):
    # Compute K value
    _ab, _dv, _un = compute_adu(gene_table, [taxons, contexts, clusters], [j_tx, j_ct, j_cl] )
    _k = Kfun(_ab, _dv, _un)
    
    
    # Store K value
    data.append( [taxons[j_tx], contexts[j_ct], clusters[j_cl], _ab, _dv, _un, _k] )
    
    # To keep track of included and missing samples    
    _id = ( gene_table["taxonomic_classification_level"] == taxons[j_tx] ) & \
           ( gene_table["Context"] ==  contexts[j_ct]) & \
           ( gene_table["cluster_id"] ==  clusters[j_cl] )
    included = included | _id


# Finishing comment
verboseprint(f"\nIncluded {sum(included)} out of {len(gene_table)} samples.")    


# Export data
# ddf = pd.DataFrame(data, columns=["Taxon","Context","Cluster","Abundance","Diversity","Univocity","k-value"])
# ddf.to_csv(out_filename, sep="\t", index=False)



