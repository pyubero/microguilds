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
from tqdm import tqdm

def verboseprint(msg):
    if VERBOSE:
        print(msg)
        
FILENAME    = 'mastertable_w_ctxt.tsv'
GENE_NAME   = 'nirs'
LEVEL_NAME  = 'Species_GTDB'
EXPORT_DATA = True 
VERBOSE = True


# ...
out_filename = f'kMatrixPerTaxon_{GENE_NAME}.csv'



def Kfun(A,C):
    return A*C



# Import mastertable
df = pd.read_csv(FILENAME, sep="\t")
verboseprint( f"Loaded mastertable with {len(df)} rows." )


# Filter by gene name
df_gene = df[ df ["gene_fun"] == GENE_NAME ]
verboseprint( f"Subtable for gene *{GENE_NAME}* has {len(df_gene)} rows." )

# Find all clusters in gene subtable
clusters   = df_gene['cluster_id'].unique()
n_clusters = len(clusters)
verboseprint( f"Found {n_clusters} clusters in gene subtable:" )
_=[ verboseprint( "\t {}".format(clustername)) for clustername in clusters]


# Find all taxons in gene subtable
taxons = df_gene[LEVEL_NAME].unique()
n_taxons= len(taxons)
verboseprint( f"There are {n_taxons} taxonomic levels according to column *{LEVEL_NAME}*.")


# Find all contexts in gene subtable
contexts = df_gene["Context"].unique()
n_contexts = len(contexts)
verboseprint( f"There are {n_contexts} contexts in gene subtable.")
_=[ verboseprint( "\t {}".format(ctxt)) for ctxt in contexts ]


# Create empty output file
if EXPORT_DATA:
    with open(out_filename, 'w+') as f: 
        f.write(f'{GENE_NAME}\n')
        f.write( ','.join( [ "%s" % _ for _ in clusters ] )+'\n' )


# Keep track of unprocessed samples and abundances
unprocessed_samples=[]


# For every taxon...
for j_taxon in tqdm(range( n_taxons)):
    
    #... create a K matrix
    K = np.zeros((3, n_clusters)) #<- the 3 is because there are 3 contexts (bathy, meso, epi)
    
    # For every clusters...
    for j_cluster in range(n_clusters):
        
        # Filter table by functional cluster and by taxon
        subtable = df_gene.copy()
        subtable = subtable[ subtable["cluster_id"] == clusters[j_cluster] ]
        subtable = subtable[ subtable[LEVEL_NAME] == taxons[j_taxon] ]

        
        # Initialize counts and abundances
        counts = np.zeros((3,)) #[0,0,0]
        abund  = np.zeros((3,)) #[0,0,0]
        
        # For every sequence...
        for _, row in subtable[subtable[LEVEL_NAME] ==taxons[j_taxon]].iterrows():
            sampleid = row["Sample_ID"]
            
            
            if sampleid in epi_list:
                counts[0] += 1
                abund[0] += row["value"]
                
            elif sampleid in meso_list:
                counts[1] += 1
                abund[1] += row["value"]
        
            elif sampleid in bathy_list:
                counts[2] += 1
                abund[2] += row["value"]
        
            else:
                unprocessed_samples.append(sampleid)

        # Compute K values for this gene>taxon>cluster
        K[:,j_cluster] = np.array([ Kfun(a, c)  for c,a in zip(counts, abund) ])

    
    # Append to output file...
    if EXPORT_DATA:
        header = f'>>{taxons[j_taxon]}\n'
        lines = []
        for _ in range( K.shape[0]) :
            lines.append( ','.join( [ '%1.4f' % k for k in K[_,:]])+'\n' )
                    
        with open(out_filename, 'a') as f: 
            f.write(header)
            _ = [ f.write(line) for line in lines ]
            
            
# Print final comments.            
if len(unprocessed_samples)==0:
    print('')
    print('All sampled were processed correctly')
else:    
    for _ in np.unique( unprocessed_samples):
        print('Sample ID: %s could not be processed.' % _ )
