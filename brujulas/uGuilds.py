# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:48:47 2022

@author: Pablo Yubero

This script reads a master table and outputs a file with K values for the 
    specified GENE_NAME per CONTEXT and LEVEL (taxonomic). Values of K are 
    computed as a function of the abundances (A) and counts of sequences (C) 
    as per Kfun(A,C).
    
    
    
* INPUTS *

FILENAME  : str :  Input of master table (as per Juan Rivas)
GENE_NAME : str :  Gene of interest to compute the K values
LEVEL_NAME: str :  Taxonomic level to generate K values, eg, Species_GTDB     
    
"""

import pandas as pd
import numpy as np
from tqdm import tqdm


FILENAME    = 'master_tab.csv'
GENE_NAME   = 'nirs'
LEVEL_NAME  = 'Species_GTDB'
EXPORT_DATA = True 
# ...
out_filename = f'kMatrixPerTaxon_{GENE_NAME}.csv'



bathy_list  = ["MP0313","MP0315","MP0530","MP0532","MP0534","MP0780","MP0782",
               "MP784", # <<---- missing 0?
               "MP0784","MP0878","MP0880","MP0882","MP1154","MP1162","MP1409",
               "MP1411","MP1519","MP1521","MP1676","MP1674","MP1847","MP1845",
               "MP2231","MP2233bis","MP2809","MP2811" ]


meso_list = ["MP0317","MP0319","MP0536","MP0538","MP0786","MP0788","MP0884",
             "MP0886","MP1164","MP1166","MP1178","MP1413","MP1415","MP1417",
             "MP1523","MP1164","MP1525","MP1677","MP1678","MP1680","MP1681",
             "MP1682","MP1849","MP1851","MP1853","MP2235","MP2237","MP2817",
             "MP2813","MP2815" ]

epi_list= ["MP2239","MP2241","MP0323", "MP2819", "MP1857","MP0311","MP1419",
           "MP1421","MP1517","MP1527","MP1529","MP0790","MP0888","MP0778",
           "MP0528","MP0540","MP1176", "MP1174", "MP2821","MP1855" , "MP1684",
           "MP0321","MP2243","MP1672"]



def Kfun(A,C):
    return A*C



# Import mastertable
df = pd.read_csv(FILENAME)


# Filter by gene name
df_gene = df[ df ["gene_fun"] == GENE_NAME ]


# Find all clusters in gene subtable
clusters   = df_gene['cluster_id'].unique()
n_clusters = len(clusters)


# Find all taxons in gene subtable
taxons = df_gene[LEVEL_NAME].unique()
n_taxons= len(taxons)


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
