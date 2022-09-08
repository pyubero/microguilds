# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""


import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


GENE_NAME   = 'amt'
K_TRANSFORM = lambda x: x
# K_TRANSFORM = lambda x: np.log10(x)




_filename = f'kMatrixPerTaxon_{GENE_NAME}.csv'



taxons = []
Kmat   = [] # Kmat will be of size ( j_taxon, j_context, j_cluster)


with open( _filename, 'r') as f:
    
    # Check if file is correct
    if f.readline().replace('\n','').lower() != GENE_NAME.lower():
        print('<< W >> Gene input does not match with gene name from file.')
        
    # Extract the names of functional clusters and the number        
    clusters = f.readline().replace('\n','').split(',')
    n_clusters = len(clusters)
    
    line = f.readline()
    
    while line!='':
    
        if line[:2] == '>>':
            #Start new matrix and add taxon name
            _K = []
            taxons.append( line.replace('\n','').replace('>',''))
            print('Starting with %s' % taxons[-1])
        
        line = f.readline()

        while line[:2] !='>>':            
            _K.append( [ float(_) for _ in  line.replace('\n','').split(',') ] )
            line = f.readline()
            
            if line=='':
                break
            
        Kmat.append( _K )
        
Kmat = K_TRANSFORM( np.array(Kmat) )
n_taxons, n_ctxts, n_clusters = Kmat.shape



MIN_CONTRIBUTION = 50.50 # in percent

# Find taxons that contribute less than X percent
K_sum_taxons = np.sum(Kmat, axis=0)
others = []
Kothers = np.zeros_like( Kmat[0,:,:] )


for j_taxon in range( n_taxons ):
    
    tx_contrib = Kmat[j_taxon,:,:] #/K_sum_taxons
    tx_contrib[ np.isnan(tx_contrib) ] = 0

    if np.all( tx_contrib < MIN_CONTRIBUTION ):
        Kothers += Kmat[ j_taxon, : , : ]
        others.append( j_taxon )
        print('%s is now considered as others' % taxons[j_taxon] )
    





















