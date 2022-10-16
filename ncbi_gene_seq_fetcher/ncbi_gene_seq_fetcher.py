# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 11:14:13 2022

@author: Pablo
"""

import requests
import pycountry
import numpy as np
import pandas as pd
from time import sleep
import xmltodict as xml2dict
from matplotlib import pyplot as plt
from tqdm import tqdm

GENE_OF_INTEREST = '16S Ribosomal rna'
FILENAME_IN = 'organismpotf.csv'
#...
GENE_OF_INTEREST = '_'.join( GENE_OF_INTEREST.split(' '))
FILENAME_OUT = 'seq_'+GENE_OF_INTEREST+'.fasta'
FILENAME_LOG = 'log_'+GENE_OF_INTEREST+'.log'
# ...
YOUR_EMAIL = 'pabloyub@protonmail.com'
TOOL = 'gene_fetcher'
DELAY = 1/2 # in seconds
# ...
URL_ELINK = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
URL_ESUMM = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
URL_EFETCH= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
URL_ESEARCH='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

def print_to_log( data ):
    with open(FILENAME_LOG , 'a+') as file:
        file.write(data+'\n')
    print(data)
    
    
# Create empty files
with open(FILENAME_OUT , 'w+') as file:    pass
with open(FILENAME_LOG , 'w+') as file:    pass

# Get list of organisms' names
organisms_names = []
with open(FILENAME_IN) as file:
    for line in file.readlines():
        organisms_names.append( line.replace('\n','') )



for jj, organism in enumerate(organisms_names):

    search_term = organism+' '+GENE_OF_INTEREST
    PARAMS_SEARCH = {
        'db' : 'nuccore', # or nucleotide
        'term' : search_term,
        'retmode' : 'json'
        }
    
    
    response = requests.get( URL_ESEARCH, params = PARAMS_SEARCH, timeout = 10)
    json = response.json()
    idlist = json['esearchresult']['idlist']
    sleep(DELAY)
    
    
    
    PARAMS_FETCH = {
        'db' : 'nuccore', # or nucleotide
        'id' : ','.join(idlist),
        'rettype' : 'fasta',
        'retmode' : 'text'
        }
    
    response = requests.get( URL_EFETCH, params = PARAMS_FETCH, timeout = 5)
    sleep(DELAY)
    
    genus = organism.split(' ')[0]
    species = organism.split(' ')[1]
    useful_id = np.nan
    species_found = []
    accesions_found=[]
    for jj, data in enumerate( response.text.split('>') ):
        
        # For every result, keep name of entry found
        if len(data)>0:
            species_found.append( ' '.join( data.split(' ')[1:-1] ) )
            accesions_found.append( data.split(' ')[0] )
            
        # Check if any entry fits perfectly the description            
        if (GENE_OF_INTEREST.lower() in data.lower()) and (genus.lower() in data.lower()) and (species.lower() in data.lower()):
            useful_id = jj
            accession_number = data.split(' ')[0]
            sequence_length  = 70*(len(data.split('\n'))-3)
    
    if np.isnan(useful_id):
        print_to_log('<W> Could not find any 16s for %s' % organism )
        print_to_log('Found instead for:')
        _=[ print_to_log('  -%s %s' % (a,b)) for a,b in zip(accesions_found, species_found) ]
    else:
        fasta = '>'+response.text.split('>')[useful_id]
        with open(FILENAME_OUT, 'a+') as file:
            file.write(fasta)
        print_to_log( f'{accession_number}\t{organism}\t{sequence_length}bp' )
    

# ##########
# Add the character > to all lines that have it missing
# with open('potf16s.fasta', 'r') as file:
#     with open('pot16s_good.fasta','w+') as file_out:
#         for line in file.readlines():
#             if ('16S' in line) and (line[0] != '>'):
#                 file_out.write( '>'+line )
#             else:
#                 file_out.write(line)
                


# ##########
# Rename tree labels


# def search_name(accession):
#     with open('pot16s.fasta','r') as file:
#         for line in file.readlines():    
#             if accession in line:
#                 return ' '.join(line.split(' ')[1:3]) +' '+accession
            
            
# with open('tree_16s.newick','r') as file:
#     line = file.readline()
# print('Line has length ', len(line) )
# leafs = line.split(':')

# for leaf in leafs:
#     if ',' in leaf:
#         old_name = leaf.split(',')[-1].replace('(','').replace(')','')
#         new_name = search_name(old_name)
#         # print(old_name, new_name)

#         # if new_name is not None:
#         line.replace( old_name, new_name)        
#         #     # print(old_name, new_name)

#         # else:
#         #     print('>W< Could not find name of %s' % old_name )
        
        
