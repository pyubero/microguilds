# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 10:55:37 2022

@author: Pablo
"""


import numpy as np
from matplotlib import pyplot as plt

# FILENAME = './data_definitive/seq_16S_ribosomal_rna.fasta'
FILENAME = './data_definitive/seq_recA_std.fasta'
# FILENAME = 'seq_recA.fasta'

FILENAME_OUT = FILENAME.split('.fasta')[0]+'_std.fasta'


headers = []
names   = []
sequences=[]

_nlines = 0
_sequence = ''
with open(FILENAME,'r') as file:
    for ii, line in enumerate(file.readlines()):
        _nlines+=1

        IS_HEADER = line[0]=='>'
        IS_EMPTY  = line[0]==''
        
        if IS_HEADER:
            headers.append(line.replace('\n','') )
            names.append( line[1:].split(' ')[0] )
            if len(_sequence)>1:
                sequences.append(_sequence)
            _sequence = ''
            
        if (not IS_HEADER) and (not IS_EMPTY):
            # print('added trozo')
            _sequence += line.replace('\n','')
sequences.append(_sequence)            

headers = np.array(headers)
names= np.array(names)
sequences = np.array(sequences)

print('Number of headers read:   %d' % len(headers))
print('Number of sequences read: %d' % len(sequences) )
print('')




new_headers = []

for header in headers:
    _accession = None
    _orgname = None
    _product = None
    
    
    if header[:4]=='>lcl':
        _accession = header.split(' ')[0][5:]
        _orgname   = header.split('orgname=')[-1].split(']')[0]
        if   'product' in header: _product = header.split('product=')[-1].split(']')[0]
        elif 'protein' in header: _product = header.split('protein=')[-1].split(']')[0]
        else:                     _product = 'unknown product'
        #...
        new_header = '>'+' '.join( [_accession, _orgname, _product])


    else:
        new_header = header
        
    new_headers.append( new_header)

new_headers = np.array( new_headers )


with open(FILENAME_OUT,'w+') as file:
    for jj in range(len(headers)):
        file.write(new_headers[jj]+'\n')
        file.write(sequences[jj]+'\n')






