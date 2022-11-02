# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 18:46:06 2022

@author: Pablo
"""

import numpy as np
from matplotlib import pyplot as plt

# FILENAME = './data_definitive/seq_16S_ribosomal_rna.fasta'
FILENAME = './data_definitive/seq_recA_std.fasta'
# FILENAME = '16s.fasta'


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


# Print names that repeat
unq_names = np.unique(names)
if len(unq_names)>len(names):
    print('---- Names that are repeated ----')    
    for name in unq_names:
        ntimes = np.argwhere((names==name))[:,0]
        if len(ntimes)>1:
            print("%s" % name)
            for idx in ntimes:
                print("    "+headers[idx])
            print('')
else:
    print('All sequences names are unique.')  
    print('')      


# Search for sequences specially longer/shorter
seqlens = np.array([ len(_) for _ in sequences])
print('---- Sequences lengths ----')
print('Min  length: %d' % np.min(seqlens))
print('Mean length: %1.1f +/- %1.1f' % (np.mean(seqlens), np.std(seqlens)))
print('Max  length: %d' % np.max(seqlens))

print('Sequences with l<1200:')
idc = np.argwhere(seqlens<800)[:,0]
if len(idc)>0:
    _=[print("    %4d\t%s" % (seqlens[ii], names[ii] ) ) for ii in idc]




plt.hist(seqlens)



