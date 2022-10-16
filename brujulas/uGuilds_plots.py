# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""


import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

# Function to convert hexadecimal colors into RGB triplets
def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


# Function to calculate Roman values
def intToRoman(num):
  
    # Storing roman values of digits from 0-9
    # when placed at different places
    m = ["", "M", "MM", "MMM"]
    c = ["", "C", "CC", "CCC", "CD", "D",
         "DC", "DCC", "DCCC", "CM "]
    x = ["", "X", "XX", "XXX", "XL", "L",
         "LX", "LXX", "LXXX", "XC"]
    i = ["", "I", "II", "III", "IV", "V",
         "VI", "VII", "VIII", "IX"]
  
    # Converting to roman
    thousands = m[num // 1000]
    hundreds = c[(num % 1000) // 100]
    tens = x[(num % 100) // 10]
    ones = i[num % 10]
  
    ans = (thousands + hundreds +
           tens + ones)
  
    return ans
    
# Function to cycle through an array in steps of equal length
# ... only used to cycle through color palettes.
def cycle_through( array, length ):
    
    N = len(array)
    cycle = np.ceil( N/length )
    
    
    x0=0
    k = 1
    new_order = [x0,]
    while len(new_order)< len(array):
        new_idx = x0 + k*cycle
        
        if new_idx<len(array):
            new_order.append( int( new_idx)     )
            #print(new_idx)
    
            k+=1
        else:
            x0 +=1
            k=0
            
    return np.array( new_order )





GENE_NAME   = 'hzsA'
MIN_PERCENTILE =1
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
            #print('Starting with %s' % taxons[-1])
        
        line = f.readline()

        while line[:2] !='>>':            
            _K.append( [ float(_) for _ in  line.replace('\n','').split(',') ] )
            line = f.readline()
            
            if line=='':
                break
            
        Kmat.append( _K )
       
        
Kmat = np.array(Kmat)
n_taxons, n_ctxts, n_clusters = Kmat.shape




# Find taxons that contribute less than X percent and aggreggate their values
# into an "others" category.
# Others will be at the *first* position in Kmat and taxon names
TAXONS_SHOWN = 35


contrib_per_taxon = np.sum(np.sum(Kmat, axis=2), axis=1)
_taxons_shown = [TAXONS_SHOWN if TAXONS_SHOWN<=n_taxons-2 else n_taxons][0]
threshold = np.sort(contrib_per_taxon)[-_taxons_shown]

not_contrib   = np.argwhere( contrib_per_taxon<threshold)[:,0]
n_not_contrib = len( not_contrib )

# ...the variable OTHERS keeps track of whether the last position is of the others class
OTHERS = n_not_contrib>0

if OTHERS:
    print(f'Found {n_not_contrib} of {n_taxons} taxons that will be aggreggated into "Others"')
    print('The complete list is:')
    _ = [ print('  %s' % taxons[idx]) for idx in not_contrib]


    Kmat2 = np.zeros( (n_taxons-n_not_contrib+1, n_ctxts, n_clusters ) )
    tx2   = ['Others',]
    h = 1 # accumulator to keep track of next valid taxon

    for i_taxon in range(n_taxons):
        if i_taxon in not_contrib:
            Kmat2[0,:,:] += Kmat[i_taxon, :,:]
        else:
            Kmat2[h,:,:] = Kmat[i_taxon, :, :]
            tx2.append( taxons[i_taxon] )
            h += 1            
    
else:
    Kmat2  = Kmat.copy()
    tx2    = taxons.copy()
    


# #################################
# ### Get "official" color list ###
# #... load colors
# df = pd.read_csv('colors_malaspina_usable.txt', sep='\t', encoding='ansi')
# #... correct taxon names
# taxons=np.array(taxons)
# taxons[np.nonzero(taxons=='c__Actinobacteria')]='p__Actinobacteriota'             
# #... assign colors
# col_w_others=[]
# for i_taxon, txn in enumerate(tx_w_others):
#     subdf = df[df["Class_GTDB"] == txn]["Color"]
#     if len(subdf)>=1:
#         color = hex_to_rgb( subdf.values[0].replace(' ','') )
#         col_w_others.append( [color[0]/255, color[1]/255, color[2]/255, 1] )
#     else:
#         print('Could not find %s in color table' % txn)
#         col_w_others.append( [1,0,0,1] )
#     # col_w_others= [ [0.5,0.5,0.5,1], ]



THRESHOLD = 0.01

#########################
######### Plots #########
y = tx2
X = np.log10( Kmat2 )
X[X<THRESHOLD] = 0

#... derived values and matrices
_ntx, _nct, _ncl = X.shape
sumX  = np.sum(X, axis=0)
maxX  = np.max( np.max( X, axis=2), axis=1)

#... general plotting properties
R_UPPER_MARGIN_REL = 1.02
BAR_WIDTH_REL = 0.3

theta = np.linspace(0.0, 2 * np.pi, _ncl, endpoint=False)
rlims = np.linspace( X.min() , R_UPPER_MARGIN_REL*sumX.max(), 4)
Dtheta = theta[1]-theta[0]
width = BAR_WIDTH_REL*Dtheta


#... display order
#I = np.argsort( maxX )[-1::-1] # to plot first those that contribute the most...
I = range( _ntx)


#... style properties and labels 
labels = [ 'F. '+intToRoman( _+1) for _ in range( _ncl ) ]
ctxt_labels= ['Epipelagic','Mesopelagic','Bathypelagic']
colors = plt.cm.tab20( np.linspace(0, 1, 20) )
mpl.rcParams['hatch.linewidth'] = 0.3  # previous pdf hatch linewidth
hatches=[ '' , '////////' ,'........']


plt.figure( figsize=(18,10), dpi=200)
for i_ctxt in range( _nct ):
    ax = plt.subplot( 1, _nct+1, i_ctxt+1, projection='polar')    
    # No mola, porque others sale on top of all, pintar others el primero pero con un zorder=9999 para que en la leyenda salga al final??

    _bottoms = np.zeros( (_ncl,) )    
    
    # First of all print the "others" data
    if OTHERS:
        ax.bar( theta, X[ 0, i_ctxt ,:], 
               width = width,
               color = [0,0,0],
               label = y[0])
        _bottoms += X[0, i_ctxt, : ]
        
        
    # For every taxon in the ordered idx list I...
    for jj, i_tx in enumerate(I):
        
        # Skip the "others" category"
        if OTHERS and i_tx==0:
            continue
        
        #
        ax.bar( theta, X[ i_tx, i_ctxt ,:], 
               width = width,
               color = colors[ jj % 20] ,
               hatch = hatches[ int(np.floor(jj/20)) ] ,
               bottom = _bottoms,
               label = y[i_tx])
        _bottoms += X[i_tx, i_ctxt, : ]



    plt.title( ctxt_labels[i_ctxt])
    ax.set_rgrids( np.ceil(rlims[1:]), labels='' )
    ax.set_thetagrids( theta*57.3, labels=labels, fontsize=8)
    ax.xaxis.grid(linewidth=0.1)
    ax.yaxis.grid(linewidth=0.2)


legend = plt.legend( fontsize=8)

ax = plt.subplot(1, _nct+1, _nct+1)
plt.axis('off')

legend.set_bbox_to_anchor( (2.1, 1.25) )

plt.savefig(f"gpattern_{GENE_NAME}.png")





