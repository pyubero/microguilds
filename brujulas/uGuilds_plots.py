# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""


import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl


def hex_to_rgb(value):
    # Function to convert hexadecimal colors into RGB triplets
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


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
    

def cycle_through( array, length ):
    # Function to cycle through an array in steps of equal length
    # ... only used to cycle through color palettes.

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

# How to nicely load our tensor?
def from_df_to_ktensor(df, column):
    # you should call before:
    # df = pd.read_csv(f"kValuesPerTaxon_{GENE_NAME}.tsv", sep="\t")
    ntaxons = len(df["Taxon"].unique())
    ncontexts = len(df["Context"].unique())
    nclusters = len(df["Cluster"].unique())
    K = df[column].to_numpy().reshape(nclusters, ntaxons, ncontexts)

    K = np.moveaxis(K, 0,2)
    return K


def score_taxon(taxon):
    _tx = taxon.split(' ')
    if (len(_tx)>1) and (_tx[1][:2]!="sp"):
        has_species = True
    else:
        has_species = False
    
    if (len(_tx)>0) and (len(_tx[0])>3) and (_tx[0][1:3]=="__") and (_tx[0][3:6]!="GCA") and (_tx[0][3:6]!="UBA"):
        has_genus = True
    else:
        has_genus = False
        
    score = len(_tx)+0.5*has_species +0.5*has_genus -1
    return score


def stackbar(ax, x, Y, *args, **kwargs):
    if "bottom" not in kwargs:
        kwargs.update( {"bottom" : np.zeros(Y.shape[1],)} )
    
    if "labels" in kwargs:
        labels = kwargs["labels"]
        kwargs.pop("labels")
        _given_labels = True
    else:
        _given_labels = False
        
    if "colors" in kwargs:
        colors = kwargs["colors"]
        kwargs.pop("colors")
        _given_colors = True
    else:
        _given_colors = False
        
    if len(Y)==0:
        return ax, kwargs["bottom"]
        
    for idx in range(Y.shape[0]):
        if _given_labels:
            kwargs.update({"label" : labels[idx]})
        
        if _given_colors:
            kwargs.update( {"color" : colors[idx,:] })
        
        # Plot
        ax.bar( x, Y[idx,:], *args, **kwargs)
        
        #... update bottoms
        new_bottoms = kwargs["bottom"] + Y[idx,:]
        kwargs.update( {"bottom" : new_bottoms} )
        
        #... remove label
        if ("label" in kwargs) and ( type(kwargs["label"])==str):
            kwargs.pop("label")
        
    return ax, kwargs["bottom"]


def do_not_overwrite_filepath(filepath):
    import os
    new_filepath = filepath
    idx = 0
    while os.path.exists(new_filepath):
        idx += 1
        split = filepath.split('.')
        fp = '.'.join(split[:-1])
        new_filepath = f"{fp}({idx}).{split[-1]}"
        
    return new_filepath

def save_figure(HFigure, filepath, overwrite=True, *args, **kwargs):
    if not overwrite:
        filepath = do_not_overwrite_filepath(filepath)
                    
    HFigure.savefig(filepath, *args, **kwargs)



# General variables
GENE_NAME = 'nirs'
LEVEL_NAME = 'Species_GTDB'
OVERWRITE = False
_filename = f'kMatrixPerTaxon_{GENE_NAME}_{LEVEL_NAME}_v2.csv'
out_filename = f"gpattern_{GENE_NAME}_{LEVEL_NAME}.png"

#... plotting options
COLOR_FILE = None # if not None, will try to load a csv file like a dict with colors per taxon
MAX_TAXONS_SHOWN = 30
K_MIN = 1 # in tpm units
TAXONOMIC_ADVANTAGE = 0 # taxons fully assigned advantage to be displayed over worse identifiable taxons

#... some styling options
DPI = 200
R_UPPER_MARGIN_REL = 1.02 # radial margin top, relative to the maximum
BAR_WIDTH_REL = 0.5
COLOR_BELOW_MIN = [0.0, 0.0, 0.0]
COLOR_UNASSIGNED= [0.3, 0.3, 0.3]
COLOR_OTHERS = [0.5, 0.5, 0.5]




df = pd.read_csv(_filename, sep='\t', header=0)
# ...
K = from_df_to_ktensor(df, column="k-value")
# ...
diversity = from_df_to_ktensor(df, column="Diversity")
taxons = from_df_to_ktensor(df, column="Taxon")[:,0,0]
contexts = from_df_to_ktensor(df, column="Context")[0,:,0]
clusters = from_df_to_ktensor(df, column="Cluster")[0,0,:]
# ...
n_taxons, n_ctxts, n_clusters = K.shape


## Find plottable DATA ##
contrib_per_taxon = np.sum(np.sum(K, axis=2), axis=1)
idx_below_min = np.argwhere(contrib_per_taxon<K_MIN)[:,0] 


# All UNASSIGNED taxons add to the class "Unassigned"
# ... we compute the taxonomic scores
tx_scores = np.array([score_taxon(tx) for tx in taxons])
idx_unassigned = np.argwhere(tx_scores==0)[:,0]


# All taxons in the top MAX_TAXONS_SHOWN add to the class "Show"
idc_left = list(set(range(n_taxons)) - set(idx_below_min) - set(idx_unassigned))
_maxshown = np.clip(MAX_TAXONS_SHOWN,0, len(idc_left))
scores = np.log10( contrib_per_taxon ) + tx_scores*TAXONOMIC_ADVANTAGE
threshold = np.sort(scores[ idc_left ])[-_maxshown]
idx_show = np.argwhere( scores[ idc_left ]>=threshold)[:,0]
idx_show = np.array( idc_left )[idx_show]


# All taxons left add to the class "Others"
idx_other = list( set(idc_left) - set(idx_show) )

  


# General plotting properties
#... style properties and labels 
clst_labels = [ 'F. '+intToRoman( _+1) for _ in range( n_clusters ) ]
# clst_labels = clusters
ctxt_labels = ['Epipelagic','Mesopelagic','Bathypelagic']
colors = plt.cm.tab20( np.linspace(0, 1, 20) )
mpl.rcParams['hatch.linewidth'] = 0.3  # previous pdf hatch linewidth
hatches=[ '' , '////////' ,'........']



#########################
######### Plots #########
# Linear version
X = K.copy()
# Log version
# X = np.log10( K )
# X[X<1] = 0



#... derived values and matrices
sumX  = np.sum(X, axis=0)    # total K per cluster per context
maxX  = np.max( np.max( X, axis=2), axis=1) # maximum contribution per taxon

theta = np.linspace(0.0, 2 * np.pi, n_clusters, endpoint=False)
rlims = np.linspace( X.min() , R_UPPER_MARGIN_REL*sumX.max(), 4)
Dtheta = theta[1]-theta[0]
width = BAR_WIDTH_REL*Dtheta


n_cols = int(np.clip(n_ctxts+1, 1, 4))
n_rows = int(np.ceil((n_ctxts+1)/ 4))
sz_cols = 4.5*n_cols
sz_rows = 6*n_rows

HFigure = plt.figure( figsize=(sz_cols, sz_rows), dpi=DPI)
bottoms = np.zeros((n_ctxts, n_clusters))

for _context in range( n_ctxts ):
    ax = plt.subplot( n_rows, n_cols, _context+1, projection='polar')    
    
    # Plot BELOW_MIN        
    ax, _bot = stackbar(ax, theta, X[idx_below_min, _context, :],
                       width = width,
                       color = COLOR_BELOW_MIN,
                       label = "Below minimum",
                       bottom = bottoms[_context, :]
                       )
    
        
    # Plot UNASSIGNED
    bottoms[_context, :] = _bot
    ax, _bot = stackbar(ax, theta, X[idx_unassigned, _context, :],
                       width = width,
                       color = COLOR_UNASSIGNED,
                       label = "Unassigned",
                       bottom = bottoms[_context, :]
                       )
    
    
    # Plot OTHERS
    bottoms[_context, :] = _bot
    ax, _bot = stackbar(ax, theta, X[idx_other, _context, :],
                       width = width,
                       color = COLOR_OTHERS,
                       label = "Others",
                       bottom = bottoms[_context, :]
                       )
    
    
    # Plot SHOWN
    bottoms[_context, :] = _bot
    ax, _bot = stackbar(ax, theta, X[idx_show, _context, :],
                       width = width,
                       colors = colors,
                       labels = list(taxons[idx_show]),
                       bottom = bottoms[_context, :]
                       )
    

    ax.set_rgrids( np.ceil(rlims[1:]), labels='' )
    ax.set_thetagrids( theta*57.3, labels=clst_labels, fontsize=8)
    ax.xaxis.grid(linewidth=0.1)
    ax.yaxis.grid(linewidth=0.2)
    plt.title( ctxt_labels[_context])

legend = plt.legend( fontsize=8)

#... place legend
ax = plt.subplot(1, n_ctxts+1, n_ctxts+1)
plt.axis('off')
legend.set_bbox_to_anchor( (2.1, 1.25) )


# Save figure
save_figure(HFigure, out_filename, overwrite=OVERWRITE, dpi=DPI )




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



