# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from tqdm import tqdm


def verboseprint(msg):
    '''Prinst a message if VERBOSE is True.'''
    if VERBOSE:
        print(msg)
      
        
def verbosebar(iterable):
    '''Generates a progress bar with tqdm if VERBOSE is True.'''
    if VERBOSE:
        return tqdm(iterable)
    else:
        return iterable
    
    
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
    
            k+=1
        else:
            x0 +=1
            k=0
            
    return np.array( new_order )


# How to nicely load our tensor?
def from_df_to_ktensor(df, data, column="k-value"):
    # you should call before:
    # df = pd.read_csv(f"kValuesPerTaxon_{GENE_NAME}.tsv", sep="\t")
    taxons, contexts, clusters = data
    ntaxons = len(taxons)
    ncontexts = len(contexts)
    nclusters = len(clusters)    
    
    Kmat = np.zeros((ntaxons, ncontexts, nclusters), dtype="object")
    idc = np.array(np.meshgrid( range(ntaxons), range(ncontexts), range(nclusters))).T.reshape(-1,3)
    for j_tx, j_ct, j_cl in verbosebar(idc):
            idx = (df["Taxon"]==taxons[j_tx]) & \
                    (df["Context"]==contexts.astype("str")[j_ct]) & \
                    (df["Cluster"]==clusters[j_cl])
            assert sum(idx)==1
            Kmat[j_tx, j_ct, j_cl] = df[idx][column].iloc[0]
         
    return Kmat

def score_taxon(taxon : str):
    '''
    Scores taxon name. 
    It returns: 
        0 if it does not have a defined species nor genus.
        1 if it only has a genus specified.
        2 if it has both, genus and species specified.
    '''
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
        
    if "hatches" in kwargs:
        hatches = kwargs["hatches"]
        kwargs.pop("hatches")
        _given_hatches = True
    else:
        _given_hatches = False
        
    if Y.shape[0] == 0:
        return ax, kwargs["bottom"]
        
    for idx in range(Y.shape[0]):
        if _given_labels:
            kwargs.update({"label" : labels[idx]})
        
        if _given_colors:
            idx_color = idx % colors.shape[0]
            kwargs.update( {"color" : colors[idx_color,:] })
            if _given_hatches:
                idx_hatch = idx // colors.shape[0]
                kwargs.update( {"hatch" : hatches[idx_hatch]})
        
        
        # Plot
        ax.bar( x, Y[idx,:], *args, **kwargs)
        
        #... update bottoms
        new_bottoms = kwargs["bottom"] + Y[idx,:]
        kwargs.update( {"bottom" : new_bottoms} )
        
        #... remove label
        if ("label" in kwargs) and ( type(kwargs["label"])==str):
            kwargs.pop("label")
        
    return ax, kwargs["bottom"]


def kmax_taxons(K, n):
    top_taxon_idc = []
    for ii in np.argsort(K.flatten())[-1::-1]:
        taxon_idx = np.argwhere(K==K.flatten()[ii])[0][0]
        if taxon_idx not in top_taxon_idc:
            top_taxon_idc.append(taxon_idx)
        if len(top_taxon_idc)>=n:
            break
    return top_taxon_idc


def kmin_taxons(K, n):
    top_taxon_idc = []
    values = K.flatten()
    for ii in np.argsort(values):
        if values[ii]==0:
            continue
        taxon_idx = np.argwhere(K==K.flatten()[ii])[0][0]
        if taxon_idx not in top_taxon_idc:
            top_taxon_idc.append(taxon_idx)
        if len(top_taxon_idc)>=n:
            break
    return top_taxon_idc


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
GENE_NAME = 'potF'
LEVEL_NAME = 'Species_GTDB'
DISPLAY_MODE = 'log10' #  ... either linear or log10
DISPLAY_KIND = "common" #  ... either common or rare
CONTRIBUTION = "summed" #  ... either single or summed
FIGSIZE = [16, 6]
OVERWRITE = True
VERBOSE = True
#...
_filename = f'kvalues_{GENE_NAME}_{LEVEL_NAME}.tsv'
out_filename = f"gpattern_{GENE_NAME}_{LEVEL_NAME}.png"

# Plotting options
MAX_TAXONS_SHOWN = 25
TAXONOMIC_ADVANTAGE = 0 # taxons fully assigned advantage over worse identifiable taxons

# Styling options
DPI = 300
N_COLORS = 20 #  <=20
R_UPPER_MARGIN_REL = 1.02 #  Radial margin top, relative to the maximum
BAR_WIDTH_REL = 0.5
COLOR_UNASSIGNED= [0.0, 0.0, 0.0, 1]
COLOR_OTHERS = [0.3, 0.3, 0.3, 1]


#######################################################
###################### LOAD DATA ######################
df = pd.read_csv(_filename, sep='\t', header=0)
# ...
taxons = df['Taxon'].unique()
contexts = np.array(["Epipelagic","Mesopelagic","Bathypelagic"])
clusters = df['Cluster'].unique()
# ...
K = from_df_to_ktensor(df, [taxons, contexts, clusters], column="k-value").astype("float")
# ...
n_taxons, n_ctxts, n_clusters = K.shape

verboseprint(f"Loaded data for:")
verboseprint(f"   {n_taxons:4d}\ttaxons")
verboseprint(f"   {n_ctxts:4d}\tcontexts")
verboseprint(f"   {n_clusters:4d}\tclusters.")



#######################################################
#################### ANALYZE DATA #####################
## Find plottable DATA ##
# contrib_per_taxon = np.sum(np.sum(K, axis=2), axis=1)

# All UNASSIGNED taxons add to the class "Unassigned"
# ... we compute the taxonomic scores
idx_unassigned = np.argwhere(np.array(taxons) == "s__")[:, 0]
taxons[idx_unassigned] = "Unassigned"


# All taxons in the top MAX_TAXONS_SHOWN add to the class "Show"
_maxshown = np.clip(MAX_TAXONS_SHOWN, 0, n_taxons)

if CONTRIBUTION == "single":
    if DISPLAY_KIND == "common":
        idx_show = kmax_taxons(K, _maxshown)   
        k_min = K[idx_show[-1],:,:].max()
        _symbol="\leq"
    elif DISPLAY_KIND == "rare":
        idx_show = kmin_taxons(K, _maxshown)   
        k_min =  min(i for i in K[idx_show[-1],:,:].flatten() if i > 0)
        _symbol = "\geq"
    else:
        assert False
        
elif CONTRIBUTION == "summed":
    if DISPLAY_KIND == "common":
        contrib_per_taxon = np.sum(np.sum(K, axis=2), axis=1)
        tx_scores = np.array([score_taxon(tx) for tx in taxons])
        scores = np.log10( contrib_per_taxon ) + tx_scores*TAXONOMIC_ADVANTAGE
        threshold = np.sort(scores)[-_maxshown]
        idx_show = np.argwhere( scores>=threshold)[:,0]
        k_min = contrib_per_taxon[idx_show].min()
        _symbol = "\leq"
    elif DISPLAY_KIND =="rare":
        print(">> ERROR << Still not implemented")
        
        
# All taxons left add to the class "Others"
idx_other = list( set(range(n_taxons)) - set(idx_show) )

# Build K matrix to be shown, by summing the Kvalues of others, and unassigned.  
taxon_labels = [f"Others, K${_symbol}${k_min:1.1f}", *taxons[idx_show]]
K_shown = np.zeros((len(idx_show)+1, n_ctxts, n_clusters) )
K_shown[0,:,:] = np.sum(K[idx_other,:,:], axis=0)
K_shown[1:,:,:]= K[idx_show,:,:]


verboseprint(f"Found {len(idx_unassigned)} unassigned taxons:")
_ = [verboseprint(f"  {tx}") for tx in taxons[idx_unassigned] ]
verboseprint("")
verboseprint(f"Will group {len(idx_other)} taxons into 'Others'.")
verboseprint(f"Showing the top {len(idx_show)} contributing taxons, with")
verboseprint(f"  K-values >= {k_min:1.3f}.")
verboseprint("")


# General plotting properties
#... style properties and labels 
# clst_labels = [ 'F. '+intToRoman( _+1) for _ in range( n_clusters ) ]
colors = plt.cm.tab20( np.linspace(0, 1, N_COLORS) )
colors = np.concatenate(([COLOR_OTHERS,], colors), axis=0)
mpl.rcParams['hatch.linewidth'] = 0.4  # previous pdf hatch linewidth
hatches=[ '' , '//////////' ,'..........']


#%%
#######################################################
###################### PLOT DATA ######################
if DISPLAY_MODE == "linear":
    X = K.copy()
    verboseprint("Displaying stacked k_values, in LINEAR mode.")
    
elif DISPLAY_MODE == "log10":
    X = np.log10( 1+K_shown )
    X[X<0] = 0
    verboseprint("Displaying stacked log10(k_values), in LOG10 mode.")


#... derived values and matrices
sumX  = np.sum(X, axis=0)    # total K per cluster per context
maxX  = np.max( np.max( X, axis=2), axis=1) # maximum contribution per taxon

theta = np.linspace(0.0, 2 * np.pi, n_clusters, endpoint=False)
# rlims = np.ceil(np.linspace(0, R_UPPER_MARGIN_REL*sumX.max(), 4))
dr = np.round(R_UPPER_MARGIN_REL*sumX.max()/4)
rlims = np.arange(0, 1+np.ceil(R_UPPER_MARGIN_REL*sumX.max()/dr)*dr, dr)
Dtheta = theta[1]-theta[0]
width = BAR_WIDTH_REL*Dtheta


n_cols = int(np.clip(n_ctxts+1, 1, 3))
n_rows = int(np.ceil((n_ctxts+1)/ 4))
FIGSIZE[1] = FIGSIZE[1] * n_rows
FIGSIZE = np.array(FIGSIZE)/2.3

HFigure = plt.figure( figsize=FIGSIZE, dpi=DPI)
bottoms = np.zeros((n_ctxts, n_clusters))

for _context in range( n_ctxts ):
    ax = plt.subplot( n_rows, n_cols, _context+1, projection='polar')    

    # Plot Others and Unassigned
    ax, _bot = stackbar(ax, theta, X[:1, _context,:],
                       width=width,
                       colors=colors[:1,:],
                       hatches=hatches,
                       labels=taxon_labels[:1])


    # Plot the rest, nicely with hatches and colors
    bottoms[_context, :] = _bot
    ax, _ = stackbar(ax, theta, X[1:, _context,:],
                       width=width,
                       colors=colors[1:],
                       hatches=hatches,
                       labels=taxon_labels[1:],
                       bottom = bottoms[_context, :])

    ax.set_rgrids( rlims, labels="" )
    ax.set_thetagrids( theta*57.3, labels=clusters, fontsize=7)
    ax.xaxis.grid(linewidth=0.1)
    ax.yaxis.grid(linewidth=0.2)
    plt.title( contexts[_context], {"fontsize" : 11}, pad=18)

# ... place legend    
plt.subplot(n_rows,n_cols, 2)
legend = plt.legend(fontsize=7, ncol=3,
                    loc="upper center", bbox_to_anchor=(0.5, -0.15))


verboseprint(f"The radial units are:")
verboseprint(f"{np.ceil(rlims)}")

# Save figure
save_figure(HFigure, out_filename, overwrite=OVERWRITE, dpi=DPI,
            bbox_inches='tight')




# #################################
# ### Get "official" color list ###
# COLOR_FILE = None # if not None, will try to load a csv file like a dict with colors per taxon
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



