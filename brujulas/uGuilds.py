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
    abundances = subtable["TPM"]
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


def export_legacy(df, filename):
    gene = df["Gene"].iloc[0]
    taxons = df["Taxon"].unique()
    contexts = df["Context"].unique()
    clusters = df['Cluster'].unique()
    new_contexts = np.array(["Epipelagic","Mesopelagic","Bathypelagic"])
    
    n_taxons= len(taxons)
    n_contexts = len(contexts)
    n_clusters = len(clusters)

    # Build Kmat
    Kmat = np.zeros((n_taxons, n_contexts, n_clusters))
    idc = np.array(np.meshgrid( range(n_taxons), range(n_contexts), range(n_clusters))).T.reshape(-1,3)
    for j_tx, j_ct, j_cl in verbosebar(idc):
            idx = (adu_table["Taxon"]==taxons[j_tx]) & \
                    (adu_table["Context"]==new_contexts.astype("str")[j_ct]) & \
                    (adu_table["Cluster"]==clusters[j_cl])
            assert sum(idx)==1
            Kmat[j_tx, j_ct, j_cl] = adu_table[idx]["k-value"]
    
    # Prepare output file
    with open(filename, 'w+') as f: 
        f.write(f'{gene}\n')
        f.write( ','.join( [ "%s" % _ for _ in clusters ] )+'\n' )
        
    # For every taxon...
    for j_taxon in verbosebar(range( n_taxons)):
        
        #... create a K matrix
        K = Kmat[j_taxon, :,:] 
        
        # For every clusters...
        for j_cluster in range(n_clusters):
            # Append to output file...
            header = f'>>{taxons[j_taxon]}\n'
            lines = []
            for _ in range( K.shape[0]) :
                lines.append( ','.join( [ '%1.4f' % k for k in K[_,:]])+'\n' )
                        
        with open(filename, 'a') as f: 
            f.write(header)
            _ = [ f.write(line) for line in lines ]


def bivariate_regression(x,y):
    '''
    Computes the bivariate regression as the angle of the covariance matrix.
    '''    
    
    x0 = np.mean(logx)
    y0 = np.mean(logy)
    C = np.cov([logx,logy])
    V = np.linalg.eig(C)[1]
    alfa = V[1][0]/V[0][0]
    beta = y0-alfa*x0
    
    return alfa, beta


FILENAME    = 'mastertable_w_ctxt.tsv'
GENE_NAME   = 'potF'
LEVEL_NAME  = 'Species_GTDB'
VERBOSE = True
EXPORT_PLOT = True
EXPORT_LEGACY = True
# ...
out_filename = f'kvalues_{GENE_NAME}_{LEVEL_NAME}.tsv'
out_plot = f"loglog_regression_{GENE_NAME}_{LEVEL_NAME}.png"


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
verboseprint("")

# Process data
idc = np.array(np.meshgrid( range(n_taxons), range(n_contexts), range(n_clusters))).T.reshape(-1,3)
included = gene_table["taxonomic_classification_level"] == np.random.rand() # to keep track of included and missing samples

adu_table = pd.DataFrame(columns=["Gene", "Taxon", "Context", "Cluster",
                                  "Abundance", "Diversity", "Univocity"])

# Compute Abundance, Diversity and Univocity
for j_tx, j_ct, j_cl in verbosebar(idc):
    # Compute K value
    _ab, _dv, _un = compute_adu(gene_table, [taxons, contexts, clusters], [j_tx, j_ct, j_cl] )
    new_row = pd.Series(
        {"Gene" : GENE_NAME,
         "Taxon":taxons[j_tx],
         "Context" : contexts[j_ct],
         "Cluster" : clusters[j_cl],
         "Abundance" : _ab,
         "Diversity" : _dv,
         "Univocity" : _un} )
    
    adu_table = pd.concat([adu_table, new_row.to_frame().T], ignore_index=True)
    
    
    # To keep track of included and missing samples    
    _id = ( gene_table["taxonomic_classification_level"] == taxons[j_tx] ) & \
           ( gene_table["Context"] ==  contexts[j_ct]) & \
           ( gene_table["cluster_id"] ==  clusters[j_cl] )
    included = included | _id

verboseprint(f"\nIncluded {sum(included)} out of {len(gene_table)} rows.")    
verboseprint(f"The sum of diversities is {sum(adu_table['Diversity'])}.")
verboseprint("")


#############
# Regression 
# ... maybe this should be refactored aka put into functions
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from matplotlib import pyplot as plt


def linear_function(x, a, b):
    return a*x + b

# Obtain data to regress    
idx = adu_table["Abundance"]>0
x = adu_table["Abundance"][idx].to_numpy().astype("float")
y = adu_table["Diversity"][idx].to_numpy().astype("float")

# Transform to loglog
ABUNDANCE_THOLD = 1e-5
logx = np.log10(ABUNDANCE_THOLD + x)
logy = np.log10(ABUNDANCE_THOLD + y)

# Linear regression
gamma, c = bivariate_regression(x,y)
r2 = r2_score(logy, logx*gamma+c)

# Print results
verboseprint("Bivariate loglog regression results:")
verboseprint(f"gamma = {gamma}")
verboseprint(f"c = {c}")
verboseprint(f"R2 = {r2}")
verboseprint("")

# Compute corrections
delta = 10**logy / np.clip( 10**linear_function(logx, gamma, c), 1, np.inf)
_delta = np.zeros(len(adu_table))
_delta[idx] = delta

# Export data
adu_table["delta"] = _delta
adu_table["k-value"] = adu_table["Abundance"]*_delta 
adu_table.to_csv(out_filename, sep="\t", index=False)
verboseprint(f"Data saved in {out_filename}.")

if EXPORT_PLOT:
    H = plt.figure( figsize=(12,4), dpi=300)
    plt.subplot(1,2,1)
    plt.scatter(logx, logy,s=10, c=logy)
    plt.plot(logx, linear_function(logx, gamma, c))
    plt.grid()
    plt.xlabel("log10 Abundance")
    plt.ylabel("log10 Diversity")
    plt.legend(["Data",f"Biv. loglog reg\n$\gamma$={gamma:.3f}\nR2={r2:.3f}"])
    
    plt.subplot(1,2,2)
    plt.scatter(x, x*delta,s=delta*10, c=logy)
    plt.grid()
    plt.colorbar(label="log10 Observed diversity")
    plt.plot(x, x, 'k', lw=0.5)
    plt.xlabel('sum(a)')
    plt.ylabel("sum(a)*delta")
    plt.xscale("log")
    plt.yscale("log")
    
    plt.suptitle(f"{GENE_NAME}, {LEVEL_NAME}")
    plt.savefig(out_plot)

if EXPORT_LEGACY:
    export_legacy(adu_table, f"legacy_kMatrixPerTaxon_{GENE_NAME}_{LEVEL_NAME}.csv")




# you should call before:


def from_df_to_ktensor(df, data, column="k-value"):
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

table = pd.read_csv(f"kvalues_{GENE_NAME}_{LEVEL_NAME}.tsv", sep="\t")

taxons = table['Taxon'].unique()
contexts = np.array(["Epipelagic","Mesopelagic","Bathypelagic"])
clusters = table["Cluster"].unique()

K = from_df_to_ktensor(table, [taxons, contexts, clusters], "k-value")
