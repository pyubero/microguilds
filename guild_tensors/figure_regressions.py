# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from matplotlib import colormaps as cm
from matplotlib import pyplot as plt
import guild_tensor_utils as gtutils
from guild_tensor_utils import verbosebar, verboseprint


FILENAME = 'mastertable.tsv'
GENE_NAMES = ["amoA", "hzsA", "nrt", "potF"]
LEVEL_NAME = 'Species_GTDB'
VERBOSE = True
EXPORT_PLOT = True
EXPORT_LEGACY = False

FIGSIZE = (12, 4)
DPI = 300
COLORS = cm.get_cmap("tab10")(np.linspace(0, 256, 200, dtype='int'))
FILENAME_OUT = "Figure_regression.png"

# Import mastertable
master_table = gtutils.check_mastertable(FILENAME, True)
verboseprint(f"Loaded mastertable with {len(master_table)} rows.")

# Validate master table

assert LEVEL_NAME in master_table.columns
assert "gene_fun" in master_table.columns
assert "cluster_id" in master_table.columns

# Standardize taxonomic column name
master_table = master_table.rename(
    columns={LEVEL_NAME: "taxonomic_classification_level"}
    )

# Create figure
H = plt.figure(figsize=FIGSIZE, dpi=DPI)


for gg, gene in enumerate(GENE_NAMES):

    # Filter by gene name
    gene_table = master_table[master_table["gene_fun"] == gene]
    verboseprint(f"Subtable for gene *{gene}* has {len(gene_table)} rows.")

    # Find clusters in gene subtable
    clusters = gene_table['cluster_id'].unique()
    n_clusters = len(clusters)

    # Find taxons in gene subtable
    taxons = gene_table["taxonomic_classification_level"].unique()
    n_taxons = len(taxons)

    # Find contexts in gene subtable
    contexts = gene_table["Context"].unique()
    n_contexts = len(contexts)

    # Process data
    idc = np.array(
        np.meshgrid(range(n_taxons), range(n_contexts), range(n_clusters))
        ).T.reshape(-1, 3)

    # Define accumulator
    adu_table = pd.DataFrame(columns=["Gene", "Taxon", "Context", "Cluster",
                                      "Abundance", "Diversity", "Univocity"])
    data = [taxons, contexts, clusters]

    for j_tx, j_ct, j_cl in verbosebar(idc):
        # Compute K value
        _ab, _dv, _un = gtutils.compute_adu(gene_table,
                                            data,
                                            [j_tx, j_ct, j_cl])
        new_row = pd.Series(
            {"Gene": gene,
             "Taxon": taxons[j_tx],
             "Context": contexts[j_ct],
             "Cluster": clusters[j_cl],
             "Abundance": _ab,
             "Diversity": _dv,
             "Univocity": _un})

        adu_table = pd.concat([adu_table, new_row.to_frame().T],
                              ignore_index=True)

    #############
    # Regression
    # ... maybe this should be refactored aka put into functions
    def linear_function(x_values, slope, offset):
        '''Returns y = slope * x + offset'''
        return slope * x_values + offset

    # Load data to regress
    idx = adu_table["Abundance"] > 0
    x = adu_table["Abundance"][idx].to_numpy().astype("float")
    y = adu_table["Diversity"][idx].to_numpy().astype("float")

    # Transform data to loglog
    ABUNDANCE_THOLD = 1e-10
    logx = np.log10(ABUNDANCE_THOLD + x)
    logy = np.log10(ABUNDANCE_THOLD + y)
    logxx = np.sort(logx)
    # Linear regression
    gamma, c = gtutils.bivariate_regression(logx, logy)
    r2 = r2_score(logy, logx * gamma + c)
    print(f"Gene: {gene} with R2={r2:0.2f}")

    label = fr"{gene}, $\gamma$={gamma:.2f}"

    ax = plt.subplot(1, 2, 1)
    ax.plot(logx, logy, '.', ms=3, color=COLORS[gg], label=label)
    ax.plot(logxx, linear_function(logxx, gamma, c), color=COLORS[gg])
    ax.grid(visible=True)
    ax.set_xlabel("log$_{10}$ Abundance")
    ax.set_ylabel("log$_{10}$ Diversity")

ax.legend(framealpha=1.0)

ax = plt.subplot(1, 2, 2)
delta = 10**logy / np.clip(10**linear_function(logx, gamma, c), 1, np.inf)
logK = np.log10(ABUNDANCE_THOLD + x*delta)

ax.plot(logx, logK, '.', ms=3, color=COLORS[gg], label=f"{gene}")
ax.plot(logx, logx, 'k', lw=0.5)
ax.grid(visible=True)
ax.set_xlabel("log$_{10}$ Abundance")
ax.set_ylabel("log$_{10}$ K")

plt.savefig(FILENAME_OUT, bbox_inches='tight')
