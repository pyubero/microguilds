# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""
import numpy as np
from matplotlib import colormaps as cm
from matplotlib import pyplot as plt
import guild_tensor_utils as gtutils
from guild_tensor_utils import verboseprint


FILENAME = 'mastertable.tsv'
GENE_NAMES = ["amoA", "hzsA", "nirK", "potF"]
LEVEL_NAME = 'Species_GTDB'
VERBOSE = True
EXPORT_PLOT = True
EXPORT_LEGACY = False

FIGSIZE = (12, 4)
DPI = 300
COLORS = cm.get_cmap("tab10")(np.linspace(0, 1, 10))
FILENAME_OUT = "Figure_regression.png"

# Import mastertable
master_table = gtutils.check_mastertable(FILENAME, True)
verboseprint(f"Loaded mastertable with {len(master_table)} rows.", VERBOSE)

# Validate master table
assert LEVEL_NAME in master_table.columns
assert "gene_fun" in master_table.columns
assert "cluster_id" in master_table.columns


def linear_function(x_values, slope, offset):
    '''Returns y = slope * x + offset'''
    return slope * x_values + offset


# Create figure
H = plt.figure(figsize=FIGSIZE, dpi=DPI)
for gg, gene in enumerate(GENE_NAMES):

    # Compute adu_table
    adu_table = gtutils.build_adu_table(master_table,
                                        gene,
                                        LEVEL_NAME,
                                        force_build=False)

    # Linear regression
    gamma, c, r2 = gtutils.compute_gamma(adu_table, VERBOSE)
    print(f"Gene: {gene} with R2={r2:0.2f}")

    # Load data to regress
    idx = adu_table["Abundance"] > 0
    x = adu_table["Abundance"][idx].to_numpy().astype("float")
    y = adu_table["Diversity"][idx].to_numpy().astype("float")

    # Transform data to loglog
    LOG_THRESHOLD = 1e-10
    logx = np.log10(LOG_THRESHOLD + x)
    logy = np.log10(LOG_THRESHOLD + y)
    logxx = np.sort(logx)

    # Plot
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
logK = np.log10(LOG_THRESHOLD + x*delta)

ax.plot(logx, logK, '.', ms=3, color=COLORS[gg], label=f"{gene}")
ax.plot(logx, logx, 'k', lw=0.5)
ax.grid(visible=True)
ax.set_xlabel("log$_{10}$ Abundance")
ax.set_ylabel("log$_{10}$ K")

plt.savefig(FILENAME_OUT, bbox_inches='tight')
