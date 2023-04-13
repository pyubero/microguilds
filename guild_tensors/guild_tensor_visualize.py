# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import guild_tensor_utils as gtutils
from guild_tensor_utils import verboseprint


# General variables
GENE_NAME = 'potF'
LEVEL_NAME = 'Species_GTDB'
CONTEXTS = np.array(["Epipelagic", "Mesopelagic", "Bathypelagic"])
# CONTEXTS = np.array(["Meso_2237", "Meso_0538", "Others"])
UNASSIGNED = ["s__", "Unspecified"]
OVERWRITE = True
VERBOSE = True

# Plotting options
MAX_TAXONS_SHOWN = 23
DISPLAY_KIND = "common"  # either common or rare
DISPLAY_MODE = 'log10'  # either linear or log10
CONTRIBUTION = "single"  # either single or summed
PROJECTION = "polar"  # either polar, or rectilinear
ORIENTATION = "horizontal"  # either horizontal or vertical
TAXONOMIC_ADVANTAGE = 0  # taxons advantage over worse identifiable taxons
THRESHOLD_LOG = 0.01

# Styling options
DPI = 300
FIGSIZE = [16, 6]
N_COLORS = 20  # <=20
R_UPPER_MARGIN_REL = 1.02  # Radial margin top, relative to the maximum
BAR_WIDTH_REL = 0.5
COLOR_OTHERS = [0.0, 0.0, 0.0, 1]
COLOR_UNASSIGNED = [0.4, 0.4, 0.4, 1]

# Other options
FILENAME_IN = f'kvalues_{GENE_NAME}_{LEVEL_NAME}.tsv'
FILENAME_OUT = f"gpattern_{GENE_NAME}_{LEVEL_NAME}.png"


# ######################################################
# ##################### LOAD DATA ######################
df = pd.read_csv(FILENAME_IN, sep='\t', header=0)

taxons = df['Taxon'].unique()
clusters = df['Cluster'].unique()
data = [taxons, CONTEXTS, clusters]
n_taxons, n_ctxts, n_clusters = len(taxons), len(CONTEXTS), len(clusters)

K = gtutils.from_df_to_ktensor(df, data, column="k-value").astype("float")

verboseprint("Loaded data for:", VERBOSE)
verboseprint(f"   {n_taxons:4d}\ttaxons", VERBOSE)
verboseprint(f"   {n_ctxts:4d}\tcontexts", VERBOSE)
verboseprint(f"   {n_clusters:4d}\tclusters.", VERBOSE)

# ######################################################
# ################### ANALYZE DATA #####################
# All UNASSIGNED taxons add to the class "Unassigned"
# ... we compute the taxonomic scores
idx_unassigned = gtutils.find_idx_of_elements(taxons, UNASSIGNED)

# All taxons in the top MAX_TAXONS_SHOWN add to the class "Show"
_maxshown = np.clip(MAX_TAXONS_SHOWN, 0, n_taxons)

if CONTRIBUTION == "single":
    if DISPLAY_KIND == "common":
        idx_show = gtutils.kmax_taxons(K, _maxshown)
        k_min = K[idx_show[-1], :, :].max()
        _symbol = r"\leq"

    elif DISPLAY_KIND == "rare":
        idx_show = gtutils.kmin_taxons(K, _maxshown)
        k_min = min(i for i in K[idx_show[-1], :, :].flatten() if i > 0)
        _symbol = r"\geq"

    else:
        assert False

elif CONTRIBUTION == "summed":
    if DISPLAY_KIND == "common":
        contrib_per_taxon = np.sum(np.sum(K, axis=2), axis=1)
        tx_scores = np.array([gtutils.score_taxon(tx) for tx in taxons])
        scores = np.log10(contrib_per_taxon) + tx_scores*TAXONOMIC_ADVANTAGE
        threshold = np.sort(scores)[-_maxshown]
        idx_show = np.argwhere(scores >= threshold)[:, 0]
        k_min = contrib_per_taxon[idx_show].min()
        _symbol = r"\leq"

    elif DISPLAY_KIND == "rare":
        print(">> ERROR << Still not implemented")


# All taxons left add to the class "Others"
idx_other = list(set(range(n_taxons)) - set(idx_show))

# Show unassigned anyway
idx_other = list(set(idx_other) - set(idx_unassigned))
idx_show = list(set(idx_show) - set(idx_unassigned))

# Build K matrix to be shown, by summing the Kvalues of others, and unassigned.
taxon_labels = [f"Others, K${_symbol}${k_min:1.1f}",
                "Unassigned",
                *taxons[idx_show]]
K_shown = np.zeros((len(idx_show)+2, n_ctxts, n_clusters))
K_shown[0, :, :] = np.sum(K[idx_other, :, :], axis=0)
K_shown[1, :, :] = np.sum(K[idx_unassigned, :, :], axis=0)
K_shown[2:, :, :] = K[idx_show, :, :]


verboseprint(f"Found {len(idx_unassigned)} unassigned taxons:", VERBOSE)
_ = [verboseprint(f"  {tx}", VERBOSE) for tx in taxons[idx_unassigned]]
verboseprint("", VERBOSE)
verboseprint(f"Will group {len(idx_other)} taxons into 'Others'.", VERBOSE)
verboseprint(f"Showing the top {len(idx_show)} contributing taxons, with",
             VERBOSE)
verboseprint(f"  K-values >= {k_min:1.3f}.", VERBOSE)
verboseprint("", VERBOSE)


# General plotting properties
# ... style properties and labels
colors = mpl.colormaps.get_cmap('tab20')(np.linspace(0, 1, N_COLORS))
colors = np.delete(colors, [14, 15], axis=0)
# colors[14] = np.array((61, 96, 76, 255))/255
# colors[15] = np.clip(colors[14] * 1.55, 0, 1)
colors = np.concatenate(([COLOR_OTHERS, COLOR_UNASSIGNED], colors), axis=0)
mpl.rcParams['hatch.linewidth'] = 0.3  # previous pdf hatch linewidth
hatches = ['', '/////////////', '.............']


# ######################################################
# ##################### PLOT DATA ######################
if DISPLAY_MODE == "linear":
    X = K_shown.copy()
    verboseprint("Displaying stacked k_values, in LINEAR mode.", VERBOSE)

elif DISPLAY_MODE == "log10":
    X = np.log10(THRESHOLD_LOG + K_shown)
    X[X < 0] = 0
    verboseprint("Displaying stacked log10(k_values), in LOG10 mode.", VERBOSE)


# ... derived values and matrices
sumX = np.sum(X, axis=0)    # total K per cluster per context
maxX = np.max(np.max(X, axis=2), axis=1)  # maximum contribution per taxon

theta = np.linspace(0.0, 2 * np.pi, n_clusters, endpoint=False)
Dtheta = theta[1]-theta[0]
dr = np.round(R_UPPER_MARGIN_REL*sumX.max()/4)
rlims = np.arange(0, 1+np.ceil(R_UPPER_MARGIN_REL*sumX.max()/dr)*dr, dr)
width = BAR_WIDTH_REL*Dtheta
n_cols, n_rows = 0, 0

if ORIENTATION == "horizontal":
    n_cols = int(np.clip(n_ctxts+1, 1, 3))
    n_rows = int(np.ceil((n_ctxts+1) / 4))
    title_padding = 18
    FIGSIZE[1] = FIGSIZE[1] * n_rows

elif ORIENTATION == "vertical":
    n_rows = int(np.clip(n_ctxts+1, 1, 3))
    n_cols = int(np.ceil((n_ctxts+1) / 4))
    title_padding = 0
    FIGSIZE[0] = FIGSIZE[0] * n_rows

FIGSIZE = np.array(FIGSIZE) / 2.3

HFigure = plt.figure(figsize=FIGSIZE, dpi=DPI)
bottoms = np.zeros((n_ctxts, n_clusters))

for _context in range(n_ctxts):
    ax = plt.subplot(n_rows, n_cols, _context+1, projection=PROJECTION)

    # Plot Others and Unassigned
    ax, _bot = gtutils.stackbar(ax, theta, X[:2, _context, :],
                                width=width,
                                colors=colors[:2, :],
                                hatches=hatches,
                                labels=taxon_labels[:2])

    # Plot the rest, nicely with hatches and colors
    bottoms[_context, :] = _bot
    ax, _ = gtutils.stackbar(ax, theta, X[2:, _context, :],
                             width=width,
                             colors=colors[2:],
                             hatches=hatches,
                             labels=taxon_labels[2:],
                             bottom=bottoms[_context, :])

    if PROJECTION == "polar":
        ax.set_rgrids(rlims, labels="")
        ax.set_thetagrids(theta*57.3, labels=clusters, fontsize=7)
        ax.tick_params(pad=-1)

    elif PROJECTION == "rectilinear":
        if _context == 0:
            ax.set_ylabel("K-value (log$_{10}$)")
        ax.set_yticks(rlims)
        ax.set_xticks(range(n_clusters), labels=clusters, fontsize=7)
        ax.grid(True)

    ax.xaxis.grid(linewidth=0.1)
    ax.yaxis.grid(linewidth=0.2)
    ax.set_axisbelow(True)
    ax.set_title(CONTEXTS[_context], {"fontsize": 11}, pad=title_padding)

# ... place legend
kwgs_legend = {}
if ORIENTATION == "horizontal":
    kwgs_legend = {
        "loc": "upper center",
        "ncol": 4,
        "bbox_to_anchor": (0.5, -0.15)
    }
elif ORIENTATION == "vertical":
    kwgs_legend = {
        "loc": "center left",
        "ncol": 1,
        "bbox_to_anchor": (1.2, 0.5)
    }

plt.subplot(n_rows, n_cols, 2)
legend = plt.legend(fontsize=6, **kwgs_legend)

# Adjust horizontal padding
plt.subplots_adjust(wspace=0.35)

# Save figure
gtutils.save_figure(HFigure, FILENAME_OUT,
                    overwrite=OVERWRITE,
                    dpi=DPI,
                    bbox_inches='tight'
                    )

# Final print
verboseprint("The radial units are:", VERBOSE)
verboseprint(f"{rlims}", VERBOSE)
