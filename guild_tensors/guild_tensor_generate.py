# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:48:47 2022

@author: Pablo Yubero

This script reads a master table and outputs a file with K values for the
specified GENE_NAME per CONTEXT and LEVEL (taxonomic). Values of K are
computed as a function of the abundances (A) and counts of sequences (C)
as per Kfun(A,C).

The mastertable should have AT LEAST the following columns:
    gene_fun : gene function by either GTDB or SQM. By specifyinh GENE_NAME you
    will effectively study a single gene function per run of the script.

    cluster_id : cluster identifier, see "aumatic sequence clustering" from
    this same project and paper. This variable can distinguish between
    implementations for the same gene, or between genes for a same route.

    Species : or some taxonomic category. Please, in case you would
    like to group taxons using a different level specify variable LEVEL_NAME.

    TPM : these are the abundances found in metagenomic samples in transcripts
    per million. If you wish to process raw_counts or some other abundance
    measure, please modify the column name in a copy of your original file.


* INPUTS *

FILENAME  : str
    Input of master table (as per Juan Rivas) in TSV <- important.

GENE_NAME : str
    Gene of interest to compute the K values

LEVEL_NAME: str
    Taxonomic level to generate K values, eg, Species_GTDB. It NEEDS to be the
    name of a column of master table.

CONTEXTS : numpy.array of str
    Please specify the contexts of interest in the desired order. If None
    they will get sorted automatically.

REGRESSION_LEVEL : str either "gene", "cluster", "context" or "context_and_cluster"
    Specifies which points are used to compute the model that relates abundance
    and diversity. It can be useful to benefit certain implementations
    (or clusters) differently for their observed abundance.

NORMALIZE_NSAMPLES : bool
    Set whether to normalize k-values by the number of samples included in each
    context. This is to correct for highly unbalanced contexts, for example,
    imagine comparing a context with a single sample, and another with 100.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import guild_tensor_utils as gtutils
from guild_tensor_utils import verboseprint


FILENAME = "mastertable.tsv"
GENE_NAME = 'potF'
LEVEL_NAME = 'Species_GTDB'
REGRESSION_LEVEL = "context_and_cluster"
CONTEXTS = np.array(["Epipelagic", "Mesopelagic", "Bathypelagic"])
NORMALIZE_NSAMPLES = False
SAMPLECODE_COLUMN = "Samples"
VERBOSE = True
EXPORT_PLOT = True
EXPORT_LEGACY = False
# ...
out_filename = f'kvalues_{GENE_NAME}_{LEVEL_NAME}.tsv'
out_plot = f"loglog_regression_{GENE_NAME}_{LEVEL_NAME}.png"

# Repair mastertable if some SQM taxonomic tags are missing
# master_table = gtutils.repair_SQMtaxonomy(FILENAME, save=True)

# Load mastertable
master_table = pd.read_csv(FILENAME, sep="\t")

# ... get rid of non-contributing rows
master_table = master_table[master_table["TPM"] > 0]
verboseprint(f"Loaded mastertable with {len(master_table)} rows.")

# Validate master table
assert LEVEL_NAME in master_table.columns
assert "gene_fun" in master_table.columns
assert GENE_NAME in master_table["gene_fun"].to_list()
assert "cluster_id" in master_table.columns

#####################
# Compute adu_table
adu_table = gtutils.build_adu_table(
    master_table,
    GENE_NAME,
    LEVEL_NAME,
    force_build=True
)

#####################
# Linear regressions
H = plt.figure("plots", figsize=(12, 4), dpi=300)
plt.subplot(1, 2, 1)

if REGRESSION_LEVEL == "gene":
    print(f"\n----- {GENE_NAME} -----")
    delta = gtutils.compute_delta(
        adu_table,
        verbose=VERBOSE,
        printfigure=H,
        printlabel=GENE_NAME
    )

elif REGRESSION_LEVEL == "cluster":
    delta = np.zeros(len(adu_table))
    for cluster in adu_table["Cluster"].unique():
        print(f"\n----- {cluster} -----")

        idx = (
            adu_table["Cluster"] == cluster
        ).to_numpy()

        delta[idx] = gtutils.compute_delta(
            adu_table[idx],
            verbose=VERBOSE,
            printfigure=H,
            printlabel=cluster
        )

elif REGRESSION_LEVEL == "context":
    delta = np.zeros(len(adu_table))
    for context in adu_table["Context"].unique():
        print(f"\n----- {context} -----")

        idx = (
            adu_table["Context"] == context
        ).to_numpy()

        delta[idx] = gtutils.compute_delta(
            adu_table[idx],
            verbose=VERBOSE,
            printfigure=H,
            printlabel=context
        )

elif REGRESSION_LEVEL == "context_and_cluster":
    delta = np.zeros(len(adu_table))
    for context in adu_table["Context"].unique():
        for cluster in adu_table["Cluster"].unique():
            print(f"\n----- {context} : {cluster} -----")

            idx_context = (adu_table["Context"] == context).to_numpy()
            idx_cluster = (adu_table["Cluster"] == cluster).to_numpy()
            idx = idx_cluster * idx_context

            delta[idx] = gtutils.compute_delta(
                adu_table[idx],
                verbose=VERBOSE,
                printfigure=H,
                printlabel=f"{context},{cluster}"
            )

else:
    raise ValueError(
        f'''Error as REGRESSION_LEVEL={REGRESSION_LEVEL} needs to be
        either 'gene', 'cluster', 'context' or "context_and_cluster'.'''
    )

adu_table["delta"] = delta


# ######################
# Normalize K-values according to the number of samples in contexts
if NORMALIZE_NSAMPLES:
    normalization = np.ones(len(adu_table))

    for ctx in CONTEXTS:
        nsamples = gtutils.compute_number_samples(master_table, ctx, SAMPLECODE_COLUMN)
        normalization[adu_table["Context"] == ctx] = nsamples

    adu_table["normalization"] = normalization
    verboseprint("k-values have been normalized.")
else:
    verboseprint("k-values have NOT been normalized.")
    adu_table["normalization"] = 1


# ####################
# Export data
adu_table["k-value"] = adu_table["Abundance"] * adu_table["delta"] / adu_table["normalization"]
adu_table.to_csv(out_filename, sep="\t", index=False)
verboseprint(f"Data saved in {out_filename}.", VERBOSE)


# ####################
# Draw plot
if EXPORT_PLOT:

    valid = adu_table["Abundance"] > 0
    x = adu_table["Abundance"][valid].to_numpy().astype("float")
    y = adu_table["Diversity"][valid].to_numpy().astype("float")

    # Transform data to loglog
    logx = np.log10(1e-10 + x)
    logy = np.log10(1e-10 + y)

    plt.figure(H)
    plt.subplot(1, 2, 1)
    plt.grid()
    # plt.legend()

    plt.subplot(1, 2, 2)
    plt.scatter(x, x * delta[valid], s=delta[valid] * 10, c=logy)
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
    _filepath = f"legacy_kMatrixPerTaxon_{GENE_NAME}_{LEVEL_NAME}.csv"
    gtutils.export_legacy(adu_table, _filepath, column="k-value",
                          contexts=CONTEXTS)
