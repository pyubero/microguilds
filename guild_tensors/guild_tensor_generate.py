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
    from this same project and paper.

    Species_GTDB : or another taxonomic category. Please, in case you would
    like to group taxons using a different level specify variable LEVEL_NAME.

    TPM : these are the abundances found in metagenomic samples in transcripts
    per million.

* INPUTS *

FILENAME  : str
    Input of master table (as per Juan Rivas) in TSV <- important.

GENE_NAME : str
    Gene of interest to compute the K values

LEVEL_NAME: str
    Taxonomic level to generate K values, eg, Species_GTDB. It NEEDS to be the
    name of a column of master table.
"""
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score
from matplotlib import pyplot as plt
import guild_tensor_utils as gtutils
from guild_tensor_utils import verbosebar, verboseprint


FILENAME = 'mastertable.tsv'
GENE_NAME = 'amt'
LEVEL_NAME = 'Species_GTDB'
VERBOSE = True
EXPORT_PLOT = True
EXPORT_LEGACY = False
# ...
out_filename = f'kvalues_{GENE_NAME}_{LEVEL_NAME}.tsv'
out_plot = f"loglog_regression_{GENE_NAME}_{LEVEL_NAME}.png"


# Import mastertable
master_table = gtutils.check_mastertable(FILENAME, True)
verboseprint(f"Loaded mastertable with {len(master_table)} rows.")

# Validate master table
assert LEVEL_NAME in master_table.columns
assert "gene_fun" in master_table.columns
assert GENE_NAME in master_table["gene_fun"].to_list()
assert "cluster_id" in master_table.columns

# Standardize taxonomic column name
master_table = master_table.rename(
    columns={LEVEL_NAME: "taxonomic_classification_level"}
    )

# Filter by gene name
gene_table = master_table[master_table["gene_fun"] == GENE_NAME]
verboseprint(f"Subtable for gene *{GENE_NAME}* has {len(gene_table)} rows.")

# Find clusters in gene subtable
clusters = gene_table['cluster_id'].unique()
n_clusters = len(clusters)

# Find taxons in gene subtable
taxons = gene_table["taxonomic_classification_level"].unique()
n_taxons = len(taxons)

# Find contexts in gene subtable
contexts = gene_table["Context"].unique()
n_contexts = len(contexts)

# Print some output
verboseprint(f"Found {n_clusters} clusters in gene subtable:", VERBOSE)
_ = [verboseprint(f"\t{clustername}", VERBOSE) for clustername in clusters]
verboseprint(f"Found {n_taxons} taxonomic levels in column *{LEVEL_NAME}*.",
             VERBOSE)
verboseprint(f"Found {n_contexts} contexts in gene subtable.", VERBOSE)
_ = [verboseprint(f"\t{ctxt}", VERBOSE) for ctxt in contexts]
verboseprint("", VERBOSE)

# Process data
idc = np.array(
    np.meshgrid(range(n_taxons), range(n_contexts), range(n_clusters))
    ).T.reshape(-1, 3)

# ... included keeps track of included and missing samples
included = gene_table["taxonomic_classification_level"] == np.random.rand()

# Define accumulator
adu_table = pd.DataFrame(columns=["Gene", "Taxon", "Context", "Cluster",
                                  "Abundance", "Diversity", "Univocity"])
data = [taxons, contexts, clusters]

for j_tx, j_ct, j_cl in verbosebar(idc):
    # Compute K value
    _ab, _dv, _un = gtutils.compute_adu(gene_table, data, [j_tx, j_ct, j_cl])
    new_row = pd.Series(
        {"Gene": GENE_NAME,
         "Taxon": taxons[j_tx],
         "Context": contexts[j_ct],
         "Cluster": clusters[j_cl],
         "Abundance": _ab,
         "Diversity": _dv,
         "Univocity": _un})

    adu_table = pd.concat([adu_table, new_row.to_frame().T], ignore_index=True)

    # To keep track of included and missing samples
    _id = (gene_table["taxonomic_classification_level"] == taxons[j_tx]) & \
          (gene_table["Context"] == contexts[j_ct]) & \
          (gene_table["cluster_id"] == clusters[j_cl])
    included = included | _id

verboseprint(f"\nIncluded {sum(included)} out of {len(gene_table)} rows.",
             VERBOSE)
verboseprint(f"The sum of diversities is {sum(adu_table['Diversity'])}.",
             VERBOSE)
verboseprint("", VERBOSE)


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

# Linear regression
gamma, c = gtutils.bivariate_regression(logx, logy)
r2 = r2_score(logy, logx * gamma + c)

# Print results
verboseprint("Bivariate loglog regression results:", VERBOSE)
verboseprint(f"gamma = {gamma}", VERBOSE)
verboseprint(f"c = {c}", VERBOSE)
verboseprint(f"R2 = {r2}", VERBOSE)
verboseprint("", VERBOSE)

# Compute correction factor delta = d_obs/d_exp
delta = 10**logy / np.clip(10**linear_function(logx, gamma, c), 1, np.inf)
_delta = np.zeros(len(adu_table))
_delta[idx] = delta

# Export data
adu_table["delta"] = _delta
adu_table["k-value"] = adu_table["Abundance"] * _delta
adu_table.to_csv(out_filename, sep="\t", index=False)
verboseprint(f"Data saved in {out_filename}.", VERBOSE)

# Draw plot
if EXPORT_PLOT:
    H = plt.figure(figsize=(12, 4), dpi=300)
    plt.subplot(1, 2, 1)
    plt.scatter(logx, logy, s=10, c=logy)
    plt.plot(logx, linear_function(logx, gamma, c))
    plt.grid()
    plt.xlabel("log10 Abundance")
    plt.ylabel("log10 Diversity")
    plt.legend(
        ["Data", fr"Biv. loglog reg\n$\gamma$={gamma:.3f}\nR2={r2:.3f}"]
        )

    plt.subplot(1, 2, 2)
    plt.scatter(x, x*delta, s=delta*10, c=logy)
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
    gtutils.export_legacy(adu_table, _filepath)
