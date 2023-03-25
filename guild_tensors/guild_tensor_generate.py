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
import numpy as np
from matplotlib import pyplot as plt
import guild_tensor_utils as gtutils
from guild_tensor_utils import verboseprint


FILENAME = 'mastertable.tsv'
GENE_NAME = 'amoA'
LEVEL_NAME = 'Species_GTDB'
VERBOSE = True
EXPORT_PLOT = True
EXPORT_LEGACY = True
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

# Compute adu_table
adu_table = gtutils.build_adu_table(master_table,
                                    GENE_NAME,
                                    LEVEL_NAME,
                                    force_build=True)


# Linear regression
def linear_function(x_values, slope, offset):
    '''Returns y = slope * x + offset'''
    return slope * x_values + offset


gamma, c, r2 = gtutils.compute_gamma(adu_table, VERBOSE)
print(f"Gene: {GENE_NAME} with R2={r2:0.2f}")

# Load data of regression
idx = adu_table["Abundance"] > 0
x = adu_table["Abundance"][idx].to_numpy().astype("float")
y = adu_table["Diversity"][idx].to_numpy().astype("float")

# Transform data to loglog
log_threshold = 1e-10
logx = np.log10(log_threshold + x)
logy = np.log10(log_threshold + y)

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
    gtutils.export_legacy(adu_table, _filepath, column="Diversity")
