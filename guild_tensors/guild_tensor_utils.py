# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:47:36 2022

@author: logslab
"""
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics import r2_score


def verboseprint(msg, verbose=True, condition=True):
    '''Prinst a message if VERBOSE is True.'''
    if verbose and condition:
        print(msg)


def verbosebar(iterable, verbose=True):
    '''Generates a progress bar with tqdm if VERBOSE is True.'''
    if verbose:
        return tqdm(iterable)
    else:
        return iterable


def compute_adu(df, inputs, loc):
    '''
    Compute A(bundance) D(iversity) and U(nivocity) from a gene master table.

    Parameters
    ----------
    df : DataFrame
        Master table at the GENE or FUNCTION levels.
    inputs : list of list of str
        Typically a list of [taxons, contexts, clusters] with the names with
        wich they appear in df.
    loc : list of int
        Typically a list [2,0,7] with the indices of the desired output. Thus,
        we would obtain data for taxons[2] in contexts[0] and clusters[7].
        This convoluted way of extracting ADU is derived from a quick linear,
        e.i. single loop, way to extract all data.

    Returns
    -------
    abundance : float
        The sum of abundances of sequences within a taxon, context and cluster.
    diversity : int
        The number of sequences within a taxon, context and cluster.
    univocity : float, between 0 and 1
        The univocity of sequences within a taxon, context and cluster.

    '''
    # Find where the interesting abundances are in the df/master table
    idc = (df["taxonomic_classification_level"] == inputs[0][loc[0]]) & \
          (df["Context"] == inputs[1][loc[1]]) & \
          (df["cluster_id"] == inputs[2][loc[2]])

    # Extract value
    subtable = df[idc]

    # Obtain abundances and relevant parameters
    abundances = subtable["TPM"]
    diversity = len(abundances)
    abundance = sum(abundances)
    univocity = 1.0
    return abundance, diversity, univocity


def export_legacy(df, filename, column="k-value"):
    '''Export data as a series of k-matrices in plain csv as in v0'''
    gene = df["Gene"].iloc[0]
    taxons = df["Taxon"].unique()
    contexts = df["Context"].unique()
    clusters = df['Cluster'].unique()
    new_contexts = np.array(["Epipelagic", "Mesopelagic", "Bathypelagic"])

    ntaxons = len(taxons)
    ncontexts = len(contexts)
    nclusters = len(clusters)

    # Build Kmat
    Kmat = np.zeros((ntaxons, ncontexts, nclusters))
    _mesh = np.meshgrid(range(ntaxons), range(ncontexts), range(nclusters))
    idc = np.array(_mesh).T.reshape(-1, 3)
    for j_tx, j_ct, j_cl in verbosebar(idc):
        idx = (df["Taxon"] == taxons[j_tx]) & \
                (df["Context"] == new_contexts.astype("str")[j_ct]) & \
                (df["Cluster"] == clusters[j_cl])
        assert sum(idx) == 1

        Kmat[j_tx, j_ct, j_cl] = df[idx][column]

    # Prepare output file
    with open(filename, 'w+', encoding="utf-8") as file:
        file.write(f'{gene}\n')
        file.write(','.join([f"{cluster}" for cluster in clusters]) + '\n')

    # For every taxon...
    for j_taxon in verbosebar(range(ntaxons)):

        # ... create a K matrix
        K = Kmat[j_taxon, :, :]

        # For every clusters...
        for _ in range(nclusters):

            # Append to output file...
            header = f'>>{taxons[j_taxon]}\n'
            lines = []
            for _ in range(K.shape[0]):
                lines.append(','.join([f'{k:1.4f}' for k in K[_, :]]) + '\n')

        with open(filename, 'a', encoding="utf-8") as f:
            f.write(header)
            _ = [f.write(line) for line in lines]


def bivariate_regression(x, y):
    '''Computes a bivariate regression as the angle of the cov matrix.'''
    C = np.cov([x, y])
    V = np.linalg.eig(C)[1]
    alfa = V[1][0]/V[0][0]
    beta = np.mean(y) - alfa * np.mean(x)
    return alfa, beta


def check_mastertable(filename: str, save=False):
    '''Checks and repairs the entries of a mastertable without a Species_SQM by
    using those given by GTDB.'''
    df = pd.read_csv(filename, sep="\t")
    _modified = False

    # Find null taxons
    idc_null = np.argwhere(df['Species_SQM'].isnull().values)[:, 0]
    if len(idc_null) == 0:
        verboseprint(f"All {len(df)} entries are apparently valid.")
    else:
        _modified = True
        verboseprint(f"Found {len(idc_null)} null classifications")
        for jj in verbosebar(idc_null):
            df.at[jj, "Domain_SQM"] = "GTDB:"+df.loc[jj, "Domain_GTDB"]
            df.at[jj, "Phylum_SQM"] = "GTDB:"+df.loc[jj, "Phylum_GTDB"]
            df.at[jj, "Class_SQM"] = "GTDB:"+df.loc[jj, "Class_GTDB"]
            df.at[jj, "Order_SQM"] = "GTDB:"+df.loc[jj, "Order_GTDB"]
            df.at[jj, "Family_SQM"] = "GTDB:"+df.loc[jj, "Family_GTDB"]
            df.at[jj, "Genus_SQM"] = "GTDB:"+df.loc[jj, "Genus_GTDB"]
            df.at[jj, "Species_SQM"] = "GTDB:"+df.loc[jj, "Species_GTDB"]

        idc_null = np.argwhere(df['Domain_SQM'].isnull().values)[:, 0]
        verboseprint(f"Still found {len(idc_null)} null classifications.",
                     len(idc_null) > 0)
    if _modified and save:
        verboseprint("Saving modified master table.")
        df.to_csv(filename, sep="\t", index=False)

    return df


def hex_to_rgb(value):
    ''' Function to convert hexadecimal colors into RGB triplets. '''
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def int_to_roman(num):
    ''' Converts an integer to the roman number system. '''
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


def cycle_through(array, length):
    ''' Returns the indices to cycle some array of size N to steps of
    a specific length.
    Typically used to cycle through color maps to increase contrast.'''
    N = len(array)
    cycle = np.ceil(N/length)

    x0 = 0
    k = 1
    new_order = [x0,]
    while len(new_order) < len(array):
        new_idx = x0 + k*cycle

        if new_idx < len(array):
            new_order.append(int(new_idx))

            k += 1
        else:
            x0 += 1
            k = 0

    return np.array(new_order)


def from_df_to_ktensor(df, data, column="k-value", verbose=True):
    ''' Creates the GUILD TENSOR, that is the K-matrices representing the
    guild. It is totally suboptimized but it is reliable.
    df :should be the gene_table, that is, a master table filtered by function.

    data : should be a list of lists of strings with the names of taxons,
    contexts and clusters in the order you wish to have in the array K.'''

    taxons, contexts, clusters = data
    ntaxons = len(taxons)
    ncontexts = len(contexts)
    nclusters = len(clusters)

    Kmat = np.zeros((ntaxons, ncontexts, nclusters), dtype="object")
    _meshgrid = np.meshgrid(range(ntaxons), range(ncontexts), range(nclusters))
    idc = np.array(_meshgrid).T.reshape(-1, 3)
    for j_tx, j_ct, j_cl in verbosebar(idc, verbose):
        idx = (df["Taxon"] == taxons[j_tx]) & \
                (df["Context"] == contexts.astype("str")[j_ct]) & \
                (df["Cluster"] == clusters[j_cl])
        assert sum(idx) == 1
        Kmat[j_tx, j_ct, j_cl] = df[idx][column].iloc[0]

    return Kmat


def score_taxon(taxon: str):
    '''
    Scores taxon name.
    It returns:
        0 if it does not have a defined species nor genus.
        1 if it only has a genus specified.
        2 if it has both, genus and species specified.
    '''
    _tx = taxon.split(' ')
    if (len(_tx) > 1) and (_tx[1][:2] != "sp"):
        has_species = True
    else:
        has_species = False

    is_UBA = (len(_tx) > 0) and (len(_tx[0] > 3)) and (_tx[0][3:6] == "UBA")
    is_GCA = (len(_tx) > 0) and (len(_tx[0] > 3)) and (_tx[0][3:6] == "GCA")
    if is_UBA or is_GCA:
        has_genus = True
    else:
        has_genus = False

    score = len(_tx) + 0.5 * has_species + 0.5 * has_genus - 1
    return score


def stackbar(ax, x, Y, *args, **kwargs):
    ''' Easily plot stacked bars with, optionally, a list of colors and
    hatches.'''

    if "bottom" not in kwargs:
        kwargs.update({"bottom": np.zeros(Y.shape[1], )})

    if "labels" in kwargs:
        labels = kwargs["labels"]
        kwargs.pop("labels")
        _given_labels = True
    else:
        labels = None
        _given_labels = False

    if "colors" in kwargs:
        colors = kwargs["colors"]
        kwargs.pop("colors")
        _given_colors = True
    else:
        colors = None
        _given_colors = False

    if "hatches" in kwargs:
        hatches = kwargs["hatches"]
        kwargs.pop("hatches")
        _given_hatches = True
    else:
        hatches = None
        _given_hatches = False

    if Y.shape[0] == 0:
        return ax, kwargs["bottom"]

    for idx in range(Y.shape[0]):
        if _given_labels:
            kwargs.update({"label": labels[idx]})

        if _given_colors:
            idx_color = idx % colors.shape[0]
            kwargs.update({"color": colors[idx_color, :]})
            if _given_hatches:
                idx_hatch = idx // colors.shape[0]
                kwargs.update({"hatch": hatches[idx_hatch]})

        # Plot
        ax.bar(x, Y[idx, :], *args, **kwargs)

        # ... update bottoms
        new_bottoms = kwargs["bottom"] + Y[idx, :]
        kwargs.update({"bottom": new_bottoms})

        # ... remove label
        if ("label" in kwargs) and (type(kwargs["label"]) == str):
            kwargs.pop("label")

    return ax, kwargs["bottom"]


def kmax_taxons(K, n):
    ''' Find the n taxons that have the largest individual contributions.
    That is, per cluster and per context.'''
    top_taxon_idc = []
    for ii in np.argsort(K.flatten())[-1::-1]:
        taxon_idx = np.argwhere(K == K.flatten()[ii])[0][0]
        if taxon_idx not in top_taxon_idc:
            top_taxon_idc.append(taxon_idx)
        if len(top_taxon_idc) >= n:
            break
    return top_taxon_idc


def kmin_taxons(K, n):
    ''' Find the n taxons that have the smallest individual contributions.
    That is, per cluster and per context.'''
    top_taxon_idc = []
    values = K.flatten()
    for ii in np.argsort(values):
        if values[ii] == 0:
            continue
        taxon_idx = np.argwhere(K == K.flatten()[ii])[0][0]
        if taxon_idx not in top_taxon_idc:
            top_taxon_idc.append(taxon_idx)
        if len(top_taxon_idc) >= n:
            break
    return top_taxon_idc


def do_not_overwrite_filepath(filepath):
    ''' Returns a filepath that will not overwrite existing files.'''
    new_filepath = filepath
    idx = 0
    while os.path.exists(new_filepath):
        idx += 1
        split = filepath.split('.')
        filepath_wo_extension = '.'.join(split[:-1])
        new_filepath = f"{filepath_wo_extension}({idx}).{split[-1]}"

    return new_filepath


def save_figure(hfigure, filepath, *args, overwrite=True, **kwargs):
    '''Save figure with overwrite protection.'''
    if not overwrite:
        filepath = do_not_overwrite_filepath(filepath)

    hfigure.savefig(filepath, *args, **kwargs)


def find_idx_of_elements(input_list, elements):
    '''Returns the indices in input list of different elements.'''
    output_idc = set([])

    for elem in elements:
        idc = np.argwhere(np.array(input_list) == elem)[:, 0]
        output_idc = output_idc.union(set(idc))
    return list(output_idc)


def build_adu_table(master_table, gene, taxonomic_level,
                    force_build=False, verbose=True):
    '''Builds an adu table WITHOUT K from a master_table.
    It loads a table if it is found in the cwd with the standard name.'''

    # Load file if it already exists
    filename_in = f'kvalues_{gene}_{taxonomic_level}.tsv'
    if (not force_build) and (os.path.isfile(filename_in)):
        adu_table = pd.read_csv(filename_in, sep='\t', header=0)

        ntaxons = len(adu_table['Taxon'].unique())
        nclusters = len(adu_table['Cluster'].unique())
        ncontexts = len(adu_table["Context"].unique())

        verboseprint("Loaded data for:", verbose)
        verboseprint(f"   {ntaxons:4d}\ttaxons", verbose)
        verboseprint(f"   {ncontexts:4d}\tcontexts", verbose)
        verboseprint(f"   {nclusters:4d}\tclusters.", verbose)
        return adu_table

    # Standardize taxonomic column name
    master_table = master_table.rename(
        columns={taxonomic_level: "taxonomic_classification_level"}
        )

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

    # Print some output
    verboseprint(f"Found {n_clusters} clusters in gene subtable:", verbose)
    _ = [verboseprint(f"\t{clustername}", verbose) for clustername in clusters]
    verboseprint(f"""Found {n_taxons} taxonomic levels
                 in column *{taxonomic_level}*.""", verbose)
    verboseprint(f"Found {n_contexts} contexts in gene subtable.", verbose)
    _ = [verboseprint(f"\t{ctxt}", verbose) for ctxt in contexts]
    verboseprint("", verbose)

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
        _ab, _dv, _un = compute_adu(gene_table,
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

    return adu_table


def compute_gamma(adu_table, verbose=True):
    '''Computes the bivariate regression to the log-transformed data
    of abundance and diversity from the adu_table to compute k-values.'''
    log_threshold = 1e-10

    # Load data to regress
    idx = adu_table["Abundance"] > 0
    x = adu_table["Abundance"][idx].to_numpy().astype("float")
    y = adu_table["Diversity"][idx].to_numpy().astype("float")

    # Transform data to loglog
    logx = np.log10(log_threshold + x)
    logy = np.log10(log_threshold + y)

    # Linear regression
    gamma, offset = bivariate_regression(logx, logy)
    r2 = r2_score(logy, logx * gamma + offset)

    # Print results
    verboseprint("Bivariate loglog regression results:", verbose)
    verboseprint(f"gamma = {gamma}", verbose)
    verboseprint(f"c = {offset}", verbose)
    verboseprint(f"R2 = {r2}", verbose)
    verboseprint("", verbose)

    return gamma, offset, r2
