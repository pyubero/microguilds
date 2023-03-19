# Generate guild compass plots
These two scripts will guide you through the extraction of meaningful data from a master table, to the creation of the associated guild tensor K, and its visualization as a series of compasses.

### Necessary starting data
Before running anything, be sure that you have a properly formatted `master_table.tsv` with at least four columns being:
- `gene_fun` includes the specific gene function associated with a given sequence,
- `clustr_id` is the functional cluster associated to the sequence,
- `Species_GTDB` is the species taxonomic level associated to the sequence,
- `TPM` is the abundance in transcripts per million of the sequence.

### Generating
Open `generate_guild_tensor.py` and go past the function definitions to where initial variables are specified. Modify FILENAME, GENE_NAME and LEVEL_NAME according to your project. Remember that the master table needs to be in TSV format.

If EXPORT_LEGACY is set to True, another file more "humanly readble" will be created in plain text to evaluate guild structure by taxon.

If EXPORT_PLOT is True, then a nice plot will be saved that includes the empirical model relating the abundance and the diversity of a sequence in panel A, and in panel B it compares sequence importance before and after favoring diversity.

It is advised that you keep VERBOSE set to True.

The script will generate a TSV table that describes the guild tensor.

### Visualization
Open `visualize_guild_tensor.py` get past function definitions and start tuning your plot. Here I will describe a few of the easily tunable parameters.

- `MAX_TAXONS_SHOWN` : As the number of taxons in a metagenomic experiment is typically large,
you can limit the number of taxons to be distinguished in the legend and with particular colors. The guild importance of taxons not shown is summed and displayed as Others. We also indicate the K value threshold that delimits shown taxons and others. It should be a positive integer.

- `DISPLAY_MODE` : As implementation importance typically spans several orders of magnitude, DISPLAY_MODE lets you choose to display stacked values of K, or of log10(K). It can be set to either "linear" or "log10".

- `CONTRIBUTION` : When selecting the N taxons that contribute the most, you can select those whose sum(k_mnl) over implementations and contexts is largest, or those with largest individual contributions in any context and any cluster. It can be set to either "summed" or "single".

- `DISPLAY_KIND` : [Untested] This lets you choose between displaying those that contribute the most or those that contribute the least. This is relevant to study rare microbiomes. It can be set to either "common" or "rare.


**IMPORTANT NOTE** Although in general sum(k_mnl) is different than sum(log10(k_mnl)) for the m-th cluster and n-th context, if MAX_TAXONS_SHOWN is sufficiently large the qualitative results are almost identical.
