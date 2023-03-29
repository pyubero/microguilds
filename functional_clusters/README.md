#  Functional clustering #
In these scripts we identify inner branches whose leaves are enriched by one, or several, environmental features. To do so, we will be using two files:

- a tree file in newick format (or any other format compatible with Biopython)
- an environmental_data file, which contains a table of taxons/species and the environmental conditions of their habitat, or for optimal growth derived from the literature.
 
Then, for every internal node of the tree, we compute the observed mean value of each environmental feature, and compare it to the distribution observed under 10^4 randomizations to obtain their z-scores (`s1_node_enrichment.py`). In this way, nodes are considered to be enriched by an environmental feature if the absolute value of the z-score is $>3$. To study the enrichment of larger branches, we focus on nodes that are significant but whose parent is not. The file `s2_tree_clustering-py` will generate a TSV file with all the information of nodes that have at least one significant environmental feature. This table can then be analyzed to manually determine the most relevant functional clusters.
 

### Tree leaf labelling
It is fundamental to have a tree whose leaves and internal nodes are properly labeled. Internal nodes can be simply labeled as `IN_X_YY` where X is the node id and YY is the bootstrap value. Leaves on the contrary need to be labelled as `XXXX_cds_YYYY_TAXONOMICTAGS` where XXXX can be a reference code, YYYY is the accession number of the sequence, an TAXONOMICTAGS are indeed, taxonomic tags as `_d_DOMAIN_p_PHYLUM_c_CLASS_o_ORDER_f_FAMILY_g_GENUS_s_SPECIES`.