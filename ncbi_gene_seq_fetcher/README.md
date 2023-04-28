# NCBI Sequence fetchers
Here I provide a few python (cli) scripts to automatically download gene sequences for **a list** of organisms.

The script `get_sequences_from_ncbi.py` retrieves refseqs, assembly or accession number (in this order of preference if some are unavailable) and then downloads and saves the sequence of interest in a fasta file. This is obtained by using NCBI's search API, and EUTILS if necessary. 

The script `get_sequences_tgtloci.py` is used when we are looking for sequences that can be accessed through the TargetedLoci project in NCBIs ftp server. For example, bacterial 16S, 5S or others. In this case the user needs to provide the TargetedLoci file, and the script will automatically extract the sequences of the CDS for the list of organisms of interest.

The script `fasta_check.py` is simply used to find whether the previous protocols have generated files with duplicate sequences, which typically raises errors in downstream applications.

The file `organisms_input.csv` is an example input file used in the reference [1].


[1] Marine ecological dynamics revealed from a microbial guild approach. https://doi.org/XXXXXXXXXXXX
