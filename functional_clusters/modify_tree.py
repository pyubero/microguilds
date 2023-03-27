from time import sleep
import requests
from Bio import Phylo
from tqdm import tqdm
import functional_clustering_utils as fcutils

# accession = ACCESSION
# db = DATABASE
# el 208 de 16S


def accession_to_taxoninfo(accession, db):
    '''Given an accession code of a protein, it retrieves the information
    on the taxon from domain to species.
    This function helps rename the leafs of some trees.'''
    def clean_string(string):
        # Remove tags
        tags = ["\n", ";", ",", "."]
        for tag in tags:
            string = string.replace(tag, " ")

        # Remove double whitespaces
        while "  " in string:
            string = string.replace("  ", " ")

        # Remove trailing white spaces
        if string[0] == " ":
            string = string[1:]

        if string[-1] == " ":
            string = string[:-1]

        return string

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params = {
              'db': db,
              'id': accession,
              'retmode': "json",
              "rettype": "json",
              'sort': 'relevance'
             }

    sleep(1/5)
    response = requests.get(url, params=params, timeout=10)

    if response.status_code > 200:
        return None

    if db == "nucleotide":
        split_word = "REFERENCE"
    elif db == "protein":
        split_word = "COMMENT"
    else:
        raise ValueError(f"Name of db {db} is invalid.")

    chunk = response.text.split("ORGANISM")[1].split(split_word)[0].split("\n")
    organism_name = clean_string(chunk[0]).split(" ")
    organism_taxonomy = clean_string(chunk[1]).split(" ")

    while len(organism_name) < 2:
        organism_name = [*organism_name, ""]
    organism_name = organism_name[:2]

    while len(organism_taxonomy) < 5:
        organism_taxonomy = [*organism_taxonomy, ""]
    organism_taxonomy = organism_taxonomy[:5]

    taxonomy = [*organism_taxonomy, *organism_name]

    tags = ["d_", "p_", "c_", "o_", "f_", "g_", "s_"]

    final = ""
    for _tag, _name in zip(tags, taxonomy):
        if len(_name) > 0:
            final = final + _tag + _name + "_"

    return final[:-1]


FILENAME_IN = "tree_rplB.newick"
FILENAME_OUT = "tree_rplB_new.newick"
TREETYPE = "newick"
DATABASE = "protein"  # either protein or nucleotide


TREE = Phylo.read(FILENAME_IN, TREETYPE)
LEAFS = TREE.get_terminals()
NODES = TREE.get_nonterminals()

print(TREE.get_nonterminals()[208].name)
print(TREE.get_terminals()[208].name)


# leaf = TREE.get_terminals()[0]
# accession = leaf.name.split(".")[0]
# db = "nucleotide"
# leaf = LEAFS[208]

for leaf in tqdm(LEAFS):
    # Uncomment for rplB and headers of the likes with cds
    OTHER = leaf.name.split("cds")[0]
    ACCESSION = leaf.name.split("cds")[-1][1:].split(".")[0]
    # Uncomment for 16S and plain accession number headers
    # OTHER = "XXX_"
    # ACCESSION = leaf.name.split(".")[0]

    TAXONINFO = accession_to_taxoninfo(ACCESSION, DATABASE)
    if TAXONINFO is None:
        print(f"Could not find info on {ACCESSION}.")

    leaf.name = OTHER + "cds_" + ACCESSION + ".1_" + TAXONINFO


for nn, node in tqdm(enumerate(NODES)):
    node.name = f"IN_{nn}_0"
    print(node.name)


Phylo.write(TREE, FILENAME_OUT, format="newick")


# Check that everything worked fine
data = fcutils.get_clade_data(FILENAME_OUT, treetype="newick")
clade_ids, clade_lfs, clade_dpt, leaf_names = data
print(clade_ids[3])
print(leaf_names[3])
