# microguilds

Welcome to `microguilds` :tada:
Microguilds is a python project that aims to simplify and automate the analysis of microbial guilds given metagenomic samples (or other types of data). Recall that an [ecological guild](https://en.wikipedia.org/wiki/Guild_(ecology)) is composed of individuals that occupy the same resource space despite exploiting it differently (methods, pathways or sequences). 

## Contents
This repo is divided into three different aspects related to microbial guild analysis:

1. `functional_clustering`: contains all scripts related with the clustering of a tree's branches given some environmental data. In this way we distinguish, in a meaningful way, between implementations of a function that have adapted to specific environmental features. [More info](https://github.com/pyubero/uGuilds/tree/main/functional_clustering).
2. `guild_tensors`: once the different functional implementations have been identified, we can quantitatively characterize the contribution of taxons and implementations to a microbial guild across different environmental contexts. These scripts help quantify, visualize and analyze microbial guilds. [More info](https://github.com/pyubero/uGuilds/tree/main/guild_tensors)
3. [[obsolete](https://github.com/pyubero/neutral-drift)] `sequence_fetcher` from NCBI's databases: the idea is that given a list of organisms in plain text file, you can automatically retrieve the sequence of any gene using either NCBI's assemblies or E-tools. [More info](https://github.com/pyubero/uGuilds/tree/main/ncbi_gene_seq_fetcher).


## Installation

To use migroguilds you don't need to install anything (for now it is not meant to be a package) but rather a collection of useful functions or cli applications. Simply clone this repo and start running our own developped applications or write your own!

- To clone the repository use: `git clone https://github.com/pyubero/microguilds.git`

## Tutorials
For a better experience, we have create some [tutorials](https://github.com/pyubero/uGuilds/tree/main/Tutorials) on how we personally use the functions available in this repo. They are somewhat deeply commented examples using mock organisms and a mock mastertable.


## Contributing and bug reports :bug:

We are always looking for people interested in helping with code development, documentation... If you wish to contribute, please first send any of us an email.

If you find any bug please do not hesitate to open an issue.

## License 
This code is offered under a GNU General Public License v3.0. 

## Citation
Please cite the following paper for academic references:

```
@article{microguilds,
   author    = "Juan Rivas and Pablo Yubero",
   title     = "Quantifying microbial guilds",
   year      = 2023
 }
```

## Credits
microguilds mainly involves Juan Rivas and Pablo Yubero (both PhD students at CNB, Madrid).

