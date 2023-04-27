import numpy as np
from Bio import Phylo
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score
import functional_clustering_utils as fcutils


def is_valid(treeleaf):
    '''Converts a leaf name to an id.'''
    species = fcutils.get_species(treeleaf.name)
    genus = fcutils.get_genus(treeleaf.name)
    if (species is None) or (genus is None):
        return False, f"{genus} {species}"

    return True, f"{genus} {species}"


GENEX = "16S"
GENEY = "potF"
FILENAME = f"data_tree_comparison_{GENEX}_{GENEY}.npz"

data = np.load(FILENAME)
nof_leafs = data["nof_leafs"]
x = data["distance_x"]
y = data["distance_y"]


# Load rplB for reference
data = np.load("data_tree_comparison_16S_rplB.npz")
_nof_leafs = data["nof_leafs"]
_x = data["distance_x"]
_y = data["distance_y"]


# ... linear regression
ii = np.isfinite(_x)
C = np.cov([_x[ii], _y[ii]])
V = np.linalg.eig(C)[1]
alfa = V[1][0]/V[0][0]
beta = np.nanmean(_y) - alfa * np.nanmean(_x)
r2 = r2_score(_y[ii], _x[ii] * alfa + beta)
residues = _y[ii] - (_x[ii] * alfa + beta)
error = np.std(residues)
print(r2)


# Identify number of identifiable bichos in each clade
tree = Phylo.read("tree_potF.newick", 'newick')
leafs = tree.get_terminals()
clades = tree.get_nonterminals()
nof_valid = []
for clade in clades:
    bool_vec = [is_valid(leaf)[0] for leaf in clade.get_terminals()]
    nof_valid.append(np.sum(bool_vec))

idx_valid = np.argwhere(np.array(nof_valid) == 2)[:, 0]
dist = [x[idx]/y[idx] for idx in idx_valid]
_ = idx_valid[np.argwhere(np.array(dist) > 3)]

for idx in [265, 635, 650, 695, 707]:
    clade = clades[idx]
    print(f"\n----- Clade #{idx} ------")
    print(f"Distance in 16S {x[idx]:.4f}")
    print(f"Distance in potF {y[idx]:.4f}")
    for leaf in clade.get_terminals():
        suc, name = is_valid(leaf)
        if suc:
            print(name)

CLUSTERIDX = [37, 257, 361, 484, 866, 869, 40]
clusternames = ["c1", "c1a", "c1b", "c1b", "c2a", "c2b", "c3"]
clustercolors = ["#d2ffff",  # c1
                 "#90f9ff",  # c1a
                 "#3fe7ff",  # c1b
                 "#3fe7ff",  # c1b
                 "#f7d5ff",  # c2a
                 "#ffb1fc",  # c2b
                 "#feffe1"]  # c3

clusterclades = []
for idx in CLUSTERIDX:
    idc = [int(cld.name.split("_")[1])
           for cld in clades[idx].get_nonterminals()]
    clusterclades.append(idc)

# Prepare some plotting options
_xx = np.array([0, 1])
style_bounds = {"lw": 1,
                "color": np.ones(3,)*0}

# Create figure
plt.figure()

# Plot potF data per cluster
for jj, idx in enumerate(CLUSTERIDX):
    plt.scatter(x[clusterclades[jj]],
                y[clusterclades[jj]],
                color=clustercolors[jj],
                marker='.',
                s=nof_leafs[clusterclades[jj]],
                edgecolors="k",
                linewidths=0.2)

# Plot rplB
plt.scatter(_x, _y,
            color=style_bounds["color"],
            alpha=0.8,
            marker='.',
            linewidths=0,
            s=_nof_leafs)

# Plot null model error lines
plt.plot(_xx, _xx * alfa + beta, **style_bounds)
plt.plot(_xx, _xx * alfa + beta + 2*error, '--', **style_bounds)
plt.plot(_xx, _xx * alfa + beta + 3*error, ':', **style_bounds)
plt.plot(_xx, _xx * alfa + beta - 2*error, '--', **style_bounds)
plt.plot(_xx, _xx * alfa + beta - 3*error, ':', **style_bounds)

# Plot clustering indices with a red edge
plt.scatter(x[CLUSTERIDX], y[CLUSTERIDX],
            color="None",
            edgecolors='r',
            s=nof_leafs[CLUSTERIDX]/4)

plt.legend([*clusternames,
            "rplB",
            "null model",
            r"$\pm2\sigma$",
            r"$\pm3\sigma$"],
           loc="center left",
           bbox_to_anchor=(1, 0.5))
plt.xlabel(f"Mean leaf distance in {GENEX}")
plt.ylabel("Mean leaf distance")
plt.ylim((0, 0.97))
plt.xlim((-0.007, 0.85))
plt.grid()
plt.tight_layout()
plt.savefig("filo_clustering-png", dpi=300)
plt.show()
