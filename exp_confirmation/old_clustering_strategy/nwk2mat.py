#!/usr/bin/env python
# Script from https://github.com/mgalardini
#


if __name__ == "__main__":
    import sys
    import pandas as pd
    import itertools
    import numpy as np
    import pickle as pkl
    from Bio import Phylo
    from tqdm import tqdm 
    
    ifile = 'ref_database_relabel.newick'
    
    
    t = Phylo.read(ifile, 'newick')
    leafs = t.get_terminals()
    nleafs= len(leafs)
    ncomb = (nleafs**2-nleafs)/2

    d = {}
    for x, y in tqdm(itertools.combinations(t.get_terminals(), 2), total=ncomb):
        v = t.distance(x, y)
        d[x.name] = d.get(x.name, {})
        d[x.name][y.name] = v
        d[y.name] = d.get(y.name, {})
        d[y.name][x.name] = v
    for x in t.get_terminals():
        d[x.name][x.name] = 0

    # Convert dictionnary to dataframe
    df = pd.DataFrame(d)

    # Convert dataframe to numpy matrix
    M = df.to_numpy()
    names = np.array( df.index.values.tolist() )

    # Sort elements such that the diagonal is all 0s
    I = []
    for ii in range(nleafs):
        idx = np.argwhere(M[ii]==0)[:,0][0]
        I.append(idx)
    I= np.array(I)
    M = M[:,I]
    names= names[I]
    
    # Check if matrix is symmetric
    assert np.allclose(A, A.T)
    
    np.savez('tree_distance.npz', M=M, names=names)
    
    
    