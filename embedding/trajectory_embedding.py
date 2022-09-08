# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 17:21:08 2022

@author: logslab
"""

import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from sklearn.manifold import TSNE, Isomap, SpectralEmbedding, LocallyLinearEmbedding, MDS #LocallyLinearEmbedding

Nworms= 10
Nsamp = 3000


e1 = [0.1, 0.5]*np.random.randn(Nsamp,2) + [0, 3]
e2 = [0.1, 0.1]*np.random.randn(Nsamp,2) + [0, 0]



WORMS = []

for _ in tqdm( range(Nworms)): 
    rr = np.random.rand()
    Ne1 = np.random.randint(0, high = Nsamp, size=int(rr*Nsamp))
    Ne2 = np.random.randint(0, high = Nsamp, size=int( (1-rr)*Nsamp))
    
    _features = np.vstack( (e1[Ne1,:] , e2[Ne2,:]) )
    H = plt.hist2d(_features[:,0], _features[:,1], bins=30)
    WORMS.append( H[0] )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    