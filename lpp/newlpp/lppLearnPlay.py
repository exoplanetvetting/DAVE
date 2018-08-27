#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 22:02:05 2018

@author: smullally
"""

import numpy as np
import matplotlib.pyplot as plt



from sklearn.datasets import make_blobs

X, y = make_blobs(1000, n_features=300, centers=4,
                  cluster_std=8, random_state=42)

fig, ax = plt.subplots(2, 2, figsize=(10, 10))
rand = np.random.RandomState(42)

for axi in ax.flat:
    i, j = rand.randint(X.shape[1], size=2)
    axi.scatter(X[:, i], X[:, j], c=y)
#%%  
from lpproj import LocalityPreservingProjection
lpp = LocalityPreservingProjection(n_components=2,n_neighbors=5)

X_2D = lpp.fit_transform(X)


plt.figure()
plt.scatter(X_2D[:, 0], X_2D[:, 1], c=y,s=3)

plt.title("Projected from 300->2 dimensions");

#%%
#now to transform a set of new ones based on that training.
newones=X[0:10] +X[0:10]*.01

new_10D=lpp.transform(newones)

plt.scatter(new_10D[:,0],new_10D[:,1],c=y[0:10],marker='x',s=12)

#%%

#%%

from sklearn.decomposition import PCA
Xpca = PCA(n_components=2).fit_transform(X)
plt.scatter(Xpca[:, 0], Xpca[:, 1], c=y);

