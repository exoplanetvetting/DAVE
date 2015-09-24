# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 16:19:41 2015

@author: fergal

"""

import numpy as np

def computeArcLength(cent_colrow, flag):
    """Compute arclength along the main eigenvector for each centroid point

    Arclength is a number which describes how much the spacecraft has
    rolled from some nominal position. The nominal position is the
    mean roll angle.

    Inspired by similar approach by Vanderburg  & Johnson (2014)

    This is done slightly differently to V&J. I don't bother computing
    arclength along the fit, I just say how far along the eigenvector
    each centroid point lies. This is noisier than their approach, but
    simpler to implement.

    Inputs:
    -------------
    cent_colrow
        (2d numpy array) Array of centroid positions. Zeroth col is column
        position, and 1st col is row position. Bad values should be stripped
        out of this array before calling compute arclenght.

    flag
        (1d np array) Array of flags showing which cadences are flagged
        as bad.

    Returns:
    ----------
    A numpy array of shape cent. The zeroth column is the value of arclength
    for each input row.
    """


    cent = cent_colrow.copy()
    #Measure eigenvectors of mean-removed centroid distribution
    c0 = nanmean(cent[~flag,0])
    r0 = nanmean(cent[~flag,1])
    cent[:,0] -= c0
    cent[:,1] -= r0
    cov = np.cov(cent[~flag, :].transpose())

    [eVal, eVec] = np.linalg.eigh(cov)

    rot = cent * 0
    rot[:,0] = cent[:,0]* eVec[0,1] + cent[:,1]*eVec[1,1]
    rot[:,1] = cent[:,0]* eVec[0,0] + cent[:,1]*eVec[0,1]
    return rot



def nanmean(a):
    """Not available in earlier versions of numpy so I implement it here"""
    idx = np.isfinite(a)
    return np.mean(a[idx])
