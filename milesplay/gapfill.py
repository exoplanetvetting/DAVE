# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 12:55:18 2016

@author: fmullall
"""

__version__ = "$Id$"
__URL__ = "$URL$"


import numpy as np
import plateau
import lsf


def gapFill(yInput, shortGapSize=4):
    """Replace Nans in a lightcurve with interpolated values

    Algorithm distinguishes between short gaps and long gaps.
    Short gaps are replaced with the average of nearby values.
    Long gaps are replaced with the fit of a cubic polynomial
    to the data near the gap.

    Gapped data is assumed to have the value of nan.

    Inputs:
    --------------
    yInput:
        (1d np array) Input flux to be filled. Assumed to be
        equally spaced.

    Optional Inputs:
    ------------------
    shortGapSize
        (int) Maximum size of a short gap. Must be at least 2.


    Returns:
    ------------
    (gapFilledFlux, indexOfGappedCadences)
    A tuple of 2 1d arrays. The first array is the gap filled
    flux, the second is a boolean array of the modified elements
    in the first array.

    """

    assert shortGapSize > 2

    y = yInput.copy()
    badIdx = np.isnan(y)
    assert not np.all(badIdx)  #Need some good data

    groups = plateau.plateau(badIdx, .5)
    for g in groups:
        lwr, upr = g
        rng = slice(lwr, upr)
        if upr < lwr + shortGapSize:
            #Short gap fill
            a1 = max(0, lwr-shortGapSize)
            a2 = min(len(y), upr+shortGapSize+1)
            y[rng] = np.nanmean( y[a1:a2])

        else:
            nPointsForFit = upr-lwr
            #long gap fill
            xSnip = np.arange(lwr-nPointsForFit, upr+nPointsForFit+1)
            xSnip -= xSnip[0]

            ySnip = y[lwr-nPointsForFit:upr+nPointsForFit+1]
            weight = np.ones_like(ySnip)

            idx = np.isnan(ySnip)
            ySnip[idx] = 6700
            weight[idx] = 1e10   #Don't fit the gap

            fObj = lsf.Lsf(xSnip, ySnip, weight, order=4)
            x = np.arange(nPointsForFit, nPointsForFit+(upr-lwr))
            y[lwr:upr] = fObj.getBestFitModel(x=x)


    return y, badIdx


def medianSmooth(y, nPointsForSmooth=48):
    assert np.all( np.isfinite(y) )

    nPoints = nPointsForSmooth

    size = len(y)
    yNew = np.zeros_like(y)
    for i in range(size):
        #This two step ensures that lwr and upr lie in the range [0,size)
        lwr = max(i-nPoints, 0)
        upr = min(lwr + 2*nPoints, size)
        lwr = upr- 2*nPoints


        offset = np.median(y[lwr:upr])
        yNew[i] = y[i] - offset

    return yNew