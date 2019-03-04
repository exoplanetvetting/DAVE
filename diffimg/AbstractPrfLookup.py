"""
Created on Sun Dec  2 14:12:41 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

import numpy as np


class AbstractPrfLookup(object):
    """Store and lookup a previously computed PRF function

    This abstract class is created in the hope that much of the functionality can
    be reused for TESS.

    To get the recorded prf, use
    getPrfForBbbox(), although this function doesn't work in the base class. See docs
    in that method for more details

    Other functions are in the works to map that PRF onto
    a given mask.

    Todo
    --------
    This class is ripe of optimization with numba
    """

    def __init__(self, path):
        self.path = path
        self.cache = dict()

        self.gridSize = None


    def abstractGetPrfForBbox(self, col, row, bboxIn, getPrfFunc, *args):
        """Get the prf for an image described by a bounding box

        This function requires as input a function to look up the PRF for a given col row
        (getPrfFunc). This function is not implemented in the base class, as it will be
        specific to each mission. The design intent is that you override getPrfForBbox()
        in the daughter class where
        you define the lookup function, then calls the parent class method.
        See KeplerPrf() for an example

        Input:
        -----------
        col, row
            (floats) Centroid for which to generate prf
        bboxIn
            (int array). Size of image to return. bbox
            is a 4 elt array of [col0, col1, row0, row1]
        getPrfFunc
            (function) Function that returns a PRF object
            This function must have the signature
            ``(np 2d array = getPrfFunc(col, row, *args)``

        Optional Inputs
        -----------------
        Any optional inputs get passed to getPrfFunc


        Returns:
        ----------
        A 2d numpy array of the computed prf at the
        given position.

        Notes:
        ------------
        If bbox is larger than the prf postage stamp returned,
        the missing values will be filled with zeros. If
        the bbox is smaller than the postage stamp, only the
        requestd pixels will be returned
        """

        bbox = np.array(bboxIn).astype(int)
        nColOut = bbox[1] - bbox[0]
        nRowOut = bbox[3] - bbox[2]
        imgOut = np.zeros( (nRowOut, nColOut) )

        #Location of origin of bbox relative to col,row.
        #This is usually zero, but need not be.
        colOffsetOut = (bbox[0] - np.floor(col)).astype(np.int)
        rowOffsetOut = (bbox[2] - np.floor(row)).astype(np.int)

        interpPrf = getPrfFunc(col, row, *args)
        nRowPrf, nColPrf = interpPrf.shape
        colOffsetPrf = -np.floor(nColPrf/2.).astype(np.int)
        rowOffsetPrf = -np.floor(nRowPrf/2.).astype(np.int)

        di = colOffsetPrf - colOffsetOut
        i0 = max(0, -di)
        i1 = min(nColOut-di , nColPrf)
        if i1 <= i0:
            raise ValueError("Central pixel column not in bounding box")
        i = np.arange(i0, i1)
        assert(np.min(i) >= 0)

        dj = rowOffsetPrf - rowOffsetOut
        j0 = max(0, -dj)
        j1 = min(nRowOut-dj, nRowPrf)
        if j1 <= j0:
            raise ValueError("Central pixel row not in bounding box")
        j = np.arange(j0, j1)
        assert(np.min(j) >= 0)

        #@TODO: figure out how to do this in one step
        for r in j:
            imgOut[r+dj, i+di] = interpPrf[r, i]

        return imgOut
