# -*- coding: utf-8 -*-

import numpy as np

__version__ = "$Id: plateau.py 2245 2015-12-15 22:53:07Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/plateau.py $"

"""
Find continuous regions in an array above some value

This code gets called in lots of places, but doesn't really belong in
any of my other modules
"""

def plateau(array, threshold):
    """Find plateaus in an array, i.e continuous regions that exceed threshold

    Given an array of numbers, return a 2d array such that
    out[:,0] marks the indices where the array crosses threshold from
    below, and out[:,1] marks the next time the array crosses that
    same threshold from below.

    Inputs:
    array       (1d numpy array)
    threshold   (float or array) If threshold is a single number, any point
                above that value is above threshold. If it's an array,
                it must have the same length as the first argument, and
                an array[i] > threshold[i] to be included as a plateau

    Returns:
    Numpy 2d array with 2 columns.


    Notes:
    To find the length of the plateaus, use
    out[:,1] - out[:,0]

    To find the length of the largest plateau, use
    np.max(out[:,1] - out[:,0])

    The algorithm fails if a value is exactly equal to the threshold.
    To guard against this, we add a very small amount to threshold
    to ensure floating point arithmetic prevents two numbers being
    exactly equal.
    """


    arr  = array.astype(np.float32)
    arr = arr - threshold + 1e-12
    arrPlus = np.roll(arr, 1)

    #Location of changes from -ve to +ve (or vice versa)
    #Last point is bogus , so we calcualte it by hand
    sgnChange = arr*arrPlus

    #Roll around can't compute sign change for zeroth elt.
    sgnChange[0] = +1
    if arr[0] > 0:
        sgnChange[0] = -1

    loc = np.where(sgnChange < 0)[0]

    if np.fmod( len(loc), 2) != 0:
        loc.resize( (len(loc)+1))
        loc[-1] = len(arr)

    if len(loc) == 0:
        return []
    return loc.reshape( (-1,2))
