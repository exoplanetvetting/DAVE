
#import sincfilter
import dave.misc.sincfilter as sincfilter
import numpy as np

def indexOfOutliers(a, threshold_sigma=4, normFreq=.45):
    """Identify indices of outlier points in a timeseries

    Inputs:
    a   (1d np array)   Input data. a is assumed to represent
                        equally spaced data

    Optional Inputs:
    threshold_sigma     How many sigma away from mean must a point
                        be to be considered an outlier
    normFreq            Input arg to sincfilter, defines how agressively
                        to smooth.

    Returns:
    1d array of booleans of length a. Value set to true if
    element is an outlier

    Description
    Data is heavily high pass filtered, with the assumption that
    residuals are mostly gaussian with a few outlying points.
    These points are id'd, and their indices returns
    """

    numPointsInFilter = np.floor(len(a)/20.)
    if numPointsInFilter % 2 == 1:
        numPointsInFilter += 1

    filt = sincfilter.highPass(np.asarray(a), normFreq, numPointsInFilter)

    rms = np.std(filt)
    idx = np.abs(filt) > threshold_sigma*rms
    return idx


