from __future__ import division, print_function
import numpy as np

"""
@TODO:
We should have one, and only one, median detrend function
Same for plateau function

Add docstrings

Should this belong to a general "lightcurve maniuplation module"?

"""




def medianDetrend(flux, binWidth):
    halfNumPoints = binWidth // 2
    medians = np.array([])
    for i in range(len(flux)):
        if i < halfNumPoints:
            medians = np.r_[medians, np.atleast_1d(np.median(flux[:i+halfNumPoints+1]))]
        elif i > len(flux) - halfNumPoints - 1:
            medians = np.r_[medians, np.atleast_1d(np.median(flux[i-halfNumPoints:]))]
        else:
            medians = np.r_[medians, np.atleast_1d(np.median(flux[i-halfNumPoints : i+halfNumPoints+1]))]
    return medians



def outlierRemoval(time, flux, precision_days=0.0205, threshold_sigma=4):
    """
    Remove single point outliers.

    Preserves consecutive outliers, and those that are evenly spaced in
    time. This protects short duration transits.

    Inputs:
    ------------
    time, flux
        (1d numpy array) Input data. Flux should have mean (or median) value
        of zero
    precision_days
        (float) Points that are evenly spaced to within this precision
        are considered periodic, and not marked as outliers

    threshold_sigma
        (float) Points more than this many sigma from zero are considered
        potential outliers.

    Returns:
    ------------
    An array of indices indicating which points are single point
    outliers. The length of this array is equal to the number of outlier
    points flagged

    TODO:
    ---------
    For consistency with rest of code, this could return an array of booleans
    """

    #Remove as much signal as possible from the data
    fluxDetrended = medianDetrend(flux, 3)

    rms = np.std(fluxDetrended)
    out1 = plateau(fluxDetrended, threshold_sigma * rms)
    out2 = plateau(-fluxDetrended, threshold_sigma * rms )
    if out1 == [] and out2 == []:
        singleOutlierIndices = []
    else:
        outliers = np.append(out1, out2).reshape(-1,2)
        # Only want groups of one outlier, since > 1 may be transit points
        singleOutlierIndices = np.sort(outliers[(outliers[:,1] - outliers[:,0] == 1)][:,0])
    # Check periodicity of outliers, with PRECISION of 0.0205 days
    # 0.0205 days = 29.52 minutes = ~length of long cadence
    outlierTimes = time[singleOutlierIndices]
    diffs = [outlierTimes[i+1] - outlierTimes[i] for i in range(0, len(outlierTimes)-1)]
    diffs = [round(d, 5) for d in diffs]
    if len(singleOutlierIndices) >= 4:
        if len(set(diffs)) == len(diffs):
            possibleTimes = np.array([])
        else:
            period = max(set(diffs), key = diffs.count) # period = most common difference
            epoch = outlierTimes[diffs.index(period)]
            possibleTimes = np.arange(epoch, outlierTimes[-1] + 0.5*period, period)
        notOutliers = []
        for i in range(len(outlierTimes)):
            if np.any((abs(possibleTimes - outlierTimes[i]) < precision_days)):
                notOutliers.append(i)
        singleOutlierIndices = np.delete(singleOutlierIndices, notOutliers)
    elif len(singleOutlierIndices) == 3:
        if abs(diffs[0] - diffs[1]) < precision_days:
            singleOutlierIndices = []
    return singleOutlierIndices


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
    exactly equal."""
    arr = np.array(array, dtype=np.float)
    arr = arr - threshold + 1e-12
    arrPlus = np.roll(arr, 1)
    # Location of changes from -ve to +ve (or vice versa)
    # Last point is bogus , so we calculate it by hand
    sgnChange = arr*arrPlus
    # Roll around can't compute sign change for zeroth elt.
    sgnChange[0] = +1
    if arr[0] > 0:
        sgnChange[0] = -1
    loc = np.where(sgnChange < 0)[0]
    if np.fmod(len(loc), 2) != 0:
        loc = np.resize(loc,(len(loc)+1))
        loc[-1] = len(arr)
    return loc
