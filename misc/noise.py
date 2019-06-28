# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 13:43:01 2016

@author: fergal

A series of metrics to quantify the noise in a lightcurve:

Includes:
x sgCdpp
x Marshall's noise estimate
o An FT based estimate of 6 hour artifact strength.
o A per thruster firing estimate of 6 hour artifact strength.

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



from scipy.signal import savgol_filter
import matplotlib.pyplot as mp
import numpy as np
#import fft
import dave.misc.fft as fft

keplerLongCadence_s = 1765.4679
keplerLongCadence_days = keplerLongCadence_s / float(86400)


def computeRollTweakAmplitude(y, nHarmonics = 3, tweakPeriod_days = .25, \
    expTime_days=None, plot=False):
    """Compute strength of roll tweak artifact in K2 data with an FT approach.

    Compute FT of lightcurve

    Optional Inputs:
    -----------------
    plot
        Show a diagnostic plot

    Returns:
    --------
    float indicating strength of correction. A value of 1 means the
    amplitude of the tweak is approx equal to the strength of all other
    signals in the FT.
    """
    if expTime_days is None:
        expTime_days = keplerLongCadence_days

    #computes FT with frequencies in cycles per days
    ft = fft.computeFft(y, expTime_days)

      #Thruster firings every 6 hours
    artifactFreq_cd = 1/tweakPeriod_days  #cycles per day

    if plot:
        mp.clf()
        mp.plot(ft[:,0], 1e6*ft[:,1], 'b-')

    metric = 0
    nPtsForMed = 50
    for i in range(1, nHarmonics+1):
        wh = np.argmin( np.fabs(ft[:,0] - i*artifactFreq_cd))

        med = np.median(ft[wh-nPtsForMed:wh+nPtsForMed, 1])
        metric += ft[wh, 1] / med

        if plot:
            mp.axvline(i*artifactFreq_cd, color='m')
    return metric / float(nHarmonics)


def computeSgCdpp_ppm(y, transitDuration_cadences=13, plot=False):
    """Estimates 6hr CDPP using Jeff Van Cleve's Savitzy-Golay technique

    An interesting estimate of the noise in a lightcurve is the scatter
    after all long term trends have been removed. This is the kernel of
    the idea behind the Combined Differential Photometric Precision (CDPP)
    metric used in classic Kepler. Jeff van Cleve devised a much simpler
    algorithm for computing CDPP using a Savitzy-Golay detrending, which
    he called Savitzy-Golay CDPP, or SG-CDPP. We implement his algorithm
    here.


    Inputs:
    ----------
    y
        (1d numpy array) normalised flux to calculate noise from. Flux
        should have a mean of zero and be in units of fractional amplitude.
        Note: Bad data in input will skew result. Some filtering of
        outliers is performed, but Nan's or Infs will not be caught.

    Optional Inputs:
    -----------------
    transitDuration_cadences
        (int) Adjust the assumed transit width, in cadences. Default is
        13, which corresponds to a 6.5 hour transit in K2
    plot
        Show a diagnostic plot

    Returns:
    ------------
    Estimated noise in parts per million.

    Notes:
    -------------
    Taken from
    svn+ssh://murzim/repo/so/trunk/Develop/jvc/common/compute_SG_noise.m
    by Jeff van Cleve
    """

    #These 3 values were chosen for the original algorithm, and we don't
    #change them here.
    window = 101
    polyorder=2
    noiseNorm = 1.40
    #Name change for consistency with original algorithm
    cadencesPerTransit = int(transitDuration_cadences)

    if cadencesPerTransit < 4:
        raise ValueError("Cadences per transit must be >= 4")

    if len(y) < window:
        raise ValueError("Can't compute CDPP for timeseries with fewer points than defined window (%i points)" %(window))

    trend = savgol_filter(y, window_length=window, polyorder=polyorder)
    detrend = y-trend

    filtered = np.ones(cadencesPerTransit)/float(cadencesPerTransit)
    smoothed = np.convolve(detrend, filtered, mode='same')

    if plot:
        mp.clf()
        mp.plot(y, 'ko')
        mp.plot(trend, 'r-')
        mp.plot(smoothed, 'g.')

    sgCdpp_ppm = noiseNorm*robustStd(smoothed, 1)*1e6
    return sgCdpp_ppm


def estimateScatterWithMarshallMethod(flux, plot=False):
    """Estimate the typical scatter in a lightcurve.

    Uses the same method as Marshall (Mullally et al 2016 submitted)

    Inputs:
    ----------
    flux
        (np 1d array). Flux to measure scatter of. Need not have
        zero mean.


    Optional Inputs:
    -----------------
    plot
        Show a diagnostic plot

    Returns:
    ------------
    (float) scatter of data in the same units as in the input ``flux``


    Notes:
    ----------
    Algorithm is reasonably sensitive to outliers. For best results
    uses outlier rejection on your lightcurve before computing scatter.
    Nan's and infs in lightcurve will propegate to the return value.
    """

    diff= np.diff(flux)

    #Remove egregious outliers. Shouldn't make much difference
    idx = sigmaClip(diff, 5)
    diff = diff[~idx]

    mean = np.mean(diff)
    mad = np.median(np.fabs(diff-mean))
    std = 1.4826*mad

    if plot:
        mp.clf()
        mp.plot(flux, 'ko')
        mp.plot(diff, 'r.')
        mp.figure(2)
        mp.clf()
        bins = np.linspace(-3000, 3000, 61)
        mp.hist(1e6*diff, bins=bins, ec="none")

        mp.xlim(-3000, 3000)
        mp.axvline(-1e6*float(std/np.sqrt(2)), color='r')
        mp.axvline(1e6*float(std/np.sqrt(2)), color='r')

    #std is the rms of the diff. std on single point
    #is 1/sqrt(2) of that value,
    return float(std/np.sqrt(2))



def singlePointDifferenceSigmaClip(a, nSigma=4, maxIter=1e4, initialClip=None):
    """Iteratively find and remove outliers in first derivative

    If a dataset can be modeled as a constant offset + noise + outliers,
    those outliers can be found and rejected with a sigma-clipping approach.
    
    If the data contains some time-varying signal, this signal must be removed
    before applying a sigma clip. This function removes the signal by applying
    a single point difference.
    
    The function computes a[i+1] - a[i], and sigma clips the result. Slowly
    varying trends will have single point differences that are dominated by noise,
    but outliers have strong first derivatives and will show up strongly in this
    metric.
    
    Inputs:
    ----------
    y
        (1d numpy array) Array to be cleaned
    nSigma
        (float) Threshold to cut at. 5 is typically a good value for
        most arrays found in practice.

    Optional Inputs:
    -------------------
    maxIter
        (int) Maximum number of iterations

    initialClip
        (1d boolean array) If an element of initialClip is set to True,
        that value is treated as a bad value in the first iteration, and
        not included in the computation of the mean and std.

    Returns:
    ------------
    1d numpy array. Where set to True, the corresponding element of y
    is an outlier.
    """
   
    #Scatter in single point difference is root 2 time larger
    #than in initial lightcurve
    threshold = nSigma/np.sqrt(2)
    
    diff1 = np.roll(a, -1) - a
    diff1[-1] = 0  #Don't trust the last value because a[-1] not necessarily equal to a
    idx1 = sigmaClip(diff1, nSigma, maxIter, initialClip)
    
    diff2 = np.roll(a, 1) - a
    diff2[0] = 0  
    idx2 = sigmaClip(diff2, nSigma, maxIter, initialClip)

    flags = idx1 & idx2
    
    #This bit of magic ensures only single point outliers are marked, 
    #not strong trends in the data. It insists that the previous point 
    #in difference time series is an outlier in the opposite direction, otherwise
    #the point is considered unflagged. This prevents marking transits as bad data.
    outlierIdx = flags
    outlierIdx &= np.roll(idx1, 1)
    outlierIdx &= (np.roll(diff1, 1) * diff1 < 0)
    return outlierIdx


def sigmaClip(y, nSigma, maxIter=1e4, initialClip=None):
    """Iteratively find and remove outliers

    Find outliers by identifiny all points more than **nSigma** from
    the mean value. The recalculate the mean and std and repeat until
    no more outliers found.

    Inputs:
    ----------
    y
        (1d numpy array) Array to be cleaned
    nSigma
        (float) Threshold to cut at. 5 is typically a good value for
        most arrays found in practice.

    Optional Inputs:
    -------------------
    maxIter
        (int) Maximum number of iterations

    initialClip
        (1d boolean array) If an element of initialClip is set to True,
        that value is treated as a bad value in the first iteration, and
        not included in the computation of the mean and std.

    Returns:
    ------------
    1d numpy array. Where set to True, the corresponding element of y
    is an outlier.
    """
    #import matplotlib.pyplot as mp
    idx = initialClip
    if initialClip is None:
        idx = np.zeros( len(y), dtype=bool)

    assert(len(idx) == len(y))

    #x = np.arange(len(y))
    #mp.plot(x, y, 'k.')

    oldNumClipped = np.sum(idx)
    for i in range(int(maxIter)):
        mean = np.nanmean(y[~idx])
        std = np.nanstd(y[~idx])

        newIdx = np.fabs(y-mean) > nSigma*std
        newIdx = np.logical_or(idx, newIdx)
        newNumClipped = np.sum(newIdx)

        #print "Iter %i: %i (%i) clipped points " \
            #%(i, newNumClipped, oldNumClipped)

        if newNumClipped == oldNumClipped:
            return newIdx

        oldNumClipped = newNumClipped
        idx = newIdx
        i+=1
    return idx


def robustMean(y, percent):
    """Compute the mean of the percent.. 100-percent percentile points

    A fast, and typically good enough estimate of the mean in the presence
    of outliers.
    """

    ySorted = np.sort( y[np.isfinite(y)] )
    num = len(ySorted)
    lwr = int( percent/100. * num)
    upr = int( (100-percent)/100. * num)

    return np.mean( ySorted[lwr:upr])


def robustStd(y, percent):
    """Compute a robust standard deviation with JVC's technique

    A fast, and typically good enough estimate of the mean in the presence
    of outliers.Cuts out 1st and 99th percentile values and computes std
    of the rest. Used by computeSgCdpp() to match the behaviour of
    Jeff van Cleve's original algorithm

    Taken from
    svn+ssh://murzim/repo/so/trunk/Develop/jvc/common/robust_std.m
    """

    ySorted = np.sort( y[np.isfinite(y)] )
    num = len(ySorted)
    lwr = int( percent/100. * num)
    upr = int( (100-percent)/100. * num)

    return np.std( ySorted[lwr:upr])

