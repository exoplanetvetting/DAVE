# -*- coding: utf-8 -*-

from __future__ import division
"""

A pure Python implementation of the BLS search. Not as fast as the
Fortran version, but no dependencies or compile problems.

Usage:
------------
Create a BlsSearch object
> pRange = [1,30]
> durList_days = np.array([1,2,4,8])/24.
> blsObj = BlsSearch(time, flux, pRange_days, durList_days)

The period range and durations are in the same units as the time array
(for K2, this is days)

Next get the properties of the strongest signal:
period, epoch, depth, duration = blsObj.getEvent()


Alternatively, you can call fbls() for an Object-Orientation-free approach.

Todo:
------------
o Add option of a condensed output to save memory, that gives an array
  of period, max(bls amplitude), epoch at max, duration at max
o Optimise min points per transit duration


Notes:
-------------
To profile code:
Add @profile decorator to functions you want to profile. From
command line run
> date; ~/bin/python ~/.local/lib/python2.7/site-packages/kernprof.py -l -v test_fbls.py > fbls.prof; date
Then edit the file fbls.prof
"""

import matplotlib.pyplot as mp
import numpy as np

import dave.fileio.kplrfits as kplrfits
import datetime




class BlsSearch(object):
    """A thin wrapper class around fbls and associated functions.
    Basically, this class stores the multitude of variables you need
    to carry around. If you don't like object orientation, call fbls()
    yourself directly"""

    def __init__(self, time, flux, periodRange, durations_days,
                 periodOverResolution=2, binsPerTransit=10, debug=False):
        """

        Initialise variables and run the bls search

        Inputs:
        ------------
        time, flux
            (1d np arrays) Times and fluxes of data. flux should have a mean
            of approximately zero.

        periodRange
            (2 element list) Minimum and maximum periods to search. Units
            the same as the time array

        durations
            (list) List of transit durations to search for. Units the same
            as the time array.


        Optional Inputs:
        -----------------
        debug
            (bool) If True, the class is initialised, but the search is
            not run. This allows you to inspect the class before it is run.
            However, the other methods do NOT check if the bls ran, so the
            they will raise exceptions if called.

        Notes:
        ---------
        This class intended for use with K2, so durations should be given
        in days. The class warns you if you seem to be input durations in hours

        """

        self.periodOverResolution = periodOverResolution
        #@TODO, just call this binsPerTransit
        self.minBinsPerTransit = binsPerTransit
        self.time = time
        self.flux = flux
        self.periodRange = periodRange

        if len(time) < 2:
            raise ValueError("Time array must have at least two elements")


        lwr, upr = periodRange
        self.periods = computePeriodList(lwr, upr, time, min(durations_days),
                                         self.periodOverResolution)

        #Allow durations to be a single value
        if not hasattr(durations_days, "__len__"):
            durations_days = [durations_days]
        assert(len(durations_days) > 0)
        self.durations = np.array(durations_days)

        if np.min(self.durations) > 1:
            print("WARN: All durations longer than 1 day. Did you input hours?")

        if lwr < np.max(self.durations):
            raise ValueError("Shortest period less than longest duration")
        if not debug:
            self.run()


    def run(self):
        """Private method. Run the bls"""
        t0 = datetime.datetime.now()
        self.blsArray = computeBlsForManyPeriods(self.time, self.flux, \
            self.durations, self.periods, self.minBinsPerTransit)

        self.elapsedTime = datetime.datetime.now() - t0

    def getEvent(self):
        """Return properites of best peak in BLS spectrum.

        This is a wrapper for  :func:`~fbls.getParamsOfIndex`
        which returns the "best", not necessarily the strongest peak
        in the BLS spectrum.

        Returns:
        A for element tuple of
        [period, epoch, depth, duration]

        """
        index = self.getIndexOfBestPeak()
        period, epoch, dur = getParamsOfIndex(self.blsArray, index, \
            self.durations, self.periods)

        epoch += np.min(self.time)

        depth = self.getTransitDepth(index)
        return period, epoch, depth, dur


    def getIndexOfBestPeak(self):
        """Thin wrapper for :func:`~findBestPeak`"""
        return findBestPeak(self.blsArray)


    def getTransitDepth(self, index):
        """Estimate transit depth based on strength of BLS spectrum
        at a given index

        Transit depth is bls strength/ sqrt(trial transit width in cadences)

        For a given location in the BLS array, this function computes the
        transit width used to compute the BLS strength at that array location,
        and then uses that to compute the depth.

        Inputs:
        ----------
        index
            (3-tuple) Location in BLS spectrum to compute depth for.
            Elements are the index of (period, duration, epoch)

        Returns:
        -----------
        Depth in fractional amplitude


        Notes:
        ------------
        There is a potential for maintenance problems here. This method assumes
        it that :func:`computeBlsForOnePeriod` is using the same algorithm.
        If that function changes, this code will be wrong"""
        blsVal = self.blsArray[index]
        minDur = np.min(self.durations)

        period = self.periods[index[0]]
        duration = self.durations[index[1]]

        nBins = self.minBinsPerTransit * period/minDur
        transitWidth = int(np.round(duration/period*nBins))

        return blsVal/np.sqrt(transitWidth)


    def compute1dBls(self):
        """Compute 1d bls spectrum

        Returns:
        ------------
        A 2d array with columns of period and strongest amplitude at
        that period.
        """
        bb = self.blsArray.max(axis=1)
        bbb = bb.max(axis=1)

        out = np.empty( (len(bbb), 2))
        out[:,0] = self.periods
        out[:,1] = bbb

        return out

    def getRuntime_sec(self):
        """Return the time taken to run the BLS"""
        return self.elapsedTime.total_seconds()


    def plot(self):
        """Plot the 2d BLS spectrum, showing the strongest peak
        across all duratons for each period/epoch pair
        """

        bb = self.blsArray.max(axis=1)
#        bb = self.blsArray[:,0,:]
        mp.imshow(bb, cmap=mp.cm.YlGnBu_r, interpolation="nearest",\
            origin="bottom", aspect="auto")
        mp.colorbar()


def fBls(time, flux, periodRange, durations):
    """
    Compute a bls spectrum of a lightcurve

    Inputs:
    ----------
    time, flux
        (1d np arrays) Times and fluxes of data. flux should have a mean
        of approximately zero.

    periodRange
        (2 element list) Minimum and maximum periods to search. Units
        the same as the time array

    durations
        (list) List of transit durations to search for. Note this
        is not necessarily in hours, but in the same units as the time array
        (which is commonly days).

    Returns:
    ----------
    Not yet decided

    Notes:
    ----------
    See bls.tex for notes on computationaly complexity.

    """

    assert(len(time) == len(flux))
    assert(np.all(np.isfinite(time)))
    assert(np.all(np.isfinite(flux)))
    assert(len(periodRange) == 2)

    if not hasattr(durations, "__len__"):
        durations = [durations]
    assert(len(durations) > 0)

    lwr, upr = periodRange
    periods = computePeriodList(lwr, upr, time, min(durations))

    bls = computeBlsForManyPeriods(time, flux, durations, periods)
    return bls, periods


def computePeriodList(Plwr, Pupr, time, duration, overres=2):
    """
    Compute the list of periods needed to comphrensively search a lightcurve
    with BLS.

    Inputs:
    -----------
    Plwr, Pupr
        (floats) Minimum and maximum periods to search

    time
        (list or array) times at which data is sampled. May be an array
        of all time values, or just a list of the first and last points.
    duration
        (float) Assumed transit duration in units of **time**
    overres
        How many times the minimum sampling should be performed

    Returns:
    -----------
    An array of periods to search at.

    Note:
    ---------
    See bls.tex in this directory for a description of this algorithm
    """


    assert(Pupr > Plwr)
    timespan = np.nanmax(time) - np.nanmin(time)
    q = duration/timespan/float(overres)

    nStep = int(np.ceil(np.log(Pupr/Plwr) / np.log(1+q)))+2
    n = np.arange(nStep)

    return Plwr*(1+q)**n


#@profile
def computeBlsForManyPeriods(t, y, durationList, periodList, minNumBinsPerTransit=10):
    """Worker function for fbls

    See fbls() for a description of inputs. Values for durationList
    and periodList should be given in the same units as t (usually days)

    Todo:
    Allow passing in weights on each data point to better account for
    varying noise.

    Do some experiments to fnd the best value for minBinsPerTransit.
    Kovacs suggests 15, I find 10 works well, I haven't tried other numbers
    """
    assert np.all(np.isfinite(t)), "Non-finite time value found"
    assert np.all(np.isfinite(y)), "Non-finite flux value found"

    durationList = np.array(durationList)
    periodList = np.array(periodList)
    if np.any(durationList <= 0):
        raise ValueError("Durations must be > 0")

    if np.any(periodList <= 0):
        raise ValueError("Periods must be > 0")

    maxNumBins = minNumBinsPerTransit * max(periodList) / min(durationList)
    maxNumBins = int(np.ceil(maxNumBins))

    out = np.zeros((len(periodList), len(durationList), int(maxNumBins)))
    dt = t - np.min(t)

    for i,period in enumerate(periodList):
        #Technically, the lightcurve should be refolded for each duration
        #But because the algorithm is only weakly sensitive to the bin width
        #above a certain value, we dramatically speed things up by
        #folding once per period
        nBins = int(minNumBinsPerTransit * period / min(durationList))

        if False:
            #Useful for debugging. This folder is much slower,
            #but well tested.
            expTime = np.median(np.diff(t))
            binned = kplrfits.foldAndBinData(dt, y, period, 0, expTime, nBins)
            binned = binned[:,1]
        else:
            binned, counts = fastFoldAndBin(dt, y, period, nBins)
            binned /= counts

#        import ipdb; ipdb.set_trace()
#        mp.clf()
#        mp.plot(binned, 'bo-')

        for j,duration in enumerate(durationList):
            #This line is wrong!
            #should be? overRed*duration/min(duration)
            #transitWidth = int(np.round(duration/period * nBins))
            transitWidth = minNumBinsPerTransit*duration/np.min(duration)

            if transitWidth > 0 and transitWidth < len(binned):
                bls = computeBlsForOnePeriod(binned, transitWidth)
#                mp.plot(bls, 'ro-')
#                mp.pause(1)
                out[i,j, :nBins] = bls
    return out

"""
def computeBlsForOnePeriod(y, transitWidth):

    if transitWidth <= 0:
        raise ValueError("Transit width must be at least 1 cadence")

    if transitWidth >= len(y):
        raise ValueError("Transit width must be shorter than length of binned lightcurve")

    template = np.ones(transitWidth)
    bls = -np.convolve(y, template, mode='same')
    assert(len(bls) == len(y))

    return bls
"""

def computeBlsForOnePeriod(y, transitWidth):

    if transitWidth <= 0:
        raise ValueError("Transit width must be at least 1 cadence")

    if transitWidth >= len(y):
        raise ValueError("Transit width must be shorter than length of binned lightcurve")

    #If y has weights, s[], r = \Sigma s_i
    r = float(transitWidth)

    #if y has weights, s is the convolution of y*s by signal
    signal = np.ones(transitWidth)
    s = np.convolve(y, signal, mode='same')
    assert(len(s) == len(y))
    bls = -s/np.sqrt(r)
    assert(len(bls) == len(y))

    return bls



#def computeTransitWidth(duration, period, nBins):
#    return int(np.round(duration/period * nBins))

#def computeNumberOfInTransitBins(minNumBinsPerTransit, period, duration):
#    tmp = minNumBinsPerTransit * period / duration
#    return int(np.round(tmp))


#@profile
def fastFoldAndBin(t, y, period, nBins):
    """Fold and bin a (possibly unevenly spaced) lightcurve.

    Most of the time required to compute a BLS spectrum is spent
    folding and refolding the data at different periods. Time spent
    optimsing this action is time well spent. This is the fastest
    function I could think of. If you can think of away to do
    the same with one call to histogram you will halve the time spent
    on the entire BLS.

    Inputs:
    -----------
    t, y
        (1d arrays) Input time and flux. Time need not be evenly spaced.
    period
        (float) Period to fold at. Units the same as **t** (usually days)
    nBins
        (int) Number of bins to fold at.

    Returns:
    -----------
    1d array. The lightcurve is folded at the input period, then binned
    into nBin points. The epoch is chosen so that time==0 is mapped to
    phi==0.

    Note:
    ---------
    When using this function in a BLS code, you'd be wise to
    choose a time offset such that time==0 happens within the dataset. Otherwise
    the phase numbers might be unreliable..
    """
    phi = np.fmod(t, period)

    eps = 1e-10  #To prevent division by zero
    counts = np.histogram(phi, nBins, normed=False)[0] + eps
    binned = np.histogram(phi, nBins, weights=y, normed=False)[0]

    return binned, counts


##@profile
#def Old_computeBlsForOnePeriod(y, transitWidth):

    #if transitWidth <= 0:
        #raise ValueError("Transit width must be at least 1 cadence")

    #if transitWidth >= len(y):
        #raise ValueError("Transit width must be shorter than length of binned lightcurve")

    ##If y has weights, s[], r = \Sigma s_i
    #r = float(transitWidth)

    ##if y has weights, s is the convolution of y*s by signal
    #signal = np.ones(transitWidth)
    #s = np.convolve(y, signal, mode='same')
    #assert(len(s) == len(y))
    #bls = -s/np.sqrt(r)
    #assert(len(bls) == len(y))

    #return bls




#def getParamsOfDetection(blsArray, indices, duration_daysList, periodList):
#    assert(blsArray.ndim == 3)
#
#    if hasattr(indices, '__len__'):
#        assert(len(indices) in [1,3])
#    else:
#        indices = np.unravel_index(indices, blsArray.shape)
#
#    maxPeriod = np.max(periodList)
#
#    period = periodList[indices[0]]
#    duration_days = duration_daysList[indices[1]]
#    epoch = maxPeriod * (indices[2]/float(blsArray.shape[2]))
#
#    return period, epoch, duration_days


#@TODO, implement local copy of median detrend
import dave.fileio.kplrfits as kf
def findBestPeak(blsArray, nPointsForSmooth=100, offset=0):
    """Find the index of the best peak in the bls spectrum.

    The best peak is not usually the strongest one. The noise floor
    of a BLS increases at long period, so this must be accounted for
    before selecting the strongest amplitude.

    Inputs:
    -------------
    blsArray
        (3d numpy array)

    Optional Inputs:
    -----------------
    nPointsForSmooth
        (int) Tuning parameter. This is the number of points near a given
        period used to determine what the typical background BLS value is.

    offset
        (int) Debugging option. Change the chosen period by this number
        of indices.

    Returns:
    ------------
    A three element tuple of indices into blsArray giving the best
    period, duration and phase of transit.
    """
    bb = blsArray.max(axis=1)
    bbb = bb.max(axis=1)

    # bbb is 1d array of strongest peak at each period.
    # We need to remove background trend before we decide which peak
    # is the best choice
    filt = kf.medianSubtract1d(bbb, nPointsForSmooth)
    iPer = np.argmax(filt)
    iPer += offset
    subArray = blsArray[iPer,:,:]
    iDur, iPhase = np.unravel_index(np.argmax(subArray), subArray.shape)

    return iPer, iDur, iPhase


def getParamsOfIndex(blsArray, indices, duration_daysList, periodList):
    """Compute period, epoch, duration and depth corresponding to a given index

    Inputs:
    blsArray
        (3d np array)  Returned by fbls()
    indices
        (tuple/list of length 3) Index into blsArray to compute
    duraiton_daysList
        (list) List of durations corresponding to the last dimension of blsArray
    periodList
        (list/array) List of periods corresponding to zeroth dimension of
        blsArray

    Returns:
    -----------
    A tuple of period, epoch, and duration

    TODO:
    ----------
    Compute depth from blsArray
    """
    assert(blsArray.ndim == 3)

    if hasattr(indices, '__len__'):
        assert(len(indices) in [1,3])
    else:
        indices = np.unravel_index(indices, blsArray.shape)

    maxPeriod = np.max(periodList)

    period = periodList[indices[0]]
    duration_days = duration_daysList[indices[1]]
    epoch = maxPeriod * (indices[2]/float(blsArray.shape[2]))

    return period, epoch, duration_days
