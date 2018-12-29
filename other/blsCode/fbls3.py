

from __future__ import division

import matplotlib.pyplot as mp
import scipy.stats as spStats
import numpy as np



"""
%This is an experimental version and should be used with caution.

TODO
x Profile code
x Test duration spacing and epoch spacing
o Write wrapper class
o Good defalt values for qMin and qMax for duration/period spacing
o Update docs for computeFbls()
o Fix failing test by doing circular convolution of some kind


I use kernprof, part of the line_profiler package to profile this code.
"""

def computePeriodListConstDuration(minPeriod_days, maxPeriod_days,
                                  timespan_days, duration_days, attenuation=.9):
    """
    Compute a list of periods to search assuming a transit duration.

    If searching a fixed set of durations at every period, use the shortest
    trial duration to compute the period list


    Inputs:
    ----------
    timespan_days
        (float) Timespan of data == max(time) - min(time)
    duration_days
        (float) Trial transit duration to be used for search.

    Optional Inputs:
    -----------------
    attenuation
        Trial duration spacing will be chosen so BLS SNR is never smaller
        than this fraction of of the maximum strengths. A between 0 and 1


    Returns
    -----------
    A 1d numpy array of periods to be searched in units of days to
    ensure any signal is detected within a fraction of the max SNR
    (assuming duration and phase measured exactly)
    """

    errMsg =     "Attenuation outside legal range"
    assert (attenuation > 0) & (attenuation < 1), errMsg

    deltaP = 2*duration_days/timespan_days * (1-attenuation)
    nPeriod = np.log(maxPeriod_days/minPeriod_days)
    nPeriod /= np.log(1 + deltaP)

    periodList = minPeriod_days * (1 + deltaP)**np.arange(nPeriod)
    return periodList


def computePeriodListForAstrophysicalDurations(minPeriod_days, maxPeriod_days,
           timespan_days, radius_m, mass_kg, attenuation=.9):
    """

    Inputs:
    ----------
    """

    timespan_sec = 86400*timespan_days

    #Duration(P) = val*P^(1/3)
    newtonG = 6.67259e-11
    val = (2*np.pi)**(2/3.)
    val *= radius_m / np.pi
    val /= (newtonG * mass_kg)**(1/3.)


    minPeriod_sec = minPeriod_days * 86400
    maxPeriod_sec = maxPeriod_days * 86400


    periodList = [minPeriod_sec]
    while periodList[-1] < maxPeriod_sec:
        period = periodList[-1]
        deltaP = 2*val/timespan_sec * (1 - attenuation)
        deltaP *= period**(4/3.)
        newPeriod = period +  deltaP
        periodList.append(newPeriod)

        print newPeriod/86400., deltaP/86400.
#        import pdb; pdb.set_trace()
    periodList_days = np.array(periodList) / 86400.
    return periodList_days





class AstrophysicalDurationsSearch():
    def __init__(self, qMin, qMax, radius_m, mass_kg, attenuation=.9):
        """
        Pick a set of astrophysically reasonable trial transit durations
        based on properties of host star.

        Proceedure is to compute duration of planet in circular orbit of
        a given period assuming zero eccentricity, then compute trial
        durations in a reasonable range of that value.

        Requires knowledge of the properties of the host star.
        TODO: Default M/R^3 is for the sun, I should change this to that
        of an A star.
        TODO: Compute reasonable defaults for qMin and qMax

        If you don't know the mass and radius of the host, use the
        default values

        Inputs:
        ------------
        qMin, qMax
            (floats) Compute trial durations in this range. qMin < 1, qMax > 1.
            For example, if qMin,qMax = .5, 2, trial durations from half to
            twice "typical" durations are computed
        radius_m, mass_kg
            (floats) Radius and mass of host star. The defaults should
            work well for most main sequence stars
        attenuation
            Trial duration spacing will be chosen so BLS SNR is never smaller
            than this fraction of of the maximum strengths. A between 0 and 1


        Notes:
        ---------
        "Typical" duration calculation based on Eqn 2 of Ofir (2014).
        Spacing based on ms.tex
        """

        assert (qMin > 0) & (qMin < 1), "qMin outside legal range"
        assert (qMax > 1), "qMax outside legal range"

        assert (attenuation > 0) & (attenuation < 1), "Attenuation outside legal range"

        newtonG = 6.67259e-11
        self.qMin = qMin
        self.qMax = qMax
        self.attenuation = attenuation
        self.radius_m = radius_m
        self.mass_kg = mass_kg

        val = (2*np.pi)**(2/3.)
        val *= radius_m / np.pi
        val /= (newtonG * mass_kg)**(1/3.)
        self.const = val


    def __call__(self, period_days):
        period_sec = period_days * 86400
        dur_sec = self.const * period_sec **(1/3.)

        tauMin = self.qMin * dur_sec
        tauMax = self.qMax * dur_sec

        stepSize = 1/self.attenuation**4
        nTrials = np.log(tauMax/tauMin) / np.log(stepSize)
        nTrials = int(np.ceil(nTrials))
        dur_sec = tauMin * (stepSize ** np.arange(nTrials))
        dur_days = dur_sec / 86400.
        return dur_days


#@profile
def computefBls(time, flux, unc, periodList, phaseOverResolution,
                durationFunc0):
    """Compute ln false alarm prob for a given lightcurve

    Inputs:
    ------------
    time, flux
        (1d np array) Lightcurve of target
    unc
        (float) Average uncertainty of flux values.

    periodList:
        (list or array) List of periods to search. Can be created, e.g by
        computePeriodList()

    phaseOverResolution
        (int) When binning the folded lightcurve, how many bins per trial transit
        duration to use. Can be computed with computePhaseOverResolution(), but
        values in the range 5-10 are typically reasonable.

    durationFunc
        (function) A function, that when called, returns the list of durations
        to search. See notes below.


    Returns:
    ---------------
    A 2d array of length **periodList**. The columns are

    0:
        Input trial transit period
    1:
        Natural log of false alarm probability of most signficant peak
        in BLS transform at this period. The row with the lowest value
        in this column gives the period of the strongest transit signal
        in the data
    2:
        The BLS amplitude of the most significant peak at this period.
        Do not use this column to decide which peak is most significant
        across different periods
    3:
        The trial transit duration corresponding to the most significant
        peak at this period
    4:
        The phase of this peak. Phase is given in units of **time** since
        the first data point.

    Notes:
    ----------
    Ofir (2014) suggests only searching a range of durations for a given
    trial period. To accomodate this approach, the durations are not passed
    as a list, but as a function (technically a *callback* function). This
    function takes one argument, the period, and returns a range of durations
    to search at that period.

    np.array = def durationFunc(period).

    This function returns an array of durations in the same units as the
    input period.

    For real world situations, a function object (like. e.g ConstantDurationSearch
    below) are much more powerful

    TODO:
    -----------
    **unc** should optionally be an array, to allow for weighting points
    by their individual flux uncertainty.
    """

    dt = time - np.min(time)

    assert len(time) == len(flux)
    assert np.all(np.isfinite(time)), "Non-finite time value found"
    assert np.all(np.isfinite(flux)), "Non-finite flux value found"
    assert np.all(np.isfinite(flux)), "Non-finite uncertainty value found"

    #If a list of durations passed in instead of a function, convert
    #that list to a function that returns the list
    if hasattr(durationFunc0, "__len__"):
        durationFunc = lambda x: durationFunc0
    else:
        durationFunc = durationFunc0

    periodList = np.array(periodList, dtype=float)

    #Initalise the output array
    blsOut = np.empty( (len(periodList), 5) )
    blsOut[:,0] = periodList

    if np.any(periodList <= 0):
        raise ValueError("Periods must be > 0")

    for i, period in enumerate(periodList):
        durationList = durationFunc(period)

        if np.any(durationList <= 0):
            raise ValueError("Durations must be > 0")

        if max(durationList) > period:
            msg = "Longest duration (%.3f) is longer than trial period (%.3f)" \
                    %(max(durationList), period)
            raise ValueError(msg)


        #Binning is the most expensive part of BLS. To save time,
        #we bin once per period, not once per period/duration pair
        nBins = phaseOverResolution * period / min(durationList)
        nBins = int(np.ceil(nBins))
        binnedFlux, counts = fastFoldAndBin(dt, flux, period, nBins)

#        mp.clf()
#        mp.plot(binnedFlux, 'k.')
#        mp.title("%.3f" %(period))
#        mp.pause(.2)
#
        bls = np.zeros( (len(durationList), nBins) )
        lnFap = np.ones_like(bls)

        #Search each duration
        for j, duration in enumerate(durationList):
            transitWidth = phaseOverResolution*duration/np.min(durationList)

            if transitWidth > 0 and transitWidth < len(binnedFlux):
                bls[j] = computeBlsForOnePeriod(binnedFlux, counts,
                                             transitWidth, unc)

                #Convert BLS amplitude to lnFAP
                nDraws = period/duration
                lnFap[j] = computeGumbelLnFap(bls[j], nDraws)


#        mp.clf()
#        mp.imshow(bls, interpolation="nearest", origin="bottom", aspect="auto")
#        mp.colorbar()

        loc = np.argmax(bls)
        durLoc, phaseLoc = np.unravel_index(loc, bls.shape)
#        print loc, durLoc, phaseLoc
        blsOut[i, 1] = lnFap[durLoc, phaseLoc]     #The false alarm prob
        blsOut[i, 2] = bls[durLoc, phaseLoc]       #The BLS amplitude
        blsOut[i, 3] = durationList[durLoc]   #The duration
        blsOut[i, 4] = period * phaseLoc / float(nBins)  #The phase

#        1/0
    return blsOut





def computeBlsForOnePeriod(y, counts, transitWidth, uncPerUnbinnedPoint):
    """Compute the BLS amplitude across all phases in a folded,
    binned lightcurve

    The BLS amplitude is defined as the sum of flux of all points
    within a **transitWidth** of a given point, divded by the square root
    of the number of summed points .

    Inputs:
    -------------
    y
        (1d np array) Flux in each bin. Each bin is assumed to have an
        equal width in time units.
    counts
        (1d np array) Number of points in each bin
    transitWidth
        (int) Width of trial pulse in units of a single bin width.
        An error is raised if transitWidth is negative, or greater than
        the number of bins.
    uncPerUnbinnedPoint
        (float) A single value representing the typical scatter in the
        unbinned photometry points.
    Returns:
    ------------
    bls
        The BLS amplitude for each bin for the given input width.


    """
    if transitWidth <= 0:
        raise ValueError("Transit width must be at least 1 cadence")

    if transitWidth >= len(y):
        raise ValueError("Transit width must be shorter than length of binned lightcurve")

    template = np.ones(transitWidth)
    bls = -np.convolve(y, template, mode='same')
    pointsPerBin = np.convolve(counts, template, mode='same')

    assert(len(bls) == len(y))

    eps = 1e-20  #To prevent division by zero
    return bls / (uncPerUnbinnedPoint * np.sqrt(pointsPerBin) + eps)



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

#@profile
def computeGumbelLnFap(sigma, nDraws):
    """Compute the false alarm probability of a signal of strength
    sigma assuming it is the largest of nDraws samples.

    Take **nDraws** samples from a gaussian distribution with mean zero
    and standard deviation of one. Take the largest sample, **sigma**. The
    probability distribution of this signal follows a Gumbel distribution.
    The probability of seeing a signal of a given strength due to noise
    alone can then be computed with this function.

    For large signals, the false alarm probabilities quickly get ludicrously
    small, and uncomputatable. This function uses the approximation
    $\ln{FAP} \\approx \ln(nDraws) + \ln( 1 - cdf)$
    for values of **sigma** > 9

    Inputs:
    ----------
    sigma
        (1d array) The strength of the signal in units of sigma.

    nDraws
        (1d array) Number of samples taken before choosing the
        given **sigma** as the largest

    Returns:
    -----------
    ln(FAP), the natural log of the false alarm probability. Shape is
    the same as sigma.

    """
    #fap = 1 - cdf**N

    lnFap = np.zeros_like(sigma)
    idx = sigma < 8

    #Low significance events
    cdf = spStats.norm.cdf(sigma[idx])
    lnFap[idx] = np.log(1 - cdf**nDraws)

    #High signif event
    lsf = spStats.norm.logsf(sigma[~idx])
    lnFap[~idx] = np.log(nDraws) + lsf

    return lnFap



def exampleDurationFunction(period):
    """An example duration function for computeFBls().
    For illustration only, shouldn't be used in production code
    """
    return np.array([1,2,3,4,5])/24.

class ConstDurationSearch():
    """A helper class for computeFBls()
    Create one of these with the list of periods you want to
    search, and pass it into the above function. Every listed duration
    will be searched at every period

    Example:
    --------------
    dur = ConstDurationSearch( np.arange(1,6)/24.)
    computeFbls(time, flux.... durationFunc=dur)

    For every period searched, the bls will be run for durations
    of 1,2,3,4 and 5 hours
    """

    def __init__(self, durList):
        """
        Inputs:
        ------------
        durList:
            (list or array) List of durations to search in units of days
        """
        self.durList = durList

    def __call__(self, period):
        return self.durList
