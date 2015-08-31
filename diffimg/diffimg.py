# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 21:08:35 2015

@author: fergal

$Id$
$URL$
"""

import numpy as np
import kplrfits

"""
Two main functions are

constructKeplerDifferenceImage uses the Kepler suitable DI algo
from Bryson13

constructK2DifferenceImage() trys to do something different for K2
"""

def constructKeplerDifferenceImage(time, cube, t0_days, duration_hrs, nDurForContinuum=1):
    """
    Generate a difference image of a TPF cube showing change in
    flux during transit.

    Inputs:
    ------------
    time
        (1d np array)  Mid times of cadences (usually in BKJD)
    cube
        (3d np array) TPF cube as created by tpf.getTargetPixelArrayFromFits()
    t0_days
        (float) Midpoint of transit. Units are the same as time
    duration_hrs
        (float) Transit duration in hours. Units are 24 times units of t0_days
    nDurForContinuum
        (float) Multiple of duration_hrs  on **either** side of the transit
        to use as the continuum

    Returns:
    -----------
    Three 2d images, the difference image, the Out of Transit Image and
    the In transit image. The diff image is calculate as oot-inTransit,
    so a transit appears as a positive feature

    Notes:
    -------
    The appearance of a transit as a positive feature mimics the behaviour
    of DV.
    """

    duration_days = duration_hrs/24.

    c1 = t0_days - duration_days*(0.5 + nDurForContinuum)
    c2 = t0_days - duration_days*0.5
    c3 = t0_days + duration_days*0.5
    c4 = t0_days + duration_days*(0.5 + nDurForContinuum)

    idx1 = (time >= c1) & (time < c2)
    idx2 = (time >= c2) & (time < c3)
    idx3 = (time >= c3) & (time < c4)

    test = np.sum(idx1) * np.sum(idx2) * np.sum(idx3)
    if test == 0:
        raise IndexError("At least one region returned no cadences")

    oot = np.sum(cube[idx1], axis=0) + np.sum(cube[idx3], axis=0)
    inTransit = np.sum(cube[idx2], axis=0)

    #Normalise
    oot /= np.sum(idx1 | idx3)
    inTransit /= np.sum(idx2)


    return oot-inTransit, oot, inTransit


def constructK2DifferenceImage(cube, indexOfCadencesInTransit, rollPhase, flags):

    #@TODO: Not all flags are bad
    badIdx = flags > 0

    ootOut = np.zeros_like( cube[0])
    diffOut = np.zeros_like(ootOut)
    nImg = 0
    for i in indexOfCadencesInTransit:
        try:
            indexBefore, indexAfter = getIndicesOfOotImages(rollPhase, i, flags)
        except ValueError:
            continue

        #TODO, look nearby for good cadences before skipping
        if badIdx[indexBefore] or badIdx[indexAfter]:
            continue

        oot = .5*(cube[indexBefore] + cube[indexAfter])
        diff = oot - cube[i]

        ootOut += oot
        diffOut += diff
        nImg += 1

    if nImg == 0:
        raise ValueError("No good Oot images created")

    return diffOut/float(nImg), ootOut/float(nImg)



def getIndicesOfOotImages(rollPhase, i, flags):
    """Compute the indices of two cadences of similar roll angle
    to the ith cadence index

    Inputs:
    --------
    rollPhase (1d np array)
        Phase of the roll angle (in arbitary units). Constructed from
        computeArcLength()[:,0]

    i (int)
        Index (not cadence number or time) of cadence for which you
        want to find the nearest cadences

    flags (1d np array)
        SAP_QUALITY flag from a lightcurve file, or the QUALITY flag
        from a TPF file


    Returns:
    -----------
    Two indices, before and after. Both lie in the range ``[0,len(rollPhase))``

    Notes:
    ---------
    * ``len(rollPhase)`` must be equal to ``len(flags)``, which is equal to the
      number of cadences in a campaign. Don't remove any bad cadences
      or the indexing will get screwed up.

    * Only works on K2 data, because it requires the Definite Thruster
      Firing flag to be set.

    * I need to check that nearest roll phase is in fact a similar roll
      phase.

    Description:
    -------------
    In Kepler, Difference images we obtained by comparing in transit
    cadences to nearby cadences. In K2, nearby isn't good enough because
    the roll of the spacecraft changes on the timescale of a transit.
    Instead we look for nearby points with similar roll angles.

    I don't create the difference image here, only find two suitable
    cadences to use as OOT points for a single in-transit cadence.

    Before constructing a difference image, check these cadences aren't
    themselves in transit.
    """
    rollPhase0 = rollPhase[i]

    dtFlag =kplrfits.SapQuality['DefiniteThruster']
    thrusterFiringIndices = np.where( flags & dtFlag )[0]
    tfi = thrusterFiringIndices
    tfiBefore = tfi[tfi < i]
    tfiAfter = tfi[tfi > i]

    lwr, upr = getThrusterFiringIntervalBefore(tfiBefore)
    if lwr is None:
        raise ValueError("No suitable interval found before requested cadence")

    indexBefore = np.argmin( np.fabs(rollPhase[lwr:upr] - rollPhase0)) + lwr

    lwr, upr = getThrusterFiringIntervalAfter(tfiAfter)
    if lwr is None:
        raise ValueError("No suitable interval found after requested cadence")

    indexAfter = np.argmin( np.fabs(rollPhase[lwr:upr] - rollPhase0)) + lwr
    return indexBefore, indexAfter



def getThrusterFiringIntervalBefore(tfIndices, minDuration = 10):
    """Returns range of points between last two thruster firings in input

    See getThrusterFiringIntervalAfter() for more details. This function
    does the same job, but returns the last two indices in the array
    meeting the criteria.
    """

    i = len(tfIndices)-1
    while i > 0:
        if tfIndices[i-1] + minDuration < tfIndices[i]:
            return tfIndices[i-1], tfIndices[i]
        i -= 1
    return None, None


def getThrusterFiringIntervalAfter(tfIndices, minDuration=10):
    """Get range of points between first two thruster firings in input

    Input:
    --------
    tfIndices (1d np array)
        A list of numbers, each number represents the index where a
        thruster firing occurs. This is not a boolean array


    Returns:
    ----------
    Values for first two numbers separated by more than minDuration

    Example:
    ---------
    self( [1,3,15,9])
    returns [3,15]
    """
    i = 0
    numTf=  len(tfIndices)
    while i < numTf:
#        print tfIndices[i], tfIndices[i+1]
        if tfIndices[i] + minDuration < tfIndices[i+1]:
            return tfIndices[i], tfIndices[i+1]
        i += 1
    return None, None



import matplotlib.pyplot as mp
import plotTpf

def plotDiffDiagnostic(diff, oot, inTransit):
    mp.clf()
    cm = mp.cm.gnuplot2
    mp.subplot(131)
    mp.title("Diff Img")
    plotTpf.plotCadence(diff, axis="relative", cmap=mp.cm.RdBu)
    mp.colorbar()
    mp.subplot(132)
    mp.title("OOT")
    plotTpf.plotCadence(oot, axis="relative", cmap=cm)
    mp.colorbar()
    mp.subplot(133)
    mp.title("In Transit")
    plotTpf.plotCadence(inTransit, axis="relative", cmap=cm)
    mp.colorbar()


"""
def example():
    import tpf

    kepid = 11904151
    t0 = 138.677668
    dur_hrs = 6

    ar = mastio.KeplerArchive()
    fits, hdr = ar.getTargetPixelFile(kepid, 1, header=True)
    time = fits['TIME']
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    res = constructKeplerDifferenceImage(time, cube, t0, dur_hrs)

    mp.figure(1)
    mp.clf()
    mp.subplot(131)
    mp.title("Diff Img")
    plotTpf.plotCadence(res[0], axis="relative")
    mp.colorbar()
    mp.subplot(132)
    mp.title("OOT")
    plotTpf.plotCadence(res[1], axis="relative")
    mp.colorbar()
    mp.subplot(133)
    mp.title("In Transit")
    plotTpf.plotCadence(res[2], axis="relative")
    mp.colorbar()

    return
    mp.figure(2)
    rng = 1.5*(dur_hrs/24.)
    idx = (time > t0-rng) & (time < t0+rng)
    print t0
    print t0 - 1.5*.125
    print 1.5*.125*48
#    plotTpf.plotTpfLc(cube[idx, 2:6, 2:6], hdr, big=True)
    plotTpf.plotTpfLc(cube[idx], hdr, big=True)
"""

