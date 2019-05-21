# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 21:08:35 2015

@author: fergal

$Id$
$URL$
"""

from __future__ import division, print_function, absolute_import

import matplotlib.pyplot as mp
import numpy as np

import dave.fileio.kplrfits as kplrfits


def constructK2DifferenceImage(cube, indexOfCadenceInTransit, rollPhase, flags):
    """Construct a difference image for a single K2 cadence

    K2's roll motion makes constructing difference images harder
    than it was in Classic Kepler. The star can not be assumed
    to be in a similar location on the focal plane in any other
    cadence. To construct difference images, this function
    interpolates between cadences of similar roll angles to
    construct an "out-of-transit" image

    Inputs:
    -------------
    cube
        (3d np array) A data cube created from a TPF file.
        See fileio.tpf.getTargetPixelArrayFromFits()
    indexOfCadenceInTransit
        (int) Specify which row of cube to construct a difference
        image for

    rollPhase
        (1d np array) An array of roll phases for each row
        of cube. len(rollPhase) == len(cube). Units of this
        array don't matter, so long as cadences with similar
        roll angles have similar values of rollPhase

    flags
        (1d array) flag values indicating bad cadences.
        Currently a non-zero value of flags indicates a bad
        cadence.


    Returns:
    ------------
    A 3 element list, [diff, oot, diagnostics]
    diff is a 2d np array corresponding to the difference image
    oot is a 2d np array corresponds to the reconstructed image
    subtracted from cube[indexOfCadence] in transit.

    diff = cube[indexOfCadenceInTransit] - oot

    diagnostics is a dictionary of details of how the image was
    constructed.

    Notes:
    ---------
    In keeping with the Classic Kepler convention, a flux decrement
    in a cadence (e.g a transit) appears as a positive feature in
    the difference image.

    Todo:
    ---------
    I should do a better job picking which cadences are bad.
    Not all flags indicate a cadence is untrustworthy.
    """


    assert(cube.shape[0] == len(rollPhase))
    i0 = indexOfCadenceInTransit

    try:
        oot, diagnostics = getInterpolatedOotImage(cube, rollPhase, flags, i0)
    except ValueError as e:
        msg = "WARN: While creating Diff Img for %i: %s" \
            %(i0, e)
        raise ValueError(msg)

    diff = oot - cube[i0]
    return diff, oot, diagnostics


def plotRollPhaseDiagnosticPlot(x, rollPhase, flags, indexOfEvent):
    if x is None:
        x = np.arange(len(rollPhase))

    if flags.dtype == bool:
        raise ValueError("Boolean flag array not accepted")

    i0 = indexOfEvent
    idx = flags == 0
    mp.plot(x[idx], rollPhase[idx], 'ko-')
    mp.axvline(i0, color='b')
    mp.axhline(rollPhase[i0], color='grey')

    dtFlag =kplrfits.SapQuality['PossibleRollTweak']
    thrusterFiringIndices = np.where( flags & dtFlag )[0]
    tfi = thrusterFiringIndices
    tfiBefore = tfi[tfi < i0]
    tfiAfter = tfi[tfi > i0]


    lwr = max(0, i0-100)
    upr = min(len(x)-1, i0+100)
    plotLwr = x[lwr]
    plotUpr = x[upr]
    mp.xlim(plotLwr, plotUpr)
    mp.ylim(-1,1)


    for w in tfi:
        if x[w] < plotLwr or x[w] > plotUpr:
            continue
        mp.axvline(x[w], color='r')

    lwr, upr = getThrusterFiringIntervalBefore(tfiBefore)
    i1, sgn, d1, d2 = getBracketingCadence(rollPhase, lwr, upr, rollPhase[i0])
    mp.plot(x[i1], rollPhase[i1], 'ro')
    mp.plot(x[i1+sgn], rollPhase[i1+sgn], 'ro')

    lwr, upr = getThrusterFiringIntervalAfter(tfiAfter)
    i1, sgn, d1, d2 = getBracketingCadence(rollPhase, lwr, upr, rollPhase[i0])
    mp.plot(x[i1], rollPhase[i1], 'ro')
    mp.plot(x[i1+sgn], rollPhase[i1+sgn], 'ro')




def getInterpolatedOotImage(cube, rollPhase, flags, i0):
    """Construct an out-of-transit image for a given cadence

    Inputs:
    -------------
    cube
        (3d np array) A data cube created from a TPF file.
        See io.tpf.getTargetPixelArrayFromFits()

    rollPhase
        (1d np array) An array of roll phases for each row
        of cube. len(rollPhase) == len(cube). Units of this
        array don't matter, so long as cadences with similar
        roll angles have similar values of rollPhase

    flags
        (1d array) flag values indicating bad cadences.
        Currently a non-zero value of flags indicates a bad
        cadence.

    i0
        (int) Specify which row of cube to construct a difference
        image for


    Returns:
    -------------
    A two element tuple:

    *  A np 2d array representing an interpolated image.
    * A dictionary containing the indices of the cadences used for interpolation. The
      keys of that dictionary are:

    rangeBefore
        (2 element tuple) The index range of between the two previous
        thruster firings
    rangeAfter
        (2 element tuple) The index range of between the two succeeding
        thruster firings
    rinBefore
        (2 element tuple) The bracking cadences interpolated to produce the
        OOT image before the event
    rinAfter
        (2 element tuple) The bracking cadences interpolated to produce the
        OOT image after the event
    """

    rollPhase0 = rollPhase[i0]
    oot = np.zeros_like(cube[0])
    diagnostics = dict()
    diagnostics['errorMsg'] = "None"

    dtFlag =kplrfits.SapQuality['DefiniteRollTweak']
    thrusterFiringIndices = np.where( flags & dtFlag )[0]
    tfi = thrusterFiringIndices
    tfiBefore = tfi[tfi < i0]
    tfiAfter = tfi[tfi > i0]

#    import pdb; pdb.set_trace()
    lwrBefore, uprBefore = getThrusterFiringIntervalBefore(tfiBefore)
    lwrAfter, uprAfter = getThrusterFiringIntervalAfter(tfiAfter)

    diagnostics['rangeBefore'] = (lwrBefore, uprBefore)
    diagnostics['rangeAfter'] = (lwrAfter, uprAfter)


    if lwrAfter is None:
        diagnostics['errorMsg'] = "No suitable range found before cadence of interest"
        return oot, diagnostics
    if lwrBefore is None:
        diagnostics['errorMsg'] = "No suitable range found after cadence of interest"
        return oot, diagnostics

    try:
        ootBefore, rinBefore = \
            getDiffFromRange(cube, lwrBefore, uprBefore, rollPhase, rollPhase0)
        diagnostics['rinBefore'] = rinBefore
    except ValueError as e:
        diagnostics['errorMsg'] = "Early diff img: %s" %(e)
        return oot, diagnostics

    #@TODO: Should I just return OOT before here?
    try:
        ootAfter, rinAfter = getDiffFromRange(cube, lwrAfter, uprAfter, \
            rollPhase, rollPhase0)
        diagnostics['rinAfter'] = rinAfter
    except ValueError as e:
        diagnostics['errorMsg'] = "Later diff img: %s" %(e)
        return oot, diagnostics

    oot = .5 * (ootBefore + ootAfter)
    return oot, diagnostics



def getDiffFromRange(cube, lwr, upr, rollPhase, rollPhase0):
    """
    Construct an interpolated difference image.

    Construct an interpolated difference image from the data
    in cube[lwr:upr] to match what would be expected if
    a cadence was observed centred on rollPhase0.

    Inputs:
    -------------
    cube
        (3d np array) A data cube created from a TPF file.
        See io.tpf.getTargetPixelArrayFromFits()

    lwr, upr
        (int) Range of cube to interpolate within.

    rollPhase
        (1d np array) An array of roll phases for each row
        of cube. len(rollPhase) == len(cube). Units of this
        array don't matter, so long as cadences with similar
        roll angles have similar values of rollPhase

    rollPhase0
        (float) the value of rollPhase to interpolate to.
        An exception is raised if rollPhase[lwr:upr] does
        not bracket rollPhase0


    Returns:
    --------------
    A two element tuple
    *  A np 2d array representing an interpolated image.
    * A tuple containing the indices of the cadences used for interpolation

    TODO:
    ------------
    Don't use cadences that are flagged in some bad way.
    """
    maxDiffBetweenAdjacentPoints = .15
    i0, sgn, d1, d2 = getBracketingCadence(rollPhase, lwr, upr, rollPhase0)


    #d1 and d2 must bracket rollPhase0. They should also
    #have similar values
    if d1*d2 < 0 and np.fabs(d2-d1) < maxDiffBetweenAdjacentPoints:
        diff = interpolateImages(rollPhase[i0], rollPhase[i0+sgn], \
            cube[i0], cube[i0+sgn], rollPhase0 )
        return diff, (i0, i0+sgn)

    raise ValueError("Can't produce difference image")


def getBracketingCadence(rollPhase, lwr, upr, rollPhase0):
    """Get bracketing cadence for a given rollphase.

    Computes i0 and sgn such that:

    * i0 is in the range [lwr, upr]
    * sgn is either +1 or -1
    * rollPhase[i0] and rollPhase[i0+sgn] bracket rollPhase0, i.e
      one value is larger than rollPhase0 and the other is lower.

    Inputs:
    ------------
    rollPhase
        (1d numpy array)  Values for rollphase
    lwr, upr
        (integers), range of values in rollphase to search
    rollPhase0
        The roll phase to bracket

    Returns:
    -----------
    i0
        (int) First bracketing cadence
    sgn
        (int) Either +1 or -1. The second bracketing cadence is i0+sgn
    d1, d2
        (int) Rollphase difference between rp[i0], rp[i0+sng] and rollPhase0
        Can be used to assess the quality of the bracket chosen.

    """
    i0 = np.argmin( np.fabs(rollPhase[lwr:upr] - rollPhase0)) + lwr

    d1 = rollPhase0 - rollPhase[i0]
    slope = rollPhase[i0+1] - rollPhase[i0-1]

    #Should I interpolate with previous or next cadence
    #If the slope is positive, and rp0 is larger than rp[i0]
    #I want to use the next cadence.
    sgn = 2*int(d1*slope > 0) - 1
    d2 = rollPhase0 - rollPhase[i0 + sgn]
    return i0, sgn, d1, d2


def interpolateImages(x0, x1, img0, img1, x):
    """Private function used by getDiffFromRange()"""
    f = (x-x0)/(x1-x0)
    return img0 + f*(img1-img0)


def getThrusterFiringIntervalBefore(tfIndices, minDuration = 10):
    """Returns range of points between last two thruster firings in input

    See getThrusterFiringIntervalAfter() for more details. This
    function does the same job, but returns the last two indices
    in the array  meeting the criteria.

    """

    nTfi = len(tfIndices)
    if nTfi == 0:
        return None, None

    if nTfi == 1:
        return 0, tfIndices[0]

    i = len(tfIndices)-1
    while i > 1:
        if tfIndices[i-1] + minDuration < tfIndices[i]:
            return tfIndices[i-1], tfIndices[i]
        i -= 1
    return None, None


def getThrusterFiringIntervalAfter(tfIndices, minDuration=10):
    """Get range of points between first two thruster firings in input

    Thruster firings tend to cluster, so we don't just want the first
    pair of firings in the array. Better is the first pair of firings that
    are separated by a minimum number of cadecnces, minDuration

    Input:
    --------
    tfIndices (1d np array)
        A list of numbers, each number represents the index where a
        thruster firing occurs. This is not a boolean array

    Optional Input:
    ---------------
    minDuration
        A pair of cadences must be separated by at least this many cadences
        to be returned.

    Returns:
    ----------
    Values for first two numbers separated by more than minDuration.



    Example:
    ---------
    ``getThrusterFiringIntervalAfter( [1,3,15,29], 10)``
    returns ``[3,15]``
    """
    numTf=  len(tfIndices)

    if numTf == 0:
        return None, None

    if numTf == 1:
        return tfIndices[0], -1

    i = 0
    while i < numTf:
        if tfIndices[i] + minDuration < tfIndices[i+1]:
            return tfIndices[i], tfIndices[i+1]
        i += 1
    return None, None




def plotDiffDiagnostic(diff, oot):
    """Only works on my machine"""
    mp.subplot(121)
#    img = np.log10(1+diff - np.min(diff[np.isfinite(diff)]))
    plotTpf.plotCadence(diff, axis="relative")
    mp.colorbar()
    mp.subplot(122)
    plotTpf.plotCadence(diff/np.sqrt(2*oot), axis="relative")
    mp.colorbar()
    mp.clim(-3.0, 3.0)



import dave.diffimg.arclen as arclen
import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf

def example():
    ar = mastio.K2Archive()

    kepid = 206103150  #A wasp planet

    fits = ar.getLongCadence(kepid, 3)
    flags = fits['SAP_QUALITY']
    cent1 = fits['MOM_CENTR1']
    cent2 = fits['MOM_CENTR2']

    fits, hdr = ar.getLongTpf(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
#    cube *= gain

    #Compute roll phase
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = arclen.computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    cadenceInTransit = 490
    diff, oot = constructK2DifferenceImage(cube,  cadenceInTransit, rollPhase, flags)
    return diff, oot
