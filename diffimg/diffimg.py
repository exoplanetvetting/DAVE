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


def constructK2DifferenceImage(cube, indexOfCadenceInTransit, rollPhase, flags, plot=False):
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
        See io.tpf.getTargetPixelArrayFromFits()
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


    Optional Inputs:
    -----------------
    plot
        (boolean) If true, create diagnostic plots

    Returns:
    ------------
    A 2 element list, [diff, oot]
    diff is a 2d np array corresponding to the difference image
    oot is a 2d np array corresponds to the reconstructed image
    subtracted from cube[indexOfCadence] in transit.

    diff = cube[indexOfCadenceInTransit] - oot


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


#    if plot:
#    idx = flags == 0
#    x = np.arange(len(rollPhase))
#        mp.figure(1)
#        mp.clf()
#        mp.plot(x[idx], rollPhase[idx], 'ko-')
#        mp.axvline(i0, color='b')
#        mp.xlim(i0-100, i0+100)
#        mp.axhline(rollPhase[i0], color='grey')

    try:
        oot, before, after = getInterpolatedOotImage(cube, rollPhase, flags, i0)
    except ValueError, e:
        msg = "WARN: No good Oot images created for %i: %s" \
            %(i0, e)
        raise ValueError(msg)

    diff = oot - cube[i0]


#    if plot:
#        mp.figure(2)
#        mp.subplot(121)
#    #    img = np.log10(1+diff - np.min(diff[np.isfinite(diff)]))
#        plotTpf.plotCadence(diff, axis="relative")
#        mp.colorbar()
#        mp.subplot(122)
#        plotTpf.plotCadence(diff/np.sqrt(2*oot), axis="relative")
#        mp.colorbar()
#        mp.clim(-3.0, 3.0)
#        mp.suptitle("Diff Img for RIN %i" %(i0))

    return diff, oot





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
    3 elements
    * A np 2d array containing the requested difference image
    * A 2 element list showing the cadence range between roll tweaks
      immediately before i0
    * A 2 element list showing the cadence range between roll tweaks
      immediately after i0
    """


    rollPhase0 = rollPhase[i0]

    dtFlag =kplrfits.SapQuality['PossibleRollTweak']
    thrusterFiringIndices = np.where( flags & dtFlag )[0]
    tfi = thrusterFiringIndices
    tfiBefore = tfi[tfi < i0]
    tfiAfter = tfi[tfi > i0]

#    import pdb; pdb.set_trace()
    lwrBefore, uprBefore = getThrusterFiringIntervalBefore(tfiBefore)

    if lwrBefore is None:
        raise ValueError("No suitable interval found before requested cadence")
    ootBefore = getDiffFromRange(cube, lwrBefore, uprBefore, rollPhase, \
        rollPhase0)


    lwrAfter, uprAfter = getThrusterFiringIntervalAfter(tfiAfter)
    if lwrAfter is None:
        raise ValueError("No suitable interval found after requested cadence")
    ootAfter = getDiffFromRange(cube, lwrAfter, uprAfter, rollPhase, rollPhase0)

    oot = .5 * (ootBefore + ootAfter)
    return oot, [lwrBefore, uprBefore], [lwrAfter, uprAfter]



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
    A np 2d array representing an interpolated image.


    TODO:
    ------------
    Don't use cadences that are flagged in some bad way.
    """

    maxDiffBetweenAdjacentPoints = .15
    i0 = np.argmin( np.fabs(rollPhase[lwr:upr] - rollPhase0)) + lwr

    d1 = rollPhase0 - rollPhase[i0]
    slope = rollPhase[i0+1] - rollPhase[i0-1]

    #Should I interpolate with previous or next cadence
    #If the slope is positive, and rp0 is larger than rp[i0]
    #I want to use the next cadence.
    sgn = 2*int(d1*slope > 0) - 1
    d2 = rollPhase0 - rollPhase[i0 + sgn]

#    x = np.arange(len(rollPhase))
#    mp.plot(x[i0], rollPhase[i0], 'ro')
#    mp.plot(x[i0+sgn], rollPhase[i0+sgn], 'ro')

    #d1 and d2 must bracket rollPhase0. They should also
    #have similar values
    if d1*d2 < 0 and np.fabs(d2-d1) < maxDiffBetweenAdjacentPoints:
        diff = interpolateImages(rollPhase[i0], rollPhase[i0+sgn], \
            cube[i0], cube[i0+sgn], rollPhase0 )
        return diff

#    import pdb; pdb.set_trace()
    raise ValueError("Can't produce difference image")


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