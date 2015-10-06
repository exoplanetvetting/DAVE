# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 08:59:40 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.centroid.plotTpf as plotTpf
import dave.centroid.centroid as centroid
import dave.diffimg.diffimg as diffimg
import dave.diffimg.arclen as arclen
import dave.centroid.prf as prf


"""
Question, am I better off fitting diff image or signif image?
"""



def example():
    k2id =  206103150
    campaign = 3

    ar = mastio.K2Archive()
    fits, hdr = ar.getLongTpf(k2id, campaign, header=True)
    hdr0 = ar.getLongTpf(k2id, campaign, ext=0)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    idx = np.isfinite(cube)
    cube[~idx] = 0  #Remove Nans

    flags = fits['QUALITY']
    ccdMod = hdr0['module']
    ccdOut = hdr0['output']

    #Compute roll phase
    llc = ar.getLongCadence(k2id, campaign)
    time= llc['TIME']
    cent1 = llc['MOM_CENTR1']
    cent2 = llc['MOM_CENTR2']
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = arclen.computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    prfObj = prf.KeplerPrf("/home/fergal/data/keplerprf")
    bbox = centroid.getBoundingBoxForImage(cube[0], hdr)

    period =  	4.1591409
    epoch = fits['time'][491]
    dur = 3.0

    out = measureDiffOffset(period, epoch, dur, time, prfObj, \
        ccdMod, ccdOut, cube, bbox, rollPhase, flags)
    return out


def measureDiffOffset(period_days, epoch_bkjd, duration_hrs, \
    time, prfObj, ccdMod, ccdOut, cube, bbox, rollPhase, flags):
    """Measure Centroid shift between intransit and difference image
    for every in-transit cadence

    Inputs:
    -----------
    period_days, epoch_bkjd, duration_hrs
        (floats) Properties of transit

    time_bkjd
        Array of times per cadence for the given campaign

    prfObj
        An object of the class prf.KeplerPrf()

    ccdMod, ccdOut
        (int) CCD module and output of image. Needed to
        create the correct PRF model

    cube
        (3d np array) A data cube created from a TPF file.
        See fileio.tpf.getTargetPixelArrayFromFits()

    bbox
        [c1, c2, r1, r2]. Define the range of columns (c1..c2)
        and rows (r1..r2)  defined by the image.
        An exception raised if the following equality not true
        img.shape = (c2-c1), (r2-r1)

    rollPhase
        (1d np array) An array of roll phases for each row
        of cube. len(rollPhase) == len(cube). Units of this
        array don't matter, so long as cadences with similar
        roll angles have similar values of rollPhase. Roll phases
        for bad cadences should be set to a bad value

    flags
        (1d array) flag values indicating bad cadences.
        Currently a non-zero value of flags indicates a bad
        cadence.

    Returns:
    -------------
    A array with 5 columns, and as many rows as there are
    in transit cadences. The columns are

    0: Relative cadence number
    1: In transit centroid column
    2: In transit centroid row
    3: Diff img centroid column
    4: Diff img centroid row

    If there is a statisically significant difference between the intransit
    and difference image centroids then the transit is most likely not
    on the target.
    """
    idx = getIndicesInTransit(period_days, epoch_bkjd, duration_hrs, time)
    wh = np.where(idx)[0]
    out = -1 * np.ones((len(wh), 5))
    for i,w in enumerate(wh):
        out[i,0] = w
        try:
            out[i, 1:] = measureInTransitAndDiffCentroidForOneImg(\
                prfObj, ccdMod, ccdOut, cube, w, bbox, rollPhase, flags, \
                hdr=None, plot=False)
        except ValueError:
            pass
        print i, len(wh)

    return out


def getIndicesInTransit(period_days, epoch_bkjd, duration_hrs, time_bkjd):
    """
    Find the cadences affected by a transit.

    Inputs::
    --------
    period_days, epoch_bkjd, duration_hrs
        (floats) Properties of transit

    time_bkjd
        Array of times per cadence for the given campaign


    Returns:
    ------------
    An array of booleans of length equal to length of time_bkjd.
    Cadences in transit are set to true, all other cadences to false
    """

    time = time_bkjd #Mneumonic
    n1 = int(np.floor((time[0] - epoch_bkjd)/period_days))
    n2 = int(np.ceil((time[-1] - epoch_bkjd)/period_days))
    dur_days = duration_hrs/24.

    inTransit = np.zeros_like(time, dtype=bool)
    for n in range(n1, n2+1):
        t0 = epoch_bkjd + n*period_days - .5*dur_days
        t1 = epoch_bkjd + n*period_days + .5*dur_days

        idx = (time >= t0) & (time <= t1)
#        import pdb; pdb.set_trace()
        inTransit |= idx

    return inTransit


def measureInTransitAndDiffCentroidForOneImg(prfObj, ccdMod, ccdOut, cube, rin, bbox, rollPhase, flags, hdr=None, plot=False):
    """Measure image centroid of in-transit and difference images

    Inputs:
    -----------
    prfObj
        An object of the class prf.KeplerPrf()


    ccdMod, ccdOut
        (int) CCD module and output of image. Needed to
        create the correct PRF model


    cube
        (3d np array) A TPF data cube as returned by
        dave.fileio.getTargetPixelArrayFromFits()

    rin
        (int) Which image to process. rin should be in the range 0..len(cube)

    bbox
        [c1, c2, r1, r2]. Define the range of columns (c1..c2)
        and rows (r1..r2)  defined by the image.
        An exception raised if the following equality not true
        img.shape = (c2-c1), (r2-r1)

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
    ---------------
    hdr
        Fits header object for TPF file. Useful if you want to plot

    plot
        (bool) Request plots.


    Returns:
    -------------
    A 4 element numpy array
    ic  In transit centroid column
    ir  In transit centroid row
    dc  Difference image centroid column
    dr  Difference image centroid row
    """
    inTrans = cube[rin]
    diff, oot= diffimg.constructK2DifferenceImage(cube, rin, rollPhase, flags)

    ootRes = centroid.fitPrfCentroidForImage(oot, ccdMod, ccdOut, bbox, prfObj)
    diffRes = centroid.fitPrfCentroidForImage(diff, ccdMod, ccdOut, bbox, prfObj)

    if plot:
        mp.subplot(121)
        plotTpf.plotCadence(inTrans, hdr)
        mp.colorbar()
        mp.subplot(122)
        plotTpf.plotCadence(diff, hdr)
        mp.colorbar()

        mp.plot(ootRes.x[0], ootRes.x[1], 'ro', ms=12, alpha=.4)
        mp.plot(diffRes.x[0], diffRes.x[1], 'r^', ms=12, alpha=.4)

    return np.array([ootRes.x[0], ootRes.x[1], diffRes.x[0], diffRes.x[1]])



