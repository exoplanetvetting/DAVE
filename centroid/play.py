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

def measureDiffOffset():
    k2id = 206103150
    campaign = 3
    rin =  489

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
    cent1 = llc['MOM_CENTR1']
    cent2 = llc['MOM_CENTR2']
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = arclen.computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    prfObj = prf.KeplerPrf("/home/fergal/data/keplerprf")
    bbox = centroid.getBoundingBoxForImage(cube[0], hdr)

    mp.clf()
    print measureInTransitAndDiffCentroidForOneImg(\
        prfObj, ccdMod, ccdOut, cube, rin, bbox, rollPhase, flags, \
        hdr, plot=True)


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



