# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 13:07:38 2015

@author: fergal

$Id$
$URL$
"""

import matplotlib.pyplot as mp
import numpy as np

import diffimg2 as diffimg
import mastio2
import kplrfits
import plotTpf
import tpf

import plotKeplerFlags as pkf


"""
Lab Notes:

The first time I picked, BKJD-2154 is during a period of minimum torque.
The diffence images according to the algo at the time look pretty good.
(newway1cadence.png).

When I look at the 8 nearby DIs (diff*.png), or the combination thereof
(newway-9cin) the results aren't so great. Most DIs look pretty good,
but the odd one looks weird.

The next step is to look at a constant star see what the performance
looks like. See wdTest.py

"""

def computeArcLength(cent_colrow, flag):
    """Compute arclength along the main eigenvector for each centroid point

    Inspired by similar approach by Vanderburg  & Johnson (2014)

    This is done slightly differently to V&J. I don't bother computing
    arclength along the fit, I just say how far along the eigenvector
    each centroid point lies. This is noisier than their approach, but
    simpler to implement.

    Inputs:
    -------------
    cent
        (2d numpy array) Array of centroid positions. Zeroth col is column
        position, and 1st col is row position. Bad values should be stripped
        out of this array before calling compute arclenght.

    Returns:
    A numpy array of shape cent. The zeroth column is the value of arclength
    for each input row.
    """


    cent = cent_colrow.copy()
    #Measure eigenvectors of mean-removed centroid distribution
    c0 = nanmean(cent[~flag,0])
    r0 = nanmean(cent[~flag,1])
    cent[:,0] -= c0
    cent[:,1] -= r0
    cov = np.cov(cent[~flag, :].transpose())

    [eVal, eVec] = np.linalg.eigh(cov)

    rot = cent * 0
    rot[:,0] = cent[:,0]* eVec[0,1] + cent[:,1]*eVec[1,1]
    rot[:,1] = cent[:,0]* eVec[0,0] + cent[:,1]*eVec[0,1]
    return rot



def nanmean(a):
    """Not available in earlier versions of numpy so I implement it here"""
    idx = np.isfinite(a)
    return np.mean(a[idx])


def plotThrusterFirings(flags, xval=None):

    if xval is None:
        xval = np.arange(len(flags))

    wh = np.where( flags & kplrfits.SapQuality['DefiniteThruster'])[0]
    for w in wh:
        mp.axvline(xval[w])



def plotDiffer(cube, hdr, index):

    diff = cube[index]- cube[index+1]
    mp.clf()

    mp.subplot(221)
    plotTpf.plotCadence(cube[index], hdr)
    mp.colorbar()

    mp.subplot(222)
    plotTpf.plotCadence(cube[index+1], hdr)
    mp.colorbar()

    mp.subplot(223)
    plotTpf.plotCadence(diff, hdr)
    mp.colorbar()

    mp.subplot(224)
    noise = np.sqrt(cube[index]+cube[index+1])
    z = diff/noise
    plotTpf.plotCadence(z, hdr)

    idx = np.isfinite(z)
    zmax = max( -np.min(z[idx]), np.max(z[idx]))
    mp.clim([-zmax, zmax])
    mp.colorbar()



def main():
    ar = mastio2.K2Archive()

    kepid = 206103150
    period = 4.16
    t0 = 2168.43 #2154.14

    fits = ar.getLongCadence(kepid, 3)
    time = fits['TIME']
    cin = fits['CADENCENO']
    flags = fits['SAP_QUALITY']
    pa = fits['SAP_FLUX']
    pdc = fits['PDCSAP_FLUX']
    cent1 = fits['MOM_CENTR1']
    cent2 = fits['MOM_CENTR2']
    badIdx = flags > 0

    i0 = np.random.randint(0, len(cin))

    fits, hdr = ar.getLongTpf(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    gain = hdr['gain']
#    cube *= gain

    #Compute roll phase
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value


    diffimg.constructK2DifferenceImage(cube,  578, rollPhase, flags)
    return
    for i in range(180, len(rollPhase)):
        mp.clf()
        try:
            diffimg.constructK2DifferenceImage(cube,  [i], rollPhase, flags)
        except ValueError:
            continue

        mp.savefig('fig-%05i.png' %(i))

def oldStuff():
    #Compute roll phase
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    #Index of transit
    mp.clf()
    time[ ~np.isfinite(time)] = -1
    i = np.argmin( np.fabs( time-t0))

    i0 = i
    indexBefore, indexAfter = diffimg.getIndicesOfOotImages(rollPhase, i0, flags)
    if True:
        mp.figure(1)
        mp.clf()
        rp0 = rollPhase[i]
        mp.axvline(time[i], color='r')
        mp.axhline(rp0, color='grey')
        mp.axvline(time[indexBefore], color='g')
        mp.axvline(time[indexAfter], color='g')
        mp.plot(time[~badIdx], rollPhase[~badIdx], 'k.')
        mp.plot(time[~badIdx], rot[~badIdx, 1], 'r.')
        plotThrusterFirings(flags,time)
#            mp.xlim([2150, 2160])
#            mp.savefig('diag%+i.png' %(i0))
        mp.xlim(2162,2172)

    fits, hdr = ar.getTargetPixelFile(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    mp.figure(2)
    mp.clf()
#    indexInTransit = np.arange(i-4, i+5)
    indexInTransit = [i0]
    diff, oot = diffimg.constructK2DifferenceImage(cube, indexInTransit, rollPhase, flags)
    diffimg.plotDiffDiagnostic(diff, oot, cube[i])
#        mp.savefig('diff%+i.png' %(i0))




















def oldStuff():
    ar = mastio2.K2Archive()

    kepid = 206103150
    period = 4.16
    t0 = 2154.14

    fits = ar.getLongCadence(kepid, 3)
    data = kplrfits.getNumpyArrayFromFitsRec(fits)
    time = fits['TIME']
    cin = fits['CADENCENO']
    flags = fits['SAP_QUALITY']
    pa = fits['SAP_FLUX']
    pdc = fits['PDCSAP_FLUX']

    print fits.dtype
    mp.clf()
    mp.plot( fits['MOM_CENTR1'], fits['MOM_CENTR2'], 'k.')
    return


    mp.figure(2)
    mp.clf()
    pkf.plotKeplerFlags(flags, time)
    mp.axvline(t0, color='r')
    mp.xlim(2153, 2156)


    mp.figure(3)
    mp.clf()
    mp.plot(cin, pdc, 'k.')
    wh = np.where( flags & kplrfits.SapQuality['DefiniteThruster'])[0]
    for w in wh:
        mp.axvline(cin[w])
    mp.xlim(100200, 103600)

#    epochs = np.arange(-10,10)
#    print epochs
#    for e in epochs:
#        t = t0 + e*period
#        mp.axvline(t, color='r')
    return


    fits, hdr = ar.getTargetPixelFile(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    diff, oot, intr = diffimg.generateDifferenceImage(time, cube, t0, 4.0)

    mp.figure(1)
    diffimg.plotDiffDiagnostic(diff, oot, intr)

    mp.figure(3)
    idx = (time > t0-1) & (time < t0+1)
    plotTpf.plotTpfLc(cube[idx, :, :], hdr, big=True)