# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 13:02:31 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import dave.fileio.kplrfits as kplrfits
import diffimg

def plotWrapper(clip):

    time = clip['serve.time']
    qFlags = clip['serve.flags']

    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    centroids = clip['diffImg.centroid_timeseries']
    rollPhase = clip['rollPhase.rollPhase']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']

#    tce = clip['eventList'][0]
#    period = tce['trapFit.period_days']
#    epoch = tce['trapFit.epoch_bkjd']
#    duration_hrs = tce['trapFit.duration_hrs']

    inTransitIndices = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, \
        duration_hrs/24., flags=flags)
    goodCentroidIndices = centroids[ centroids[:,1]>1, 0].asarray().astype(int)


#    mp.figure(1)
#    multiPanelPlotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, \
#        inTransitIndices, goodCentroidIndices, qFlags)

    mp.figure(2)
    plotCentroidOffsets(centroids)


def plotCentroidOffsets(centroids):

    idx =centroids[:,1] > 0
    ootCol = centroids[idx, 'diff_col']
    ootRow = centroids[idx, 'diff_row']

    #itr => in transit
    itrCol = centroids[idx, 'intr_col']
    itrRow = centroids[idx, 'intr_row']

    diffC = ootCol - itrCol
    diffR = ootRow - itrRow

    meanDiffCol = np.mean(diffC)
    rmsDiffCol = np.std(diffC - meanDiffCol)
    meanDiffRow = np.mean(diffR)
    rmsDiffRow = np.std(diffR - meanDiffRow)

    #Compute, offset, unc and std
    unc = np.hypot(rmsDiffCol, rmsDiffRow) / np.sqrt(np.sum(idx))

    mp.clf()
    mp.plot(diffC, diffR, 'bo', mec="none", ms=8)
    mp.plot([0], [0], 'ko', ms=10)
    mp.axhline(0, lw=1, color='#222222')
    mp.axvline(0, lw=1, color='#222222')

    ax = mp.gca()
    mp.plot(meanDiffCol, meanDiffRow, 'ro', mec="none", ms=10)
    for i in range(1,4):
        radius = i*unc
        circle = mp.Circle((meanDiffCol, meanDiffRow), radius, color='r', alpha=.2)
        ax.add_artist(circle)

    mp.xlabel(r"$\Delta$ Column")
    mp.ylabel(r"$\Delta$ Row")


def multiPanelPlotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, \
    inTransitIndices, goodCentroidIndices, qFlags):

    fig = mp.gcf()
    fig.set_size_inches(11, 8.5)

    nPanel = 3
    start = np.min(time[~flags])
    deltaT = np.max(time[~flags]) - start
    deltaT /= float(nPanel)

    for i in range(nPanel):
        ax = mp.subplot(2*nPanel, 1, 2*i+1)
        plotTimeseries(time, flux, flags, inTransitIndices, \
            goodCentroidIndices, qFlags)
        mp.ylabel("Frac Flux")

        mp.subplot(2*nPanel, 1, 2*i+2, sharex=ax)
        plotTimeseries(time, rollPhase, flags, inTransitIndices, \
            goodCentroidIndices, qFlags)
        mp.ylim(-1.5,1.5)
        mp.ylabel("Roll Phase")

        mp.xlim(start + i*deltaT, start + (i+1)*deltaT)
        mp.xlabel("Time (BKJD)")


#def plotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, \
#    inTransitIndices, centroids, qFlags):
#
#    idx =
#
#
##    for i0 in np.where(idx)[0]:
##        i1, i2 = diffimg.getIndicesOfOotImages(rollPhase, i0, flags)
##        indexBefore.append(i1)
##        indexAfter.append(i2)
#
#
#    print centroids[:10]
#
#    mp.clf()
#    ax = mp.subplot(211)
#    plotThrusterFirings(qFlags, time, color='#888888', lw=.4)
#    mp.plot(time[~flags], flux[~flags], 'ko', ms=4, alpha=.8)
#    mp.plot(time[inTransitIndices], flux[inTransitIndices], 'rs')
#    mp.plot(time[goodCentroidIndices], flux[goodCentroidIndices], 'bo')
#
#    mp.subplot(212, sharex=ax)
#    plotThrusterFirings(qFlags, time, color='#888888', lw=.4)
#    mp.plot(time[~flags], rollPhase[~flags], 'ko', ms=4, alpha=.8)
#    mp.plot(time[inTransitIndices], rollPhase[inTransitIndices], 'ro')
#    mp.plot(time[goodCentroidIndices], rollPhase[goodCentroidIndices], 'bo')
#
##    mp.xlim(2144,2156)

def plotTimeseries(time, y, flags, inTransitIndices, goodCentroidIndices, qFlags):
    plotThrusterFirings(qFlags, time, color='#888888', lw=.4)
    mp.plot(time[~flags], y[~flags], 'ko', ms=4, alpha=.8)
    mp.plot(time[inTransitIndices], y[inTransitIndices], 'rs')
    mp.plot(time[goodCentroidIndices], y[goodCentroidIndices], 'bo')


def plotThrusterFirings(qualityFlags, xval=None, **kwargs):

    flags = qualityFlags.astype(np.int32)  #Mneumonic
    if xval is None:
        xval = np.arange(len(flags))

    wh = np.where( flags & kplrfits.SapQuality['DefiniteRollTweak'])[0]
#    wh = np.where( flags & kplrfits.SapQuality['PossibleRollTweak'])[0]

    for w in wh:
        mp.axvline(xval[w], **kwargs)
