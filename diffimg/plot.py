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
import dave.misc.covar as covar


def plotWrapper(clip):
    """Wrapper function for difference image centroid diagnostic plots

    Call this function from the exporter.

    Inputs:
    -----------
    clip
        A clipboard object. Should have the following keys: serve, detrend, diffImg, rollPhase, trapFit

    Returns:
    -----------
    Two figure handles. The zeroth figure handle shows the flux and rollPhase
    plot, the first shows the centroid offset plot

    Outputs:
    ------------
    Two figures are produced.

    """
    time = clip['serve.time']
    qFlags = clip['serve.flags']

    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    centroids = clip['diffImg.centroid_timeseries']
    rollPhase = clip['rollPhase.rollPhase']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']
    epic = clip['value']

#    tce = clip['eventList'][0]
#    period = tce['trapFit.period_days']
#    epoch = tce['trapFit.epoch_bkjd']
#    duration_hrs = tce['trapFit.duration_hrs']

    inTransitIndices = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, \
        duration_hrs/24., flags=flags)
    goodCentroidIndices = centroids[ centroids[:,1]>1, 0].asarray().astype(int)


    f1 = mp.figure(1)
    mp.clf()
    multiPanelPlotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, \
        inTransitIndices, goodCentroidIndices, qFlags)

    f2 = mp.figure(2)
    mp.clf()
    titleStr = plotCentroidOffsets(centroids)
    titleStr = "EPIC: %i  %s" %(epic, titleStr)
    mp.title(titleStr)

    return f1, f2


def plotCentroidOffsets(centroids):
    """Plot the centroid offsets in column and row

    Inputs:
    ----------
    centroids:
        (Nca) As stored in ``clip['diffImg.centroid_timeseries']``


    Returns:
    ----------
    A string containing the probability the measured centroids are consistent
    with no offset

    Output:
    ------------
    A plot is added to the current figure.
    """
    idx =centroids[:,1] > 0
    ootCol = centroids[idx, 'diff_col']
    ootRow = centroids[idx, 'diff_row']

    #itr => in transit
    itrCol = centroids[idx, 'intr_col']
    itrRow = centroids[idx, 'intr_row']

    diffC = ootCol - itrCol
    diffR = ootRow - itrRow

    mp.plot(diffC, diffR, 'bo', mec="none")
    mp.plot(0,0, 'ko', ms=12)
    mp.axhline(0, color='k', lw=.5)
    mp.axvline(0, color='k', lw=.5)

    covar.plotErrorEllipse(diffC, diffR)

    mp.xlabel(r"$\Delta$ Column")
    mp.ylabel(r"$\Delta$ Row")

    probOffset, chiSq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
    titleStr = "Prob. On Target: %.1e: $\chi^2$: %.3f" %(1-probOffset, chiSq)
    return titleStr


def multiPanelPlotDiffImgCentroidsDiagnostic(time, flux, flags, rollPhase, \
    inTransitIndices, goodCentroidIndices, qFlags):
    """Plot the flux and roll phase time series.

    A multi panel plot with three panels for the flux time series,
    and three for the roll phase time series. Cadences in transit are marked
    with blue or red squares. Blue indicates a difference image was successfully
    created, red indicatese no difference image.

    Each time series is spread over three panels to better show the details.
    Panels showing common sections of the timeseries are highlighted
    with the same background colour.

    Inputs:
    ----------
    time
        (1d numpy array) X axis for plot

    flux
        (1d numpy array) Flux time series

    flags
        (1d numpy boolean array) Cadences where flags is set are considered
        bad and not plotted.

    rollPhase
        (1d numpy array) Roll phase time series

    inTransitIndices
        (1d numpy boolean array) Which cadences are in transit

    goodCentroidIndices
        (1d numpy boolean array) Which cadences have centroids measured.


    Returns:
    ----------
    **None**

    Output:
    ------------
    A plot is added to the current figure.
    """
    fig = mp.gcf()
    fig.set_size_inches(11, 8.5)

    nPanel = 3
    start = np.min(time[~flags])
    deltaT = np.max(time[~flags]) - start
    deltaT /= float(nPanel)

    colour = ["#FFDDDD", "#DDFFDD", "#DDDDFF"]

    for i in range(nPanel):
        ax = mp.subplot(2*nPanel, 1, 2*i+1, axisbg=colour[i])
        plotTimeseries(time, 1e3*flux, flags, inTransitIndices, \
            goodCentroidIndices, qFlags)
        mp.ylabel("Frac Flux (ppk)")

        mp.subplot(2*nPanel, 1, 2*i+2, sharex=ax, axisbg=colour[i])
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
    """Plot a time series

    Also marks cadences with thruster firings with vertical lines on the plot,
    and the location of cadences in transit, and cadences with measured centroids
    """
    plotThrusterFirings(qFlags, time, color='#888888', lw=.4)
    mp.plot(time[~flags], y[~flags], 'ko', ms=2, alpha=.8)
    mp.plot(time[inTransitIndices], y[inTransitIndices], 'rs')
    mp.plot(time[goodCentroidIndices], y[goodCentroidIndices], 'bo')


def plotThrusterFirings(qualityFlags, xval=None, **kwargs):
    """Mark cadences with thruster firings

    Inputs:
    ----------
    qualityFlags
        (1d numpy array) Quality flags as found in lightcurves and TPF
        files. This is not a boolean array, but an array of flag masks.
        Must include a bit for thruster firings.

    Optional Inputs:
    -----------------
    xval
        (1d array) Values of x axis to plot (eg the time, or cadence number)
        Defeault is 1.. len(qualityFlags)

    All other optional arguments passed to matplotlibs axvline

    Returns:
    ----------
    **None**

    Output:
    ----------
    A series of vertical lines are added to the current plotx
    """

    flags = qualityFlags.astype(np.int32)  #Mneumonic
    if xval is None:
        xval = np.arange(len(flags))

    wh = np.where( flags & kplrfits.SapQuality['DefiniteRollTweak'])[0]

    for w in wh:
        mp.axvline(xval[w], **kwargs)