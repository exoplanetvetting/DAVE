# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 20:14:24 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import dave.fileio.kplrfits as kplrfits
import dave.pipeline.clipboard as dpc
import dave.pipeline.gather as gather
import dave.pipeline.pipeline as pl

def wedgePlot(clipList):
    """Create a wedge plot by plotting period and epoch of all events.

    Based on similar diagnostic plot used by SOC pipeline.

    Inputs:
    -----------
    clipList
        (list) list of filenames of clips to process

    Returns:
    -----------
    A figure handle for the plot created.

    Notes:
    ----------
    Each target is represented by a semi-transparent grey dot
    indicating the period and epoch of that target. Those targets
    with the value of disposition.isCandidate are marked in red.
    """

    epic, vals = gather.gatherFunction(clipList, getPeriodEpochDuration)

    period = np.array(map(lambda x: x[0], vals))
    epoch = np.array(map(lambda x: x[1], vals))
    isCand = np.array(map(lambda x: x[3], vals))

    mp.clf()
    mp.plot(period, epoch, 'ko', alpha=.4, label="All Targets")
    mp.plot(period[isCand], epoch[isCand], 'ro', \
        ms=10, label="Candidates")

    mp.xlabel("Period (days)")
    mp.ylabel("Epoch (BKJD)")
    mp.legend(loc=0)
    return mp.gcf()

def skyLinePlot(clipList):
    """A plot of which cadences contribute to the most transits

    Based on similar plot created by Jessie Christiansen for the
    SOC pipeline.

    Inputs:
    -----------
    clipList
        (list) list of filenames of clips to process

    """

    epic, vals = gather.gatherFunction(clipList, getPeriodEpochDuration)

    clip = dpc.loadClipboard(clipList[0])
    clip = pl.serveTask(clip)
    time = clip['serve.time']
    flags= clip['detrend.flags']

    period = np.array(map(lambda x: x[0], vals))
    epoch = np.array(map(lambda x: x[1], vals))
    duration_days = np.array(map(lambda x: x[2], vals)) / 24.
    isCand = np.array(map(lambda x: x[3], vals))

    skyLine = time*0
    candSkyLine = time*0
    for i in range(len(period)):
        idx = kplrfits.markTransitCadences(time, period[i], epoch[i], \
            duration_days[i], flags=flags)
        skyLine[idx] += 1

        if isCand[i]:
            candSkyLine[idx] += 1

    mp.clf()
    mp.step(time[~flags], skyLine[~flags], 'b-', lw=2, \
        label="All targets")
    mp.step(time[~flags], candSkyLine[~flags], 'r-', lw=2, \
        label="Candidates")

    mp.xlabel("Time (BKJD)")
    mp.ylabel("Number of Transits on Cadence")
    return mp.gcf()

def getPeriodEpochDuration(clip):
    p = clip['trapFit.period_days']
    e = clip['trapFit.epoch_bkjd']
    d = clip['trapFit.duration_hrs']
    isCand = clip['disposition.isCandidate']

    return p, e, d, isCand