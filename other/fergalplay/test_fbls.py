# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:54:08 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import shelve
import fbls


def loadData():
    sh = shelve.open('wasp47b.clip')
    clip = sh['clip']

    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    time = time[~flags]
    flux = flux[~flags]

    return time, flux


import dave.trapezoidFit.trapfit as dtf
def makeTestData():
    t = np.linspace(0, 60, 60*48)
    period = 2.41
    epoch = 1.1
    y = dtf.trapezoid_model_onemodel(t, period, epoch, 100, 6, 1, 15).modellc

    y /= np.mean(y)
    y-=1
    return t,y



def test_longDurations():
    """Test that bls is not computed if transit duration is too long"""
    time,flux = loadData()

    periodList = [10.]
    duration_daysList = [11]

    blsArray = fbls.computeBlsForManyPeriods(time, flux, duration_daysList, periodList)


    assert(blsArray.shape[:2] == (1,1))
    assert(np.all(blsArray==0))


def test_periodGeneration():
    periods = fbls.computePeriodList(1,5, [2220, 2140], 2)

    assert(periods[0] < 1.01)
    assert(periods[-1] >= 5)


def test_fbls1(plot=False):
    time, flux = makeTestData()

    duration_daysList = np.array([2,4,6,8]).astype(float) / 24.
    periodList = [2., 2.3, 2.4, 2.41, 2.42, 2.5, 3.0]
    blsArray = fbls.computeBlsForManyPeriods(time, flux, duration_daysList, periodList)

    index = np.argmax(blsArray)
    per, epc, dur=  fbls.getParamsOfIndex(blsArray, index,\
        duration_daysList, periodList)

    if plot:
        print per, epc+time[0], dur*24
        makePlots(time, flux, periodList, blsArray, per, epc)

    assert(per> 2.3 and per < 2.42)
    assert(int(dur*24) == 6)
    assert(epc >1.09 and epc < 1.11)

    return blsArray



def test_fbls2(plot=False):
    time,flux = loadData()

    duration_daysList = np.array([2,4,6,8]).astype(float) / 24.
    blsArray, periods = fbls.fBls(time, flux, [1,5], duration_daysList)

#    print periods[0], periods[-1]
    index = list(fbls.findBestPeak(blsArray))
    per, epc, dur=  fbls.getParamsOfIndex(blsArray, index,\
        duration_daysList, periods)

    if plot:
        print index
        print per, epc+time[0], dur*24
        makePlots(time, flux, periods, blsArray, per, epc)

    assert(per> 4.1 and per < 4.2)
    assert(int(dur*24) == 4)
#    print epc
    assert(epc >1.6 and epc < 1.7)

    return blsArray, periods, duration_daysList





def makePlots(time, flux, periods, blsArray, per, epc):
        mp.figure(1)
        mp.clf()
        mp.plot(time, flux, 'bo', mec="none")

        for i in range(17):
            mp.axvline(epc+np.min(time)+ i*per, color='r')

        mp.figure(2)
        mp.clf();
        bb = blsArray.max(axis=1)
        mp.imshow(bb, cmap=mp.cm.YlGnBu_r, interpolation="nearest", origin="bottom", aspect="auto")
#        mp.clim(vmin=0)
        mp.colorbar()

        mp.figure(3)
        mp.clf()
        bbb = bb.max(axis=1)
        mp.plot(periods, bbb, 'b-')
        mp.xlabel("Period")




if __name__ == "__main__":
    test_fbls2()