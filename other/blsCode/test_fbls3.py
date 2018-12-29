# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 11:57:25 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np
import fbls3 as fbls


def test_blsStrength():

    size = 100
    width = 10
    depth = 1
    y = np.zeros(size)
    counts = np.ones(size)

    y[40:40+width] -= depth

    bls = fbls.computeBlsForOnePeriod(y, counts, width, 1)

    expectedBls = (depth*width) / (np.sqrt(width))

    assert( np.fabs( np.max(bls) - expectedBls) < .01*expectedBls)
    mp.plot(bls, 'ko')
    mp.axhline(expectedBls, color='r')


def test_durationAttenuation():
    attenuation = .9
    timespan_days = 90
    period = 5
    size = timespan_days*50
    sigma = 1
    depth = 20

    rSun = 6.96e8
    mSun = 2e30

    x = np.linspace(0, timespan_days, size)
    y = sigma * np.random.randn(size)

    #Compute a plausible duration for transit
    durObj = fbls.AstrophysicalDurationsSearch(.333, 5, rSun, mSun, attenuation)
    duration_days = np.mean(durObj(period))
    duration_cadences = int(duration_days/float(timespan_days) * size)

    #Add in some transits
    phase = 25
    for i in np.arange(0, size, size*period/float(timespan_days)):
        y[i+phase:i+phase+duration_cadences] -= depth

    periodList = [5]
    blsArray = fbls.computefBls(x, y,sigma, periodList, 10,
                   durObj)

    n0 = duration_cadences * timespan_days/period
    expectedMax = (depth*n0) / (sigma*np.sqrt(n0))
    assert( blsArray[0,2] > attenuation*expectedMax)


def test_phaseAttenuationPass():
    phaseAttenuation(63)


def test_phaseAttenuationFail():
    phaseAttenuation(29.117)


def phaseAttenuation(period):
    #29.771 causes this test to fail because the period is folded to split
    #the transit over the end of the folded array
#    period = 10 + np.random.rand()*90
    size = 10000
    sigma = 1
    depth = 20
    duration = 10


    x = np.linspace(0, size, size)
    y = sigma * np.random.randn(size)

    phase = 25
    for i in np.arange(0, size, period):
        y[i+phase:i+phase+duration] -= depth


    n0 = duration * size/period
    expectedMaxBls = (depth*n0) / (sigma*np.sqrt(n0))


    isPass = True
    mp.clf()
    mp.axhline(expectedMaxBls)
    for phaseOverResolution in np.arange(1, 8):
        nBins = phaseOverResolution * period / duration
        nBins = int(np.ceil(nBins))
        binnedFlux, counts = fbls.fastFoldAndBin(x, y, period, nBins)

        mp.plot(binnedFlux, 'k.-')

        transitWidth = phaseOverResolution
        expectedAttenuation = 1 - 1/ float(2*phaseOverResolution)

        if transitWidth > 0 and transitWidth < len(binnedFlux):
            bls = fbls.computeBlsForOnePeriod(binnedFlux, counts,
                                             transitWidth, sigma)

            obs = np.max(bls)
            expected = expectedMaxBls * expectedAttenuation
            if obs < expected:
                isPass = False

    assert(isPass)

def test_periodAttenuation():

    attenuation = .9
    timespan_days = 90
    period = 5
    size = timespan_days*50
    sigma = 1
    depth = 20

    rSun = 6.96e8
    mSun = 2e30

    x = np.linspace(0, timespan_days, size)
    y = sigma * np.random.randn(size)
    durObj = fbls.AstrophysicalDurationsSearch(.99, 1.01, rSun, mSun, attenuation)
    duration_days = durObj(period)[0]
    duration_cadences = int(duration_days/float(timespan_days) * size)
#    duration_cadences = 4

    #Add in some transits
    phase = 25
    for i in np.arange(0, size, size*period/float(timespan_days)):
        y[i+phase:i+phase+duration_cadences] -= depth

    periodList = fbls.computePeriodListForAstrophysicalDurations(4,6, \
                    timespan_days, rSun, mSun, attenuation)

    blsArray = fbls.computefBls(x, y,sigma, periodList, 10,
                   [duration_days])

    n0 = duration_cadences * timespan_days/period
    expectedMax = (depth*n0) / (sigma*np.sqrt(n0))
    assert( np.max(blsArray[:,2]) > attenuation*expectedMax)

#    mp.clf()
#    mp.plot(blsArray[:,0], blsArray[:,2])
#    mp.axhline(expectedMax, color='r')


def test_strengthVPhaseOverres():
    """Test that measured strength at a given period is not
    sensitive to the choice of phase over resolution

    """

    transitWidth = lambda x: [10]   #Function that always returns 10
    period = 200
    size = 1000
    sigma = 1

    x = np.arange(size)
    y = sigma * np.random.randn(size)

    #If the BLS is under resolved we get an underestimated signal strength.
    #So we start at a phase over resolution of 16 for the purpose of this
    #test. Your actual needs may tolerate a lower overresolution.
    overresList = 2** np.arange(4,10)
    signal = np.zeros_like(overresList)

    for i in range(len(overresList)):
        blsArray = fbls.computefBls(x, y,sigma, [period], overresList[i],
                       transitWidth)
        signal[i] = np.min(blsArray[:,1])

    signal -= np.mean(signal)

    assert np.max( np.fabs(signal) < 1e-3)



def test_strengthVDuration():
    """Test that measured signal strength is strongest at the
    true duration of the transit
    """

    period = 200
    size = 1000
    sigma = 1
    overres = 10
    injectedWidth = 10

    for kk in range(30):
        x = np.arange(size)
        y = sigma * np.random.randn(size)

        #Add a small signal
        y[400:400+injectedWidth] -= 20
        trialWidths = np.arange(2, 100, 4)
        signal = np.zeros_like(trialWidths)

        for i, tw in enumerate(trialWidths):
            durFunc = lambda x: [tw]
            blsArray = fbls.computefBls(x, y,sigma, [period], overres,
                       durFunc)
            signal[i] = np.min(blsArray[:,1])

        wh = np.argmin( np.fabs(trialWidths - injectedWidth) )
        assert(signal[wh] <  .95*np.min(signal))


def test_strengthVTimespan():
    """Test that signal strength is independent of timespan"""

    size = 10000
    sigma = 1
    overres = 10
    width = 10
    #BLS strength should be constant with timespan, but FAP should
    #decrease for a signal of a given strength when timespan increases
    #because the number of draws increases

    x = np.arange(size)
    y = sigma * np.random.randn(size)

    #Add a small signal
    y[10:10+width] -= 2

    for iter in range(100):
        nList = [30, 100, 300, 1000, 3000, 10000]
        signal = np.zeros_like(nList)
        for i,n in enumerate(nList):
            durFunc = lambda x: [width]
            blsArray = fbls.computefBls(x[:n], y[:n], sigma, [n], overres,
                       durFunc)

            signal[i] = np.max(blsArray[:,2])
        assert np.all( np.fabs(signal - np.mean(signal) ) < .001)


def test_strengthVScatter():
    """False alarm prob should be independent of the amount of scatter,
    all else being equal"""
    period = [200]
    size = 1000
    sigma = 1
    overres = 10
    width = [6]

    x = np.arange(size)
    noise = np.random.randn(size)

    mp.clf()
    sigmaList = [1,2,4,8]
    signal = np.zeros_like(sigmaList)

    for i,sigma in enumerate(sigmaList):
        y = sigma * noise
        blsArray = fbls.computefBls(x, y,sigma, [period], overres,
                                    width)
        signal[i] = np.min(blsArray[:,1])


    signal /= np.mean(signal) - 1
    assert np.all( np.fabs(signal) < .01)





def plot_computeGumbelLnFap():
    """Show continuous behaviour around the break between analytic
    and approximate methods (around sigma==8), and non-linear behaviour
    at low sigma and high values of nDraw"""

    mp.clf()
    sigma = np.linspace(1, 15, 100)
    for nDraws in [1e1, 1e2, 1e4, 1e6]:
        lnFap = fbls.computeGumbelLnFap(sigma, nDraws)
        print lnFap
        mp.plot(sigma, lnFap, '.-')







def plotPeriodTrend2():
    overres = 1
    transitWidth = 3
    size = 10000
    sigma = 1

    periodList = np.linspace(3, 600, 50)

    mp.figure(2)
    mp.clf()

    i = 0
    while True:
        blsArray = computeBlsOfNoise(size, sigma, transitWidth, periodList,
                                     overres)
        #Plot
        mp.figure(1)
        mp.clf()
        bins = np.linspace(-20, 0, 30)
        mp.hist(blsArray, bins=bins)
        nPeriod = len(periodList)
        thres1 = np.log(1/float(nPeriod))
        thres2 = np.log(1/float(size))

        mp.axvline(thres1, color='r')
        mp.axvline(thres2, color='g')
        mp.pause(.1)

        n1 = np.sum(blsArray < thres1)
        n2 = np.sum(blsArray < thres2)
#        print n1, n2

        mp.figure(2)
        print i, n1
        mp.plot(i, n1, 'ro')
        mp.plot(i, n2, 'go', ms=10, mec="none")
        i += 1
        mp.pause(.1)


def computeBlsOfNoise(size, sigma, trialTransitWidth, periodList, phaseOverres):

    #Compute lnFAP
    x = np.arange(size)
    y = sigma * np.random.randn(size)
    lnFapArray = fbls.computefBls(x, y, sigma, periodList,
                       [trialTransitWidth], phaseOverres)

    idx= lnFapArray > 1e3
    lnFapArray[idx] = 0
    #Compute BLS spectrum
    nPeriod = len(periodList)
    blsArray = np.zeros( nPeriod )
    for j in range(nPeriod):
        blsArray[j] = np.min( lnFapArray[j,:,:])

    return blsArray


