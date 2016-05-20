
from __future__ import division
import matplotlib.pyplot as mp
import numpy as np
import mastio2 as mastio

import noise

def main():
    
    mp.figure(1)
    mp.clf()
    ax = mp.subplot(211)
    #ax = mp.subplot(211)
    #x, y = testData()
    x, y, flag = jupiter1()
    cin = np.arange(len(x))
        
    #import pdb; pdb.set_trace()
    idx = noise.singlePointDifferenceSigmaClip(y, 4, initialClip=flag)
    mp.plot(cin[~flag],y[~flag], 'k.')
    mp.plot(cin[idx], y[idx], 'ro', ms=10)
    mp.plot(cin[flag], y[flag], 'gs')
    mp.axvline(499, color='b')
    
    mp.subplot(212, sharex=ax)
    mp.plot(cin, y - np.roll(y, -1), 'ko-')
    mp.axvline(499, color='b')
    
    outliers = np.where(np.bitwise_xor(idx, flag))[0]
    outliers = np.bitwise_xor(idx, flag).astype(int)
    
    mp.figure(2)
    mp.clf()
    mp.plot(np.convolve(outliers, outliers, mode='same'))
    
    
def testData():    
    period = 100.
    x = np.arange(1000)
    y = 4 + np.sin(2*np.pi/period*x) + .1*np.random.randn(len(x))
    
    y[400] += .4
    y[500] -= 1
    
    y[100:102] -= 4
    return x, y


def jupiter1():
    return loadMast(211418729, 5)


def jupiter2():
    return loadMast(210775710, 4)

def loadMast(kepid, campaign):    
    ar = mastio.K2Archive()
    fits = ar.getLongCadence(kepid, campaign)

    time = fits['TIME']
    flux = fits['PDCSAP_FLUX']
    flag = ~(np.isfinite(time) & np.isfinite(flux))
    
    return time, flux, flag