
"""
Created on Wed Nov 28 22:08:56 2018

Run the SWEET test to check for out of transit variations in the flux at the period of the 
the claimed transit.

BLS and BLS-like searches often mis-fire on variability. A periodic signal in the time-series
gets flagged as a transit. SWEET tests for this by fitting a sine wave to the out-of-transit
flux and searching for a statistically signifant amplitude.

The SWEET test was introduced in Thompson et al. (2017), with J. "Wiggles" Coughlin and F. Mullally as 
the original authors. This code may not exactly reproduce that code, as it's reproduced from memory


@author: fergal
"""

from __future__ import print_function
from __future__ import division

from dave.fileio.kplrfits import markTransitCadences
from pdb import set_trace as debug
import dave.misc.lsf  as lsf 
import numpy as np


def runSweetTest(time, flux, period, epoch, duration, threshold_sigma=3):
    result = computeSweetMetrics(time, flux, period, epoch, duration)

    msg = []
    if result[0, -1] > threshold_sigma:
        msg.append("WARN: SWEET test finds signal at HALF transit period")
    if result[1, -1] > threshold_sigma:
        msg.append("WARN: SWEET test finds signal at the transit period")
    if result[2, -1] > threshold_sigma:
        msg.append("WARN: SWEET test finds signal at TWICE the transit period")

    if len(msg) == 0:
        msg = ["OK: SWEET finds no out-of-transit variability at transit period"]
        
    out = dict()
    out['msg'] = "\n".join(msg)
    out['amp'] = result
    return out
    

def computeSweetMetrics(time, flux, period, epoch, duration):
    """
    period, epoch and duration all in same units (e.g days)
    """
    
    assert len(time) == len(flux)
    
    out = []
    idx = markTransitCadences(time, period, epoch, duration)
    for per in [period/2., period, 2*period]:
        phase = np.fmod(time - epoch + per, per)
    
        amp, ampUnc = SweetFitOotFlux(phase[~idx], flux[~idx])
        out.append( [amp, ampUnc, amp/ampUnc] )
        
    return np.array(out)

import matplotlib.pyplot as plt

def SweetFitOotFlux(phase, flux):
    period = np.max(phase)
    fObj = lsf.Lsf(phase, flux, None, 2, lsf.sine, period=period)
    
    #plt.figure()
    #plt.clf()
    #plt.plot(phase, flux, 'k.')
    #plt.plot(phase, fObj.getBestFitModel(), 'r.')

    amp, phase, ampUnc, phaseUnc = lsf.computeAmplitudeAndPhaseWithUnc(fObj)
    return amp, ampUnc
