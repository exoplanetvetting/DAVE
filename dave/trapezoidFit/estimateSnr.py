# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:48:49 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import dave.trapezoidFit.trapfit as tf
import dave.misc.outliers as outliers
import numpy as np
import dave.fileio.kplrfits as kplrfits
import dave.misc.noise as noise

def getSnrOfTransit(time_days, flux_frac, unc, flags, period_days, phase_bkjd, \
    duration_hrs, depth_frac):
    """

    Inputs:
    ------------
    flux_frac
        (1d np array) Flux in fractional amplitude. The mean of this array
        should be zero for sane data. The trapezoid fit takes data with
        a mean of 1, the conversion is done within this function
    """

    idx = np.isfinite(time_days) & (np.isfinite(flux_frac))
    idx = idx & ~flags
    ioblk = tf.trapezoid_fit(time_days[idx], 1+flux_frac[idx], unc[idx], \
                  period_days, phase_bkjd, duration_hrs, \
                  1e6*depth_frac, fitTrialN=13, fitRegion=10.0, \
                  errorScale=1.0, debugLevel=0, \
                  sampleN=15)

    #Taken from trapfit.py around lines 434
    out = dict()
    out['period_days'] = period_days
    out['epoch_bkjd'] = ioblk.timezpt + ioblk.bestphysvals[0]
    out['duration_hrs'] = 24* ioblk.bestphysvals[2]
    out['ingress_hrs'] = out['duration_hrs'] * ioblk.bestphysvals[3]
    out['depth_frac'] = ioblk.bestphysvals[1]

    #compute modelat all input time values
    subSampleN= 15
    time_days[~np.isfinite(time_days)] = 0  #Hide the Nans from one_model
    assert(np.all(np.isfinite(time_days)))
    ioBlock = tf.trapezoid_model_onemodel(time_days, period_days, \
                out['epoch_bkjd'], 1e6*out['depth_frac'], out['duration_hrs'], \
                out['ingress_hrs'], subSampleN)
    
    out['bestFitModel'] = ioBlock.modellc - 1  #Want mean of zero
    out['snr'] = estimateSnr(time_days, flux_frac, flags, out['period_days'], \
                    out['epoch_bkjd'], out['duration_hrs'], out['depth_frac'])

    #out['bestFitModel'] = time_days*0
    #out['snr'] = -1
    return out



def estimateSnr(time, flux, flags, period_days, epoch_bkjd, \
    duration_hrs, depth_frac, nDurForClip=2):
    """Estimate the SNR of a transit.

    SNR is defined as (transit depth) / (rms scatter). Transit depth
    is an input parameter, snr is calculated with the Marshall method

    Inputs:
    -------------
    time, flux, flags
        (np 1d arrays) arrays of time flux and flag values. All flag
        values > 0 are treated as though they indicate bad data.

    period_days, epoch_bkjd, duration_hrs, depth_frac
        (floats) Parameters of transit

    Optional Inputs:
    ----------------
    nDurForClip
        Points within nDurForClip*duration_hrs around each transit are
        excluded from the estimate of noise.

    """

    dur_days = duration_hrs / 24.
    idx = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, \
        dur_days, nDurForClip, flags=flags)


    if np.all(idx):
        msg = "All cadences seem to be in or near transit: "
        msg += "Period %.1f Duration %.2f hrs" %(period_days, duration_hrs)
        raise ValueError(msg)

    idx |= flags > 0  #Remove data flagged as bad
    idx |= ~np.isfinite(time)  #Or otherwise NaN
    idx |= ~np.isfinite(flux)  #Or otherwise NaN
    idx |= outliers.indexOfOutliers(flux)    #Remove outliers

    #No good cadences for some reason.
    if np.all(idx):
        raise ValueError("No good cadences found for noise estimate. Check transit duration")

    assert( np.all(np.isfinite(flux[~idx])))

    expTime_days = np.median(np.diff(time[~idx]))
    duration_cadences = dur_days/expTime_days

    #Duration must be at least 4 cadences of sgCdpp will crash
    duration_cadences = max(duration_cadences, 4)
    rms = noise.computeSgCdpp_ppm(flux[~idx], duration_cadences)*1e-6

    idx = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, \
        dur_days, 1, flags=flags)
    nCadenceInTransit = np.sum(idx)

    return depth_frac/rms * np.sqrt(nCadenceInTransit)


def estimateScatterWithMarshallMethod(flux):
    """Estimate the typical scatter in a lightcurve.

    Uses the same method as Marshall (Mullally et al 2015 submitted)

    Inputs:
    ----------
    flux
        (np 1d array). Flux to measure scatter of


    Returns:
    ------------
    (float) scatter of data in the same units as in the input ``flux``


    Notes:
    ----------
    Algorithm is reasonably sensitive to outliers. For best results
    uses outlier rejection on your lightcurve before computing scatter.
    """

    diff= np.diff(flux)
    mean = np.mean(diff)
    mad = np.median(np.fabs(diff-mean))
    std = 1.4826*mad

    #std is the rms of the diff. std on single point
    #is 1/sqrt(2) of that value,
    return std/np.sqrt(2)