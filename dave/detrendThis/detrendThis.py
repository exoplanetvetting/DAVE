from __future__ import division, print_function
import logging
import numpy as np
# from martinsff import martinsff
import extract_lc
extract_lc = reload(extract_lc)
from astropy.stats import median_absolute_deviation as MAD


def twostepreplace(arr1, arr2, idx1, idx2, zpt):

    midx = np.arange(len(arr1))[idx1][zpt:][idx2]
    arr1[midx] = arr2

    return arr1


def detrendThat(time, flux, xpos, ypos, ferr=None, qflags=None,
                inpflag=None,
                ap=4.0,
                cadstep=300):

    """
    code to run self flat field on K2 data
    
    Keyword arguments:
    time -- time time array
    flux -- the array of brightnesses
    xpos -- the x-pixel position
    ypos -- the y-pixel position
    ferr -- flux error, will be assumed to be uniform if None
    qflags -- data quality flags
    inpflag -- flags to use to explude data, non-zero values are removed
    ap -- a random number
    """

    if ferr is None:
        ferr = np.ones_like(time)

    if inpflag is None:
        inpflag = np.ones_like(time)

    if qflags is None:
        qflags = np.zeros_like(time)

    # we are going to remove all the bad values
    # we should never pass out the time
    # this should never trigger
    badmask = np.isfinite(flux) * (inpflag > 0.)

    time = np.copy(time[badmask])
    flux = np.copy(flux[badmask])
    ferr = np.copy(ferr[badmask])
    xpos = np.copy(xpos[badmask])
    ypos = np.copy(ypos[badmask])

    flatlc = extract_lc.medfilt(time, flux, window=3)

    n_chunks = 6
    cadstep = np.int(np.floor(len(time) / n_chunks))  # 600
    zpt = len(time) % cadstep
    if zpt == cadstep:
        zpt = 0
    logging.info("%d pts %d chunks step=%d zpt=%d ap=%f ",
                 len(time), n_chunks, cadstep, zpt, ap)

    outflux, correction, thr_cad = extract_lc.run_C0_detrend(
        time, flatlc, xpos, ypos, cadstep=cadstep, skip=None)

    not_thr = ~thr_cad
    assert len(outflux) == len(correction)
    assert len(correction[not_thr]) == len(time[zpt:][not_thr])
    corflux = (flux[zpt:][not_thr] /
               np.median(flux[zpt:][not_thr]) /
               correction[not_thr])

    corflatflux = (flatlc[zpt:][not_thr] /
                   np.median(flatlc[zpt:][not_thr]) /
                   correction[not_thr])

# The 1.4826 and *4 factors make this similar to a 4-sigma cut.
# this will never trigger
    mad_cut = 1.4826*MAD(corflatflux-1.)*4
    keep = np.abs(corflatflux-1.) < mad_cut
    
    newflags = np.zeros(len(flux), dtype=bool)
    newflags[zpt:][not_thr] = 1

    return corflux, corflatflux, correction[not_thr], newflags





