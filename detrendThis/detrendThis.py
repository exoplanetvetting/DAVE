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


def detrendThat(time, flux, xpos, ypos, ferr=None, qflags=None, inpflag=None,
                ap=4.0):

    if ferr is None:
        ferr = np.ones_like(time)

    if inpflag is None:
        inpflag = np.ones_like(time)

    if qflags is None:
        qflags = np.zeros_like(time)

    # we are going to remove all the bad values
    # we should never pass out the time
    badmask = np.isfinite(flux) * (inpflag > 0.)

    time = time[badmask]
    flux = flux[badmask]
    ferr = ferr[badmask]
    xpos = xpos[badmask]
    ypos = ypos[badmask]

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
    corflux = (flux[zpt:][not_thr] /
               np.median(flux[zpt:][not_thr]) /
               correction[not_thr])

    corflatflux = (flatlc[zpt:][not_thr] /
                   np.median(flatlc[zpt:][not_thr]) /
                   correction[not_thr])

# The 1.4826 and *4 factors make this similar to a 4-sigma cut.
    mad_cut = 1.4826*MAD(corflatflux-1.)*4
    keep = np.abs(corflatflux-1.) < mad_cut

    outcorflux = np.ones(len(badmask)) * np.nan
    outcorflatflux = np.ones(len(badmask)) * np.nan
    outcorrection = np.ones(len(badmask)) * np.nan

    outcorflux = twostepreplace(outcorflux, corflux, badmask, not_thr, zpt)
    outcorflatflux = twostepreplace(outcorflatflux, corflatflux, badmask,
                                    not_thr, zpt)
    outcorrection = twostepreplace(outcorrection, correction, badmask,
                                   not_thr, zpt)

    return outcorflux, outcorflatflux, outcorrection





