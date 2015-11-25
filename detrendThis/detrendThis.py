from __future__ import division, print_function
import logging, os
import numpy as np
import matplotlib.pyplot as plt
from martinsff import martinsff
import extract_lc
extract_lc = reload(extract_lc)
from astropy.stats import median_absolute_deviation as MAD
import astropy
from astropy.table import Table
import glob
from numpy import transpose
import astropy.io.ascii as ascii

def detrendThat(time, flux, xpos, ypos, ferr=None, qflags=None, inpflag=None, ap=4.0):

    if ferr is None:
        ferr = np.ones_like(time)

    if inpflag is None:
        inpflag = np.ones_like(time)

    if qflags is None:
        qflags = np.zeros_like(time)

    
    #we are going to remove all the bad values
    # we should never pass out the time
    badmask = np.isfinite(flux) * (inpflag > 0.)

    time = time[badmask]
    flux = flux[badmask]
    ferr = ferr[badmask]
    xpos = xpos[badmask]
    ypos = ypos[badmask]

    flatlc = extract_lc.medfilt(time,flux,window=3)

    n_chunks = 6
    cadstep = np.int(np.floor(len(time) / n_chunks)) #600
    zpt = len(time) % cadstep
    if zpt==cadstep:
        zpt = 0
    logging.info("%d pts %d chunks step=%d zpt=%d ap=%f ", 
                 len(time), n_chunks, cadstep, zpt, ap)

    outflux, correction, thr_cad = extract_lc.run_C0_detrend(
        time, flatlc, xpos, ypos, cadstep=cadstep, skip=None)
    
    not_thr = ~thr_cad
    corflux = (flux[zpt:][not_thr]/
        np.median(flux[zpt:][not_thr])/
        correction[not_thr])

    corflatflux = (flatlc[zpt:][not_thr]/
        np.median(flatlc[zpt:][not_thr])/
        correction[not_thr])

# The 1.4826 and *4 factors make this similar to a 4-sigma cut.    
    mad_cut = 1.4826*MAD(corflatflux-1.)*4
    keep = np.abs(corflatflux-1.) < mad_cut


    outcorflux = np.zeros(len(badmask)) * np.nan
    outcorflatflux = np.zeros(len(badmask)) * np.nan
    outcorrection = np.zeros(len(badmask)) * np.nan

    outcorflux[badmask * not_thr] = corflux
    outcorflatflux[badmask * not_thr] = corflatflux
    outcorrection[badmask * not_thr] = correction

    return outcorflux, outcorflatflux, outcorrection





