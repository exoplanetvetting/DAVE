# -*- coding: utf-8 -*-
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from martinsff import martinsff
import extract_lc
#from cdpp.py import twohr_cdpp
#from cdpp.py import sixhr_cdpp
from astropy.stats import median_absolute_deviation as MAD

def K2_Detrend(fn,):    
    DANCe = fn.split('/')[-1].split('_')[0]    
    time, flux, xbar, ybar = np.genfromtxt(fn, unpack = True)    
    m1 = np.isfinite(flux)
    
    time = time[m1]
    flux = flux[m1]
    xbar = xbar[m1]
    ybar = ybar[m1]
    
    flatlc = extract_lc.medfilt(time,flux,window=3)
    zpt = len(time)%300
    
    outflux, correction, thr_cad = extract_lc.run_C0_detrend(
    time, flatlc, xbar, ybar, cadstep=300, skip=1828)
    
    not_thr = ~thr_cad
    corflux = (flux[zpt:][not_thr]/
        np.median(flux[zpt:][not_thr])/
        correction[not_thr])

    corflatflux = (flatlc[zpt:][not_thr]/
        np.median(flatlc[zpt:][not_thr])/
        correction[not_thr])
    
    mad_cut = 1.4826*MAD(corflatflux-1.)*4
    keep = np.abs(corflatflux-1.) < mad_cut
    plt.plot(time[zpt:][not_thr],corflatflux,c='r')
    plt.plot(time[zpt:][not_thr][keep],corflatflux[keep],c='b')

    #CDPP Calculations
    #thr = twohr_cdpp(flux)
    #shr = sixhr_cdpp(flux)

    plt.plot(time,flatlc)
    plt.figure(figsize=[15,6])
    plt.plot(time,flux/np.median(flux),'bo',markersize=1)
    plt.plot(time[zpt:][not_thr][keep],corflux[keep],marker='.')
    plt.xlabel('Time [d]')
    plt.ylabel('Normalized Flux')
    #plt.title(str(DANCe) + ' 6.5_Hr_CDPP_' + str(shr) + ' | 2.5_Hr_CCPD_' + str(thr))
    #plt.savefig(fn.split('.')[-1] + '.png')
    plt.show()
    
    #plt.close('all')