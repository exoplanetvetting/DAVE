# -*- coding: utf-8 -*-
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from martinsff import martinsff
import extract_lc
extract_lc = reload(extract_lc)
from astropy.stats import median_absolute_deviation as MAD
import astropy
from astropy.table import Table

def K2_DetrendRev2(fn,):    
    #DANCe = fn.split('/')[-1].split('_')[0]    
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
    #plt.plot(time[zpt:][not_thr],corflatflux,c='r')
    #plt.plot(time[zpt:][not_thr][keep],corflatflux[keep],c='b')

    #CDPP Calculations
    #thr = twohr_cdpp(flux)
    #shr = sixhr_cdpp(flux)

    #plt.plot(time,flatlc)
    #plt.figure(figsize=[15,6])
    #plt.plot(time,flux/np.median(flux),'bo',markersize=1)
    #plt.plot(time[zpt:][not_thr][keep],corflux[keep],marker='.')
    #plt.xlabel('Time [d]')
    #plt.ylabel('Normalized Flux')
    #plt.title(str(DANCe) + ' 6.5_Hr_CDPP_' + str(shr) + ' | 2.5_Hr_CCPD_' + str(thr))
    #plt.savefig(fn.split('.')[-1] + '.png')
    #plt.show()
    
    #plt.close('all')
    compiled_table = astropy.table.Table()
    compiled_table['time'] = time[zpt:][not_thr] #('time', time[zpt:][not_thr])
    compiled_table['corflatflux'] = corflatflux
    compiled_table['corflux'] = corflux
    compiled_table['keep'] = keep
    compiled_table['flux'] = flux[zpt:][not_thr]
    compiled_table['xbar'] = xbar[zpt:][not_thr]
    compiled_table['ybar'] = ybar[zpt:][not_thr]
    compiled_table['flatflux'] = flatlc[zpt:][not_thr]
        
    DANCe = fn.split('/')[-1].split('.')[-1]
    np.savetxt('/Users/bryanmann/Documents/NASA_Kepler_2.0/All_Light_Curves/Lightcurves_RADec_ap2.0_BJM/Detrended_Data/' + DANCe, compiled_tbl, fmt='ascii', delimiter=',')
    #print(compiled_table)
    #compiled_table.close

def DetrenderRev2(dat_list):
    x = Table.read(dat_list, format='ascii')
    for row in x:
        F = str(row[0]) #Table.read(str(row[0]), format='ascii')
        K2_DetrendRev2(F)