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

def K2_DetrendRev3(fn,):   
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

    # Adds the detrended data to an ascii table.
    compiled_table = astropy.table.Table()
    compiled_table['time'] = time[zpt:][not_thr]
    #compiled_table['corflatflux'] = corflatflux
    compiled_table['corflux'] = corflux
    #compiled_table['keep'] = keep
    #compiled_table['flux'] = flux[zpt:][not_thr]
    #compiled_table['flatflux'] = flatlc[zpt:][not_thr]
    compiled_table['xbar'] = xbar[zpt:][not_thr]
    compiled_table['ybar'] = ybar[zpt:][not_thr]
    
    # Generates the DANCe # for the file title.
    DANCe = fn.split('/')[-1].split('.')[0]
    # Saves the detrended data to a file. Need to update the path to the correct folder.
    np.savetxt('/Users/bryanmann/Documents/NASA_Kepler_2.0/All_Light_Curves/Lightcurves_RADec_ap2.0_BJM/Detrended_Data/Detrended_ACF_Data/Raw_Detrended_Data' + DANCe + '.0.dat', compiled_table, delimiter=',')
    
def DetrenderRev3(dat_list):
    # dat_list should be an ascii file with all of the files you want to detrend. This list can be created
    # using the latest revision of the dat_list_compiler script included in this python folder. You may need
    # to comment out certain lines for the desired function.
    x = Table.read(dat_list, format='ascii')
    for row in x:
        F = str(row[0])
        K2_DetrendRev3(F)