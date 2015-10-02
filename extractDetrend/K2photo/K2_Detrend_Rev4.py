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
import glob
from numpy import transpose

def K2_DetrendRev4(fn,):
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
    # To use this in conjunction with the ACF included in this package, you must
    # comment out corflatflux, keep, flux, and flatflux, then run the script on 
    # your data set.
    compiled_table = astropy.table.Table()
    compiled_table['time'] = time[zpt:][not_thr]
    compiled_table['corflatflux'] = corflatflux
    compiled_table['corflux'] = corflux
    compiled_table['keep'] = keep
    compiled_table['flux'] = flux[zpt:][not_thr]
    compiled_table['flatflux'] = flatlc[zpt:][not_thr]
    compiled_table['xbar'] = xbar[zpt:][not_thr]
    compiled_table['ybar'] = ybar[zpt:][not_thr]
    
    # Generates the DANCe # for the file title.
    DANCe = fn.split('/')[-1].split('.')[0]
    # Saves the detrended data to a file. Need to update the path to the correct folder.
    np.savetxt('/Users/bryanmann/Documents/NASA_Kepler_2.0/All_Light_Curves/Lightcurves_RADec_ap2.0_BJM/Detrended_Data/Detrended_ACF_Data/Raw_Detrended_Data_' + DANCe + '.0.dat', compiled_table, delimiter=',')
    
def DetrenderRev4(file_pathway):
    # file_pathway should be the directory that contains all of the .dat files you 
    # wish to detrend. You may need to comment out or lightly edit certain lines 
    # in this code. This is the command you runt to detrend the data!
    x = glob.glob("%s/*.dat" % file_pathway)
    y = transpose(x)
    z = Table.read(y, format='ascii')
    for row in z:
        F = str(row[0])
        K2_DetrendRev4(F)