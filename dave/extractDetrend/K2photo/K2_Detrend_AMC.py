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
import astropy.io.ascii as ascii

def K2_DetrendRev4(fn,):
    time, flux, xbar, ybar = np.genfromtxt(fn, unpack = True, skip_header=1)    
    m1 = np.isfinite(flux)
    
    time = time[m1]
    flux = flux[m1]
    xbar = xbar[m1]
    ybar = ybar[m1]
    
    flatlc = extract_lc.medfilt(time,flux,window=3)
#    zpt = len(time)%300
    zpt = len(time)%327
    
    outflux, correction, thr_cad = extract_lc.run_C0_detrend(
#    time, flatlc, xbar, ybar, cadstep=300, skip=1828)
    time, flatlc, xbar, ybar, cadstep=327, skip=1828)
    
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
#    DANCe = fn.split('/')[-1].split('.')[0]
    substr = fn.split('/')[-1]
    end = substr.find('.dat')
    DANCe = substr[:end]

    newtable = {'Dates': time[zpt:][not_thr], 'Flux': flux[zpt:][not_thr], 'Corrflux': corflux, 'Xpos': xbar[zpt:][not_thr], 'Ypos': ybar[zpt:][not_thr]}
    ascii.write(newtable, '/Users/acody/Data/K2/Field_0/M35/Lightcurves_RADec_ap3.0_v4_AMCdetrend/'+DANCe + '_detrended.dat', names=['Dates','Flux', 'Corrflux','Xpos','Ypos'])

# Create some plots
    plt.clf()
    plt.subplot(211)
    plt.plot(time[zpt:][not_thr], flux[zpt:][not_thr]/np.median(flux[zpt:][not_thr]), 'bo', markersize=2)
    plt.xlabel('Time [d]')
    plt.ylabel('Flux/Median flux')
    plt.ylim((np.median(flux[zpt:][not_thr]/np.median(flux[zpt:][not_thr]))-4.5*np.std(flux[zpt:][not_thr]/np.median(flux[zpt:][not_thr])),np.median(flux[zpt:][not_thr]/np.median(flux[zpt:][not_thr]))+
      4.5*np.std(flux[zpt:][not_thr]/np.median(flux[zpt:][not_thr]))))
    plt.title = DANCe

    plt.subplot(212)
    plt.plot(time[zpt:][not_thr], corflux/np.median(corflux), 'bo', markersize=2)
    plt.xlabel('Time [d]')
    plt.ylabel('Flux/Median flux')
    plt.ylim((np.median(corflux/np.median(corflux))-4.5*np.std(corflux/np.median(corflux)),np.median(corflux/np.median(corflux))+4.5*np.std(corflux/np.median(corflux))))
 
    plt.savefig('/Users/acody/Data/K2/Field_0/M35/Lightcurves_RADec_ap3.0_v4_AMCdetrend/'+DANCe + '_detrended.png')
    
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
