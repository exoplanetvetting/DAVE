#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 13:16:57 2018

@author: smullally
"""


from astropy.io import fits
from astropy.time import Time
import numpy as np

    

def writeDetrendSeriesFitsFile(clip):
    """
    Given a clipboard shelf write out the detrended light curve and the time
    to the first extension.
    The resulting FITS file mimics the Kepler DV time series files.
    """
    
    #General information we want to capture in this file 
    #That is stored in the clipboard
    
    epic = int(clip['value'])
    campaign = clip['config']['campaign']
    
    outputFileName= "%s/%s/epic%s-c%02u_%s_davefit.fits" % \
           (clip['config']['onepageBasename'],str(epic),str(epic),int(campaign),str(clip['config']['detrendType']))    
    
    fl = clip['detrend']['flags']
    time = clip['extract']['time']
 
    raw = clip['extract.rawLightcurve']
    flux = clip['detrend.flux_frac']
    trapsnr=clip['trapFit.snr']
    trapper=clip['trapFit.period_days']
    trapdur=clip['trapFit.duration_hrs']
    trapdepth=clip['trapFit.depth_frac']*1.0e6
    epoch=clip['trapFit.epoch_bkjd']
    trapModel=clip['trapFit.bestFitModel']
    spantime=[2454833+clip['extract.time'][0],2454833+clip['extract.time'][-1]]
    t=Time(spantime,format='jd',scale='tdb')
    
    #Need to check this phase calculation against how DV does it.
    #I want phase in the units of the period with transit
    #at 0.25 phase.
    phase = (time-epoch+trapper) % (trapper)
    
    
    c1 = fits.Column(name='TIME', array=time, format='D')
    c2 = fits.Column(name='PHASE',array=phase,format='D')
    c3 = fits.Column(name='LC_DETREND', array=flux, format='E')
    c4 = fits.Column(name='MODEL_INIT', array=trapModel, format='E')
    
    #Gather Primary Header Info
    phdr=fits.Header()
    phdr['MISSION'] = ('K2','Mission name')
    phdr['TELESCOP'] = ('Kepler', 'telescope')
    phdr['OBSMODE'] = ('long cadence', 'mode of observations')
    phdr['INSTRUME'] = ('Kepler Photometer','instrument')
    phdr['OBJECT'] = ('EPIC %u' % epic, 'string version of target id')
    phdr['EPICID'] = (epic, 'unique EPIC target identifier')
    phdr['CAMPAIGN'] = (campaign, 'K2 Campaign number')
    phdr['LCSOURC'] = (clip['config']['detrendType'], 'Lightcurve detrending source')
    
    hdr=fits.Header()
    hdr['TRAPPER'] = trapper
    hdr.comments['TRAPPER'] = "Trapezoidal Fit Period [days]"
    hdr['TRAPTZER'] = epoch
    hdr.comments['TRAPTZER'] = "Trapezoidal Fit Transit Time [BKJD]"
    hdr['TRAPDEPT'] = trapdepth
    hdr.comments['TRAPDEPT'] = "Trapezoidal Fit Depth [ppm]"
    hdr['TRAPSNR'] = (clip['trapFit.snr'],'transit signal-to-noise ratio')
    hdr['TRAPDUR'] = (clip['trapFit.duration_hrs'],'Trapeoidal Fit Duration [hrs]')
    hdr['TIMESYS'] = ('TDB', 'time system is barycentric JD')   
    hdr['BJDREFI'] = (2454833, 'integer part of BJD reference date')
    hdr['BJDREFF'] = (0, 'fraction of the day BJD reference date')
    hdr['TSTART'] = (clip.extract.time[0], 'observation start time in BJD-BJDREF')
    hdr['TEND'] = (clip.extract.time[-1], 'observation stop time in BJD-BJDREF')
    hdr['DATE-OBS'] = (t[0].isot,'TSTART as UTC calendar date')
    hdr['DATE-END'] = (t[1].isot,'TSTOP as UTC calendar date')
    

    fitstable=fits.BinTableHDU.from_columns([c1,c3,c2,c4],header=hdr,name="Signal1")
    
    primary = fits.PrimaryHDU(header=phdr)
    
    hdus=fits.HDUList(hdus=[primary,fitstable])
    
    hdus.writeto(outputFileName,overwrite=True)
    