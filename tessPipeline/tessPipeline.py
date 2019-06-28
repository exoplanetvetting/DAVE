#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Create tasks to run a vetting pipeline on tess data.
Created on Tue Nov 27 21:23:14 2018

@author: smullally

"""

from __future__ import division

from pdb import set_trace as debug

# from dave.lpp.lppTransform import computeLPPTransitMetric # save for later
import dave.tessPipeline.tessfunc as tessfunc
import dave.pipeline.clipboard as clipboard
import dave.vetting.ModShift as ModShift
from dave.pipeline.task import task
import dave.fileio.tpf as tpf
from astropy.io import fits
from scipy.signal import savgol_filter
import dave.vetting.RoboVet as RoboVet
import dave.misc.covar as covar
import os
import glob
import pandas as pd

import numpy as np

######################################

def runOne(config, returnClip=False):
    """Run the pipeline on a single target.

    Inputs:
    ------------
    k2id
        (int) Epic id of target to run on.

    config
        (dict) Dictionary of configuration parameters as created by, e.g
        loadMyConfiguration()

    Returns:
    ---------
    A clipboard containing the results.

    Notes:
    ---------
    Don't edit this function. The pipeline can recover gracefully from
    errors in any individual task, but an error in this function will crash
    the pipeline
    """

    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config


    #Check that all the tasks are properly defined
    print( "Checking tasks exist")
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        print( "Running %s" %(t))
        f = eval(t)
        clip = f(clip)

    if returnClip:
        return clip

######################################

@task
def serveTask(clip):

    source_ = clip['config.detrendType']

    sector = clip['config.sector']
    tic = clip['config.tic']
    planNum = clip['config.planetNum']
    localPath = clip['config.dvtLocalPath']

    if source_ == "tess_2min":

        dvt, hdr, tpf_, hdr_tpf = tessfunc.serve(sector, tic, planNum, localPath, source_)


#        print(hdr)
#        xxxx
#        print(hdr_tpf)
#        xxxx

        cube = tpf.getTargetPixelArrayFromFits(tpf_, hdr_tpf)

#        print(hdr)
#        print(dvt.shape)
#        print(cube.shape)
#        xxxx

        out = dict()
        out['time'] = dvt['TIME']
        out['cube'] = cube
        out['tpfHeader'] = hdr_tpf
        try:
            out['detrendFlux'] = dvt['LC_DETREND']
            flux_tmp = dvt['LC_DETREND']
            out['modelFlux'] = dvt['MODEL_INIT']
        except:
            out['detrendFlux'] = (dvt['PDCSAP_FLUX']/np.nanmedian(dvt['PDCSAP_FLUX'])) - 1.
            flux_tmp = (dvt['PDCSAP_FLUX']/np.nanmedian(dvt['PDCSAP_FLUX'])) - 1.
            out['modelFlux'] = (dvt['PDCSAP_FLUX']/np.nanmedian(dvt['PDCSAP_FLUX'])) - 1.

        out['flags'] = np.isnan(flux_tmp)

        par = dict()
        par['orbitalPeriod_days'] = clip['config.period']#hdr['TPERIOD']#
        par['epoch_btjd'] = clip['config.tepoch']#hdr['TEPOCH']#
        par['transitDepth_ppm'] = clip['config.tdepth']#hdr['TDEPTH']#
        par['transitDuration_hrs'] = clip['config.tdur']#hdr['TDUR']#

        clip['serve'] = out
        clip['serve.param'] = par

    #Enforce contract
        clip['serve.time']
        clip['serve.cube']
        clip['serve.detrendFlux']
        clip['serve.flags']
        clip['serve.modelFlux']
        clip['serve.param.orbitalPeriod_days']
        clip['serve.param.epoch_btjd']
        clip['serve.param.transitDepth_ppm']
        clip['serve.param.transitDuration_hrs']

        out = dict()
        out['flux_frac'] = flux_tmp
        out['flags'] = np.isnan(flux_tmp)
        clip['detrend'] = out
        clip['detrend.flux_frac']
        clip['detrend.flags']

        out = dict()
        out['period'] = clip['config.period']#hdr['TPERIOD']
        out['epoch'] = clip['config.tepoch']#hdr['TEPOCH']
        clip['bls'] = out
        clip['bls.period']
        clip['bls.epoch']

        out = dict()
        out['time'] = dvt['TIME']
        out['rawLightcurve'] = flux_tmp
        clip['extract'] = out
        clip['extract.time']
        clip['extract.rawLightcurve']


    elif source_ == "tess_FFI":

        lc, hdr, tpf_, hdr_tpf = tessfunc.serve(sector, tic, planNum, localPath, source_)
#        lc, hdr, tpf_, hdr_tpf = tessfunc. serve_TESS_FFI(sector, tic, planNum, localPath)

        cube = tpf_ #tpf.getTargetPixelArrayFromFits(tpf_, hdr_tpf)

        out = dict()
        out['time'] = lc['TIME']
        out['cube'] = cube
        out['tpfHeader'] = hdr_tpf
        out['detrendFlux'] = lc['PDCSAP_FLUX']
        out['flags'] = np.isnan(lc['PDCSAP_FLUX'])

        par = dict()
        par['orbitalPeriod_days'] = clip['config.period']
        par['epoch_btjd'] = clip['config.tepoch']
        par['transitDepth_ppm'] = clip['config.tdepth']
        par['transitDuration_hrs'] = clip['config.tdur']

        clip['serve'] = out
        clip['serve.param'] = par

    #Enforce contract
        clip['serve.time']
        clip['serve.cube']
        clip['serve.detrendFlux']
        clip['serve.flags']
        clip['serve.param.orbitalPeriod_days']
        clip['serve.param.epoch_btjd']
        clip['serve.param.transitDepth_ppm']
        clip['serve.param.transitDuration_hrs']

        out = dict()
        out['flux_frac'] = lc['PDCSAP_FLUX']
        out['flags'] = np.isnan(lc['PDCSAP_FLUX'])
        clip['detrend'] = out
        clip['detrend.flux_frac']
        clip['detrend.flags']

        out = dict()
        out['period'] = clip['config.period']#hdr['TPERIOD']
        out['epoch'] = clip['config.tepoch']#hdr['TEPOCH']
        clip['bls'] = out
        clip['bls.period']
        clip['bls.epoch']

        out = dict()
        out['time'] = lc['TIME']
        out['rawLightcurve'] = lc['PDCSAP_FLUX']
        clip['extract'] = out
        clip['extract.time']
        clip['extract.rawLightcurve']


    elif source_ == "eleanor":

        prefix, suffix = "hlsp_eleanor_tess_ffi_tic", "_s01_tess_v0.1.8_lc.fits"

        fits_fn = "%s%i%s" %(prefix, int(tic), suffix)
        fits_fn = os.path.join(localPath, fits_fn)

        hdu = fits.open(fits_fn)
#	print(hdu[1].columns)
#	header_ = hdu[0].header#hdu[1].header

        out = dict()
        out['tpfHeader'] = hdu[0].header
        out['cube'] = hdu[1].data['TPF']
        out['time'] = hdu[1].data['TIME']# - min(hdu[1].data['TIME'])
        out['detrendFlux'] = (hdu[1].data['CORR_FLUX']/np.nanmedian(hdu[1].data['CORR_FLUX'])) - 1.
        flags_ = np.isfinite(hdu[1].data['CORR_FLUX'])
        flags_[hdu[1].data['QUALITY'] == 0] = False
        out['flags'] = flags_

        par = dict()
        par['orbitalPeriod_days'] = clip['config.period']
        par['epoch_btjd'] = clip['config.tepoch']
        par['transitDepth_ppm'] = clip['config.tdepth']
        par['transitDuration_hrs'] = clip['config.tdur']

        clip['serve'] = out
        clip['serve.param'] = par

    	#Enforce contract
        clip['serve.time']
        clip['serve.cube']
        clip['serve.detrendFlux']
        clip['serve.flags']
        clip['serve.param.orbitalPeriod_days']
        clip['serve.param.epoch_btjd']
        clip['serve.param.transitDepth_ppm']
        clip['serve.param.transitDuration_hrs']

        out = dict()
        out['period'] = clip['config.period']
        out['epoch'] = clip['config.tepoch']
        clip['bls'] = out
        clip['bls.period']
        clip['bls.epoch']

        out = dict()
        out['time'] = hdu[1].data['TIME']# - min(hdu[1].data['TIME'])
        out['rawLightcurve'] = (hdu[1].data['RAW_FLUX']/np.nanmedian(hdu[1].data['RAW_FLUX'])) - 1.
        clip['extract'] = out
        clip['extract.time']
        clip['extract.rawLightcurve']


    elif source_ == "ryan":

# GET DATA
        def get_data_(tic, sector):
            file_path_ = "/Users/vkostov/Desktop/Ideas_etc/TESS/tesseract/ryan_LCs/data/"

            for fn_ in glob.glob(file_path_ + str(tic) + "_sector0" + str(sector) + "*.lc"):
                data = pd.read_csv(fn_, sep = ' ')
                data = np.array(data)
                time = data[:,0]
                mag_raw = data[:,1]
                flux = 10**(-0.4*data[:,1])

            return time, flux
#
#
#DETREND DATA
        def detrend_data_(time, flux, sector):
            x = time
            y = flux/np.nanmedian(flux)
            mu = np.median(y)
            y = (y / mu - 1)# * 1e3

# Identify outliers
            m = np.ones(len(y), dtype=bool)
            for i in range(10):
                y_prime = np.interp(x, x[m], y[m])
                smooth = savgol_filter(y_prime, 121, polyorder=3)
                resid = y - smooth
                sigma = np.sqrt(np.mean(resid**2))
                m0 = np.abs(resid) < 3*sigma
                if m.sum() == m0.sum():
                    m = m0
                    break
                m = m0

# Only discard positive outliers
            if sector < 3:
                m = (resid < 3*sigma) & (resid > -5*sigma)
            else:
                m = (resid < 2*sigma) & (resid > -5*sigma)
# REMOVE DATA AROUND DATA GAPS/DISCONTINUITIES!!!!
            if sector < 3:
                m[(x - np.min(x) >= 12.25) & (x - np.min(x) <= 15.3)] = False
            elif sector == 3:
                m[(x - np.min(x) >= 11.5) & (x - np.min(x) <= 14.5)] = False
                m[(x - np.min(x) <= 2.25)] = False
                m[(x - np.min(x) > 26.5)] = False
            elif sector == 4:
                m[(x - np.min(x) >= 7.) & (x - np.min(x) <= 15.)] = False
                m[(x - np.min(x) >= 25.)] = False
#
# Shift the data so that it starts at t=0. This tends to make the fit
# better behaved since t0 covaries with period.
#    print(min(x), min(x[m]))
            x_ref = np.min(x[m]) if sector != 3 else np.min(x)
#    x -= x_ref

# Make sure that the data type is consistent
            x_bak_, y_bak_ = x, y
            x = np.ascontiguousarray(x[m], dtype=np.float64)
            y = np.ascontiguousarray(y[m], dtype=np.float64)
            smooth = np.ascontiguousarray(smooth[m], dtype=np.float64)

            x_detrend_ = x
            y_detrend_ = y - smooth

            return x_detrend_, y_detrend_, x_ref, x_bak_, y_bak_, smooth, m
#
#
# Get all sectors
        sector = [1,2,3,4]

        time, flux_raw, time_detrend_, flux_detrend_ = 0., 0., 0., 0.

        for jj in range(len(sector)):
            try:
                if sector[jj] == 1:
                    time_raw, flux_raw = get_data_(tic, sector[jj])
                    time = time_raw
                    time_detrend_, flux_detrend_, t0, x_bak_, y_bak_, smooth, m = detrend_data_(time, flux_raw, sector[jj])
                    time_detrend_ = time_detrend_ + 0.*t0
                else:
                    time_, flux_ = get_data_(tic, sector[jj])
                    time_tmp_, flux_tmp_, t0_, x_tmp_, y_tmp_, smooth, m_tmp = detrend_data_(time_, flux_, sector[jj])
                    time_detrend_ = np.hstack((time_detrend_,time_tmp_+0.*t0_))
                    flux_detrend_ = np.hstack((flux_detrend_, flux_tmp_))
            except:
                do_ = 'nothing'

        time_detrend_, flux_detrend_ = time_detrend_[1:], flux_detrend_[1:]
#	time_detrend_ -= min(time_detrend_)

        out = dict()
#	out['tpfHeader'] = hdu[0].header
        out['cube'] = np.zeros((10,5,5))#hdu[1].data['TPF']
        out['time'] = time_detrend_
        out['detrendFlux'] = flux_detrend_
        flags_ = np.zeros(len(flux_detrend_), dtype = bool)
        out['flags'] = flags_

        par = dict()
        par['orbitalPeriod_days'] = clip['config.period']
        par['epoch_btjd'] = clip['config.tepoch']
        par['transitDepth_ppm'] = clip['config.tdepth']
        par['transitDuration_hrs'] = clip['config.tdur']

        clip['serve'] = out
        clip['serve.param'] = par

    	#Enforce contract
        clip['serve.time']
        clip['serve.cube']
        clip['serve.detrendFlux']
        clip['serve.flags']
        clip['serve.param.orbitalPeriod_days']
        clip['serve.param.epoch_btjd']
        clip['serve.param.transitDepth_ppm']
        clip['serve.param.transitDuration_hrs']

        out = dict()
        out['period'] = clip['config.period']
        out['epoch'] = clip['config.tepoch']
        clip['bls'] = out
        clip['bls.period']
        clip['bls.epoch']

        out = dict()
        out['time'] = time_detrend_
        out['rawLightcurve'] = flux_detrend_
        clip['extract'] = out
        clip['extract.time']
        clip['extract.rawLightcurve']

        out = dict()
        out['time'] = time_detrend_
        out['flux_frac'] = flux_detrend_
        out['flags'] = np.zeros(len(flux_detrend_), dtype = bool)

        clip['detrend'] = out
        clip['detrend.flux_frac']
        clip['detrend.flags']

    return clip

######################################

@task
def detrendTask(clip):

    out = dict()

    x = clip['serve.time']
    y = clip['serve.detrendFlux']/np.nanmedian(clip['serve.detrendFlux'])
    flags = clip['serve.flags']

    mu = np.median(y)
    y = (y / mu - 1)

# Identify outliers
    m = np.ones(len(y), dtype=bool)
    for i in range(10):
        y_prime = np.interp(x, x[m], y[m])
        smooth = savgol_filter(y_prime, 121, polyorder=3)
        resid = y - smooth
        sigma = np.sqrt(np.mean(resid**2))
        m0 = np.abs(resid) < 3*sigma
        if m.sum() == m0.sum():
            m = m0
            break
        m = m0

# Only discard positive outliers
    m = (resid < 3*sigma) & (resid > -5*sigma)

# REMOVE DATA AROUND DATA GAPS/DISCONTINUITIES!!!!
    m[(x - np.min(x) >= 12.5) & (x - np.min(x) <= 15.5)] = False
#
# Shift the data so that it starts at t=0. This tends to make the fit
# better behaved since t0 covaries with period.
    x_ref = np.min(x[m]) if sector != 3 else np.min(x)

# Make sure that the data type is consistent
    x_bak_, y_bak_ = x, y
    x = np.ascontiguousarray(x[m], dtype=np.float64)
    y = np.ascontiguousarray(y[m], dtype=np.float64)
    smooth = np.ascontiguousarray(smooth[m], dtype=np.float64)

    x_detrend_ = x
    y_detrend_ = y - smooth

    detrendFlux = y - smooth

    out = dict()
    out['time'] = time_days
    out['flux_frac'] = detrendFlux
    out['flags'] = ~m

    clip['detrend'] = out
    clip['detrend.flux_frac']
    clip['detrend.flags']

    return clip
######################################

@task
def blsTask(clip):
    from astropy.stats import BoxLeastSquares

    time_days = clip['serve.time']
    flux_norm = clip['serve.detrendFlux'] if clip['config.detrendType'] == "tess_2min" else clip['detrend.flux_frac']
    flags = clip['serve.flags'] if clip['config.detrendType'] == "tess_2min" else clip['detrend.flags']
    if clip['config.detrendType'] != "ryan":
    	flux_norm[flags] = 0

    min_period_ = 0.5
    max_period_ = int((0.8*(max(time_days) - min(time_days))))
    period_grid = np.exp(np.linspace(np.log(min_period_), np.log(max_period_), 50000))
    durations_ = 0.05+0.05*np.linspace(0,4,5)

    time_days_ref_ = min(time_days)
    time_days -= time_days_ref_
#    print(np.nanmean(flux_norm), max(flux_norm), min(flux_norm))
#    print(min(time_days))
#    xxxxx

    bls = BoxLeastSquares(time_days, flux_norm)
    bls_power = bls.power(period_grid, durations_, oversample=20)

    index = np.argmax(bls_power.power)
    period = bls_power.period[index]
    epoch = bls_power.transit_time[index]
    depth = bls_power.depth[index]
    duration_hrs = 24.*bls_power.duration[index]
    bls_model = bls.model(time_days, period, duration_hrs/24., epoch)

    out = clipboard.Clipboard()
    out['period'] = period#clip['serve.param.orbitalPeriod_days']#
    out['epoch'] = epoch + time_days_ref_#clip['serve.param.epoch_btjd']#
    out['duration'] = duration_hrs#clip['serve.param.transitDuration_hrs']#
    out['model'] = bls_model
    out['depth'] = depth * 1e6#clip['serve.param.transitDepth_ppm']/1e6#

    clip['bls'] = out

    ##Enforce contract
    clip['bls.period']
    clip['bls.epoch']
    clip['bls.duration']
    clip['bls.model']
    clip['bls.depth']

    return clip
######################################

import dave.trapezoidFit.estimateSnr as tf
@task
def trapezoidFitTask(clip):

    time_days = clip['serve.time']
    flux_norm = clip['serve.detrendFlux'] if clip['config.detrendType'] == "tess_2min" else clip['detrend.flux_frac']
    flags = clip['serve.flags'] if clip['config.detrendType'] == "tess_2min" else clip['detrend.flags']

    period_days = clip['bls.period']#clip['serve.param.orbitalPeriod_days']#
    duration_hrs = clip['bls.duration']#clip['serve.param.transitDuration_hrs']#
    phase_bkjd = clip['bls.epoch']#clip['serve.param.epoch_btjd']#
    depth_frac = clip['bls.depth']#clip['serve.param.transitDepth_ppm']/1e6#

#    print(period_days, duration_hrs, phase_bkjd, depth_frac)

    #We don't know these values.
    unc = np.ones_like(flux_norm)
    unc[flags] = 1e99
    if clip['config.detrendType'] != "ryan":
    	flux_norm[flags] = 0

    assert(np.all(np.isfinite(time_days[~flags])))
    assert(np.all(np.isfinite(flux_norm[~flags])))

#    print(time_days[0:10], clip['serve.param.epoch_btjd'])
#    print(flux_norm[0:10])
#    print(flags[0:10])

#    time_days -= min(time_days)

    out = tf.getSnrOfTransit(time_days, flux_norm, unc, flags, \
        period_days, phase_bkjd, duration_hrs, depth_frac)

    assert(len(time_days) == len(out['bestFitModel']))

    clip['trapFit'] = out
    clip['trapFit.period_days']
    clip['trapFit.epoch_bkjd']
    clip['trapFit.duration_hrs']
    clip['trapFit.ingress_hrs']
    clip['trapFit.depth_frac']
    clip['trapFit.bestFitModel']
    clip['trapFit.snr']

#    print(clip['trapFit'])
#    print()
#    print(clip['bls'])

    return clip
######################################

@task
def modshiftTask(clip):

    time = clip['serve.time']

    if clip['config.detrendType'] == "tess_2min":
        flux = clip['serve.detrendFlux']
        model = clip['serve.modelFlux']
        model = clip['trapFit.bestFitModel']
    elif clip['config.detrendType'] == "eleanor":
        flux = clip['detrend.flux_frac']
        model = clip['trapFit.bestFitModel']
    elif clip['config.detrendType'] == "ryan":
        flux = clip['serve.detrendFlux']
        model = clip['trapFit.bestFitModel']

    period_days = clip['trapFit.period_days']#clip['serve.param.orbitalPeriod_days']
    epoch_btjd = clip['trapFit.epoch_bkjd']#clip['serve.param.epoch_btjd']

    #Filter out nans
    idx = np.isnan(time) | np.isnan(flux) | np.isnan(model)
    time = time[~idx]
    flux = flux[~idx]
    model = model[~idx]

    tic = clip['config.tic']
    basePath = clip['config.modshiftBasename']
    ticStr = "%016i" %(tic)
    basename = tessfunc.getOutputBasename(basePath, ticStr)

    # Name that will go in title of modshift plot
    objectname = "TIC %012i" % (tic)

    modplotint = 1  # Change to 0 or anything besides 1 to not have modshift produce plot
    plotname = "%s-%02i-%04i" % (basename, np.round(period_days*10), np.round(epoch_btjd))

    out = ModShift.runModShift(time, flux, model, plotname, objectname, period_days, epoch_btjd, modplotint)
    clip['modshift'] = out

    #Enforce contract
    clip['modshift.mod_Fred']
    clip['modshift.mod_ph_pri']
    clip['modshift.mod_secdepth']
    clip['modshift.mod_sig_pri']

    return clip
######################################

#from dave.lpp.loadLppData import MapInfo # save for later
@task
def lppMetricTask(clip):

    class clipToLppInputClass(object):

        def __init__(self, clip):
            """
            create a TCE class from the clipboard info
            """
            self.time=clip['serve.time']
            self.tzero=clip['serve.param.epoch_btjd']
            self.dur=clip['serve.param.transitDuration_hrs']
            self.period=clip['serve.param.orbitalPeriod_days']
            self.flux=clip['serve.detrendFlux']
            self.mes = 10.

    data=clipToLppInputClass(clip)
    mapInfoFile = clip['config.lppMapFile']
    mapInfoObj = MapInfo(mapInfoFile)
    normTLpp,rawTLpp,transformedTransit=computeLPPTransitMetric(data,mapInfoObj)

    out = dict()
    out['TLpp'] = normTLpp
    out['TLpp_raw'] = rawTLpp

    clip['lpp'] = out

    #Enforce contract
    clip['lpp.TLpp']

    return clip
######################################

from dave.tessPipeline.sweet import runSweetTest
@task
def sweetTask(clip):
    time = clip['serve.time']
    flux = clip['serve.detrendFlux'] if clip['config.detrendType'] == "tess_2min" else clip['detrend.flux_frac']
    period_days = clip['serve.param.orbitalPeriod_days']
    epoch_btjd = clip['serve.param.epoch_btjd']
    duration_hrs = clip['serve.param.transitDuration_hrs']

    idx = np.isnan(time) | np.isnan(flux)
    time = time[~idx]
    flux = flux[~idx]

    duration_days = duration_hrs / 24.
    result = runSweetTest(time, flux, period_days, epoch_btjd, duration_days)

    clip['sweet'] = result

    #Enforce contract
    clip['sweet.msg']
    clip['sweet.amp']

    return clip
######################################

#Won't work until serve loads up a TPF file
from dave.tessPipeline.pertransitcentroids import measurePerTransitCentroids
@task
def centroidsTask(clip):

    time = clip['serve.time']
    cube = clip['serve.cube']
    period_days = clip['trapFit.period_days']#clip['serve.param.orbitalPeriod_days']
    epoch_btjd = clip['trapFit.epoch_bkjd']#clip['serve.param.epoch_btjd']
    duration_hrs = clip['trapFit.duration_hrs']#clip['serve.param.transitDuration_hrs']

    duration_days = duration_hrs / 24.

#    print(period_days, epoch_btjd, duration_days)
#    xxxx

    res = measurePerTransitCentroids(time, cube, period_days, epoch_btjd, duration_days, plotFilePattern=None)

    res['method'] = "Fast Gaussian PSF fitting"
    clip['diffImgCentroids'] = res

    #Enforce contract
    clip['diffImgCentroids.results']
    return clip
######################################

@task
def vbkPsfCentroidsTask(clip):

    def raw_moment(data, iord, jord):
        nrows, ncols = data.shape
        y, x = np.mgrid[:nrows, :ncols]
        data = data * x**iord * y**jord
        return data.sum()

    def intertial_axis(data):
        data_sum = data.sum()
        m10 = raw_moment(data, 1, 0)
        m01 = raw_moment(data, 0, 1)
        x_bar = m10 / data_sum
        y_bar = m01 / data_sum
        u11 = (raw_moment(data, 1, 1) - x_bar * m01) / data_sum
        u20 = (raw_moment(data, 2, 0) - x_bar * m10) / data_sum
        u02 = (raw_moment(data, 0, 2) - y_bar * m01) / data_sum
        cov = np.array([[u20, u11], [u11, u02]])
        return x_bar, y_bar, cov

    def make_lines(eigvals, eigvecs, mean, i):
        std = np.sqrt(eigvals[i])
        vec = 3 * std * eigvecs[:,i] / np.hypot(*eigvecs[:,i])
        x, y = np.vstack((mean-vec, mean, mean+vec)).T
        return x, y

    np.set_printoptions(precision = 2, suppress = True)

    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    flags = clip['detrend.flags']

    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    duration_hrs = clip['trapFit.duration_hrs']
    epic = clip['value']

    cube = clip['serve.cube']
#    cube_BKG = clip['serve.cube_BKG']

    hdr_ = clip['serve.tpfHeader']
#    print hdr_

    if clip['config.detrendType'] == "tess":

    	col_zero_, row_zero_ = int(hdr_['1CRV4P']), int(hdr_['2CRV4P'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['1CRPX4']), row_zero_ + int(hdr_['2CRPX4'])

    if clip['config.detrendType'] == "eleanor":

    	col_zero_, row_zero_ = int(hdr_['CRPIX1']), int(hdr_['CRPIX2'])
    	epic_Col, epic_Row = col_zero_ + int(hdr_['TPF_W']), row_zero_ + int(hdr_['TPF_H'])


    inTransitIndices = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, duration_hrs/24., flags=flags)

    oot_cadence_ = np.where((inTransitIndices == False) & (flags == False))
    oot_cadence = np.asarray(oot_cadence_).flatten()

    itr_cadence_ = np.where(inTransitIndices == True)
    itr_cadence = np.asarray(itr_cadence_).flatten()
#
#
# GET IN-TRANSIT, BEFORE TRANSIT, AND AFTER TRANSIT CADENCES ONLY
#
#
    transit_number_ = 1
    transits_ = [itr_cadence[0]]
    tmp_idx_, = np.where(oot_cadence < transits_[0])
    no_transits_ = oot_cadence[tmp_idx_[-1]]

    for ii in range(1,len(itr_cadence)):
        if itr_cadence[ii] - itr_cadence[ii-1] > 10:
            transit_number_ += 1
            transits_ = np.hstack((transits_ , itr_cadence[ii-1], itr_cadence[ii]))
    transits_ = np.hstack((transits_ , itr_cadence[-1]))
#
#
# CALCULATE CENTROIDS
#
#
    use_before_after_images_per_transit_ = 'y'

    if use_before_after_images_per_transit_ == 'y':

        itrCol, itrRow, itr_cov = [], [], []
        ootCol, ootRow, oot_cov = [], [], []
        diffCol, diffRow, diff_cov = [], [], []

        cube_back_ = cube
        if clip['config.detrendType'] == "tess":
            cube = cube[:,epic_Col-col_zero_-3:epic_Col-col_zero_+3, epic_Row-row_zero_-3:epic_Row-row_zero_+3]
        if clip['config.detrendType'] == "eleanor":
            cube = cube[:,epic_Col-col_zero_-10:epic_Col-col_zero_+10, epic_Row-row_zero_-10:epic_Row-row_zero_+10]

        ss_ = cube.shape
        itr_mean_cube_ = np.zeros((transit_number_, ss_[1], ss_[2]))
        oot_mean_cube_ = np.zeros((transit_number_, ss_[1], ss_[2]))
        diff_mean_cube_ = np.zeros((transit_number_, ss_[1], ss_[2]))

        for ii in range(transit_number_):

            number_of_cadences_in_transit_ = transits_[2*ii+1] - transits_[2*ii]

            idx_in_transit_ = np.linspace(transits_[2*ii], transits_[2*ii+1], int(number_of_cadences_in_transit_+1))
            idx_in_transit = [int(aa) for aa in idx_in_transit_]

            idx_before_, = np.where(oot_cadence < transits_[2*ii])
            idx_before = oot_cadence[idx_before_[-1-number_of_cadences_in_transit_:]]

            idx_after_, = np.where(oot_cadence > transits_[2*ii+1])
            idx_after = oot_cadence[idx_after_[0:number_of_cadences_in_transit_+1]]

            itr_mean_img_by_transit_ = np.nanmean(cube[idx_in_transit,:,:], axis = 0)
            before_tr_mean_img_by_transit_ = np.nanmean(cube[idx_before,:,:], axis = 0)
            after_tr_mean_img_by_transit_ = np.nanmean(cube[idx_after,:,:], axis = 0)

            oot_mean_img_by_transit_ = 0.5*(before_tr_mean_img_by_transit_ + after_tr_mean_img_by_transit_)

            diff_mean_img_by_transit_ = oot_mean_img_by_transit_ - itr_mean_img_by_transit_
#		diff_mean_img_by_transit_ = diff_mean_img_by_transit_ - np.min(diff_mean_img_by_transit_).flatten()

            itrCol_by_transit_, itrRow_by_transit_, itr_cov_by_transit_ = intertial_axis(itr_mean_img_by_transit_)
            ootCol_by_transit_, ootRow_by_transit_, oot_cov_by_transit_ = intertial_axis(oot_mean_img_by_transit_)
            diffCol_by_transit_, diffRow_by_transit_, diff_cov_by_transit_ = intertial_axis(diff_mean_img_by_transit_)

            itrCol, itrRow = np.hstack((itrCol, itrCol_by_transit_)), np.hstack((itrRow, itrRow_by_transit_))
            ootCol, ootRow = np.hstack((ootCol, ootCol_by_transit_)), np.hstack((ootRow, ootRow_by_transit_))
            diffCol, diffRow = np.hstack((diffCol, diffCol_by_transit_)), np.hstack((diffRow, diffRow_by_transit_))

            itr_mean_cube_[ii,:,:] = itr_mean_img_by_transit_
            oot_mean_cube_[ii,:,:] = oot_mean_img_by_transit_
            diff_mean_cube_[ii,:,:] = diff_mean_img_by_transit_

#		else:
#			print 'BAD transit =', str(ii)
#			continue

        itr_mean_img_ = np.nanmean(itr_mean_cube_, axis = 0)
        oot_mean_img_ = np.nanmean(oot_mean_cube_, axis = 0)
        diff_mean_img_ = oot_mean_img_ - itr_mean_img_#np.nanmedian(diff_mean_cube_, axis = 0)

#	print itrRow + 1*row_zero_, itrCol + 1*col_zero_
#    	print ootRow + 1*row_zero_, ootCol + 1*col_zero_
#    	print diffRow + 1*row_zero_, diffCol + 1*col_zero_
#	print epic_Row, epic_Col
#	print row_zero_, col_zero_
#	xxxxx

    else:
        itr_mean_img_ = np.nanmean(cube[itr_cadence_], axis = 0)
        oot_mean_img_ = np.nanmean(cube[oot_cadence_], axis = 0)
        diff_mean_img_ = oot_mean_img_ - itr_mean_img_

        itrCol, itrRow, itr_cov = intertial_axis(itr_mean_img_)
        ootCol, ootRow, oot_cov = intertial_axis(oot_mean_img_)
        diffCol, diffRow, diff_cov = intertial_axis(diff_mean_img_)
#    	itrCol, itrRow = 1*col_zero_ + itrCol, 1*row_zero_ + itrRow
#    	ootCol, ootRow = 1*col_zero_ + ootCol, 1*row_zero_ + ootRow
#    	diffCol, diffRow = 1*col_zero_ + diffCol, 1*row_zero_ + diffRow

        print(itrRow + 1*row_zero_, itrCol + 1*col_zero_)
        print(ootRow + 1*row_zero_, ootCol + 1*col_zero_)
        print(diffRow + 1*row_zero_, diffCol + 1*col_zero_)

    return clip

@task
def dispositionTask(clip):
    """Decide whether an event is a planet candidate or not

    """
    out = clipboard.Clipboard(isSignificantEvent=True, isCandidate=True, reasonForFail="None")

# CENTROID VET
    minProbForFail = 0.99
    centVet = {'Warning':"None"}
    try:
    	centroids = clip['diffImgCentroids.results']

    	ootCol_prf, ootRow_prf = np.mean([centroids[:,0],centroids[:,4]], axis = 0), np.mean([centroids[:,1],centroids[:,5]], axis = 0)
    	diffCol_prf, diffRow_prf = centroids[:,2], centroids[:,3]
    	diffC, diffR = (ootCol_prf - diffCol_prf), (ootRow_prf - diffRow_prf)

    	prob, chisq = covar.computeProbabilityOfObservedOffset(diffC, diffR)
    except ValueError as e:
        centVet['Warning'] = "Probability not computed: %s" %(e)
        prob = 0
        chisq = 0

    centVet['probabilityOfOffset'] = prob
    centVet['chiSquaredOfOffset'] = chisq
    centVet['numTransitsWithCentroids'] = int( np.sum(centroids[:,0] > 0))

    centVet['isCentroidFail'] = False
    if np.isfinite(prob):
        if prob > minProbForFail:
            centVet['isCentroidFail'] = True
            out['isCandidate'] = False
            out['reasonForFail'] = "Centroid offset probability is %.1e" %(prob)

    out['centroidVet'] = centVet

# FLUX VET
    try:
        fluxVet = RoboVet.roboVet(clip.modshift)
        assert(fluxVet['disp'] in ["candidate", "false positive"])

    except:
        fluxVet ={}
        fluxVet['disp'] = 'candidate'
        fluxVet['not_trans_like'] = 0
        fluxVet['sig_sec'] = 0
        fluxVet['comments'] = 'NO_MODSHIFT'

    out['fluxVet'] = fluxVet

    lpp_th = 4  #Threshold for LPP
    if clip.lpp.TLpp > lpp_th:
        fluxVet['disp'] = 'false positive'
        fluxVet['comments']= fluxVet['comments'] + "-LPP_TOO_HIGH"
        fluxVet['not_trans_like']=1

    sweet_th = 3.5
    if (clip.sweet.amp[0,-1] > sweet_th) | \
       (clip.sweet.amp[1,-1] >sweet_th) | \
       (clip.sweet.amp[2,-1] > sweet_th):
           fluxVet['disp']='false positive'
           fluxVet['comments']= fluxVet['comments'] + "-SWEET_FAIL"
           fluxVet['not_trans_like']=1

    if fluxVet['disp'] == "false positive":
        out['isCandidate'] = False
        out['reasonForFail'] = fluxVet['comments']

        if fluxVet['not_trans_like'] > 0:
            out['isSignificantEvent'] = False

#    clip['fluxVet']= fluxVet
#    clip['disposition']= fluxVet
    clip['disposition'] = out
    clip['disposition.isSignificantEvent']
    clip['disposition.isCandidate']
    clip['disposition.reasonForFail']

    return clip
