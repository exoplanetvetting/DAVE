#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Create tasks to run a vetting pipeline on tess data.
Created on Tue Nov 27 21:23:14 2018

@author: smullally

"""

from __future__ import division

from pdb import set_trace as debug

from dave.lpp.newlpp.lppTransform import computeLPPTransitMetric
import dave.tessPipeline.tessfunc as tessfunc
import dave.pipeline.clipboard as clipboard
import dave.vetting.ModShift as ModShift
from dave.pipeline.task import task

import numpy as np

 

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


@task
def serveTask(clip):
    
    sector = clip['config.sector']
    tic = clip['config.tic']
    planNum = clip['config.planetNum']
    localPath = clip['config.dvtLocalPath']
    
    dvt, hdr = tessfunc.serve(sector, tic, planNum, localPath)

    out = dict()
    out['time'] = dvt['TIME']
    out['detrendFlux'] = dvt['LC_DETREND']
    out['modelFlux'] = dvt['MODEL_INIT']
    
    par = dict()
    par['orbitalPeriod_days'] = hdr['TPERIOD']
    par['epoch_btjd'] = hdr['TEPOCH']
    par['transitDepth_ppm'] = hdr['TDEPTH']
    par['transitDuration_hrs'] = hdr['TDUR']
    
    clip['serve'] = out
    clip['serve.param'] = par
 
    #Enforce contract
    clip['serve.time']
    clip['serve.detrendFlux']
    clip['serve.modelFlux']
    clip['serve.param.orbitalPeriod_days']
    clip['serve.param.epoch_btjd']
    clip['serve.param.transitDepth_ppm']
    clip['serve.param.transitDuration_hrs']     

    return clip




from dave.lpp.newlpp.loadLppData import MapInfo
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


@task
def modshiftTask(clip):
    
    time = clip['serve.time']
    flux = clip['serve.detrendFlux']
    model = clip['serve.modelFlux']
    period_days = clip['serve.param.orbitalPeriod_days']
    epoch_btjd = clip['serve.param.epoch_btjd']

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

    modplotint = 0  # Change to 0 or anything besides 1 to not have modshift produce plot
    plotname = "%s-%02i-%04i" % (basename, np.round(period_days*10), 
                                    np.round(epoch_btjd))

    out = ModShift.runModShift(time, flux, model, \
        plotname, objectname, period_days, epoch_btjd, modplotint)
    clip['modshift'] = out

    #Enforce contract
    clip['modshift.mod_Fred']
    clip['modshift.mod_ph_pri']
    clip['modshift.mod_secdepth']
    clip['modshift.mod_sig_pri']
    return clip


from dave.tessPipeline.sweet import runSweetTest
@task 
def sweetTask(clip):
    time = clip['serve.time']
    flux = clip['serve.detrendFlux']
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


#Won't work until serve loads up a TPF file
from dave.tessPipeline.pertransitcentroids import measurePerTransitCentroids
@task
def centroidsTask(clip):
    
    time = clip['serve.time']
    cube = clip['serve.cube']
    period_days = clip['serve.param.orbitalPeriod_days']
    epoch_btjd = clip['serve.param.epoch_btjd']
    duration_hrs = clip['serve.param.transitDuration_hrs']     
    
    duration_days = duration_hrs / 24.
    res = measurePerTransitCentroids(time, cube, period_days, epoch_btjd, 
                                     duration_days, plotFilePattern=None)

    res['method'] = "Fast Gaussian PSF fitting"
    clip['diffImgCentroids'] = res
    
    #Enforce contract
    clip['diffImgCentroids.results']
    return clip
