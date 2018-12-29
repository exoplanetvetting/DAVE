# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 06:28:16 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import numpy as np

import dave.pipeline.clipboard as clipboard
import dave.pipeline.pipeline as pl
import dave.pipeline.task as task
import dave.fileio.nca as nca

import dave.lpp.calcLPPoctave as lpp

import os


def runOne(k2id, config):

    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = k2id
    #clip['dataStorePath'] = config['dataStorePath']

    #Check that all the tasks are properly defined
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        f = eval(t)
        clip = f(clip)

    return clip


def loadMultiEventConfig():
    cfg = dict()
    cfg['debug'] = True
    cfg['campaign'] = 3
    cfg['timeout_sec'] = 120
    cfg['nPointsForMedianSmooth'] = 2*48

    cfg['maxNumEvents'] = 10
    cfg['blsMinPeriod'] = 0.5
    cfg['blsMaxPeriod'] = 30

    #Vetting parameters
    cfg['minSnrForDetection'] = 10.
    cfg['maxLppForTransit'] = 0.03
    cfg['minCentroidSignifForFp'] = 3.0


    #The LPP mapping file is used by LPP to define the regions of param
    #space where the transits cluster.
    path = lpp.getLppDir()
    cfg['lppMapFilePath'] = \
        os.path.join(path, "octave/maps/mapQ1Q17DR24-DVMed6084.mat")
    cfg['modshiftBasename'] = \
        os.path.join(os.environ['HOME'],"daveOutput","modshift")
    cfg['onepageBasename'] = \
        os.path.join(os.environ['HOME'],"daveOutput","onepage")
    #Location of the place all the light curves and TPF files are stored
    cfg['dataStorePath'] = os.path.join(os.environ['HOME'],".mastio/k2")

    #Location of the model PRF fits files.
    cfg['prfPath'] = os.path.join(os.environ['HOME'], ".mastio/keplerprf")

    #My front end
    tasks = """pl.checkDirExistTask pl.serveTask pl.extractLightcurveTask
        pl.computeCentroidsTask pl.rollPhaseTask pl.cotrendDataTask
        pl.detrendDataTask singleEventSearchTask pl.saveOnError""".split()
    cfg['taskList'] = tasks

    cfg['singleEventTaskList'] = "blsTask trapezoidFitTask vetTask".split()

    cfg['clipSavePath'] = "./clips"
    cfg['keysToIgnoreWhenSaving'] = ["serve"]
    return cfg


@task.task
def multiEventSearchTask(clip):


    maxNumEvents= clip['config.maxNumEvents']

    clip['eventList'] = []
    for i in range(maxNumEvents):
        subClip = searchForEvent(clip)
        clip['eventList'].append(subClip)

        if 'exception' in subClip:
            clip['exception'] = subClip['exception']
            clip['backtrace'] = subClip['backtrace']
            return clip

        if subClip['disposition.isSignificantEvent'] == False:
            return clip

    #Iteration limit reached.
    return clip


@task.task
def singleEventSearchTask(clip):


    clip['eventList'] = []
    subClip = searchForEvent(clip)

    if 'exception' in subClip.keys():
        clip['exception'] = subClip['exception']
        clip['backtrace'] = subClip['backtrace']

    clip['eventList'].append(subClip)
    return clip


import dave.fileio.kplrfits as kplrfits
def searchForEvent(clip):
    subClip = clip.shallowCopy()

    originalKeyList = subClip.keys()

    #Set the flags attribute of the new subclip
    #Problem with this code is it closely tied to the behaviour
    #of multiEventSearchTask
    try:
        tmp = clip.eventList[-1]
        flags = tmp['flags']
    except (IndexError, KeyError):
        flags = clip['detrend.flags']
    subClip['flags'] = flags


    #@TODO List of tasks to run should be config param
    subClip = pl.placeholderBls(subClip)
    subClip = trapezoidFitTask(subClip)
    subClip = modshiftTask(subClip)
    subClip = measureDiffImgCentroidsTask(subClip)
    subClip = dispositionTask(subClip)

    newKeys = list(set(subClip.keys()) - set(originalKeyList))
    out = clipboard.Clipboard(__meta__=subClip['__meta__'])
    for k in newKeys:
        out[k] = subClip[k]


    #Mark all locations for this event as data not to be used.
    time = subClip['serve.time']
    period_days = subClip['trapFit.period_days']
    epoch_bkjd = subClip['trapFit.epoch_bkjd']
    duration_days = subClip['trapFit.duration_hrs'] / 24.

#    assert(np.all(np.isfinite(time[~flags])))
#    assert(np.any(flags))
    idx = kplrfits.markTransitCadences(time, period_days, epoch_bkjd, \
        duration_days, numberOfDurations=2, flags=flags)

    out['flags'] = flags | idx

    return out








import dave.trapezoidFit.estimateSnr as tf
@task.task
def trapezoidFitTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['flags']
    period_days = clip['bls.period']
    duration_hrs = clip['bls.duration_hrs']
    phase_bkjd = clip['bls.epoch']  #Check this what BLS returns
    depth_frac = clip['bls.depth']

    #We don't know these values.
    unc = np.ones_like(flux_norm)
    unc[flags] = 1e99
    flux_norm[flags] = 0


    assert(np.all(np.isfinite(time_days[~flags])))
    assert(np.all(np.isfinite(flux_norm[~flags])))
    out = tf.getSnrOfTransit(time_days, flux_norm,\
        unc, flags, \
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
    return clip


import dave.trapezoidFit.trapfit as trapFit
import dave.vetting.ModShift as ModShift
@task.task
def modshiftTask(clip):

    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    fl = clip['flags']

    epic = clip['value']
    basename = clip['config.modshiftBasename'] + "%010i" %(epic)
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    dur_hrs =  clip['trapFit.duration_hrs']
    ingress_hrs = clip['trapFit.ingress_hrs']
    depth_ppm = 1e6*clip['trapFit.depth_frac']

    subSampleN= 15
    ioBlock = trapFit.trapezoid_model_onemodel(time[~fl], period_days, \
        epoch_bkjd, depth_ppm, dur_hrs, \
        ingress_hrs, subSampleN)
    model = ioBlock.modellc -1   #Want mean of zero
#    model *= -1  #Invert for testing

    basename = "%s-%010i" %(basename, epic)
    
    modplotint=1  # Change to 0 or anything besides 1 to not have modshift produce plot
    
    out = ModShift.runModShift(time[~fl], flux[~fl], model, basename, \
        "OBJECTNAME", period_days, epoch_bkjd, modplotint)

    clip['modshift'] = out

    #I don't know which values are important, so I can't enfornce contract yet
    return clip


import dave.diffimg.centroid as cent
import dave.diffimg.prf as prf
@task.task
def measureDiffImgCentroidsTask(clip):

    #Measuring centroids requires a lot of input params
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']  #Check this what BLS returns
    duration_hrs = clip['trapFit.duration_hrs']

    cube = clip['serve.cube']
    cube[ ~np.isfinite(cube) ] = 0
    tpfHeader0 = clip['serve.tpfHeader0']
    tpfHeader = clip['serve.tpfHeader']
    ccdMod = tpfHeader0['MODULE']
    ccdOut = tpfHeader0['OUTPUT']
    bbox = cent.getBoundingBoxForImage(cube[0], tpfHeader)
    rollPhase = clip['rollPhase.rollPhase']
    prfPath = clip['config.prfPath']
    prfObj = prf.KeplerPrf(prfPath)

    time_days = clip['serve.time']
    flags = clip['flags']

#    import pdb; pdb.set_trace()
    out,log = cent.measureDiffOffset(period_days, epoch_bkjd, duration_hrs, \
        time_days, prfObj, ccdMod, ccdOut, cube, bbox, rollPhase, flags)

    #Set column names
    out = nca.Nca(out)
    out.setLookup(1, "rin intr_col intr_row diff_col diff_row".split())

    clip['diffImg'] = {'centroid_timeseries':out, 'log':log}

    clip['diffImg.centroid_timeseries']
    return clip


import dave.vetting.RoboVet as RoboVet
@task.task
def dispositionTask(clip):
    """Decide whether an event is a planet candidate or not

    TODO:
    Much of this should be parcelled off into a function
    """

    #Thresholds
    snrThreshold = clip['config.minSnrForDetection']
#    lppThreshold = clip['config.maxLppForTransit']
    offsetThreshold_sigma = clip['config.minCentroidSignifForFp']

    #Data on which to make a decision
    snr = clip['trapFit.snr']
    modshiftDict = clip['modshift']
    centroidArray = clip['diffImg.centroid_timeseries']

    out = clipboard.Clipboard(isSignificantEvent=True, isCandidate=True, \
        reasonForFail="None")
    if snr < snrThreshold:
        out['isSignificantEvent'] = False
        out['isCandidate'] = False
        out['reasonForFail'] = "SNR (%.1f) below threshold %.1f" \
            %(snr, snrThreshold)
        return out


    #Parse modshift results
    fluxVetDict = RoboVet.roboVet(modshiftDict)
    out['fluxVet'] = fluxVetDict
    assert(fluxVetDict['disp'] in ["candidate", "false positive"])

    if fluxVetDict['disp'] == "false positive":
        out['isCandidate'] = False
        out['reasonForFail'] = fluxVetDict['comments']
        return out

    #Compute centroid offset and significance
    result = cent.measureOffsetInTimeseries(centroidArray)
    out['centroidVet'] = result
    signif = result['signif']
    offset = result['offset']

    if signif > offsetThreshold_sigma:
        out['isCandidate'] = False
        out['reasonForFail'] = "Centroid offset of %.2f (%.1f sigma) detected" \
            %( offset, signif)
        return out

    clip['disposition'] = out

    #Enforce contract
    clip['disposition.isSignificantEvent']
    clip['disposition.isCandidate']
    clip['disposition.reasonForFail']
    return clip
