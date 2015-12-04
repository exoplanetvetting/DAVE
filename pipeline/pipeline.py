# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:54:57 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"




import numpy as np


import dave.fileio.kplrfits as kplrfits
import dave.lpp.calcLPPoctave as lpp
import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.fileio.nca as nca
import clipboard
import task
import os


def runOne(k2id, config):
    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = k2id

    #Check that all the tasks are properly defined
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        f = eval(t)
        clip = f(clip)

    return clip



def loadDefaultConfig():
    cfg = dict()
    cfg['debug'] = True
    cfg['campaign'] = 3
    cfg['nPointsForMedianSmooth'] = 2*48
    cfg['blsMinPeriod'] = 0.5
    cfg['blsMaxPeriod'] = 30

    #The LPP mapping file is used by LPP to define the regions of param
    #space where the transits cluster.
    path = lpp.getLppDir()
    cfg['lppMapFilePath'] = os.path.join(path, "octave/maps/mapQ1Q17DR24-DVMed6084.mat")
    cfg['modshiftBasename'] = "./modshift"

    #Location of the model PRF fits files.
    cfg['prfPath'] = os.path.join(os.environ['HOME'], ".mastio/keplerprf")

    #My front end
    tasks = """serveTask extractLightcurveTask
        computeCentroidsTask rollPhaseTask cotrendDataTask detrendDataTask
        blsTask trapezoidFitTask lppMetricTask
        measureDiffImgCentroidsTask saveClip""".split()
    cfg['taskList'] = tasks


    cfg['clipSavePath'] = "./clips"
    cfg['keysToIgnoreWhenSaving'] = ["serve"]
    return cfg






@task.task
def serveTask(clip):
    k2id = clip['value']
    campaign = clip['config.campaign']

    clip['serve'] = loadTpfAndLc(k2id, campaign)

    #Enforce contract. (Make sure expected keys are in place)
    clip['serve.time']
    clip['serve.cube']
    clip['serve.socData']
    clip['serve.tpfHeader']
    return clip


@task.task
def extractLightcurveTask(clip):
    data = clip['serve.socData']

    flagValues = clip.get('serve.flags', data[:, 'SAP_QUALITY'])
    flux = data[:, 'SAP_FLUX']

    #Convert flags to a boolean.
    mask = kplrfits.getMaskForBadK2Data()
    flags = (flagValues & mask).astype(bool)

    #Flag bad values
    flags[~np.isfinite(flux)] = True
    flags[flux<1] = True

    #Placeholder. Use the SOC PA data for the lightcurve
    out = dict()
    out['rawLightcurve'] = flux
    clip['extract'] = out
    clip['extract.source'] = "SOC PA Pipeline"
    clip['extract.flags'] = flags

    #Enforce contract
    clip['extract.rawLightcurve']
    return clip


@task.task
def cotrendDataTask(clip):
    """Produce a cotrended lightcurve in units of fractional amplitude"""

    data = clip['serve.socData']
    flags = clip['extract.flags']
    flux = data[:, 'PDCSAP_FLUX']

    flags |= ~np.isfinite(flux)

    #Remove dc offset
    dcOffset = np.median( flux[~flags])
    flux = (flux/ dcOffset) - 1
    clip['cotrend'] = {'flux_frac': flux}
    clip['cotrend.dcOffset'] = dcOffset
    clip['cotrend.flags'] = flags
    clip['cotrend.source'] = "SOC PDC Pipeline"

    #Enforce contract
    clip['cotrend.flux_frac']
    return clip


@task.task
def detrendDataTask(clip):
    flux = clip['cotrend.flux_frac']
    flags = clip['cotrend.flags']

    nPoints = clip['config.nPointsForMedianSmooth']

    #When you detrend, you must do something about the gaps and bad values.
    #This is the simplest possible thing. Replace all bad/missing data with
    #zeros. This is a placehold. Bad data inside a transit is replaced with
    #a zero, which is not what you want.
    flux[flags] = 0

    #Do a simple detrend.
    detrend = kplrfits.medianSubtract1d(flux, nPoints)
    clip['detrend'] = dict()
    clip['detrend.flux_frac'] = detrend
    clip['detrend.flags'] = flags
    clip['detrend.source'] = "Simple Median detrend"

    assert(detrend is not None)
    return clip


@task.task
def computeCentroidsTask(clip):
    data = clip['serve.socData']

    cent_colrow = np.empty( (len(data), 2))
    cent_colrow[:,0] = data[:, 'MOM_CENTR1']
    cent_colrow[:,1] = data[:, 'MOM_CENTR2']
    clip['centroids'] = {'cent_colrow': cent_colrow}
    clip['centroids.source'] = "SOC PA Pipeline"

    #Enforce contract
    clip['centroids.cent_colrow']
    return clip


import dave.diffimg.arclen as arclen
@task.task
def rollPhaseTask(clip):

    centColRow = clip['centroids.cent_colrow']
    flags = clip['extract.flags']
    rot = arclen.computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    clip['rollPhase'] = {'rollPhase':rollPhase}
    return clip


import dave.blsCode.bls_ktwo as bls
@task.task
def blsTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['detrend.flags']
    minPeriod = clip['config.blsMinPeriod']
    maxPeriod = clip['config.blsMaxPeriod']

    #Zero out the bad data. This crashes BLS
#    flux_norm[flags] = 0
#    assert(np.all( np.isfinite(flux_norm)))

    idx = flags == 0
    period, epoch, duration, depth, bls_search_periods, convolved_bls = \
        bls.doSearch(time_days[idx], flux_norm[idx], minPeriod, maxPeriod)

    out = clipboard.Clipboard()
    out['period'] = period
    out['epoch'] = epoch
    out['duration_hrs'] = duration * 24
    out['depth'] = depth
    out['bls_search_periods'] = bls_search_periods
    out['convolved_bls'] = convolved_bls
    clip['bls'] = out

    #Enforce contract
    clip['bls.period']
    clip['bls.epoch']
    clip['bls.duration_hrs']
    return clip


@task.task
def placeholderBls(clip):
    """Debugging code. Returns the ephemeris of the largest event in
    K2Id 206103150
    """
    out = clipboard.Clipboard()
    out['period'] = 4.15892
    out['epoch'] = 2145.76
    out['duration_hrs'] = 1.94443
    out['depth'] = .01112825

    clip['bls'] = out
    return clip


@task.task
def lppMetricTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']+1
    period_days = clip['bls.period']
    duration_hrs = clip['bls.duration_hrs']
    phase_bkjd = clip['bls.epoch']  #Check this what BLS returns
    mapFile = clip['config.lppMapFilePath']

    #Place holder, use Susan's version when it shows up.
    TLpp, Y, binnedFlux = lpp.fergalVersion(time_days, flux_norm, mapFile,\
        period_days, duration_hrs, phase_bkjd)

    out = dict()
    out['TLpp'] = TLpp
    out['binnedFlux'] = binnedFlux

    clip['lpp'] = out

    #Enforce contract
    clip['lpp.TLpp']
    
    return clip


import dave.trapezoidFit.estimateSnr as tf
@task.task
def trapezoidFitTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['detrend.flags']
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
    out = tf.getSnrOfTransit(time_days[~flags], flux_norm[~flags],\
        unc[~flags], flags[~flags], \
        period_days, phase_bkjd, duration_hrs, depth_frac)
    clip['trapFit'] = out

    #compute modelat all input time values
    subSampleN= 15
    ioBlock = trapFit.trapezoid_model_onemodel(time_days, period_days, \
        out['epoch_bkjd'], 1e6*out['depth_frac'], out['duration_hrs'], \
        out['ingress_hrs'], subSampleN)
    clip['trapFit.bestFitModel'] = ioBlock.modellc - 1  #Want mean of zero

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
    fl = clip['detrend.flags']

    basename = clip['config.modshiftBasename']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    dur_hrs =  clip['trapFit.duration_hrs']
    ingress_hrs = clip['trapFit.ingress_hrs']
    depth_ppm = 1e6*clip['trapFit.depth_frac']

    subSampleN= 15
    ioBlock = trapFit.trapezoid_model_onemodel(time[~fl], period_days, \
        epoch_bkjd, depth_ppm, dur_hrs, \
        ingress_hrs, subSampleN)
    model = ioBlock.modellc   #Want mean of zero


    out = ModShift.runModShift(time[~fl], flux[~fl], model, basename, \
        period_days, epoch_bkjd)

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
    flags = clip['serve.flags']

#    import pdb; pdb.set_trace()
    out,log = cent.measureDiffOffset(period_days, epoch_bkjd, duration_hrs, \
        time_days, prfObj, ccdMod, ccdOut, cube, bbox, rollPhase, flags)

    #Set column names
    out = nca.Nca(out)
    out.setLookup(1, "rin intr_col intr_row diff_col diff_row".split())

    clip['diffImg'] = {'centroid_timeseries':out, 'log':log}

    clip['diffImg.centroid_timeseries']
    return clip



#@task.task
#def plotDiagnosticsTask(clip):
#
#    time = clip['serve.time']
#    rawLc = clip['extract.rawLightCurve']
#    cotrendLc = clip['cotrend.cotrendedLightcurve']
#    detrendLc = clip['detrend.detrendedLightcurve']
#
#    plotting.plotDiagnosticLightcurves(time, rawLc, cotrendLc, detrendLc)
#    return clip


def saveOnError(clip):
    """Note this is not a task, because it should run even if
    an exception is raised"""
    if 'exception' in clip.keys():
        saveClip(clip)
    return clip


import shelve
def saveClip(clip):
    value = clip['value']
    campaign = clip['config.campaign']
    path = clip.get('config.clipSavePath', ".")

    #The problem with this is how do I which tasks to run
    #when I restore?
    keysToSkip = clip.get('config.keysToIgnoreWhenSaving', [])

    fn = "clip-%09i-%02i.shelf" %(value, campaign)
    sh = shelve.open(os.path.join(path, fn))
    for k in clip.keys():
        if k in keysToSkip:
            sh[k] = "Clip not saved"
        else:
            sh[k] = clip[k]
    sh.close()
    return clip

def loadTpfAndLc(k2id, campaign):
    ar = mastio.K2Archive()

    out = dict()
    fits, hdr = ar.getLongTpf(k2id, campaign, header=True)
    hdr0 = ar.getLongTpf(k2id, campaign, ext=0)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    out['cube'] = cube
    out['tpfHeader'] = hdr
    out['tpfHeader0'] = hdr0

    fits, hdr2 = ar.getLongCadence(k2id, campaign, header=True)
    data = kplrfits.getNumpyArrayFromFitsRec(fits)
    lookup = """ TIME TIMECORR CADENCENO
                 SAP_FLUX SAP_FLUX_ERR SAP_BKG SAP_BKG_ERR
                 PDCSAP_FLUX PDCSAP_FLUX_ERR SAP_QUALITY
                 PSF_CENTR1 PSF_CENTR1_ERR PSF_CENTR2 PSF_CENTR2_ERR
                 MOM_CENTR1 MOM_CENTR1_ERR MOM_CENTR2 MOM_CENTR2_ERR
                 POS_CORR1 POS_CORR2""".split()
    data = nca.Nca(data)
    data.setLookup(1, lookup)
    out['socData'] = data
    out['time'] = fits['TIME']
    out['flags'] = fits['SAP_QUALITY']
    return out


if __name__ == "__main__":
    cfg = loadDefaultConfig()
    runOne(206103150, cfg)


