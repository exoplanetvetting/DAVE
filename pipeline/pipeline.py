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
import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.fileio.nca as nca
import clipboard
import task



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

    #My front end
    tasks = """serveTask extractLightcurveTask
        computeCentroidsTask cotrendDataTask detrendDataTask
        placeholderBls trapezoidFitTask   """.split()

    #My transit finder and triage
    """ placeholderBls trapezoidFitTask  computeLppMetricTask"""
    cfg['taskList'] = tasks
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

    flags = clip.get('serve.flags', data[:, 'SAP_QUALITY'])
    flux = data[:, 'SAP_FLUX']

    #Convert flags to a boolean.
    flags = flags > 0

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



import dave.blsCode.bls_ktwo as bls
@task.task
def runBlsTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['detrend.flags']
    minPeriod = clip['config.blsMinPeriod']
    maxPeriod = clip['config.blsMaxPeriod']

    #Zero out the bad data. This crashes BLS
#    flux_norm[flags] = 0
#    assert(np.all( np.isfinite(flux_norm)))

    out = clipboard.Clipboard()
#    import pdb; pdb.set_trace()
    period, epoch, duration, depth, bls_search_periods, convolved_bls = \
        bls.doSearch(time_days, flux_norm, minPeriod, maxPeriod)

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

#import dave.lpp.calcLPPoctave as lpp
#@task.task
#def computeLppMetricTask(clip):
#    time_days = clip['serve.time']
#    flux_norm = clip['detrend.flux']
#    period_days = clip['bls.period']
#    duration_hrs = clip['bls.duration']
#    phase_bkjd = clip['bls.epoch']  #Check this what BLS returns
#    mapFile = clip['config.lppMapFilePath']
#
#    #Place holder, use Susan's version when it shows up.
#    TLpp, Y, binnedFlux = lpp.fergalVersion(time_days, flux_norm, mapFile,\
#        period_days, duration_hrs, phase_bkjd)
#
#    out = dict()
#    out['TLpp'] = TLpp
#
#    clip['lpp'] = out
#
#    #Enforce contract
#    clip['lpp.TLpp']
#    return clip
#
#
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

    print depth_frac
    print period_days

    #We don't know these values.
    unc = np.ones_like(flux_norm)

    assert(np.all(np.isfinite(time_days[~flags])))
    assert(np.all(np.isfinite(flux_norm[~flags])))

    out = tf.getSnrOfTransit(time_days[~flags], flux_norm[~flags],\
        unc[~flags], flags[~flags], \
        period_days, phase_bkjd, duration_hrs, depth_frac)

    clip['trapFit'] = out


    clip['trapFit.period_days']
    clip['trapFit.epoch_bkjd']
    clip['trapFit.duration_hrs']
    clip['trapFit.depth_frac']
    clip['trapFit.snr']
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

    #The problem with this is how do I which tasks to run
    #when I restore?
    keysToSkip = clip.get('config.keysToIgnoreWhenSaving', [])

    fn = "clip-%09i-%02i.shelf" %(value, campaign)
    sh = shelve.open(fn)
    for k in clip.keys():
        if k in keysToSkip:
            sh[k] = "Clip not saved"
        else:
            sh[k] = clip[k]
    sh.close()


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


