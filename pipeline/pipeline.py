# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:54:57 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import multiprocessing.pool as pool
import multiprocessing
import contextlib

import numpy as np

import dave.fileio.kplrfits as kplrfits
import dave.fileio.mastio as mastio
import task
import dave.fileio.tpf as tpf
import dave.fileio.nca as nca


def runList(k2idList):

    cfg = loadDefaultConfig()

    #Tell the task module about the tasks.
    for t in cfg['taskList']:
        f = eval(t)
        task.__dict__[f] = f

    print task.__dict__.keys()


    parallel = cfg.get('debug', False)

    f = lambda x: task.runOne(x, cfg)

    f(k2idList[0])
#    out = map(f, k2idList)
#    p = pool.Pool(count)
#    count = multiprocessing.cpu_count() - 1
#    with contextlib.closing(pool.Pool(count)) as p:
#        out = task.parmap.map(task.runOne, k2idList, cfg, pool=p, parallel=parallel)

    return out


def loadDefaultConfig():
    cfg = dict()
    cfg['debug'] = True


    tasks = """serveTask extractLightcurveTask
        computeCentroidsTask cotrendDataTask""".split()
    cfg['taskList'] = tasks
    return cfg






@task.task
def serveTask(clip):
    k2id = clip['value']
    campaign = clip['campaign']

    clip['serve'] = loadTpfAndLc(k2id, campaign)

    #Enforce contract. (Make sure expected keys are in place)
    clip['serve.cube']
    clip['serve.socData']
    clip['serve.tpfHeader']
    return clip


@task.task
def extractLightcurveTask(clip):
    data = clip['serve.data']

    #Placeholder. Use the SOC PA data for the lightcurve
    out = dict()
    out['rawLightCurve'] = data[:, 'SAP_FLUX']
    clip['extract'] = out
    clip['extract.source'] = "SOC PA Pipeline"

    #Enforce contract
    clip['extract.rawLightcurve']
    return clip


@task.task
def cotrendDataTask(clip):
    data = clip['serve.data']

    clip['cotrend'] = {'cotrendedLightcurve': data[:, 'PDCSAP_FLUX']}
    clip['cotrend.source'] = "SOC PDC Pipeline"

    #Enforce contract
    clip['cotrend.cotrendedLightcurve']
    return clip



import dave.blsCode.bls_ktwo as bls
@task.task
def runBlsTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['cotrend.cotrendedLightcurve']
    minPeriod = clip['config.blsMinPeriod']
    maxPeriod = clip['config.blsMaxPeriod']


    out = clipboard.Clipboard()
    period, epoch, duration, depth, bls_search_periods, convolved_bls = \
        bls.doSearch(time_days, flux_norm, minPeriod, maxPeriod)

    out['period'] = period
    out['epoch'] = epoch
    out['duration'] = duration
    out['depth'] = depth
    out['bls_search_periods'] = bls_search_periods
    out['convolved_bls'] = convolved_bls
    clip['bls'] = out

    #Enforce contract
    clip['bls.period']
    clip['bls.epoch']
    clip['bls.duration']
    return clip


@task.task
def computeCentroidsTask(clip):
    data = clip['serve.data']

    cent_colrow = np.empty( (len(data), 2))
    cent_colrow[:,0] = data[:, 'MOM_CENTR1']
    cent_colrow[:,1] = data[:, 'MOM_CENTR2']
    clip['centroids'] = {'cent_colrow': cent_colrow}
    clip['centroids.source'] = "SOC PA Pipeline"

    #Enforce contract
    clip['centroids.cent_colrow']
    return clip


import dave.lpp.calcLPPoctave as lpp
@task.task
def computeLppMetric(clip):
    time_days = clip['serve.time']
    flux_norm = clip['cotrend.cotrendedLightcurve']
    period_days = clip['bls.period']
    duration_hrs = clip['bls.duration']
    phase_bkjd = clip['bls.epoch']  #Check this what BLS returns
    mapFile = clip['config.lppMapFilePath']

    #Place holder, use Susan's version when it shows up.
    TLpp, Y, binnedFlux = lpp.fergalVersion(time_days, flux_norm, mapFile,\
        period_days, duration_hrs, phase_bkjd)

    out = dict()
    out['TLpp'] = TLpp

    clip['lpp'] = out

    #Enforce contract
    clip['lpp.TLpp']
    return clip

def loadTpfAndLc(k2id, campaign):
    ar = mastio.K2Archive()

    out = dict()
    fits, hdr = ar.loadLongTpf(k2id, campaign, header=True)
    hdr0 = ar.loadLongTpf(k2id, campaign, ext=0)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    out['cube'] = cube
    out['tpfHeader'] = hdr
    out['tpfHeader0'] = hdr0

    fits, hdr2 = ar.loadLongCadence(k2id, campaign, header=True)
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
    return out




runList([206103150])