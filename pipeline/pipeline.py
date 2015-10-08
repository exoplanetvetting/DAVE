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
import clipboard
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

    f = lambda x: runOne(x, cfg)

    f(k2idList[0])
#    out = map(f, k2idList)
#    p = pool.Pool(count)
#    count = multiprocessing.cpu_count() - 1
#    with contextlib.closing(pool.Pool(count)) as p:
#        out = task.parmap.map(task.runOne, k2idList, cfg, pool=p, parallel=parallel)

    return out


def runOne(value, config):
    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = value

    for t in taskList:
        f = eval(t)
        clip = f(clip)

    return clip


def loadDefaultConfig():
    cfg = clipboard.Clipboard()
    cfg['debug'] = True
    cfg['campaign'] = 3

    cfg['blsMinPeriod'] = 0.5
    cfg['blsMaxPeriod'] = 30

    tasks = """serveTask extractLightcurveTask
        computeCentroidsTask cotrendDataTask""".split()
    cfg['taskList'] = tasks
    return cfg






@task.task
def serveTask(clip):
    k2id = clip['value']
    campaign = clip['config.campaign']

    clip['serve'] = loadTpfAndLc(k2id, campaign)

    #Enforce contract. (Make sure expected keys are in place)
    clip['serve.cube']
    clip['serve.socData']
    clip['serve.tpfHeader']
    return clip


@task.task
def extractLightcurveTask(clip):
    data = clip['serve.socData']

    #Placeholder. Use the SOC PA data for the lightcurve
    out = dict()
    out['rawLightCurve'] = data[:, 'SAP_FLUX']
    out['source'] = "SOC PA Pipeline"
    clip['extract'] = out

    #Enforce contract
    clip['extract.rawLightCurve']
    return clip


@task.task
def cotrendDataTask(clip):
    data = clip['serve.socData']

    clip['cotrend'] = {'cotrendedLightcurve': data[:, 'PDCSAP_FLUX']}
    clip['cotrend.source'] = "SOC PDC Pipeline"

    #Enforce contract
    clip['cotrend.cotrendedLightcurve']
    return clip


@task.task
def detrendDataTask(clip):
    time = clip['serve.time']
    cotrendFlux = clip['cotrend.cotrendedLightcurve']

    #Stub function
    #Copy the cotrended data without alteration
    out = dict()
    out['detrendedLightcurve'] = cotrendFlux
    clip['detrend'] =

    #Enforce contract
    out['d
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


import dave.blsCode.bls_generic as bls
@task.task
def runBlsTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['cotrend.cotrendedLightcurve']
    minPeriod = clip['config.blsMinPeriod']
    maxPeriod = clip['config.blsMaxPeriod']

    clip['bls'] = bls.doSearch(time, flux, minPeriod, maxPeriod)

    #Enforce contract
    clip['bls.period']
    clip['bls.epoch']
    clip['bls.duration']

def loadTpfAndLc(k2id, campaign):
    ar = mastio.K2Archive()

    out = dict()
    fits, hdr = ar.getLongTpf(k2id, campaign, header=True)
    hdr0 = ar.getLongTpf(k2id, campaign, ext=0)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    out['time'] = fits['TIME']
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
    return out




#runList([206103150])
