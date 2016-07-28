# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 15:58:28 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import numpy as np
import dave.pipeline.clipboard as dpc

import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf


#Bad values in input detrendings are replaced with this value
BAD_FLUX_VALUE = np.nan


#def loadMultipleDetrendings(cfg):
#
#
#    epic = cfg['value']  #Is this right
#    campaign = cfg['campaign']
#    dataStorePath = cfg['dataStorePath']
#    detrendTypes = cfg['detrendTypes']

def loadMultipleDetrendings(epic, campaign, dataStorePath, detrendTypes):
    #Definition of the different kinds of detrendings and how to
    #parse the output of their archive class
    dTypeDict = dict()
    dTypeDict['PDC'] = (mastio.K2Archive(), pdcParser, "PDC")
    dTypeDict['AGP'] = (mastio.K2SCArchive(), agpParser, "K2SC")
    dTypeDict['EVEREST'] = (mastio.EverestArchive(), everestParser, "Everest")
    dTypeDict['SFF'] = (mastio.VanderburgArchive(), sffParser,
                        "Vanderburg SFF")

    out = dpc.Clipboard()

    #Load the TPF data cube
    ar = mastio.K2Archive(dataStorePath)
    fits, hdr = ar.getLongTpf(epic, campaign, header=True, mmap=False)
    hdr0 = ar.getLongTpf(epic, campaign, ext=0, mmap=False)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    out['time'] = fits['TIME']
    out['cube'] = cube
    out['tpfHeader'] = hdr
    out['tpfHeader0'] = hdr0

    #Load PA data
    fits = ar.getLongCadence(epic, campaign)
    out['rawFlux'] = fits['SAP_FLUX']


    #Load lightcurves from a specific detrending, and replace
    #the pdc time series with the new detrending
    nDetrend = 0
    for i, dType in enumerate(detrendTypes):
        key = dType.upper()
        if key not in dTypeDict:
            raise IOError("Unrecognised detrended %s" %(key))

        ar, parser, label = dTypeDict[key]

        data = ar.getLongCadence(epic, campaign)
        flux = parser(out['time'], data)

        typeName = "type%i" %(i+1)
        out[typeName] = key
        fluxName = "flux%i" %(i+1)
        out[fluxName] = flux
        labelName = "label%i" %(i+1)
        out[labelName] = label

        nDetrend += 1


    out['numDetrendings'] = nDetrend

    #Enforce contract
    out['flux1']
    out['label1']
    out['cube']

    return out


#
# Parsers. These functions take the return value from the detrenders
# archive object, and extract the flux column. They also adjust the
# output so len(output) == len(time). Missing values from the detrended
# flux are filled with BAD_FLUX_VALUE
#

def pdcParser(time, fits):
    """Extract the detrended flux for PDC data"""
    flux = fits['PDCSAP_FLUX']
    idx = np.isnan(flux)
    flux[idx] = BAD_FLUX_VALUE
    return flux

def agpParser(time, agpFits):
    """Extract the detrended flux for Suzanne Aigrain's K2SC data"""
    agpTime = agpFits['time']
    agpFlux = agpFits['flux']

    flux = np.zeros_like(time)
    idx = mapTime2ToIndexOfTime1(time, agpTime)
    flux[idx] = agpFlux
    flux[~idx] = BAD_FLUX_VALUE
    return flux


def everestParser(time, fits):
    flux = fits['FLUX']
    assert len(flux) == len(time)
    idx = np.isnan(flux)
    flux[idx] = BAD_FLUX_VALUE
    return flux


def sffParser(time, data):
    sffTime = data[:,0]
    sffFlux = data[:, 1]

    flux = np.zeros_like(time)
    idx = mapTime2ToIndexOfTime1(time, sffTime)
    flux[idx] = sffFlux
    flux[~idx] = BAD_FLUX_VALUE
    return flux






def mapTime2ToIndexOfTime1(time1, time2):
    """Compute the indices of time1 that correspond closest to values of time2

    Some K2 detrended lightcurves (like PDC) contain elements for every
    cadence, others (like sff) only contain elements for good data.
    This function enables you to figure out which rows of time1 correspond
    to the same times in time2, so you can place place the first and second
    data sets in contiguous arrays

    Typical Usage:
    ----------------
    A typical usage would look like
    ::
        idx = mapTime2ToIndexOfTime1(data1[:, 'TIME'], data2[:, 'TIME'])
        data1[idx, 'FLUX2'] = data2[:, 'FLUX']
        data1[~idx, 'FLUX2'] = BAD_VALUE

    Inputs:
    ------------
    time1:
        (1d np array) Array of times to map into
    time2
        (1d np array) Array of times to map from. Typically ``len(time2) < len(time1)``

    Returns:
    -----------
    A boolean array of length ``time1``

   """

    t1 = np.atleast_2d(time1)
    t2 = np.atleast_2d(time2)
    dt = t1 - t2.transpose()
    dt = np.fabs(dt)
    assert dt.ndim == 2
    dt = np.nanmin(dt, axis=0)
    idx = dt < 1e-8  #Two times must agree very well to be accepted
    assert len(idx) == len(time1)
    return idx
