#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 20:32:12 2018
Functions to correctly fold and bin a light curve.
Calculate the lpp metric: transform to lower dimensions, knn

Depends on class from reading in a previously created LPP metric Map

Depends on reading in the light curve to data structure.

input is a class called data 
data contains
data.time (days)
data.tzero (day)
data.dur (hours)
data.period (days)
data.flux (normalized to 0)

After foldBinLightCurve it contains
data.binned
After transform it contains
data.lpp_transform


@author: smullally
"""
from __future__ import division
import numpy as np
import scipy.signal as signal
from sklearn.neighbors import NearestNeighbors


def computeLPPTransitMetric(data,mapInfo):
    """
    This function takes a data class with light curve info
    and the mapInfo with information about the mapping to use.
    It then returns a lpp metric value.
    """
    
    binFlux, binPhase=foldBinLightCurve(data,mapInfo.ntrfr,mapInfo.npts)
    
    #Dimensionality Reduction and knn parts
    rawTLpp=computeRawLPPTransitMetric(binFlux,mapInfo)
    
    #Normalize by Period Dependence
    normTLpp=periodNormalLPPTransitMetric(rawTLPP,mapInfo)
    
    return normTLpp
    
    
    

def runningMedian(t,y,dt,runt):
    """
    Take a running median of size dt
    Return values at times given in runt
    """
    newy=np.zeros(len(y))
    newt=np.zeros(len(y))
    
    srt = np.argsort(t)
    newt = t[srt]
    newy = y[srt]

    runy=[]
    for i in range(len(runt)):      
        tmp=[]
        for j in range(len(newt)):     
            if (newt[j] >= (runt[i]-dt)) and (newt[j] <= (runt[i]+dt)):
                tmp.append(newy[j])
                
        if np.isnan(np.nanmedian(np.array(tmp))) :
            runy.append(0)
        else:
            runy.append(np.nanmedian(np.array(tmp)))
    
    return(list(runt),runy)

                
    

def foldBinLightCurve (data, ntrfr, npts):
    """
    Fold and bin light curve for input to LPP metric calculation
    
    data contains time, tzero, dur, priod and flux (around zero)
    
    ntrfr -- number of transit fraction for binning around transit ~1.5
    npts -- number of points in the final binning.
    
    """

    #Create phase light curve
    phaselc =np.mod((data.time-(data.tzero-0.5*data.period))/data.period,1)
    flux=data.flux
    #Determine the fraction of the time the planet transits the star.
    #Insist that ntrfr * transit fraction
    if ~np.isnan(data.dur) & (data.dur >0):
        transit_dur = data.dur
    else:
        transit_dur = 0.2 * data.period/24.
    
    transit_fr=transit_dur/24./data.period
    if (transit_fr * ntrfr) > 0.5 :
        transit_fr = 0.5/ntrfr
        
    #Specify the out of transit (a) and the in transit regions
    endfr = .03
    a = np.concatenate((np.arange(endfr,(0.5-endfr),1/npts) , \
                        np.arange((0.5+endfr),(1-endfr),1/npts)), axis=None)
    b =np.concatenate((np.arange((0.5-ntrfr*transit_fr),\
                                (0.5+ntrfr*transit_fr),(4*ntrfr*transit_fr)/npts)),axis=None)
    
    [runta,runya] = runningMedian(phaselc,flux,1.5/npts,a)
    [runtb,runyb] = runningMedian(phaselc,flux,(5*ntrfr*transit_fr)/npts,b)

    #Combine the two sets of bins
    runymess=np.array(runya + runyb)
    runtmess = np.array(runta + runtb)

    srt=np.argsort(runtmess)
    runy=runymess[srt]
    runt=runtmess[srt]
    
    #Scale the flux by the depth so everything has the same depth.
    #Catch or dividing by zero is to not scale.
    scale = -1*np.min(runyb)
    if scale != 0:
        scaledFlux=runy/scale
    else:
        scaledFlux=runy
    
    binnedFlux=scaledFlux
    phasebins=runt
    
    return binnedFlux,phasebins


def computeRawLPPTransitMetric(binFlux,mapInfo):
    """
    Perform the matrix transformation with LPP
    Do the knn test to get a raw LPP transit metric number.
    """
    
    
    
    