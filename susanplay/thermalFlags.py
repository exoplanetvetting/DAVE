# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 16:58:43 2016

@author: smullall
"""


import dave.susanplay.mainSusan as mS
import dave.pipeline.pipeline as pipe
import numpy as np
import dave.pipeline.clipboard as c

def countThermFlags(clip):

    thermal=dict()
    #Create the light curves.
    clip['config']['dataStorePath']='/home/smullall/Science/datastore'

    clip=pipe.serveTask(clip)
    clip=pipe.trapezoidFitTask(clip)  
    
    #Get just the interesting flags
    thruster=2**20;
    safemode=2**1;
    desat=2**5
    isbad=np.bitwise_and(clip.serve.flags,thruster+safemode+desat) != 0
    thermal['isBad'] = isbad

    time=clip.serve.time
    period=clip.trapFit.period_days
    epoch=clip.trapFit.epoch_bkjd
    #phi = np.fmod(time-epoch + .25*period, period)
    phiorig = (time-epoch + .25*period) % period
    phi = phiorig[np.isfinite(phiorig)] 
    dur=clip.trapFit.duration_hrs/24;
    
    #Calculate the phase range of the folded transit model.
    phi1=0.25*period - 0.5*dur;
    phi2=0.25*period + 0.5*dur;
    thermal['phimin']=phi1;
    thermal['phimax']=phi2;
    thermal['inTransCadTot'] = len( phi[[phi>phi1] and [phi<phi2]] )    
    
    #How many isbads exist in that phase range.
    phiisbad=phiorig[isbad]
    countBad=0
    for (i,v) in enumerate(phiisbad):
        if v > phi1 and v< phi2:
            countBad=countBad+1
            
    thermal['inTransCadBad'] = countBad
    thermal['numTrans']=np.floor((time[-1]-time[0])/period)
    clip['thermal']=thermal
    
    
    return clip


def getFractionThermal(clip):
    clip=countThermFlags(clip)
    if clip.thermal.numTrans >= 2:
        frac=clip.thermal.inTransCadBad/clip.thermal.numTrans
    else:
        frac=0;
    return frac;