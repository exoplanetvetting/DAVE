# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 09:36:25 2016

@author: smullall
"""

import dave.susanplay.mainSusan as mS
import dave.pipeline.pipeline as pipe
import dave.pipeline.main as main
import numpy as np
#import dave.pipeline.plotting as pp
import matplotlib.pyplot as plt
import dave.susanplay.sueplotting as sp
import dave.plot.multipage as mp
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
    
#Code to test how many thruster firings land on the transits.
#
epic='212275451'

clip=c.loadClipboard('/soc/nfs/so-nfs/dave/c6-v2/%s/c%s-06.clip' % (epic,epic))

clip=countThermFlags(clip)
print clip.thermal

#%%
isbad=clip.thermal.isBad
time=clip.serve.time
period=clip.trapFit.period_days
epoch=clip.trapFit.epoch_bkjd
bestmodel=clip.trapFit.bestFitModel
phi = (time-epoch + .25*period) % period
dur=clip.trapFit.duration_hrs/24;

plt.figure()
plt.plot(phi[isbad],isbad[isbad]*-1*clip.trapFit.depth_frac,'ro')
plt.plot(phi,bestmodel,'b.')
plt.xlim([clip.thermal.phimin-dur,clip.thermal.phimax+dur])
plt.ylim([-1.1*clip.trapFit.depth_frac,0.05*clip.trapFit.depth_frac])
plt.title(epic)

#plt.figure()
#plt.plot(time,isbad,'rx')*-1*clip.trapFit.depth_frac
#plt.plot(time,bestmodel,'b-.')
#plt.title(epic)

#%%

import dave.pipeline.gather as g
cliploc='/soc/nfs/so-nfs/dave/c6-v2/'
sfile='/soc/nfs/so-nfs/dave/c6-v2/c6Candidates.txt'
clipar=np.loadtxt(sfile,dtype='str',usecols=[0],skiprows=0)
cliplist=list()
for v in clipar:
    cliplist.append("%s/%s/c%s-06.clip" % (cliploc,v,v))
#%%
import dave.susanplay.thermalFlags as tf

fun=tf.getFractionThermal

(epics,values)=g.gatherFunctionB(cliplist[0:20],fun)