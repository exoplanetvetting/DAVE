# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 13:31:28 2016

@author: smullall
"""

"""
Code to Perform transit injection into K2 light curves
"""

import numpy as np
import dave.trapezoidFit.trapfit as tf
import dave.pipeline.pipeline as dpp
import dave.pipeline.fergalmain as dpf



def injectTransitClip(clip):
    """
    The clipboard needs to contain
       clip.inject:  period_days,epoch_bkjd,depth_ppm,duration_hrs,ingress_hrs
       clip.serve: needs to contain the serve and cotrend keys
       
       the output clip has altered the clip.cotrend.flux_frac with
       transits specified in clip.inject
    """
    time_days=clip.serve.time    
    inj=clip.inject    
    influx=clip.cotrend.flux_frac
    
    
    subSampleN= 15
    time_days[~np.isfinite(time_days)] = 0  #Hide the Nans from one_model
    assert(np.all(np.isfinite(time_days)))
    ioBlock = tf.trapezoid_model_onemodel(time_days, inj['period_days'], \
                inj['epoch_bkjd'], inj['depth_ppm'], inj['duration_hrs'], \
                inj['ingress_hrs'], subSampleN)
    
    injModel = ioBlock.modellc - 1  #Want mean of zero

    clip.cotrend['flux_frac'] = influx + injModel
    
    clip.cotrend['source'] = "%s %s" % (clip.cotrend.source, "Injected")
    
    return clip


def findInjTransits(clip):
    """
    Now that the light curve has a transit in it, run the rest of the pipeline.
    dpp.lppMetricTask 
        dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.saveClip 
    """
    
    clip=dpp.detrendDataTask(clip)
    clip=dpp.fblsTask(clip)
    clip=dpp.trapezoidFitTask(clip)
    clip=dpp.dispositionTask(clip)
    clip=dpp.lppMetricTask(clip)
    clip=dpp.saveClip(clip)
    
    return clip

def loadInjConfig():
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    tasks = """dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveTask
            dpp.cotrendDataTask""".split() #dpp.runExporterTask
        
    #sfftasks ="""dpp.checkDirExistTask dpp.serveTask dpp.extractLightcurveFromTpfTask
    #    dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendSffDataTask
    #    dpp.detrendDataTask dpp.fblsTask dpp.trapezoidFitTask dpp.lppMetricTask 
    #    dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
    #    dpp.saveClip """.split()   
        
    cfg['taskList'] = tasks
    
    #cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
#    cfg['taskList'].insert(9, "dpp.lppMetricTask") #Not in parallel
   # cfg['taskList'] = cfg['taskList'][:8]  #Limited set
    
    cfg['minSnrForDetection'] = 5
    cfg['blsMinPeriod'] = 0.5
    cfg['blsMaxPeriod'] = 30
    cfg['maxLppForTransit']= 0.007093282347242814

    cfg['keysToIgnoreWhenSaving'] = "serveTask"
    davePath = "/soc/nfs/so-nfs/dave/injection/"
    cfg['modshiftBasename'] =  davePath
    cfg['onepageBasename']  = davePath
    cfg['clipSavePath'] = davePath
    cfg['dataStorePath'] = '/home/smullall/Science/DAVE/datastore/data/'
    cfg['prfPath'] = '/home/smullall/Science/DAVE/datastore/prf/'
    return cfg

def testInjResults(clip):
    """
    Look at what was found and determine if it matches the input.
    """
    pdiff=1 - (clip.inject.period_days / clip.trapFit.period_days)  
  
    
    if pdiff < .01:
        clip.inject['wasFound'] = True
    else:
        clip.inject['wasFound'] = False
        
    return clip

#---------
#import pandas as pd
#C=pd.DataFrame(data=[],index)
import pdb
outfile='/soc/nfs/so-nfs/dave/injection'
fid=
targetlistfile='/home/smullall/Science/DAVE/dave/susanplay/forinjections.txt'
periodRange=(1,30)
depthRange=(3000,50000)
epochoffset=10  #days
N=10

data=np.loadtxt(targetlistfile,usecols=(0,1), delimiter=",",comments='#')

targetchoice=np.random.randint(0,len(data),size=N)

cfg=loadInjConfig()

for n in targetchoice:
    epic=data[n,0]
    campaign=data[n,1]
    inject={}
    inject['period_days']=np.random.rand(1)*(periodRange[1]-periodRange[0])
    inject['depth_ppm']=np.random.rand(1)*(depthRange[1]-depthRange[0])
    inject['duration_hrs']=6.0
    inject['ingress_hrs']=1.0
    
    
    cfg['campaign']=campaign
    clip=dpf.runOne(epic,cfg,True)
    clip['inject']=inject
    clip.inject['epoch_bkjd']=clip.serve.time[0]+epochoffset
    #pdb.set_trace()        
    
    clip=injectTransitClip(clip)
    
    clip=findInjTransits(clip)
    
    clip=testInjResults(clip)
    
    print clip.inject
    print clip.inject['wasFound'], clip.trapFit.period_days, clip.trapFit.depth_frac*1e6 
        


    