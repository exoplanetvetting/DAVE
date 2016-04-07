# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 16:32:28 2016

@author: smullall
"""


import numpy as np
import dave.pipeline.fergalmain as dpf

#%%

def loadSoConfig():
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 7
    tasks = """dpp.checkDirExistTask dpp.serveLocalTask dpp.extractLightcurveTask
        dpp.computeCentroidsTask dpp.rollPhaseTask dpp.cotrendDataTask
        dpp.detrendDataTask dpp.fblsTask dpp.trapezoidFitTask dpp.lppMetricTask 
        dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.saveClip """.split() #dpp.runExporterTask
        
    sfftasks ="""dpp.checkDirExistTask dpp.serveLocalTask dpp.extractLightcurveTask
        dpp.computeCentroidsTask dpp.rollPhaseTask 
        dpp.extractLightcurveFromTpfTask dpp.cotrendSffDataTask
        dpp.detrendDataTask dpp.fblsTask dpp.trapezoidFitTask dpp.lppMetricTask 
        dpp.modshiftTask dpp.measureDiffImgCentroidsTask dpp.dispositionTask
        dpp.saveClip """.split()   
        
    cfg['taskList'] = sfftasks
    
    #cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
#    cfg['taskList'].insert(9, "dpp.lppMetricTask") #Not in parallel
   # cfg['taskList'] = cfg['taskList'][:8]  #Limited set
    
    cfg['minSnrForDetection'] = 5
    cfg['blsMinPeriod'] = 0.7
    cfg['blsMaxPeriod'] = 40
    cfg['maxLppForTransit']= 0.007093282347242814

    cfg['keysToIgnoreWhenSaving'] = "serveTask"
    davePath = "/soc/nfs/so-nfs/dave/c7-sff/"
    cfg['modshiftBasename'] =  davePath
    cfg['onepageBasename']  = davePath
    cfg['clipSavePath'] = davePath
    cfg['dataStorePath'] = '/home/smullall/Science/DAVE/datastore/data/'
    cfg['prfPath'] = '/home/smullall/Science/DAVE/datastore/prf/'
    return cfg


#%%
epicList = np.loadtxt("/home/smullall/Science/DAVE/C7/epiclist", usecols=(0,), delimiter=",")
#print epicList
cfg = loadSoConfig()
#clip=dpf.runOne(epicList[4],cfg,True)
dpf.runAll(dpf.runOne, epicList[200:1000], cfg)