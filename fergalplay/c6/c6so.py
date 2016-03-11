# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:05:29 2016

@author: fergal

Run DAVE on the SO C6 proposal. Do your analysis in a separate module


Things to try:
Don't return clip
Don't invoke shelve (writeClip)
Explicit del of fits objects.
Explict copy of time and flux
Adjust parmap to close processes sooner
"""

__version__ = "$Id$"
__URL__ = "$URL$"


import numpy as np

from dave.pipeline import gather
import dave.pipeline.fergalmain as dpf

from glob import glob
import os

def loadSoConfig():
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 6
    #cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
#    cfg['taskList'].insert(9, "dpp.lppMetricTask") #Not in parallel
    cfg['taskList'] = cfg['taskList'][:11]  #Limited set

    print cfg['taskList']
    
    cfg['minSnrForDetection'] = 5
    cfg['blsMinPeriod'] = 1
    cfg['blsMaxPeriod'] = 40

    davePath = os.path.join(os.environ['HOME'],"daveOutput","memLeak")
    cfg['modshiftBasename'] =  davePath
    cfg['onepageBasename']  = davePath
    cfg['clipSavePath'] = davePath

    return cfg


def main():
    epicList = np.loadtxt("GO6086.txt", usecols=(0,), delimiter=",")

    cfg = loadSoConfig()
    dpf.runAll(dpf.runOne, epicList[:], cfg)
#    return dpf.runOne(epicList[0], cfg)

#    for f in epicList[:100]:
#        dpf.runOne(f, cfg)


def reRun():
    epic, exception = reportExceptions()
    cfg = loadSoConfig()
    dpf.runAll(dpf.runOne, epic, cfg)


def partialRerun():
    epicList = np.loadtxt("GO6086.txt", usecols=(0,), delimiter=",")

    pattern = "/home/fergal/daveOutput/c6so/*/*.clip"
    runList = glob(pattern)
    exceptionList = reportExceptions()[0]

    set1 = set(runList) - set(exceptionList)
    set2 = set(epicList) - set1
    runList = list(set2)

    cfg = loadSoConfig()
    dpf.runAll(dpf.runOne, runList, cfg)


def reportExceptions():
    pattern = "/home/fergal/daveOutput/c6so/*/*.clip"
    fList  = np.array(glob(pattern))
    print len(fList)

#    idx = np.random.choice(len(fList), 100, replace=False)
#    exceptions = gather.gatherValue(fList[idx], 'exception')
    exceptions = gather.gatherValue(fList, 'exception')

    return exceptions


if __name__ == "__main__":
    main()
#    partialRerun()
#    reRun()