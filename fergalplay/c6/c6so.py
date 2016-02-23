# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:05:29 2016

@author: fergal

Run DAVE on the SO C6 proposal. Do your analysis in a separate module

"""

__version__ = "$Id$"
__URL__ = "$URL$"


import numpy as np

import dave.pipeline.fergalmain as dpf
import os

def loadSoConfig():
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 6
    cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
#    cfg['taskList'].insert(9, "dpp.lppMetricTask") #Not in parallel


    cfg['minSnrForDetection'] = 5
    cfg['blsMinPeriod'] = 1
    cfg['blsMaxPeriod'] = 40

    davePath = os.path.join(os.environ['HOME'],"daveOutput","c6so")
    cfg['modshiftBasename'] =  davePath
    cfg['onepageBasename']  = davePath
    cfg['clipSavePath'] = davePath

    return cfg


def main():
    epicList = np.loadtxt("GO6086.txt", usecols=(0,), delimiter=",")

    cfg = loadSoConfig()
    dpf.runAll(dpf.runOne, epicList[10:], cfg)
#    dpf.runOne(epicList[:1], cfg)




if __name__ == "__main__":
    main()