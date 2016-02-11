# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:05:29 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"


import dave.pipeline.clipboard as dpc
import dave.pipeline.fergalmain as dpf
import matplotlib.pyplot as mp
import numpy as np

from glob import glob
import comparebls as cb


def loadKeesConfig():
    cfg = dpf.loadMyConfiguration()
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 5
    cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
    cfg['taskList'][7] = "cb.fblsTask"

    cfg['minSnrForDetection'] = 3
    return cfg


def main():
    epicList = np.loadtxt("kees-c5.txt", usecols=(0,), delimiter=",")
    epicList = epicList[:100]


    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 5
    cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
    cfg['searchTaskList'][0] = "cb.fblsTask"

    dpf.cb = cb
    for epic in epicList[:]:
        print "Running %s" %(epic)
        clip = dpf.runOne( int(epic), cfg)

#        try:
#            time = clip['serve.time']
#            flux = clip['detrend.flux_frac']
#            flags = clip['detrend.flags']
#
#            tce = clip['eventList'][0]
#            isReal = tce['disposition.isSignificantEvent']
#            mp.clf()
#            mp.plot(time[~flags], flux[~flags], 'ko')
#            mp.title("%s: %i" %(epic, isReal))
#            mp.pause(.01)
#            mp.savefig("wd%s.png" %(epic))
#        except Exception, e:
#            print "Exception raised: %s %s" %(epic, e)




def examine():
    clipList = glob("clips/*.shelf")


    for filename in clipList:
        clip = dpc.loadClipboard(filename)
        print filename,
        if 'exception' in clip.keys():
            print "Error: ", clip['exception']
        else:
            try:
                print clip['disposition.isCandidate'], \
                    clip['disposition.reasonForFail']
            except KeyError, e:
                print "Error: %s" %(e)