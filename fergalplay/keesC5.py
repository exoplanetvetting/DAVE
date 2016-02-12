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
import numpy as np

from glob import glob


def loadKeesConfig():
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 5
    cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips

    cfg['minSnrForDetection'] = 3
    return cfg


def main():
    epicList = np.loadtxt("kees-c5.txt", usecols=(0,), delimiter=",")
    epicList = epicList[:10]


    cfg = loadKeesConfig()

    for epic in epicList[:]:
        print "Running %s" %(epic),
        clip = dpf.runOne( int(epic), cfg)
        time_sec = np.sum(clip.__meta__.store.values())
        print "... took %.1f sec" %(time_sec)



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