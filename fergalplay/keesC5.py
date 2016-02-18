# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:05:29 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"


import matplotlib.pyplot as mp
import numpy as np

import dave.pipeline.clipboard as dpc
import dave.pipeline.fergalmain as dpf
import dave.pipeline.exporter as exporter
import numpy as np

from glob import glob


def loadKeesConfig():
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 5
    cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
#    cfg['taskList'].insert(9, "dpp.lppMetricTask") #Not in parallel
    cfg['clipSavePath'] = "./clips"

    cfg['minSnrForDetection'] = 3
    return cfg


def main():
    epicList = np.loadtxt("kees-c5.txt", usecols=(0,), delimiter=",")
    epicList = epicList[:]

    cfg = loadKeesConfig()
    dpf.runAll(dpf.runOne, epicList, cfg)


#    for epic in epicList:
#        print "Running %s" %(epic),
#        clip = dpf.runOne( int(epic), cfg)
#        time_sec = np.sum(clip.__meta__.store.values())
#        print "... took %.1f sec" %(time_sec)



def examine():
    clipList = glob("clips/*.clip")

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


def getPcList():
    clipList = glob("clips/*.clip")

    pcList = []
    for filename in clipList:
        clip = dpc.loadClipboard(filename)
        try:
            isCand = clip['disposition.isCandidate']
            if isCand:
                pcList.append(filename)
        except KeyError:
            print filename, clip.exception
    return pcList

def writeTable(pcList):
    for fn in pcList:
        clip = dpc.loadClipboard(fn)
        exporter.exporterTask(clip)


import dave.pipeline.pipeline as pl
import dave.pipeline.plotting as dpp
import comparebls as cb
def display(pcList):
    for i, pc in enumerate(sorted(pcList)):
        print i, pc
        clip = dpc.loadClipboard(pc)
        clip = pl.serveTask(clip)

        mp.figure(1)
        mp.clf()
        dpp.plotData(clip)

        mp.figure(2)
        mp.clf()
        cb.compareFits(clip)

        mp.pause(1)
        raw_input("Press ENTER to continue")

if __name__ == "__main__":
    main()