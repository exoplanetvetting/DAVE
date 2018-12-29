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
import dave.pipeline.gather as gather

from glob import glob
from join import join

def loadSoConfig():
    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['campaign'] = 6
    cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
    cfg['clipSavePath'] = "./c5BlsClips"


    cfg['minSnrForDetection'] = 3
    cfg['blsMinPeriod'] = 0.1
    cfg['blsMaxPeriod'] = 4
    return cfg


def main():
    epicList = np.loadtxt("kees-c5.txt", usecols=(0,), delimiter=",")
    epicList = epicList[:]

    cfg = loadKeesConfig()
    dpf.runAll(dpf.runOne, epicList[:], cfg)


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
        displayOne(clip)

        mp.pause(1)
        raw_input("Press ENTER to continue")


def displayOne(clip):
    mp.figure(1)
    mp.clf()
    dpp.plotData(clip)

    mp.figure(2)
    mp.clf()
    cb.compareFits(clip)


def blsSummaryPlot(clipList, num=None):
    mags = np.loadtxt("kees-c5.mags", delimiter="|")
    nPad = mags.shape[1]

#    epic, blsArray = gather.gatherValue(clipList, 'bls.convolved_bls')
    epic, blsArray = gather.gatherFunction(clipList, getBls)
    #Strip out occasional bls spectrum with non standard length
    lengths = np.array(map(lambda x: len(x), blsArray))
    typicalLength = int(np.median(lengths))
    idx = lengths == typicalLength

    epic = np.array(epic)[idx]
    blsArray = np.array(blsArray)[idx]
    obj = np.column_stack([epic, blsArray])


    obj2 = join(mags, 0, None, obj, 0, None, dtype=object)
    nPad += 1
#    print obj2.shape

    magCol = 3
    idx = np.argsort(obj2[:, magCol])
    obj2 = obj2[idx]

    mag = obj2[:,magCol]
    blsArray = np.vstack(obj2[:, -1])
    print blsArray.shape
#    return obj2

    mp.clf()
    mp.imshow(blsArray, interpolation="nearest", origin="bottom",\
        aspect="auto", cmap=mp.cm.YlGnBu_r)
    mp.colorbar()
    return blsArray


import dave.fileio.kplrfits as kf
def getBls(clip):
    bls = clip['bls.convolved_bls']
    filt = kf.medianSubtract1d(bls, 100)
    return filt



if __name__ == "__main__":
    main()