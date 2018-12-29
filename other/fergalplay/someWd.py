# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:05:29 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"


import dave.pipeline.fergalmain as dpf
import matplotlib.pyplot as mp
import numpy as np

from glob import glob
import shelve


def main():
    epicList = """
    206169444
    206181269
    206197016
    206210269
    206216038
    206227856
    206233343
    206238453
    206242757
    206261621
    206270755
    206276634
    206284230
    206302487
    206335190
    206342949
    206364301""".split()

    cfg = dpf.loadMyConfiguration()
    cfg['debug'] = False
    cfg['taskList'][-1] = "dpp.saveClip"  #Save all clips
#    cfg['searchTaskList'][0] = "placeholderBls"

    for epic in epicList[:]:
        print "Running %s" %(epic)
        clip = dpf.runOne( int(epic), cfg)

        try:
            time = clip['serve.time']
            flux = clip['detrend.flux_frac']
            flags = clip['detrend.flags']

            tce = clip['eventList'][0]
            isReal = tce['disposition.isSignificantEvent']
            mp.clf()
            mp.plot(time[~flags], flux[~flags], 'ko')
            mp.title("%s: %i" %(epic, isReal))
            mp.pause(.01)
            mp.savefig("wd%s.png" %(epic))
        except Exception, e:
            print "Exception raised: %s %s" %(epic, e)




def examine():
    clipList = glob("clips/*.shelf")


    for filename in clipList:
        clip = shelve.open(filename)
        print filename,
        if 'exception' in clip.keys():
            print clip['exception']
        else:
            tce = clip['eventList'][0]
            print tce['disposition.isSignificantEvent']