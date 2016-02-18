# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 21:00:33 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import dave.lpp.calcLPPoctave as lpp


import dave.pipeline.fergalmain as fm
import dave.pipeline.clipboard as dpc
import dave.pipeline.task as task

from oct2py import Oct2Py
import scipy.io as sio
from glob import glob
import os

def expt1():
    """Confirm that a load() function in octave will hang in parallel
    mode.

    Simple test case. Run a task that includes an octave load()
    call.
    """

    #Set debug to False to run in parallel and create crash
    taskList = ["testTask"]
    cfg = {'debug':False, 'taskList':taskList}

    fm.testTask = testTask
#    clip = fm.runOne(1234, cfg)
#    print clip

    #This causes a crash
    idList = range(10)
    fm.runAll(fm.runOne, idList, cfg)


def expt2():
    """See if code works when mat loading is done by scipy"""
    #Set debug to False to run in parallel and create crash
    taskList = ["loadMapTask"]
    cfg = {'debug':False, 'taskList':taskList}

    fm.loadMapTask = loadMapTask

#    fm.runOne(1234, cfg)
    idList = range(10)
    fm.runAll(fm.runOne, idList, cfg)



def expt3():
    clipList = glob("../clips/*.clip")

    epics = np.array(map(lambda x: x[-17:-8], clipList))
    print epics[0]

    #Set debug to False to run in parallel and create crash
    taskList = ["loadTask", "dpp.serveTask", "dpp.lppMetricTask"]
    cfg = {'debug':True, 'taskList':taskList}
    cfg['debug'] = False

    fm.loadTask = loadTask
#    fm.runOne(epics[0], cfg)
    fm.runAll(fm.runOne, epics[:9], cfg)

def loadTask(clip):
    epic = clip['value']
    debug = clip['config.debug']
    fn = "../clips/c%s-05.clip" %(epic)

    clip = dpc.loadClipboard(fn)
    clip['config.debug'] = debug
    return clip


@task.task
def testTask(clip):

    path = lpp.getLppDir()
    octave = Oct2Py()

    octave.addpath(path)
    octave.addpath(os.getcwd())

    mapFile = os.path.join(path,
       "octave/maps/mapQ1Q17DR24-DVMed6084.mat")
    total = octave.testFunc(mapFile)

    clip['lpp'] = total
    print clip['value'], total
    return clip



@task.task
def loadMapTask(clip):
    value = clip.value
    path = lpp.getLppDir()
    mapFile = os.path.join(path,
       "octave/maps/mapQ1Q17DR24-DVMed6084.mat")

    #This isn't sufficient. I have to massage this to
    #become acceptable to oct2py as a matlab struct of arrays
    print "Cp1 %i" %value
    lppObj = sio.loadmat(mapFile, squeeze_me=True, struct_as_record=False)['map']
    print "Cp2 %i" %value

    struct = dict()
    struct['nDim'] = lppObj.nDim
    struct['knn'] = lppObj.knn
    struct['knnGood'] = lppObj.knnGood
    struct['mapped'] = lppObj.mapped
    clip['lppObj'] = struct

    print "Cp3 %i" %value

    octave = Oct2Py()
    path = lpp.getLppDir()
    octave.addpath(path)
#    octave.addpath(os.getcwd())

#    mapping = octave.testFunc2(struct)
#    print "mapping is %s" %(str(mapping))
#    clip['mapping'] = mapping

    return clip


#def loadLppMatAsStruct(mapFile):
#    """Load the LPP mat file in a format that can be passed to octave
#
#    This is a work around specific to the fact that octave's load() command
#    hangs when run in parallel (probably due to a race condition).
#
#    You probably don't want to reuse this function in any other code.
#
#    Inputs:
#    ----------
#    mapFile
#        (string) Full path to LPP's .mat file you want to load
#
#    Returns:
#    ------------
#    A dict. Pass this dict to octave and it will be seen as a struct.
#    """
#    lppObj = sio.loadmat(mapFile, squeeze_me=True, struct_as_record=False)['map']
#    struct = dict()
#    struct['nDim'] = lppObj.nDim
#    struct['knn'] = lppObj.knn
#    struct['knnGood'] = lppObj.knnGood
#    struct['mapped'] = lppObj.mapped
