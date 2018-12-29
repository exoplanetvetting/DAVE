# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:24:17 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import numpy as np

import dave.pipeline.task as task
import dave.pipeline.pipeline as pl

from multiprocessing import pool
import multiprocessing
import contextlib
import parmap

def main():
    data = np.loadtxt("./exampleTargets/C3/K2C3cat.txt", usecols=(0,))

    cfg = pl.loadDefaultConfig()

    taskList = cfg['taskList']
#    for i in range(len(taskList)):
#        taskList[i] = "pl.%s" %(taskList[i])
#    cfg['taskList'] = taskList

    cfg['taskList'] = taskList[:10]
    print cfg['taskList']

    count = multiprocessing.cpu_count() - 1
    p = pool.Pool(count)
    print count

    cfg['debug'] = False
    parallel = cfg.get('debug', False)
    parallel= False

    #Pool doesn't release threads even when it runs to completion.
    #Problem not related to exceptions being raised
    with contextlib.closing(pool.Pool(count)) as p:
        out = parmap.map(pl.runOne, data[1:3], cfg, parallel=parallel)
    p.join()
    p.close()

    return out