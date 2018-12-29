# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:10:56 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np


"""
This is the task to run the fortran bls code. It's not used
in the normal pipeline, but I keep it here in case we ever need
to test against it.
"""

import dave.blsCode.bls_ktwo as bls
@task.task
def blsTask(clip):
    time_days = clip['serve.time']
    flux_norm = clip['detrend.flux_frac']
    flags = clip['detrend.flags']
    minPeriod = clip['config.blsMinPeriod']
    maxPeriod = clip['config.blsMaxPeriod']

    #Zero out the bad data. This crashes BLS
#    flux_norm[flags] = 0
#    assert(np.all( np.isfinite(flux_norm)))

    idx = flags == 0
    period, epoch, duration, depth, bls_search_periods, convolved_bls = \
        bls.doSearch(time_days[idx], flux_norm[idx], minPeriod, maxPeriod)

    out = clipboard.Clipboard()
    out['period'] = period
    out['epoch'] = epoch
    out['duration_hrs'] = duration * 24
    out['depth'] = depth
    out['bls_search_periods'] = bls_search_periods
    out['convolved_bls'] = convolved_bls
    clip['bls'] = out

    #Enforce contract
    clip['bls.period']
    clip['bls.epoch']
    clip['bls.duration_hrs']
    return clip
