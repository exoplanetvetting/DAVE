# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:35:57 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np
import os

import dave.pipeline.clipboard as dpc
import loadMultipleDetrendings as lmd

def main():

    cfg = dpc.Clipboard()
    cfg['value']  = 211816003
    cfg['campaign'] = 5
    cfg['dataStorePath'] = os.path.join(os.environ['HOME'], ",mastio")
    cfg['detrendTypes'] = ["PDC", "Everest", "Agp", "sff"]

    #If this doesn't crash, then all detrendings were loaded.
    return lmd.loadMultipleDetrendings(cfg)