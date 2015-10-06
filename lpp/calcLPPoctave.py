# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 21:03:09 2015

@author: sthomp
"""

import numpy as np
import matplotlib.pyplot as plt
from oct2py import octave
#import time as timer

#t0=timer.time()
#mapfile='/home/sthomp/DAVE/origLPP/maps/mapQ1Q17DR24-DVMed6084.mat'


def calcLPPone(time,flux,mapFile,period,duration,phase):
    """
    Calculate the LPP transit metric given a time, flux (detrended)
    
    inputs
    ----------
    time : Time array in days
        array
    period : in days    
        float
    duration : in hours
        float
    phase : in days
        .float
    This runs octave code.
    
    outputs
    -------
    Tlpp : LPP transit metric value
    
    binneFlux :  The sorted, folded, binned flux values input to LPP
    """
        
    octave.addpath('/home/sthomp/DAVE/origLPP/transitLike/')
    octave.addpath('/home/sthomp/DAVE/origLPP/createLightCurves/')
    octave.addpath('/home/sthomp/DAVE/origLPP/drtoolbox/')
    octave.addpath('/home/sthomp/DAVE/origLPP/drtoolbox/techniques/')
    octave.addpath('/home/sthomp/DAVE/origLPP/stats/')
    octave.addpath('/home/sthomp/DAVE/origLPP/createMatrix/')
    
    Tlpp, Y, binnedFlux = octave.calcLPPMetricLCarray(time,flux,period,duration,phase,mapFile)


    return Tlpp , binnedFlux
    
