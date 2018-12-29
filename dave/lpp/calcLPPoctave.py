# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 21:03:09 2015

@author: sthomp
"""

import numpy as np
import matplotlib.pyplot as plt
from oct2py import Oct2Py
import os
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

    binnedFlux :  The sorted, folded, binned flux values input to LPP
    """

    octave = Oct2Py()
    octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/transitLike')
    octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/createLightCurves/')
    octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/drtoolbox/')
    octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/drtoolbox/techniques/')
    #octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/drtoolbox')

    Tlpp, Y, binnedFlux = octave.calcLPPMetricLCarray(time,flux,period,duration,phase,mapFile)


    return Tlpp , binnedFlux



def fergalVersion(time, flux, mapFile, period, duration, phase):
    path = getLppDir()

    #Create a new instance for each time we run LPP. The
    #oct2py.octave is not threadsafe and will crash when run in
    #parallel. Oct2Py() won't

    with Oct2Py() as octave:
        octave.addpath(path)
        octave.addpath(path + "/octave/transitLike")
        octave.addpath(path + "/octave/createLightCurves/")
        octave.addpath(path + "/octave/drtoolbox/")
        octave.addpath(path + "/octave/drtoolbox/techniques")

        Tlpp, Y, binnedFlux = octave.calcLPPMetricLCarray(\
            time,flux,period,duration,phase,mapFile)


    return Tlpp.copy(), Y.copy(), binnedFlux.copy()



def getLppDir():
    """Get the path where LPP stores its .m files"""
    pathSep = "/"
    path = os.path.realpath(__file__)
    path = pathSep.join(path.split(pathSep)[:-1])
    return path
