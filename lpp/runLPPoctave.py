# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 21:10:16 2015

@author: sthomp
"""

import numpy as np
import matplotlib.pyplot as plt
from oct2py import octave
import time as timer
import calcLPPoctave as lpp


t0=timer.time()
mapfile='/home/sthomp/DAVE/origLPP/maps/mapQ1Q17DR24-DVMed6084.mat'

octave.addpath('/home/sthomp/DAVE/origLPP/transitLike/')
octave.addpath('/home/sthomp/DAVE/origLPP/createLightCurves/')
octave.addpath('/home/sthomp/DAVE/origLPP/drtoolbox/')
octave.addpath('/home/sthomp/DAVE/origLPP/drtoolbox/techniques/')
octave.addpath('/home/sthomp/DAVE/origLPP/stats/')
octave.addpath('/home/sthomp/DAVE/origLPP/createMatrix/')

time=np.arange(1,80,.04167)
flux=10*(np.sin(2*np.pi*time/20))**4;

period=10
phase=15;
dur=15;

plt.figure();
plt.plot(time,flux,'.')

Tlpp, Y, binnedFlux = octave.calcLPPMetricLCarray(time,flux,period,dur,phase,mapfile)

plt.figure()
plt.subplot(211)
plt.plot(binnedFlux)
plt.subplot(212)
plt.plot(Y,'r.')

print Tlpp

t1=timer.time()
t=t1-t0;
print t