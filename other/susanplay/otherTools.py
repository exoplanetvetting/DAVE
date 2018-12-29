# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 13:46:03 2015

@author: sthomp
"""

import numpy as np
import nufft as nu
#import matplotlib.pyplot as plt
#from scipy import interpolate
#from scipy import signal

#Fit a sine wave at a series of harmonics.


def dofft(time,flux, over):
    """Take a fourier transform using nufft
    """
    
    dt=(time[3]-time[1])/2.0;
    n=len(time)
    endf=1/(dt*2.0); #Nyquist Frequency
    step=endf/(n*over);
    
    freq=np.arange(0,endf,step)
    freq.size
    
    fft=nu.nufft3(time,flux,freq)
    amp=2*np.sqrt(fft.real**2 + fft.imag**2)
    
      
    return(freq,amp)

