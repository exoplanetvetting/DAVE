# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 21:48:33 2015

@author: sthomp
"""

import matplotlib.pyplot as plt
import dave.pipeline.plotting as pp
import numpy as np


def summaryPlot(output):
    """
    Plot the data and some info from the clipboard
    """
    epicid=str(output['value'])
    blsdur=output['bls.duration_hrs']
    trapsnr=output['trapFit.snr']
    trapper=output['trapFit.period_days']
    logTlpp = np.log10(output['lpp.TLpp'])
    
    plt.clf()
    plt.subplot(2,2,(1,2))
    pp.plotData(output)
    titlewords="%s dur=%.2f h"  %(epicid, blsdur)
    plt.title(titlewords)
    plt.subplot(223)
    pp.plotFolded(output)
    titlewords="dur=%.2f h P=%.2f d SNR=%.1f " % (blsdur,trapper ,trapsnr)
    plt.title(titlewords)
    plt.xlim((0,trapper))
    
    plt.subplot(224)
    plt.cla()
    plt.plot(np.arange(1,142),output['lpp.binnedFlux'],'bo')
    lab="logLPP=%.2f" % (logTlpp)
    plt.xlim((-.9,141.9))
    plt.title(lab)
    
def indivPlot(output):
    """Plot individual transits."""
    
    
def setlimits(frac):
    