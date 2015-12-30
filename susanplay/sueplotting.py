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
    
    plt.figtext(.25,.95,output['vet.fluxVet.disp'])
    
def indivPlot(clip,ndur):
    """Plot individual transits
    only plot the first six if that many exist.
    Plot ndur around the center of the epochs
    ndur = 2 is typical.
    """
    factor=1000; #multiplicative factor for flux
    tminus=2000;
    nt=30;
    time_days = clip['serve.time'] - tminus
    flux_zero = clip['detrend.flux_frac']*factor
    epoch = clip['trapFit.epoch_bkjd'] - tminus
    period = clip['trapFit.period_days']
    dur = clip['trapFit.duration_hrs']
    depth=clip['bls.depth']*factor;
    flags=clip['detrend.flags'];
    trapfit=clip['trapFit.bestFitModel']*factor;
    epic = str(clip['value'])
    
    time=time_days[~flags]
    flux=flux_zero[~flags]
    model=trapfit
    
    #Find number of transits available in the light curve.    
    ntransit = (time[len(time)-1]-time[0])/period;
    print np.floor(ntransit)
    ind=np.arange(-1*nt,nt)
    
    #Get time of first transit
    posBegEpochs = np.arange(-1*nt,nt) * period + (epoch - (ndur/2)*dur/24) 
    posEndEpochs = np.arange(-1*nt,nt) * period + (epoch + (ndur/2)*dur/24)
    
    #Which beginging epochs occur before the last time?
    #And which ending epochs occur after the first time?
    #I want those for which both are true.
    want = (posBegEpochs < time[len(time)-1]) & (posEndEpochs > time[0]);
    #print want
    
    plt.clf()
    #Plot first four
    c=1;
    for i in ind[want]:
        idx=np.where(ind==i)
        if c<=4:
            plt.subplot(2,3,c)
            plt.plot(time,flux,'k.')
            plt.xlim(posBegEpochs[idx[0]],posEndEpochs[idx[0]])
            plt.ylim(-1.8*depth,0.8*depth)
            ax=plt.gca()
            plt.text(.1,.1, str(c),transform=ax.transAxes)
            c=c+1;
    #Plot odd and even light curves
    plt.subplot(2,3,5)
    #plt.plot(time_days,trapfit,'r--')
    
    for i in ind[want]:
        if np.mod(i,2)==1:
            plt.plot(time-i*period,flux,'k.')
    
    plt.xlim(posBegEpochs[nt],posEndEpochs[nt])
    plt.ylim(-1.8*depth,0.8*depth)
    plt.title('Odd')
    
    plt.subplot(2,3,6)
    
    for i in ind[want]:
        if np.mod(i,2)==0:
            plt.plot(time-i*period,flux,'k.')
    plt.xlim(posBegEpochs[nt],posEndEpochs[nt])
    plt.ylim(-1.8*depth,0.8*depth)
    plt.title('Even')
        
    
    plt.figtext(.48,.95,epic,size=14)
    plt.figtext(.2,.95,clip['vet.fluxVet.disp'])
    #Plot the even and odd ones
    
    
    
    
    
    
    
    
#def setlimits(frac):
    