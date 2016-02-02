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
    
    try:
        plt.figtext(.15,.96,output['vet.fluxVet.disp'],color='r',fontsize=14)
        plt.figtext(.12,.93,output['vet.reasonForFail'],color='r',fontsize=14)
    except:
        pass
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
    depth=np.abs(clip['bls.depth'])*factor;
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
        if c<=6:
            plt.subplot(3,3,c)
            plt.plot(time,flux,'k.')
            plt.xlim(posBegEpochs[idx[0]],posEndEpochs[idx[0]])
            plt.ylim(-2*depth,0.9*depth)
            ax=plt.gca()
            plt.text(.1,.1, str(i),transform=ax.transAxes,color='m',fontsize=13)
            c=c+1;
    #Plot odd and even light curves
    plt.subplot(3,3,7)
    #plt.plot(time_days,trapfit,'r--')
    
    mid=  posBegEpochs[nt]+(posEndEpochs[nt] - posBegEpochs[nt])/2 
    print mid
    
    for i in ind[want]:
        if np.mod(i,2)==1:
            plt.plot(time-i*period,flux,'k.')
    
    plt.xlim(posBegEpochs[nt],posEndEpochs[nt])
    plt.ylim(-2*depth,0.9*depth)
    plt.plot((mid-0.5*dur/24,mid+0.5*dur/24),(-1*depth,-1*depth),'c-',linewidth=3) 
    ax=plt.gca()
    plt.text(.1,.1, 'Odd',transform=ax.transAxes,color='m',fontsize=13)
    
    plt.subplot(3,3,8)
    
    for i in ind[want]:
        if np.mod(i,2)==0:
            plt.plot(time-i*period,flux,'k.')
    plt.xlim(posBegEpochs[nt],posEndEpochs[nt])
    plt.ylim(-2*depth,0.9*depth)
    plt.plot((mid-0.5*dur/24,mid+0.5*dur/24),(-1*depth,-1*depth),'c-',linewidth=3)
    ax=plt.gca()
    plt.text(.1,.1, 'Even',transform=ax.transAxes,color='m',fontsize=13)


    plt.subplot(3,3,9) 
    for i in ind[want]:
        if np.mod(i,1)==0:
            plt.plot(time-i*period/2,flux,'k.')
    plt.xlim(posBegEpochs[nt],posEndEpochs[nt])
    plt.ylim(-2*depth,0.9*depth)
    plt.plot((mid-0.5*dur/24,mid+0.5*dur/24),(-1*depth,-1*depth),'c-',linewidth=3)
    ax=plt.gca()
    plt.text(.1,.1, 'half Period',transform=ax.transAxes,color='m',fontsize=13)
        
    
    plt.figtext(.48,.95,epic,size=14)
    
    try:
        
        plt.figtext(.12,.95,clip['vet.fluxVet.disp'],color='r',fontsize=13)
        plt.figtext(.8,.95,str(period),fontsize=13)
    #Plot the even and odd ones
    except:
        pass
    
    
    
    
    
    
    
#def setlimits(frac):
    