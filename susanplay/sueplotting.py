# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 21:48:33 2015

@author: sthomp
"""

import matplotlib.pyplot as plt
import dave.pipeline.plotting as pp
import numpy as np


def summaryPlot1(output):
    """
    Plot the data and some info from the clipboard
    output is a clipboard.
    """
    epicid=str(output['value'])
    blsper=output['bls.period']
    trapsnr=output['trapFit.snr']
    trapper=output['trapFit.period_days']
    trapdur=output['trapFit.period_days']
    logTlpp = np.log10(output['lpp.TLpp'])
    
    plt.clf()
    plt.subplot(2,2,(1,2))
    pp.plotData(output)
    titlewords="%s P=%.2f"  %(epicid,blsper)
    plt.title(titlewords)
    plt.subplot(223)
    pp.plotFolded(output)
    titlewords="P=%.2f d dur=%.2f h SNR=%.1f " % (trapper, trapdur,trapsnr)
    plt.title(titlewords)
    plt.xlim((0,trapper))
    
    plt.subplot(224)
    plt.cla()
    try:
        plt.plot(np.arange(1,142),output['lpp.binnedFlux'],'bo')
        lab="logLPP=%.2f" % (logTlpp)
        plt.xlim((-.9,141.9))
        plt.title(lab)
    except:
        pass
    
    try:
        plt.figtext(.15,.96,output['disposition.fluxVet.disp'],color='r',fontsize=14)
        plt.figtext(.12,.93,output['disposition.reasonForFail'],color='r',fontsize=10)
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
    depth=np.abs(clip['trapFit.depth_frac'])*factor;
    flags=clip['detrend.flags'];
    trapfit=clip['trapFit.bestFitModel']*factor;
    epic = str(clip['value'])
    
    time=time_days[~flags]
    flux=flux_zero[~flags]
    #model=trapfit
    
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
        plt.figtext(.11,.98,'Individual Transits',color='k',fontsize=12)
        plt.figtext(.11,.96,clip['disposition.fluxVet.disp'],color='r',fontsize=13)
        plt.figtext(.13,.91,clip['disposition.reasonForFail'],color='r',fontsize=10)
        plt.figtext(.8,.95,('P=%f d' % (period)),fontsize=13)
    #Plot the even and odd ones
    except:
        pass
    
def blsPlot(clip):
    """Plot the bls spectrum and info about the BLS search"""
    
    
    periods=clip.bls.bls_search_periods;
    bls=1e6*clip.bls.convolved_bls;
    highper=clip.bls.period;
    highdur=clip.bls.duration_hrs    
    
    
    plt.subplot(211)
    plt.plot(periods,bls,'k')
    plt.plot(highper,max(bls),'rv',ms=9,alpha=.6)
    plt.xscale('log')
    plt.xlim(min(periods),1.05)
    plt.title(('%u' %(clip.value)))
    plt.ylabel('bls (ppm)')    
    
    plt.subplot(212)
    plt.plot(periods,bls,'k')
    plt.plot(highper,max(bls),'rv',ms=9,alpha=.6)
    plt.xscale('log')
    plt.xlim(.95,max(periods))
    plt.xlabel('log10 (Period-days)')
    plt.ylabel('bls (ppm)')
    

    plt.figtext(.11,.97,'BLS Spectrum',color='k',fontsize=12)
    plt.figtext(.11,.94,clip['disposition.fluxVet.disp'],color='r',fontsize=12)
    plt.figtext(.65,.92,('blsP=%.2f vs. trP=%.2f days' % (clip.bls.period,clip.trapFit.period_days)),fontsize=13)
    #Mark the largest point.
    
    
    
    
#def setlimits(frac):
    