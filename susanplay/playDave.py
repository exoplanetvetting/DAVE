# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:52:03 2015

@author: sthompson
"""

#Specify an EPIC
EPICnum = 206103150;
camp = 3;

import mastio
import kplrfits
import matplotlib.pyplot as plt
import numpy as np
import blsCode.bls_ktwo as k2bls
path='/Users/sthompson/data/K2/'

ar = mastio.K2Archive(path)
fits = ar.getLongCadence(EPICnum, camp)
tpf = ar.getLongTpf(EPICnum, camp)

data = kplrfits.getNumpyArrayFromFitsRec(fits)

time=data[:,0]
pdc=data[:,7]
col=7
#plt.plot(time,pdc)
nPoints=10
detrend = kplrfits.medianSubtract(data, nPoints,fCol=col)
plt.figure(1)
plt.cla()
plt.plot(time,detrend[:,col],'.')

#Quick Clean
std=np.nanstd(detrend[:,col])
m=np.nanmean(detrend[:,col])
nstd=30
want=(detrend[:,col] > m-nstd*std) & (detrend[:,col] < m+nstd*std) 

t=time[want]
dcflux=detrend[want,col]
plt.figure(2)
plt.cla()
plt.plot(t,dcflux,'g.')
plt.title('Clean and Detrended')
#%%

numlook=[0,1,2]
newLC=dcflux
periods=np.zeros(len(numlook))
depths=np.zeros(len(numlook))
durs=np.zeros(len(numlook))
for i,v in enumerate(numlook) :
    
    lcflux=newLC
    minPeriod=.4
    maxPeriod=20
    period, epoch, duration, depth, bls_search_periods, convolved_bls = \
            k2bls.doSearch(t, lcflux, minPeriod, maxPeriod)
    periods[i]=period
    depths[i]=depth
    durs[i]=duration
    
    
    plt.figure(3)
    plt.cla()
    plt.plot(bls_search_periods,convolved_bls,'-')
    plt.title(str(EPICnum))
    
    phases=np.mod((t-epoch),period)+period/2
    plt.figure(4)
    plt.cla()
    plt.plot(phases,dcflux,'r.')
    plt.plot(phases+period,dcflux,'m.')
    plt.title(str(period))
#%
    import trapfit as tf
    stdev=np.std(dcflux/1e6)
    error=stdev*np.ones(len(dcflux))
    trapFitout=tf.trapezoid_fit(t, (lcflux/1e6)+1, error, \
                      period, epoch, duration*24, -1*depth, \
                      fitTrialN=10, fitRegion=4.0, errorScale=1.0, debugLevel=0,
                      sampleN=15, showFitInterval=30)
    #%
    plt.figure()
    plt.cla()
    plt.subplot(311)
    plt.plot(t,lcflux/1e6+1,'ro')
    plt.plot(t,trapFitout.modellc,'b--')
    plt.title(str(period))
   
    
    plt.subplot(312)
    plt.plot(phases,lcflux/1e6+1,'ro')
    plt.plot(phases,trapFitout.modellc,'b.')
    plt.plot(phases+period,lcflux/1e6+1,'ro')
    plt.plot(phases+period,trapFitout.modellc,'b.')
    
    diffflux=lcflux/1e6+1 - trapFitout.modellc
    plt.subplot(313)
    plt.plot(phases,diffflux,'go')
    newLC=diffflux*1e6