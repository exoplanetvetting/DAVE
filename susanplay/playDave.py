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


import blsCode.bls_ktwo as k2bls
minPeriod=.33
maxPeriod=30
period, epoch, duration, depth, bls_search_periods, convolved_bls = \
        k2bls.doSearch(t, dcflux, minPeriod, maxPeriod)

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