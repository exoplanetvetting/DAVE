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
#%%

t0=timer.time()
mapfile='/home/sthomp/DAVE/origLPP/maps/mapQ1Q17DR24-DVMed6084.mat'

octave.addpath('/home/sthomp/DAVE/origLPP/transitLike/')
octave.addpath('/home/sthomp/DAVE/origLPP/createLightCurves/')
octave.addpath('/home/sthomp/DAVE/origLPP/drtoolbox/')
octave.addpath('/home/sthomp/DAVE/origLPP/drtoolbox/techniques/')
octave.addpath('/home/sthomp/DAVE/origLPP/stats/')
octave.addpath('/home/sthomp/DAVE/origLPP/createMatrix/')
#%%
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

#%%

#Read in the c2 vanderberg list
octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/transitLike')
octave.addpath('/home/sthomp/DAVE/dave/octave/createLightCurves/')
octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/drtoolbox')
octave.addpath('/home/sthomp/DAVE/dave/lpp/octave/drtoolbox/techniques/')
    
c2File='/home/sthomp/DAVE/vanTCEs/c2/c2tce.txt'
mapfile='/home/sthomp/DAVE/dave/lpp/octave/maps/mapQ1Q17DR24-DVMed6084.mat'
fileout='/home/sthomp/DAVE/vanTCEs/results/c2tce_lpp-shortP.txt'

data=np.loadtxt(c2File,dtype='float',skiprows=1,delimiter=',')
np.shape(data)
percol=1
phcol=2
durcol=4
epiccol=0

dataloc='/home/sthomp/DAVE/vanTCEs/c2/ep'

fp=open(fileout,'w')

for (i,v) in enumerate(np.arange(0,2510)):
    filename=dataloc + str(int(data[i,epiccol])) + 'flat.csv'
    #print filename
    lc=np.loadtxt(filename,dtype=float,skiprows=0,delimiter=',',usecols=[0,1])
    time=lc[:,0]
    flux=lc[:,1]
    period=data[i,percol]
    phase=data[i,phcol]
    duration=data[i,durcol]*24;
    
    if period < 10:
        try:
            Tlpp,binnedFlux = lpp.calcLPPone(time,flux,mapfile,period,duration,phase)
    #    Tlpp, Y, binnedFlux = octave.calcLPPMetricLCarray(time,flux,period,duration,phase,mapfile)
        except:
            Tlpp=0.0999;
        
        out="%s %f %f %f  %f\n" % (str(int(data[i,epiccol])), period, phase, duration, Tlpp)
        print out 
        fp.write(out)
    #"%s %f %f %f  %f\n", str(int(data[i,epiccol])), period, phase, duration    plt.pause(.5), Tlpp)
    
#    if Tlpp < .005:
        tit="Tlpp=%s  period=%s" % (Tlpp,period)
        plt.figure(1)
        plt.clf()
        plt.plot(binnedFlux,'ro-')
        plt.title(tit)
        plt.pause(.1)

    

fp.close()
    
#%%
outdata=np.loadtxt(fileout,dtype='float',delimiter=None)

#Read in Vandenbergs list of candidtes and junk
candfile='/home/sthomp/DAVE/vanTCEs/c2/C2candidates-r3.csv'
cdata=np.loadtxt(candfile,dtype='string', delimiter=',', usecols=[0,1])

epic=cdata[:,0]
vet=cdata[:,1]

cand=np.zeros(len(outdata[:,0]))
for (i,v) in enumerate(outdata[:,0]):
    ind=np.where(str(int(v)) == epic)
    
    if vet[ind] == ' C':
        cand[i]= 1
    
#%%
want=cand==1
plt.figure()
plt.clf()
plt.plot(np.log10(outdata[:,1]),np.log10(outdata[:,4]),'o')
plt.plot(np.log10(outdata[want,1]), np.log10(outdata[want,4]),'r.',ms=12)
plt.plot([-1,1],[-2.2,-2.2],'g--')