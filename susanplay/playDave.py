# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:52:03 2015

@author: sthompson
"""



import mastio
import kplrfits
import matplotlib.pyplot as plt
import numpy as np
import blsCode.bls_ktwo as k2bls


def FindPlanets(EPICnum,camp,datapath,outpath,nplanets):

    ar = mastio.K2Archive(datapath)
    fits = ar.getLongCadence(EPICnum, camp) 
    data = kplrfits.getNumpyArrayFromFitsRec(fits)
    
    time=data[:,0]
    #pdc=data[:,7]
    col=7
    #plt.plot(time,pdc)
    nPoints=12
    detrend = kplrfits.medianSubtract(data, nPoints,fCol=col)
    #plt.figure(1)
    #plt.cla()
    #plt.plot(time,detrend[:,col],'.')
    
    #Quick Clean
    std=np.nanstd(detrend[:,col])
    m=np.nanmean(detrend[:,col])
    nstd=4
    want= (detrend[:,col] < m+nstd*std) 
    
    t=time[want]
    dcflux=detrend[want,col]
    plt.figure(2)
    plt.cla()
    plt.plot(t,dcflux,'g.')
    plt.title('Clean and Detrended')
    
    numlook=np.arange(0,nplanets)
    newLC=dcflux
    periods=np.zeros(len(numlook))
    depths=np.zeros(len(numlook))
    durs=np.zeros(len(numlook))
    epochs=np.zeros(len(numlook))
    chisqs=np.zeros(len(numlook))
    good=np.zeros(len(numlook))-99999.99
    
    for i,v in enumerate(numlook) :
        
        lcflux=newLC
        minPeriod=.5
        maxPeriod=25        
        period, epoch, duration, depth, bls_search_periods, convolved_bls = \
                k2bls.doSearch(t, lcflux, minPeriod, maxPeriod)


        
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
        chisqs[i]=trapFitout.chi2min
        periods[i]=period
        depths[i]=depth
        durs[i]=duration
        epochs[i]=epoch
        if (-1*depth > 2*stdev):        
            good[i]=6999.99
        
        #%
        plt.figure(6)
        plt.clf()
        plt.subplot(311)
        plt.plot(t,lcflux/1e6+1,'ro',ms=4)
        plt.plot(t,trapFitout.modellc,'b--')
        titlestr="EPIC%u  P=%.3f" % (EPICnum,period)
        plt.title(titlestr)
       
        
        plt.subplot(312)
        plt.plot(phases,lcflux/1e6+1,'ro',ms=3)
        plt.plot(phases,trapFitout.modellc,'b.',ms=2)
        plt.plot(phases+period,lcflux/1e6+1,'ro',ms=3)
        plt.plot(phases+period,trapFitout.modellc,'b.',ms=2)
        plt.xlim((period-0.33*period,period+0.8*period))
        
        diffflux=lcflux/1e6+1 - trapFitout.modellc
        plt.subplot(313)
        plt.plot(phases,diffflux,'go',ms=4)
        newLC=diffflux*1e6
        filename="%sEPIC%u-c%u-%i.png" % (outpath,EPICnum,camp,i)
        plt.savefig(filename,format="png")
        
    return periods,epochs,durs,depths,chisqs,good
#------------------------
    
#Specify an EPIC
EPICnum = 206103150;
camp = 3;
datapath='/Users/sthompson/data/K2/'
outpath='/Users/sthompson/K2/DAVE/Play/c3planetFind/'
nplanets=2

infile='/Users/sthompson/K2/DAVE/Play/c3planetFind/lowCdppList345.csv'
#infile='/Users/sthompson/K2/DAVE/Play/c3planetFind/tmp'
inlist=np.loadtxt(infile,dtype='float',usecols=(0,1),delimiter=',')

outfile='/Users/sthompson/K2/DAVE/Play/c3planetFind/lowCdppResults.dat'

fid=open(outfile,'w')

for i,v in enumerate(inlist[:,0]):
    
    EPICnum=int(inlist[i,0])
    camp=int(inlist[i,1])    
    per,eps,durs,deps,chis,good = FindPlanets(EPICnum,camp,datapath,outpath,nplanets)

    for i, v in enumerate(per):
        outstr="%u  %f  %f %f %f  %f   %f\n" % (EPICnum, per[i],eps[i],-1*deps[i],durs[i],chis[i],good[i])
        fid.write(outstr)

fid.close()
