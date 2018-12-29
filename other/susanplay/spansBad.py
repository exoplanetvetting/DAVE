# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 20:25:37 2016

@author: sthomp
"""

import numpy as np

#Try to locate cadences of transit ingress or egress
#Then check if any of those cadences are flagged with thruster firings.

def flagIngressEgress(clip):
    """
    Take a clip and return an array of flags
    That give the ingress and egress times.
    """
    
    epoch=clip.trapFit.epoch_bkjd;
    period=clip.trapFit.period_days;
    ingress=clip.trapFit.ingress_hrs/(24.0);
    duration=clip.trapFit.duration_hrs/(24.0);
    qflags=clip.serve.flags
    time=clip.serve.time;

    thruster=2**20;
    safemode=1;

    t=np.bitwise_and(qflags,thruster)/thruster
    s=np.bitwise_and(qflags,safemode)/safemode    
    inflags=t | s
    
    t1=epoch-0.5*duration;
    t4=epoch+0.5*duration;
    
    nstart=np.ceil((time[0]-epoch)/period)
    nend=np.ceil((time[-1]-epoch)/period)
    ntransit=nend-nstart;
    n=np.linspace(nstart,nend-1,nend-nstart)
    
    count=0
    for i in n:
        start=t1+i*period;
        stop=start+ingress;
        inflagged=checkTimeSpan(start,stop,time,inflags)
        
        stop=t4+i*period;
        start=stop-ingress;
        egflagged=checkTimeSpan(start,stop,time,inflags) 
        print start,stop
        print egflagged,inflagged
                
        if (egflagged+inflagged)>0:
            count=count+1
    
        
    return (count,ntransit)

def checkTimeSpan(start,stop,time,flags):
    """
    return number of flags between the start and stop times.
    flags should be true false for flags you care about.
    """
    
    want=(time>=start) & (time<=stop);
    nflags=np.sum(flags[want]==1)
    print flags[want]
    
    return nflags
    