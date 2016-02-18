# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:32:54 2016

@author: sthomp
"""

import dave.susanplay.mainSusan as mS
#import dave.pipeline.pipeline as pipe
import numpy as np
#import dave.pipeline.plotting as pp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import dave.pipeline.pipeline as pipe

vers="20150212";
infile='/home/sthomp/DAVE/playK2/ch42_bls.txt'
outfile='/home/sthomp/DAVE/playK2/c6question%s.txt' %(vers)
#outcand='/home/sthomp/DAVE/playK2/k2_list_cand.txt'


cfg = mS.loadMyConfiguration()
cfg['debug'] = False
cfg['blsMaxPeriod'] = 2.2
cfg['blsMinPeriod'] = .3
cfg['taskList']="""dpp.serveTask dpp.extractLightcurveTask dpp.cotrendDataTask dpp.detrendDataTask dpp.blsTask""".split()
cfg['clipSavePath']="/home/sthomp/daveOutput/bls"
#cfg['prfPath']='morejunk/junk';

indata=np.loadtxt(infile,dtype='float',delimiter=',',comments='#',usecols=[0,1,2])
want = indata[:,1]==6
rdata=indata[want,:]
idx=np.argsort(rdata[:,2])
rundata=rdata[idx,:]

#%%

periods=np.linspace(0.32,2.1,num=25000)
bls=np.zeros((410,25000))
fid=open(outfile,'a')
span=[55,410]
n=0
for i,v in enumerate(rundata[span[0]:span[1],0]):    
    epicid=np.int(v)
    acti=i+span[0];
    cfg['campaign'] = int(rundata[acti,1])
    print epicid,cfg['campaign']

    clip=mS.runOne(epicid, cfg)
    
    #plt.plot(clip.bls.bls_search_periods,clip.bls.convolved_bls,'r')
    
    try:
        f=interp1d(clip.bls.bls_search_periods,clip.bls.convolved_bls)
        ibls=f(periods)
        #bls[n,:]=ibls;
        clip.bls['ibls']=ibls;
        clip.bls['iperiods']=periods
        clip.bls['kepmag']=rundata[acti,2]
        clip=pipe.saveClip(clip)
        n=n+1;
    except AttributeError:
        pass

#np.save('/home/sthomp/DAVE/playK2/ch42-bls.npy', bls)    
#%%
#Code to gather up the bls spectra into a 2d array belongs here.
import dave.pipeline.clipboard as cb
listdir='/home/sthomp/daveOutput/bls/list.txt'
clipnames=np.loadtxt(listdir,dtype='string')
bls=np.zeros((len(clipnames),25000))
kepmag=np.zeros((len(clipnames),1))
for i,v in enumerate(clipnames):
    c=cb.loadClipboard(v)
    bls[i,:]=c.bls.ibls;
    kepmag[i]=c.bls.kepmag;
    
idx=np.argsort(kepmag,axis=0)

#%%
    
#X, Y = np.meshgrid(periods,np.linspace(0,499,500))
Z=bls[idx[:,0],:]/np.max(bls[idx[:,0],:],axis=1, keepdims=True)
plt.clf()
plt.set_cmap('gist_earth')
isf=np.isnan(Z)
Z[isf]=0
plt.imshow(Z,extent=(0.32,2.1,len(clipnames),0),aspect="auto",interpolation='nearest')
plt.title('Channel 42  BLS Spectra')
plt.xlabel('Period (days)')
plt.savefig('ch42blsSpectrum.png')