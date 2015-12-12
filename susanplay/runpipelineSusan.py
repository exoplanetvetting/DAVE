# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:20:23 2015

@author: sthomp
"""

import dave.pipeline.pipeline as pipe
import numpy as np
import dave.pipeline.plotting as pp
import matplotlib.pyplot as plt

infile='/home/sthomp/DAVE/playK2/k2_go3049.txt'
#infile='/home/sthomp/DAVE/playK2/k2_short.txt'
outfile='/home/sthomp/DAVE/playK2/k2_go3049_bad.txt'
fid=open(outfile,'a')
#%%

cfg = pipe.loadDefaultConfig()
cfg['debug'] = False
cfg['modshiftBasename']='/home/sthomp/daveOut/jj/modshift';
cfg['prfPath']='morejunk/junk';

data=np.loadtxt(infile,dtype='float',delimiter=',',comments='#')
#%%
for i,v in enumerate(data[72:74,0]):
    
    epicid=np.int(v)
    print epicid
    output=pipe.runOne(epicid, cfg)
        
    try:
        print output['exception']
        rep="%u --bad -- %s\n" % (epicid, output['exception'])
        fid.write(rep)
        print ' ! !  EXCEPTION EXCEPTION ! ! ! ! '
        oneout="/home/sthomp/DAVE/playK2/epic%u.exc" % epicid
        fida=open(oneout,'w')
        fida.write(output['backtrace'])
        fida.close()
        
    except KeyError:
        print epicid
        rep='%u  -- good \n' % epicid
        fid.write(rep)
        
        plt.figure(1)
        plt.clf()
        plt.subplot(211)
        pp.plotData(output)
        titlewords="%s dur=%.2f h"  %(str(epicid), output['bls.duration_hrs'])
        plt.title(titlewords)
        plt.subplot(212)
        pp.plotFolded(output)
        titlewords="dur=%.2f h P=%.2f d SNR=%.1f logLPP=%.2f" % (output['trapFit.duration_hrs'], output['trapFit.period_days'],output['trapFit.snr'], np.log10(output['lpp.TLpp']))
        plt.title(titlewords)
        outfig="fig%s.png" % str(epicid)
        plt.savefig(outfig)
        
        plt.figure(2)
        plt.clf()
        plt.plot(np.arange(1,142),output['lpp.binnedFlux'],'bo')
        plt.pause(.2)
        
    
fid.close()