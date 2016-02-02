# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:20:23 2015

@author: sthomp
"""

import dave.pipeline.pipeline as pipe
import numpy as np
import dave.pipeline.plotting as pp
import matplotlib.pyplot as plt
import dave.susanplay.sueplotting as sp

#infile='/home/sthomp/DAVE/playK2/k2_go3049.txt'
infile='/home/sthomp/DAVE/playK2/k2_list.txt'
outfile='/home/sthomp/DAVE/playK2/k2_list_bad.txt'
#outcand='/home/sthomp/DAVE/playK2/k2_list_cand.txt'
fid=open(outfile,'a')
#%%

cfg = pipe.loadDefaultConfig()
cfg['debug'] = False
cfg['modshiftBasename']='/home/sthomp/daveOutput/k2list/vet';
#cfg['prfPath']='morejunk/junk';

data=np.loadtxt(infile,dtype='float',delimiter=',',comments='#')
#%%
for i,v in enumerate(data[0:58,0]):    
    epicid=np.int(v)
    print epicid
#    try:
    output=pipe.runOne(epicid, cfg)
#    except KeyError:
#        print "Silly Jeff error!"
#        continue
        
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
        
        plt.figure(1)
        sp.summaryPlot(output)
        outfig="%sfig%s.png" % (cfg['modshiftBasename'],str(epicid))
        plt.savefig(outfig)
        plt.figure(2)
        sp.indivPlot(output,6)
        outfig="%sind%s.png" % (cfg['modshiftBasename'],str(epicid))
        plt.savefig(outfig)
        plt.pause(.1)
        info = "%u   %s   %s\n" % (epicid, output['vet.fluxVet.disp'],output['vet.reasonForFail'])
        fid.write(info)
        
        
    
    
fid.close()