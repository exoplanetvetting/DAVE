# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:20:23 2015

@author: sthomp
"""

import dave.susanplay.mainSusan as mS
import dave.pipeline.pipeline as pipe
import dave.pipeline.main as main
import numpy as np
#import dave.pipeline.plotting as pp
import matplotlib.pyplot as plt
import dave.susanplay.sueplotting as sp
import dave.plot.multipage as mp
import os
    
#%%   
#infile='/home/sthomp/DAVE/playK2/k2_go3049.txt'
vers="c6v1";
infile='/home/sthomp/DAVE/dave/fergalplay/c6/GO6086.txt'
outfile='/home/sthomp/daveOutput/c6/%s.txt' %(vers)
#outcand='/home/sthomp/DAVE/playK2/k2_list_cand.txt'
#fid=open(outfile,'a')

cfg = mS.loadMyConfiguration()

davePath = os.path.join(os.environ['HOME'],"daveOutput","c6")
cfg['modshiftBasename'] =  davePath
cfg['onepageBasename']  = davePath
cfg['clipSavePath'] = davePath
cfg['timeoout_sec'] = 200

cfg['stellarPar']=('Rad','Teff','Mass','logg')    
cfg['stellarFile']='/home/sthomp/DAVE/dave/etc/k2EpicCatalogStellarTable5.txt'

#cfg['modshiftBasename']='/home/stho
#cfg['prfPath']='morejunk/junk';

indata=np.loadtxt(infile,dtype='float',delimiter=',',comments='#',usecols=[0,3])
data=indata
#%%
span=[498,500]
for i,v in enumerate(data[span[0]:span[1],0]):    
    epicid=np.int(v)
    acti=i+span[0]
    cfg['campaign'] = 6
    print epicid,cfg['campaign']

    clip=mS.runOne(epicid, cfg)
     
    try:
        print clip['exception']
        rep="%u --bad -- %s\n" % (epicid, clip['exception'])
        print ' ! !  EXCEPTION EXCEPTION ! ! ! ! '
        oneout="/home/sthomp/daveOutput/c6/epic%s.exc" % epicid
        #fida=open(oneout,'w')
        #fida.write(clip['backtrace'])
        #fida.close()
        
    except KeyError:
        #print "%u   %s" % (epicid,clip.disposition.isCandidate)
        pass
    
    
