# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:20:23 2015

@author: sthomp
"""

import dave.pipeline.pipeline as pipe
import numpy as np


infile='/home/sthomp/DAVE/playK2/k2_search.txt'
#infile='/home/sthomp/DAVE/playK2/k2_short.txt'
outfile='/home/sthomp/DAVE/playK2/k2_search_bad.txt'
fid=open(outfile,'a')
#%%

cfg = pipe.loadDefaultConfig()
cfg['debug'] = False

data=np.loadtxt(infile,dtype='float',delimiter=',',comments='#')

for i,v in enumerate(data[50:114,0]):
    
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
        
    
fid.close()