# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:06:56 2016

@author: smullall
"""

import dave.susanplay.mainSusan as mS
import numpy as np
import matplotlib.pyplot as plt
import dave.pipeline.gather as gr
import dave.pipeline.exporter as ex


cliploc='/soc/nfs/so-nfs/dave/c6-v2/'

listfile="%s/%s" % (cliploc,'cliplist')

clipar=np.loadtxt(listfile,dtype='str')

cliplist=list()

for v in clipar:
    cliplist.append("%s/%s" % (cliploc,v))

#%%    
func=ex.writeTableLine
(epics,table)=gr.gatherFunctionB(cliplist, func)

outtable='/soc/nfs/so-nfs/dave/c6-v2/c6results.txt'
fid=open(outtable,'w')

for e in table:
    fid.write("%s\n" % e)

fid.close()
#%%
func=ex.writeCandidates
(epics,table)=gr.gatherFunctionB(cliplist,func)

outtable='/soc/nfs/so-nfs/dave/c6-v2/c6Candidates.txt'
fid=open(outtable,'w')

for e in table:
    if e != "":
        fid.write("%s\n" % e)
fid.close()
#%%
sfile='/soc/nfs/so-nfs/dave/c6-v2/c6Candidates.txt'
clipar=np.loadtxt(sfile,dtype='str',usecols=[0],skiprows=0)
cliplist=list()
for v in clipar:
    cliplist.append("%s/%s/c%s-06.clip" % (cliploc,v,v))
#%%

func=ex.createOutputs
(epics,out)=gr.gatherFunctionB(cliplist[2:],func)

