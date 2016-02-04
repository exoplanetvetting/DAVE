# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:20:23 2015

@author: sthomp
"""

import dave.susanplay.mainSusan as mS
#import dave.pipeline.pipeline as pipe
import numpy as np
#import dave.pipeline.plotting as pp
import matplotlib.pyplot as plt
import dave.susanplay.sueplotting as sp
import dave.plot.multipage as mp


def createExportString(clip, delimiter=" ", badValue="nan"):
    """Create a line of text for the exporter
    
    Inputs:
    ------------
    clip
        A clipboard object
    
    Optional Inputs:
    -----------------
    delimiter:
        (string) The character, or set of characters to separate elements
       
    badValue
        (string) String to be output when a value isn't present
    Returns:
    -----------
    Two strings. The first is the text to be exported. The second is the 
    list of keys that were exported
    """
    keysForExport = (   ('value' , '%10i'), \
                        ('trapFit.period_days', '%7.3f'), \
                        ('trapFit.epoch_bkjd', '%12.6f'), \
                        ('trapFit.duration_hrs', '%7.3f'), \
                        ('trapFit.snr', '%6.2f'), \
                        ('disposition.isSignificantEvent', '%1i'), \
                        ('disposition.isCandidate', '%i'), \
                        ('disposition.reasonForFail', '%s'), \
                    )
                    
    hdr = []
    text = []
    
    for tup in keysForExport:
        key, fmt = tup
        hdr.append(key)
        try:
            text.append( fmt % (clip[key]))
        except KeyError:
            text.append(badValue)
            
    text = delimiter.join(text)
    hdr = delimiter.join(hdr)
    return text, hdr
    
#%%   
#infile='/home/sthomp/DAVE/playK2/k2_go3049.txt'
vers="20150128";
infile='/home/sthomp/DAVE/playK2/KEES_2016_01_13.txt'
outfile='/home/sthomp/DAVE/playK2/KEES_2016_01_13_out%s.txt' %(vers)
#outcand='/home/sthomp/DAVE/playK2/k2_list_cand.txt'
fid=open(outfile,'a')
#%%

cfg = mS.loadMyConfiguration()
cfg['debug'] = False
cfg['modshiftBasename']='/home/sthomp/DAVE/playK2/Kees/k';
#cfg['prfPath']='morejunk/junk';

indata=np.loadtxt(infile,dtype='float',delimiter=None,comments='#',usecols=[0,5])
want=(indata[:,1]>=3) & (indata[:,1] <=5)
data=indata[want,:]
#%%
span=[55,100]
for i,v in enumerate(data[span[0]:span[1],0]):    
    epicid=np.int(v)
    acti=i+span[0]
    cfg['campaign'] = int(data[acti,1])
    print epicid,cfg['campaign']

    output=mS.runOne(epicid, cfg)
     
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
        
        try:
        
#            plt.figure(1)
#            sp.summaryPlot(output)
#            outfig="%sfig%s.png" % (cfg['modshiftBasename'],str(epicid))
#            plt.savefig(outfig)
#            plt.figure(2)
#            sp.indivPlot(output,6)
#            outfig="%sind%s.png" % (cfg['modshiftBasename'],str(epicid))
#            plt.savefig(outfig)
#            plt.pause(.1)
            if output.disposition.isCandidate == 1:
                disp="CANDIDATE"
            else:
                disp="FALSE POSITIVE"
            output.disposition.finaldisp=disp
            
            period=output.bls.period
            #info = "%u  %6.2f  %s   %s\n" % (epicid, period, disp, output['disposition.reasonForFail'])
            line, hdr= createExportString(output,delimiter="|")        
            fid.write("%s  |%s\n" % (line, disp))
            
            outfile="%s%s-%smp.pdf" % (cfg['modshiftBasename'],str(epicid),vers)
            info="version %s\n%u\n%s\n%s" % (vers,epicid,disp,output.disposition.reasonForFail)
            mp.plot_all_multipages(outfile,output,info)
            print outfile
            
        except KeyError, IOError:
            print epic
            print "Could not print outfile"
        
fid.close()


