# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:48:23 2016

@author: smullall
"""

import dave.pipeline.multiPagePlot as mpp
import dave.pipeline.pipeline as dpp
import dave.stellar.readStellarTable as stel
import os
#%%
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
def createOutputs(clip,logTableFileName):
    """
    Read in a clip and create outputs.
    and appends results to a table (logTableFile)
    """
  
    intext=logTableFileName
    outfile="%s/%s/epic%s-mp.pdf" % (clip.config['modshiftBasename'],clip.value,clip.value)
    
    dpp.plotTask(clip)    
    
    if clip.disposition.isSignificantEvent:
        clip=stel.addStellarToClip(clip)
        clip=stel.estimatePlanetProp(clip)
        
        mpp.plot_multipages(outfile,clip,intext)
        
        file2="%s/%s/%s-modshift.pdf" % (clip.config.modshiftBasename,clip.value,clip.value)
        file3="%s/%s/%s-onepage.pdf" % (clip.config.onepageBasename,clip.value,clip.value)
    
        cmd="pdftk %s %s %s output %s/%s/%s-all.pdf" % (outfile,file2,file3,clip.config['modshiftBasename'],clip.value,clip.value)
        os.system(cmd)     
    
    
    (outtxt,hdr)=createExportString(clip, delimiter=" ", badValue="nan")
    
    fid=open(logTableFileName,'a')
    fid.write("%s\n" % outtxt)
    fid.close()
    
    

    
   
    
    
    
    