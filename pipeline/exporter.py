# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:48:23 2016

@author: smullall
"""

import dave.pipeline.multiPagePlot as mpp

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
def createOutputs(clip,logTableFile):
    """
    Read in a clip and create outputs.
    and appends results to a table (logTableFile)
    """
    
    intext=logTableFile
    outfile="%s/%s.pdf" % (clip.config['modshiftBasename'],clip.value)
    
    mpp.plot_multipages(outfile,clip,intext)
    
    
    
    
    