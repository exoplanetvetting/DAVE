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
                        ('stellar.Teff', '%4.1f'),\
                        ('stellar.Rad', '%5.3f'),\
                        ('stellar.dis', '%5.1f'),\
                        ('planet.rad_earth', '%5.2f'),\
                        ('planet.sma_au','%6.4f'),\
                        ('disposition.isSignificantEvent', ' %1i'), \
                        ('disposition.isCandidate', ' %i'), \
                        ('disposition.reasonForFail', ' %s'), \
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
def createOutputs(clip):
    """
    Read in a clip and create outputs.
    and appends results to a table (logTableFile)
    """
    #Some of this needs to not be hardwired here.  This is ugly.
    clip['config']['exportLoc']='/soc/nfs/so-nfs/dave/c6'
    clip['config']['onepageBasename']=clip['config']['exportLoc']
    clip['config']['dataStorePath']='/external_disk/K2/data'
    epic=str(int(clip.value))
    
    clip=dpp.serveTask(clip)    
    
    
    dpp.plotTask(clip)  
    
    cmd=""
    
    try:
        if clip.disposition.isSignificantEvent:
            clip['config']['stellarPar']=['Mass','Rad','Teff','dis','rho','prov','logg'] 
            clip['config']['stellarFile']='/home/smullall/Science/DAVE/dave/etc/k2EpicCatalogStellarTable5.txt' 

            clip=stel.addStellarToClip(clip)
            clip=stel.estimatePlanetProp(clip)
            
            outfile="%s/%s/epic%s-mp.pdf" % (clip.config['exportLoc'],epic,epic)
            mpp.plot_multipages(outfile,clip,clip.config.clipSavePath)
            
            file2="%s/%s/%s-modshift.pdf" % (clip.config.exportLoc,epic,epic)
            file3="%s/%s/%s-onepage.pdf" % (clip.config.onepageBasename,epic,epic)
        
            cmd="gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s/%s/%s-all.pdf %s  %s" % (clip.config.exportLoc,epic,epic, outfile,file2)
            #cmd="pdftk %s %s %s output %s/%s/%s-all.pdf" % (outfile,file2,file3,clip.config['exportLoc'],epic,epic)
            os.system(cmd)     
            
    except (KeyError,AttributeError),e:
        cmd="None"
        print epic, e
    
    return cmd
    
def writeTableLine(clip):
    
    clip['config']['stellarPar']=['Mass','Rad','Teff','dis','rho','prov','logg'] 
    clip['config']['stellarFile']='/home/smullall/Science/DAVE/dave/etc/k2EpicCatalogStellarTable5.txt' 
    try:
        if clip.disposition.isSignificantEvent:
            clip=stel.addStellarToClip(clip)
            clip=stel.estimatePlanetProp(clip)
    except (KeyError,AttributeError):
        pass
    
    outtxt,hdr = createExportString(clip, delimiter=" ", badValue="nan")

    return outtxt
    
   
    
def writeCandidates(clip):
    
    clip['config']['stellarPar']=['Mass','Rad','Teff','dis','rho','prov','logg'] 
    clip['config']['stellarFile']='/home/smullall/Science/DAVE/dave/etc/k2EpicCatalogStellarTable5.txt' 
    try:
        if clip.disposition.isSignificantEvent:
            clip=stel.addStellarToClip(clip)
            clip=stel.estimatePlanetProp(clip)
            outtxt,hdr = createExportString(clip, delimiter=" ", badValue="nan")
        else:
            outtxt=""
    except (KeyError,AttributeError):
        outtxt=""
    
    return outtxt
    
    
    