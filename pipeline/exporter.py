# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:48:23 2016

@author: smullall
"""

import dave.pipeline.multiPagePlot as mpp
import dave.pipeline.pipeline as dpp
import dave.stellar.readStellarTable as stel
import os
from pdb import set_trace as db
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
    keysForExport = (   ('value' , '%09i'), \
                        ('trapFit.period_days', '%7.3f'), \
                        ('trapFit.epoch_bkjd', '%12.6f'), \
                        ('trapFit.duration_hrs', '%7.3f'), \
                        ('trapFit.snr', '%6.2f'), \
                        ('stellar.Teff', '%6.1f'),\
                        ('stellar.Rad', '%5.3f'),\
                        ('stellar.dis', '%6.1f'),\
                        ('planet.rad_earth', '%6.2f'),\
                        ('planet.sma_au','%6.3f'),\
                        ('config.detrendType','%8s'),\
                        ('disposition.fluxVet.not_trans_like', '\t%1i'), \
                        ('disposition.fluxVet.sig_sec', ' %1i'),\
                        ('disposition.centroidVet.isCentroidFail', ' %1i'),\
                        ('disposition.isCandidate', ' %1i'), \
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
    clip['config']['exportLoc']='/soc/nfs/so-nfs/dave/c7-pdc/'
    clip['config']['onepageBasename']=clip['config']['exportLoc']
    clip['config']['dataStorePath']='/home/smullall/Science/datastore'
    epic=str(int(clip.value))
    
    #print clip.serve

    try:
        clip.serve.time
    except AttributeError:
        clip=dpp.serveLocalTask(clip)    
        print('hi serve')        
        print(clip.serve)
    try:
        print(clip.exception)
    except AttributeError:
            pass
    #dpp.plotTask(clip)  
    
    cmd=""

    try:
        if clip.disposition.isSignificantEvent:
            clip['config']['stellarPar']=['Mass','Rad','Teff','dis','rho','prov','logg'] 
            #clip['config']['stellarFile']='/home/smullall/Science/DAVE/dave/etc/k2EpicCatalogStellarTable5.txt' 

            clip=stel.addStellarToClip(clip)
            clip=stel.estimatePlanetProp(clip)
            
            outfile="%s/%s/epic%s-mp.pdf" % (clip.config['exportLoc'],epic,epic)
            mpp.plot_multipages(outfile, clip, outfile)
            
            #fig = plt.figure(1, figsize=figuresize, dpi=dotperinch)
            #pp.summaryPlot1(clip)
            
            #file2="%s/%s/%s-modshift.pdf" % (clip.config.exportLoc,epic,epic)
            #file3="%s/%s/%s-onepage.pdf" % (clip.config.onepageBasename,epic,epic)
            cmd=''
            #cmd="gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s/%s/%s-all.pdf %s  %s" % (clip.config.exportLoc,epic,epic, outfile,file2)
            #cmd="pdftk %s %s %s output %s/%s/%s-all.pdf" % (outfile,file2,file3,clip.config['exportLoc'],epic,epic)
            print(cmd)
            os.system(cmd)     
            
    except (KeyError,AttributeError,TypeError) as e:
        cmd="None"
        print(epic, e)
        print("No Exports")
    
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
        if clip.disposition.isCandidate:
            clip=stel.addStellarToClip(clip)
            clip=stel.estimatePlanetProp(clip)
            outtxt,hdr = createExportString(clip, delimiter=" ", badValue="nan")
        else:
            outtxt="none"
    except (KeyError,AttributeError) as e:
        outtxt=e
        print(clip.value,e)
    
    return outtxt
    
    
    
