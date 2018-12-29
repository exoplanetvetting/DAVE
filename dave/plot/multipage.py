# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:55:11 2016

@author: sthomp

A multi plot creator to create Vetting type documents for our 
detected signals.
This code is depricated, see pipeline.multiPagePlot
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import dave.pipeline.plotting as pp
import dave.diffimg.plot as dip 
 
def plot_all_multipages(outfile,clip,intext):
    """Take a clipboard, clip, and create plots
    put these plots all into one multi paged document
    specifieed by outfile
    """
    
    dotperinch=300    
    figuresize=(10,8)
    # The PDF document
    pdf_pages = PdfPages(outfile)
      # Create a figure instance (ie. a new page) 
    #pdf_pages.attach_note(('KIC %u   [%u]' % (clip.value, clip.disposition.isCandidate)),positionRect=[100,200,10,400])
    fig =  plt.figure(figsize=figuresize, dpi=dotperinch)  
    plt.figtext(0.5,0.5,intext,color='r',fontsize=15)
    fig.patch.set_visible(False)
    plt.gca().axis('off')
    pdf_pages.savefig(fig)
    plt.close()
    
    fig = plt.figure(figsize=figuresize, dpi=dotperinch)
    # Plot whatever you wish to plot
    pp.summaryPlot1(clip)
    # Done with the page
    pdf_pages.savefig(fig)
    plt.close()
 
    fig = plt.figure(figsize=figuresize, dpi=dotperinch)
    pp.indivTransitPlot(clip,7)    
    pdf_pages.savefig(fig)
    plt.close()

    try:
        fig = plt.figure(figsize=figuresize, dpi=dotperinch)
        pp.blsPlot(clip)    
        pdf_pages.savefig(fig)
        plt.close()
    except AttributeError:
        pass

    try:    
        fig =  plt.figure(figsize=figuresize, dpi=dotperinch)
        plt.figtext(0.2,0.35,clip.disposition,color='b',fontsize=14)
        plt.title('Disposition Information in Clipboard')
        pdf_pages.savefig(fig)
        plt.close()
    except:
        pass
    try:    
        #Plot centroid plots
        (fig1,fig2)=dip.plotWrapper(clip)
    except:
        fig1=plt.plot()
        fig2=plt.plot()
        
    pdf_pages.savefig(fig2)
    pdf_pages.savefig(fig1)
    plt.close()
    plt.close()

    
    # Write the PDF document to the disk
    pdf_pages.close()
    
    plt.close()

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
