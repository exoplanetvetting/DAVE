# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:55:11 2016

@author: sthomp

A multi plot creator to create Vetting type documents for our 
detected signals.
"""
import matplotlib as mpl
mpl.use('Agg')  #Use Agg backend to prevent X server use when headless
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import dave.pipeline.plotting as pp
import dave.diffimg.plot_vbk as dip 
import dave.stellar.readStellarTable as stel
from pdb import set_trace as bp
 
def plot_multipages(outfile,clip,intext,figtype="png"):
    """Create a line of text for the exporter
    
    Inputs:
    ------------
    outfile
        A name for the pdf file to be written to
    clip
        A clipboard object
    intext
        A string containing information to be put on the first page.
        
    
    Takes a clipboard, clip, and create plots
    put these plots all into one multi paged document
    specifieed by outfile
    """

    plt.ioff()
    dotperinch=120    
    figuresize=(11,8)
    # The PDF document
    pdfout=outfile + '.pdf'
    pdf_pages = PdfPages(pdfout)
      # Create a figure instance (ie. a new page) 
    #pdf_pages.attach_note(('KIC %u   [%u]' % (clip.value, clip.disposition.isCandidate)),positionRect=[100,200,10,400])
    
    fig =  plt.figure(1, figsize=figuresize, dpi=dotperinch)  
   
    plt.figtext(0.2,0.85,intext,color='r',fontsize=15)
    plt.figtext(0.15,0.2,clip.disposition,color='b',fontsize=14)
    plt.title('Disposition Information for EPIC %u' % clip['value'])
    
    if clip.disposition.isSignificantEvent:   
         steltxt=clip.get('stellar',defaultValue='No Stellar Available')
         plt.figtext(0.7,0.3,steltxt)
         plantxt=clip.get('planet',defaultValue='No Planet Param. Avail.')
         plt.figtext(0.7,0.2,plantxt)

    fig.patch.set_visible(False)
    plt.gca().axis('off')
    plt.savefig(pdf_pages,format='pdf')
    outname="%s-%s.%s" %(outfile,"info",figtype)
    
    plt.savefig(outname)
    plt.close()
    
    fig = plt.figure(1, figsize=figuresize, dpi=dotperinch)
    # Plot whatever you wish to plot
    pp.summaryPlot1(clip)
    # Done with the page
    pdf_pages.savefig()
    outname="%s-%s.%s" %(outfile,"summary",figtype)
    fig.savefig(outname)
    plt.close()

    fig = plt.figure(figsize=figuresize, dpi=dotperinch)
    # Plot whatever you wish to plot
    pp.plotData(clip, nPanel=3)
    pdf_pages.savefig()
    outname="%s-%s.%s" % (outfile,"data",figtype)
    fig.savefig(outname)
    plt.close()

 
    fig = plt.figure(figsize=figuresize, dpi=dotperinch)
    pp.indivTransitPlot(clip,6)    
    pdf_pages.savefig()
    outname="%s-%s.%s" %(outfile,"indiv",figtype)
    plt.savefig(outname)
    plt.close()

    try:
        fig = plt.figure(figsize=figuresize, dpi=dotperinch)
        pp.blsPlot(clip)    
        pdf_pages.savefig()
        outname="%s-%s.%s" % (outfile,"bls",figtype)
        plt.savefig(outname)
        plt.close()
    except AttributeError:
        pass

    
    try:    
        #Plot centroid plots
        (fig2,fig1)=dip.plotWrapper(clip)
    except ValueError:
        fig1=plt.figure()
        fig2=plt.figure()

    fig1.set_size_inches(figuresize)
    pdf_pages.savefig(fig1, dpi=dotperinch)
    outname="%s-%s.%s" % (outfile,"cent1",figtype)
    fig1.savefig(outname)
    fig2.set_size_inches(figuresize)        
    pdf_pages.savefig(fig2, dpi=dotperinch)
    outname="%s-%s.%s" % (outfile,"cent2",figtype)
    fig2.savefig(outname)
    plt.close()


    fig=plt.figure(figsize=figuresize,dpi=dotperinch)
    pp.lppDiagnostic(clip)
    pdf_pages.savefig()
    outname="%s-%s.%s" % (outfile,"lpp",figtype)
    plt.savefig(outname)
    plt.close()    
        
    
    # Write the PDF document to the disk
    pdf_pages.close()
    
    plt.close()

