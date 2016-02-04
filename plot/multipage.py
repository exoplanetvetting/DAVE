# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:55:11 2016

@author: sthomp

A multi plot creator to create Vetting type documents for our 
detected signals.
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import sueplotting as sp
 
 
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
    pdf_pages.savefig(fig)
    plt.close()
    
    fig = plt.figure(figsize=figuresize, dpi=dotperinch)
    # Plot whatever you wish to plot
    sp.summaryPlot(clip)
    # Done with the page
    pdf_pages.savefig(fig)
    plt.close()
 
    fig = plt.figure(figsize=figuresize, dpi=dotperinch)
    sp.indivPlot(clip,6)    
    pdf_pages.savefig(fig)

    fig = plt.figure(figsize=figuresize, dpi=dotperinch)
    sp.blsPlot(clip)    
    pdf_pages.savefig(fig)
    plt.close()

    fig =  plt.figure(figsize=figuresize, dpi=dotperinch)  
    plt.figtext(0.2,0.35,clip.disposition,color='b',fontsize=14)
    plt.title('Disposition Information in Clipboard')
    pdf_pages.savefig(fig)
    plt.close()
    
    # Write the PDF document to the disk
    pdf_pages.close()
    
    plt.close()
    