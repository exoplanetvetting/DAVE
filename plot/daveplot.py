# Code to make diagnostic plots 

import sys
import Gnuplot, Gnuplot.funcutils
import numpy
import math

def onepage(basename,time,rawflux,detrendflux,modelflux,period,epoch):

    """Make diagnostic plots

    Inputs:
    -------------
    basename
        Base name of the output pdf
    time
        The full time array
    fullflux
        The full raw flux time series 
    detrendflux
        The detrended flux time series
    modelflux
        The fitted model flux time series
    period
        The period of the system in days.
    epoch
        The epoch of the system in days.
        

    Returns:
    -------------
    Nothing at present. Maybe an error code?

    Output:
    ----------
    A pdf of plots
    """


    

    g = Gnuplot.Gnuplot(persist=1)
    g("set terminal pdfcairo font 'Droid Sans Mono,6'")
    g("set output '"+basename+".pdf'")
    g("set xlabel 'Time (BKJD)'")
    g("set format y '%8.1E'")
    g("set ylabel 'Raw Flux'")
    g("set multiplot layout 3,1")
    #g("set 
    g.plot(Gnuplot.Data(time,rawflux, with_='points pt 7 lc 7 ps 0.05'))
    g("set ylabel 'Detrended Flux'")
    g.plot(Gnuplot.Data(time,modelflux, with_='lines lt 1 lc 1 lw 1'),Gnuplot.Data(time,detrendflux, with_='points pt 7 lc 7 ps 0.05'))
    
    
    # Plot phased data
    phase = (time-epoch)/period - numpy.around((time-epoch)/period)
    p = phase.argsort()
    phase = phase[p]
    detrendflux = detrendflux[p]
    modelflux = modelflux[p]
    
    g("set xlabel 'Phase'")
    g.plot(Gnuplot.Data(phase,modelflux, with_='lines lt 1 lc 1 lw 3'),Gnuplot.Data(phase,detrendflux, with_='points pt 7 lc 7 ps 0.05'))
    g("unset multiplot")



    
    