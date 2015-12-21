# Code to make diagnostic plots 

import sys
import Gnuplot, Gnuplot.funcutils
import numpy

def onepage(basename,time,rawflux,detrendflux,modelflux):

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
    model
        The fitted model

    Returns:
    -------------
    Nothing at present. Maybe an error code?

    Output:
    ----------
    A pdf of plots
    """

    # Remove NAN's
    #badvals = ~numpy.isnan(rawflux)
    #time        =        time[badvals]
    #rawflux     =     rawflux[badvals]
    #detrendflux = detrendflux[badvals]
    #modelflux   =   modelflux[badvals]
    
    #time        =        time[~flags]
    #rawflux     =     rawflux[~flags]
    #detrendflux = detrendflux[~flags]
    #modelflux   =   modelflux[~flags]

    #print(len(rawflux))

    g = Gnuplot.Gnuplot(persist=1)
    g("set terminal pdfcairo font ',6'")
    g("set output '"+basename+".pdf'")
    g("set xlabel 'Time'")
    g("set format y '%4.2E'")
    g("set ylabel 'Raw Flux'")
    g("set multiplot layout 3,1")
    #g("set 
    g.plot(Gnuplot.Data(time,rawflux, with_='points pt 7 lc 7 ps 0.05'))
    g("set ylabel 'Detrended Flux'")
    g.plot(Gnuplot.Data(time,detrendflux, with_='points pt 7 lc 7 ps 0.5'),Gnuplot.Data(time,modelflux, with_='lines lt 1 lc 1 lw 1')  )
    g("unset multiplot")



    
    