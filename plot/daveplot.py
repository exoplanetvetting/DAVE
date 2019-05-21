# Code to make diagnostic plots

import sys
import PyGnuplot as Gnuplot#Gnuplot, Gnuplot.funcutils
import numpy
import math

def onepage(pdfname, epicname, time, rawflux, detrendflux, modelflux, \
    period, epoch, duration):

    """Make diagnostic plots

    Inputs:
    -------------
    basename
        Base name of the output pdf file.
    time
        The full time array.
    fullflux
        The full raw flux time series.
    detrendflux
        The detrended flux time series.
    modelflux
        The fitted model flux time series.
    period
        The period of the system in days.
    epoch
        The epoch of the system in days.
    duration
        The duration of the transit in days.


    Returns:
    -------------
    Nothing at present. Maybe an error code?

    Output:
    ----------
    A pdf of plots
    """

    # Convert detrended and model fluxs to ppm
    detrendflux = detrendflux*1E6
    modelflux = modelflux*1E6


    # Prepare phased light curves
    phase = (time-epoch)/period - numpy.around((time-epoch)/period)
    p = phase.argsort()
    phase = phase[p]
    detrendfluxphased = detrendflux[p]
    modelfluxphased = modelflux[p]

    # Take the phased model and repeat it for each cycle so that the full time-series model light curve is full phased model repeated, and not just the model value for each individual point
    tstart = time[0]
    tend   = time[len(time)-1]
    fullphase = numpy.empty(0)
    fullmod   = numpy.empty(0)
    fphase = (tstart-epoch)/period - numpy.around((tstart-epoch)/period)  # Phase of the first data point
    for t in numpy.arange(tstart,tend,period):
        for p in phase:
            fullphase = numpy.append(fullphase,t + p*period - fphase*period)
        fullmod   = numpy.append(fullmod,modelfluxphased+ (tstart-epoch)/period)

    # Extend phase array to repeat two cycles so I can plot from -0.25 to 1.25 later on
    phase = numpy.append(phase,phase+1)
    detrendfluxphased = numpy.append(detrendfluxphased,detrendfluxphased)
    modelfluxphased = numpy.append(modelfluxphased,modelfluxphased)

    # Figure out odd and even phased
    dobperiod = 2*period
    oddphase = (time-epoch)/dobperiod - numpy.around((time-epoch)/dobperiod)
    p = oddphase.argsort()
    oddphase = oddphase[p]
    oddphase = 2*oddphase
    odddetrendfluxphased = detrendflux[p]
    oddmodelfluxphased = modelflux[p]

    evnphase = (time-epoch+period)/dobperiod - numpy.around((time-epoch+period)/dobperiod)
    p = evnphase.argsort()
    evnphase = evnphase[p]
    evnphase = 2*evnphase
    evndetrendfluxphased = detrendflux[p]
    evnmodelfluxphased = modelflux[p]



    # Start plotting
    g = Gnuplot.Gnuplot(persist=1)
    g("set terminal pdfcairo size 11in,8.5in font 'Droid Sans Mono,12'")
    g("set output '"+pdfname+"-onepage.pdf'")
    g("set xlabel 'Time (BKJD)'")
    g("set format y '%7.1E'")
    g("set ylabel 'Raw Flux'")
    g("set multiplot layout 4,1 title 'EPIC " + str(epicname) + ", P = " + str(period)[:8] + " days, E = " + str(epoch)[:11] + " BKJD' font ',22'")

    # Raw Flux
    g("set xrange [] writeback")
    g.plot(Gnuplot.Data(time,rawflux, with_='points pt 7 lc 7 ps 0.25'))

    # Detrended full time-series will full model repeated
    g("set format y '%7.0f'")
    g("set ylabel 'Detrended Flux (ppm)'")
    g("set xrange restore")
    g.plot(Gnuplot.Data(fullphase,fullmod, with_='lines lt 1 lc 1 lw 3'),Gnuplot.Data(time,detrendflux, with_='points pt 7 lc 7 ps 0.25'))

    # Phased full light curve
    g("set xlabel 'Phase'")
    g("set xrange [-0.25 to 1.25]")
    g.plot(Gnuplot.Data(phase,modelfluxphased, with_='lines lt 1 lc 1 lw 10'),Gnuplot.Data(phase,detrendfluxphased, with_='points pt 7 lc 7 ps 0.25'))

    # Phased Zoomed Primary
    g("set size 0.36,0.25")
    g("set origin 0.3,0.0")
    g("set label 1 'All' at graph 0.5, graph 0.9 center")
    g("set xrange [" + str(-2.5*duration/period) + " to " + str(2.5*duration/period) + "]")
    g("set yrange [] writeback")
    g("set y2range [] writeback")
    g("set ylabel ' '")
    g("set ytics ('       ' 0)")
    g.plot(Gnuplot.Data(phase,modelfluxphased, with_='lines lt 1 lc 1 lw 10'),Gnuplot.Data(phase,detrendfluxphased, with_='points pt 7 lc 7 ps 0.25'))

    # Odd transits
    g("set size 0.36,0.25")
    g("set origin 0.0,0.0")
    g("set label 1 'Odd' at graph 0.5, graph 0.9 center")
    g("set xrange [" + str(-2.5*duration/period) + " to " + str(2.5*duration/period) + "]")
    g("set yrange restore")
    g("set ytics auto")
    g("set ylabel 'Detrended Flux (ppm)'")
    g("set format y '%7.0f'")
    g("unset y2label")
    g.plot(Gnuplot.Data(oddphase,oddmodelfluxphased, with_='lines lt 1 lc 1 lw 10'),Gnuplot.Data(oddphase,odddetrendfluxphased, with_='points pt 7 lc 7 ps 0.25'))

    # Even transits
    g("set size 0.36,0.25")
    g("set origin 0.64,0.0")
    g("set label 1 'Even' at graph 0.5, graph 0.9 center")
    g("unset ytics")
    g("unset ylabel")
    g("set y2label 'Detrended Flux (ppm)'")
    g("set y2tics mirror")
    g("set yrange restore")
    g("set y2range restore")
    g("set format y '%7.0f'")
    g("set xrange [" + str(-2.5*duration/period) + " to " + str(2.5*duration/period) + "]")
    g.plot(Gnuplot.Data(evnphase,evnmodelfluxphased, with_='lines lt 1 lc 1 lw 10'),Gnuplot.Data(evnphase,evndetrendfluxphased, with_='points pt 7 lc 7 ps 0.25'))


    g("unset multiplot")




