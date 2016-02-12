
import matplotlib.pyplot as plt
import numpy as np

import dave.fileio.kplrfits as kplrfits

def plotDiagnosticLightcurves(time, rawLc, cotrendLc, detrendLc, path="."):

    plt.clf()
    plt.gcf().set_size_inches((8,10))
    plt.subplot(311)
    plt.plot(time, rawLc, 'k.')
    plt.ylabel("Raw Lightcurve")

    plt.subplot(312)
    plt.plot(time, cotrendLc, 'k.')
    plt.ylabel("Cotrended Lightcurve")

    plt.subplot(313)
    plt.plot(time, detrendLc, 'k.')
    plt.xlabel('Time (BKJD)')
    plt.ylabel("Detrended Lightcurve")


def plotData(clip, nPanel=3):
    """Plot pa and pdc lightcurves"""
    epic = clip['value']
    campaign = clip['config.campaign']
    fl = clip['detrend.flags']
    time = clip['serve.time']
    raw = clip['extract.rawLightcurve'] / 1000.
    flux = clip['detrend.flux_frac']

    markTransits = True
    try:
        per = clip['trapFit.period_days']
        epc = clip['trapFit.epoch_bkjd']
        dur_days = clip['trapFit.duration_hrs']/24.
    except KeyError, e:
        print "WARN: %s" %(e)
        markTransits = False

    fig = plt.gcf()
    fig.set_size_inches(11, 8.5)

    colour = ["#FFF8F8", "#F8FFF8", "#F4F4FF"]
    start = np.min(time[~fl])
    deltaT = np.max(time[~fl]) - start
    deltaT /= float(nPanel)

    rawRange = np.percentile(raw[~fl], [1,99])

    for i in range(nPanel):
        ax = plt.subplot(2*nPanel, 1, 2*i+1, axisbg=colour[i])
        plt.plot(time[~fl], raw[~fl], 'ko', ms=2, alpha=.8)
        plt.ylim(rawRange)
        plt.ylabel("Raw flux/1000")

        if markTransits:
            plotTransitRegions(time[~fl], per, epc, dur_days)

        plt.subplot(2*nPanel, 1, 2*i+2, sharex=ax, axisbg=colour[i])
        #Plotting bad data cadences turned off
#        plt.plot(time[fl], 0*time[fl], 'mo', ms=8, mec="none")
        plt.plot(time[~fl], flux[~fl], 'ko', ms=2, alpha=.8)
        plt.ylabel("Detrended")


        if markTransits:
            plotTransitRegions(time[~fl], per, epc, dur_days)

        plt.xlim(start + i*deltaT, start + (i+1)*deltaT)
        plt.xlabel("Time (BKJD)")

    plt.suptitle("EPIC %i    Campaign %i" %(epic, campaign))

def plotTransitRegions(time, period_days, epoch_bkjd, duration_days, **kwargs):
    tmin = np.min(time)
    tmax = np.max(time)

    n1 = int(np.floor( (tmin-epoch_bkjd)/period_days) )
    n2 = int(np.ceil(  (tmax-epoch_bkjd)/period_days) )
    print n1, n2
    color = kwargs.pop('color','#888888')
    alpha = kwargs.pop('alpha', 1)

    for n in range(n1, n2+1):
        t0 = epoch_bkjd + n*period_days
        lwr = t0 - .5*duration_days
        upr = t0 + .5*duration_days
#        print lwr, upr
        plt.axvspan(lwr, upr, color=color, alpha=alpha)

def plotFolded(clip, doublePeriod = False):
    fl = clip['detrend.flags']
    time = clip['serve.time']

#    tce = clip['eventList'][0]
    tce = clip
    flux = clip['detrend.flux_frac']
    period = tce['trapFit.period_days']
    epoch = tce['trapFit.epoch_bkjd']
#    model = tce['trapFit.bestFitModel']

    if doublePeriod:
        period *= 2
#    phi = np.fmod(time- epoch + .25*period, period)
    if epoch > time[0]:
        diff = epoch-time[0]
        epoch -= period*np.ceil(diff/period)
        print "Reducing epoch"
        print

    plt.clf()
    print epoch, time[0]
    phi = np.fmod(time-epoch + .25*period, period)
#    phi = np.fmod(time, period)
    plt.plot(phi[~fl], 1e6*flux[~fl], 'ko', ms=4)
    plt.plot(period+ phi[~fl], 1e6*flux[~fl], 'o', color='#888888', ms=4, mec="none")

#    x = phi[~fl]
#    y = 1e6*model[~fl]
#    idx = np.argsort(x)
#
#    plt.plot(x[idx], y[idx], 'r-')
#    plt.plot(period+x[idx], y[idx], 'r-')
#    plt.ylabel("Fractional Amplitude (ppm)")
#    plt.xlabel("Phase (days)")


