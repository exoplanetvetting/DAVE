
import matplotlib.pyplot as plt
import numpy as np


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


def plotData(clip):
    fl = clip['detrend.flags']
    time = clip['serve.time']
    raw = clip['extract.rawLightcurve']
    flux = clip['detrend.flux_frac']

    raw = raw/np.median(raw[~fl]) - 1

    plt.plot(time[~fl], raw[~fl], 'r-')
    plt.plot(time[~fl], flux[~fl], 'k.')

def plotFolded(clip):
    fl = clip['detrend.flags']
    time = clip['serve.time']
    flux = clip['detrend.flux_frac']
    period = clip['trapFit.period_days']

    phi = np.fmod(time- .25*period, period)
    plt.plot(phi[~fl], flux[~fl], 'k.')

    #Fudge to account for the fact that bls doesn't return the *right* model
    model = clip['trapFit.bestfitModel']
    x = phi[~fl]
    n = min(len(x), len(model))
    plt.plot(x[:n], model[:n], 'r.')

