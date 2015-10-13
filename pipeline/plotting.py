
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

