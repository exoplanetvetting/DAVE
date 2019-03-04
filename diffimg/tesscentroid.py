# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Thu Nov  8 16:19:30 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np

from dave.fileio.plateau import plateau
import dave.fileio.kplrfits as kplrfits
import dave.fileio.pyfits as pyfits
import dave.fileio.tpf as ktpf
#import kepler.apj as apj
import psffit


def show():
    path = '/Users/vkostov/Desktop/Ideas_etc/DAVE_test/data/hlsp_tess-data-alerts_tess_phot_00307210830-s02_tess_v1_tp.fits'
    fits, hdr = pyfits.getdata(path, header=True)

    cube = ktpf.getTargetPixelArrayFromFits(fits, hdr)

    for i in range(1000, 1002):
        plt.clf()
        mn = np.fabs(np.min(cube[i,:,:])) + 1
        plt.imshow(np.log10(cube[i,:,:]), origin="bottom")
#        plt.imshow(cube[i,:,:], origin="bottom", cmap=plt.cm.bone)
#        plt.clim(-20, 100)
        plt.colorbar()
        plt.title(i)
        plt.pause(1)

def tic_307210830_02_01():
    """First TCE on TIC 307210830 in second sector """
    tic = 307210830
    sector = 2

    period_days = 3.69061
    epoch_btjd = 1356.2038
    duration_days = 1.2676/24.

    main(tic, sector, period_days, epoch_btjd, duration_days)


def tic_307210830_02_03():
    """First TCE on TIC 307210830 in second sector """
    tic = 307210830
    sector = 2

    period_days = 2.25301
    epoch_btjd = 1355.2867
    duration_days = 1.0185/24.

    outpattern = "t%11i-s%02i-c03" %(tic, sector)

    main(tic, sector, period_days, epoch_btjd, duration_days, outpattern)



def main(tic, sector, period_days, epoch_btjd, duration_days, outpattern):

    path = '/Users/vkostov/Desktop/Ideas_etc/DAVE_test/data/hlsp_tess-data-alerts_tess_phot_%011i-s%02i_tess_v1_tp.fits'
    path = path %(tic, sector)
    fits, hdr = pyfits.getdata(path, header=True)
    cube = ktpf.getTargetPixelArrayFromFits(fits, hdr)
    cube = cube[:, 3:9, 2:8]

    time = fits['TIME']
    isnan = np.isnan(time)
    time = time[~isnan]
    cube = cube[~isnan]

    transits = getIngressEgressCadences(time, period_days, epoch_btjd, duration_days)
    with open('%s.cent.txt' %(outpattern), 'w') as fp:
        for i in range(len(transits)):
            print("Transit %i" %(i))
            cin = transits[i]
            res = measureCentroidShift(cube, cin, True)

            plt.suptitle('%s-trans%02i' %(outpattern, i))
            plt.savefig('%s-trans%02i.png' %(outpattern, i))

            pattern = "%.6f " * len(res)
            pattern = pattern + "\n"
            fp.write( pattern % tuple(res))

def plotCentroids(fn):
    plt.clf()
#    apj.pre()
    data = np.loadtxt(fn)
    plt.plot(data[:,0], data[:,1], 'ko', label="Before")
    plt.plot(data[:,2], data[:,3], 'ro', label="Difference")
    plt.plot(data[:,4], data[:,5], 'co', label="After")

    for i in range(len(data)):
        plt.text(data[i,2], data[i,3], ' %i' %(i))

    plt.xlabel("Column")
    plt.ylabel("Row")
    plt.title(fn)
#    apj.post()
#    apj.pgid()

def measureCentroidShift(cube, cin, plot=True):
    before, after, diff = generateDiffImg(cube, cin, plot=plot)
    plt.pause(.01)

    print("Before...")
    guess = pickInitialGuess(before)
    beforeSoln = psffit.fitPrf(before, psffit.gaussianWithConstantSkyPrf, guess)

    print("Diff...")
    guess = pickInitialGuess(diff)
    diffSoln = psffit.fitPrf(diff, psffit.gaussianWithConstantSkyPrf, guess)

    print("After...")
    guess = pickInitialGuess(after)
    afterSoln = psffit.fitPrf(after, psffit.gaussianWithConstantSkyPrf, guess)

    if not np.all( map(lambda x: x.success, [beforeSoln, diffSoln, afterSoln]) ):
        print("WARN: Not all fits converged for [%i, %i]" %(cin[0], cin[1]))

    out = []
    out.extend(beforeSoln.x[:2])
    out.extend(diffSoln.x[:2])
    out.extend(afterSoln.x[:2])
    return out

def pickInitialGuess(img):

    r0, c0 = np.unravel_index( np.argmax(img), img.shape)

    guess = [c0+.5, r0+.5, .5, np.max(img), np.median(img)]
    return guess

def getIngressEgressCadences(time, period_days, epoch_btjd, duration_days):
    assert np.all(np.isfinite(time))

    idx = kplrfits.markTransitCadences(time, period_days, epoch_btjd, duration_days)
    transits = np.array(plateau(idx, .5))

    return transits


def generateDiffImg(cube, transits, plot=False):
    """Generate a difference image.

    Also generates an image for each the $n$ cadedences before and after the transit,
    where $n$ is the number of cadences of the transit itself

    Inputs
    ------------
    cube
        (np 3 array) Datacube of postage stamps
    transits
        (2-tuples) Indices of the first and last cadence

    Optional Inputs
    -----------------
    plot
        (Bool) If true, generate a diagnostic plot


    Returns
    -------------
    Three 2d images,

    before
        The sum of the n cadences before transit (where n is the number of in-transit cadences
    after
        The sum of the n cadences after transit
    diff
        The difference between the flux in-transit and the average of the flux before and after

    Notes
    ---------
    When there is image motion, the before and after images won't be identical, and the difference
    image will show distinct departures from the ideal prf.
    """

    dur  = transits[1] - transits[0]
    s0, s1 = transits - dur
    e0, e1 = transits + dur

    before = cube[s0:s1].sum(axis=0)
    during = cube[transits[0]:transits[1]].sum(axis=0)
    after = cube[e0:e1].sum(axis=0)

    diff = .5 * (before + after) - during
#    diff = before - during
#    diff = after - before

    if plot:
        plt.clf()
        plt.subplot(221)
        plt.imshow(before, origin='bottom')
        plt.title("Before")
        plt.colorbar()

        plt.subplot(222)
        plt.imshow(after, origin='bottom')
        plt.title("After")
        plt.colorbar()

        plt.subplot(223)
        plt.imshow(after - before, origin='bottom', cmap=plt.cm.RdYlBu_r)
        plt.title("After - Before")
        plt.colorbar()

        plt.subplot(224)
        plt.imshow(diff, origin='bottom', cmap=plt.cm.RdYlBu_r)
        plt.title("Diff")
        plt.colorbar()

    return before, after, diff
