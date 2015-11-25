from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits as pyfits
import glob

from astropy.stats.funcs import median_absolute_deviation as MAD
from scipy.ndimage import label
from photo_test import raw_moment, intertial_axis, plot_bars

from kepfit import lsqclip
from martinsff import martinsff

"""
this is the code to make light curves for Debra and the planet hunters
I'm going to write it so it will work on individual files
"""

def bg_sub(fla):
    """
    subtract the background from a series of images
    by assuming the aperture is large enough to be
    predominantly background
    """
    for i in xrange(np.shape(fla)[0]):
        fla[i,:,:] = fla[i,:,:] - np.nanmedian(fla[i,:,:])
    return fla

def run_C0_data_extract(fn, second_half_only=True,
    qual_cut=False, return_qual=False, toss_resat=True,
    bg_cut=5, skip_cads=None, skip=None):

    #for C0, lets just ise the good data post safe mode
    if second_half_only and skip_cads == None:
        skip = 1976
    elif skip is not None:
        skip = skip_cads
    else:
        skip = 0

    with pyfits.open(fn) as f:
        time = f[1].data['TIME'][skip:] - f[1].data['TIME'][0]
        fluxarr = f[1].data['FLUX'][skip:]
        quality = f[1].data['QUALITY'][skip:]

    if qual_cut:
        time = time[quality == 0]
        fluxarr = fluxarr[quality == 0,:,:]
    elif toss_resat:
        # data the cadences where there is a wheel
        # resetuation event
        time = time[quality != 32800]
        fluxarr = fluxarr[quality != 32800,:,:]

    # fix dodgy data
    # the first C0 data release included zeros
    # this will be changed later but we need this
    # fix for now
    fluxarr[fluxarr == 0] = np.nan

    #subtract background
    flux_b = bg_sub(fluxarr)

    # create a median image to calculate where
    # the pixels to use are: is this good?
    # this needs to be investiagated.
    flatim = np.nanmedian(flux_b,axis=0)

    #fixd pixels that are X MAD above the median
    vals = flatim[np.isfinite(flatim)].flatten()
    mad_cut = 1.4826 * MAD(vals) * bg_cut

    region = np.where(flatim > mad_cut,1,0)
    lab = label(region)[0]

    #find the central pixel
    imshape = np.shape(flatim)
    centralpix = [1+imshape[0] // 2,1+imshape[1] // 2]

    #find brightest pix within 9x9 of central pix
    centflatim = flatim[centralpix[0]-4:centralpix[0]+4,
                centralpix[1]-4:centralpix[1]+4]
    flatimfix = np.where(np.isfinite(centflatim),centflatim,0)
    brightestpix = np.unravel_index(flatimfix.argmax(), centflatim.shape)
    bpixy, bpixx = brightestpix

    regnum = lab[centralpix[0]-4+bpixy,centralpix[1]-4+bpixx]
    if regnum == 0:
        print('WARNING, no star was found in light curve, \
            {} light curve will be junk!'.format(fn))

    lc = np.zeros_like(time)
    xbar = np.zeros_like(time)
    ybar = np.zeros_like(time)

    #there is a loop that performs the aperture photometry
    #lets also calcualte the moments of the image

    #make a rectangular aperture for the moments thing
    ymin = np.min(np.where(lab == regnum)[0])
    ymax = np.max(np.where(lab == regnum)[0])
    xmin = np.min(np.where(lab == regnum)[1])
    xmax = np.max(np.where(lab == regnum)[1])
    #print(np.where(lab == regnum))
    momlims = [ymin,ymax+1,xmin,xmax+1]

    for i,fl in enumerate(fluxarr):
        lc[i] = np.sum(fl[lab == regnum])
        momim = fl[momlims[0]:momlims[1],
                    momlims[2]:momlims[3]]
        momim[~np.isfinite(momim)] == 0.0
        xbar[i], ybar[i], cov = intertial_axis(momim)

    if return_qual:
        return None
    else:
        return (time,lc, xbar /
            np.mean(xbar), ybar / np.mean(xbar), regnum)


def run_C0_detrend(time, lc, xbar, ybar, skip, cadstep=200):
    #some reasonable defaults
    npoly_cxcy = 3 # Order of ploynomial fit to target centroids
    sigma_cxcy = 10.0 # Sigma-clipping threshold for fit to target centroids [sigma]
    npoly_ardx = 5 # Order of ploynomial fit for thruster firing detection
    npoly_dsdt = 1 # Order of ploynomial fit for thruster firing detection
#    sigma_dsdt = 1.5 # Sigma-clipping threshold for thruster firing detection [sigma]
    sigma_dsdt = 3.0 # Sigma-clipping threshold for thruster firing detection [sigma]
    npoly_arfl = 3 # Order of ploynomial for for arclength-flux calibration
    sigma_arfl = 3.0 # Sigma-clipping threshold for arclength-flux calibration [sigma]

    #normalize lc
    lc = (lc / np.median(lc))
    #fudge some yerr
    yerr = np.ones_like(lc) * MAD(lc) * 1.4826

    # we're going to step through the light curve
    # i don't like the data at the start but love it at
    # the end so we're going to use
    npt = len(time)
    outflux = np.array([])
    outcorr = np.array([])
    thr_cad = np.array([],dtype=bool)
    for t in np.arange(npt%cadstep,npt,cadstep):
        tran = slice(t,t+cadstep)

        intime = time[tran]
        indata = lc[tran]
        centr1 = xbar[tran]
        centr2 = ybar[tran]

        outf, outc, th_c = martinsff(intime,indata,
                centr1,centr2,
                npoly_cxcy,sigma_cxcy,npoly_ardx,
                npoly_dsdt,sigma_dsdt,npoly_arfl,
                sigma_arfl,2,'crapfile',
                0)
        outflux = np.r_[outflux,outf]
        outcorr = np.r_[outcorr,outc]
        thr_cad = np.r_[thr_cad,th_c]

    return outflux,outcorr, thr_cad

import untrendy
def medfilt(time,flux,window=10.0):
    """
    window in days
    """
    flux = flux / untrendy.median(
        time,flux,dt=window)
    return flux


if __name__ == '__main__':
    ## test on a single file

    filename = glob.glob(
        '/Volumes/K2_ENG/C0_all/ktwo202126847-c00_lpd-targ.fits.gz')[0]
    time,lc,xbar,ybar,regnum = run_C0_data_extract(filename)

    mask = np.ones(len(lc),dtype=bool)
    mask[range(11,len(lc),12)] = False
    time,lc,xbar,ybar = time[mask], lc[mask], xbar[mask], ybar[mask]

    #lets flatten the light curve
    flatlc = medfilt(time,lc,window=1.5)
    zpt = len(time)%400.

    outflux, correction, thr_cad = run_C0_detrend(
        time,flatlc,xbar,ybar,cadstep=400)

    not_thr = ~thr_cad

    corflux = (lc[zpt:][not_thr]/
                np.median(lc[zpt:][not_thr])/
                correction[not_thr])
    corflatflux = (flatlc[zpt:][not_thr]/
                np.median(flatlc[zpt:][not_thr])/
                correction[not_thr])


    outname = filename.split(
        '/')[-1].split(
        'pd-targ.fits.gz')[0]
    np.savetxt(
        '/Users/tom/Projects/Debra_data/data/{}.txt'.format(
            outname),
        np.array([time[zpt:][not_thr],
                    corflux,corflatflux,
                    correction[not_thr]]).T)


    plt.figure()
    plt.plot(time[zpt:][not_thr],
        lc[zpt:][not_thr]/np.median(lc[zpt:][not_thr])
        /correction[not_thr])
    plt.scatter(time[zpt:],
        flatlc[zpt:],s=2,color='b')

















