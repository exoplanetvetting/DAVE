#Tools for manipulating target pixel files
#Tom Barclay
#Fergal Mullally
#Knicole Colon


#Light Curve Extraction from TPFs (Optimal Apertures)
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import matplotlib
import dave.fileio.mastio as mastio

from astropy.io import fits as pyfits
from photo_test import raw_moment, intertial_axis, plot_bars

from astropy.stats.funcs import median_absolute_deviation as MAD
from scipy.ndimage import label


def getTargetPixelArrayFromFits(fits, hdr, column='FLUX'):
    nRows = int(hdr['NAXIS2'])
    #All imgs are the same size, so it doesn't matter which TDIM
    #I pick here.
    shape = eval(hdr['TDIM5'])
    #The reverse puts the columns on the x-axis and row
    #on the y-axis when you call mp.imshow()
    shape = tuple(reversed(shape))

    tpfArray = np.empty((nRows, shape[0], shape[1]))

    for i, cadence in enumerate(fits[column]):
        tpfArray[i,:,:] = cadence.reshape(shape)

    return tpfArray


def lightCurveFromTpfCube(cube, mask=None):
    """Extract a lightcurve from a TPF cube.

    Inputs:
    cube    (3d numpy array) As returned by getTargetPixelArrayFromFits()
    mask    (2d numpy array) Which pixels to use in extraction. Default
            is all finite pixels

    Returns
    A 1d numpy array

    Notes:
    Pixels marked Nan (i.e outside the aperture) are ignored regardless
    of the mask value.

    Default mask includes all finite pixels.

    To extract only the optimal aperture use
    mask = np.where(mask == 3, 1, 0)
    lightCurveFromTpfCube(cube, mask)
    """

    if mask is None:
        mask = np.ones( cube[0,:,:].shape)

    n, nR, nC = cube.shape
    out = np.zeros(n)
    for i in range(n):
        img = cube[i, :, :]
        idx = np.bitwise_and(mask > 0, np.isfinite(img))
        out[i] = np.sum(img[idx])

    return out


def bg_sub(fla):
    """
    subtract the background from a series of images
    by assuming the aperture is large enough to be
    predominantly background
    """
    for i in xrange(np.shape(fla)[0]):
        fla[i,:,:] = fla[i,:,:] - np.nanmedian(fla[i,:,:])
    return fla


def optimalAperture(t_time, t_fluxarr, t_quality, qual_cut=False, return_qual=False, toss_resat=False, bg_cut=5, skip=0):
    """
    This routine determines an optimal apertures and outputs the flux (i.e. a light curve) from a TPF.


    Inputs:
    ------------
    t_time = 1D array of 'TIME' from TPF
    t_fluxarr = 1D array of 'FLUX' from TPF
    t_quality = 1D array of 'QUALITY' from TPF
    qual_cut = exclude cadences with a non-zero quality flag; this is False by default
    return_qual = if True then nothing is returned; this is True by default
    toss_resat = exclude cadences where there is a wheel resaturation event; this is True by default
    bg_cut = threshold to find pixels that are bg_cut * MAD above the median
    skip = index of first cadence that should be used in the time series

    Outputs:
    ------------
    time = 1D array of time in BKJD
    lc = 1D array of flux measured in optimal aperture
    xbar = 1D array of x-coordinate of target centroid
    ybar = 1D array of y-coordinate of target centroid
    regnum = integer value of brightest pixel
    lab = 2D array identifying pixels used in aperture

    Usage:
    ------------
    tpf,tpf_hdr = ar.getLongTpf(k2id, campaign, header=True)

    tpf_time = tpf['TIME']
    tpf_flux = tpf['FLUX']
    tpf_quality = tpf['QUALITY']

    time,lc,xbar,ybar,regnum,lab = optimalAperture(tpf_time, tpf_flux, tpf_quality, qual_cut=False, return_qual=False, toss_resat=True, bg_cut=5, skip=0)
    """

    time = t_time[skip:]
    fluxarr = t_fluxarr[skip:]
    quality = t_quality[skip:]

    if qual_cut:
        time = time[quality == 0]
        fluxarr = fluxarr[quality == 0,:,:]
    elif toss_resat:
        # cadences where there is a wheel resaturation event
        time = time[quality != 32800]
        fluxarr = fluxarr[quality != 32800,:,:]

    #remove any nans
    fluxarr[fluxarr == 0] = np.nan

    #subtract background
    flux_b = bg_sub(fluxarr)

    # create a median image to calculate where the pixels to use are
    flatim = np.nanmedian(flux_b,axis=0)

    #find pixels that are X MAD above the median
    vals = flatim[np.isfinite(flatim)].flatten()
    mad_cut = 1.4826 * MAD(vals) * bg_cut
    
    flatim[np.isnan(flatim)] = 0.
    region = np.where(flatim > mad_cut,1,0)
    lab = label(region)[0]

    #find the central pixel
    imshape = np.shape(flatim)
    centralpix = [1+imshape[0] // 2,1+imshape[1] // 2]

    #find brightest pix within 9x9 of central pix
    #this assumes target is at center of postage stamp which I think is ok
    centflatim = flatim[centralpix[0]-4:centralpix[0]+4,
                centralpix[1]-4:centralpix[1]+4]
    flatimfix = np.where(np.isfinite(centflatim),centflatim,0)
    brightestpix = np.unravel_index(flatimfix.argmax(), centflatim.shape)
    bpixy, bpixx = brightestpix

    #use all pixels in the postage stamp that are X MAD above the median
    #this identifies location of brightest pixel only
    regnum = lab[centralpix[0]-4+bpixy,centralpix[1]-4+bpixx]

    lc = np.zeros_like(time)
    xbar = np.zeros_like(time)
    ybar = np.zeros_like(time)

    #make a rectangular aperture for the moments thing
    ymin = np.min(np.where(lab == regnum)[0])
    ymax = np.max(np.where(lab == regnum)[0])
    xmin = np.min(np.where(lab == regnum)[1])
    xmax = np.max(np.where(lab == regnum)[1])

    momlims = [ymin,ymax+1,xmin,xmax+1]

    #loop that performs the aperture photometry
    for i,fl in enumerate(flux_b):
        lc[i] = np.sum(fl[lab == regnum])
        #lc[i] = np.sum(fl[np.where(lab == 1)]
        momim = fl[momlims[0]:momlims[1],
                    momlims[2]:momlims[3]]
        momim[~np.isfinite(momim)] == 0.0
        xbar[i], ybar[i], cov = intertial_axis(momim)


    if return_qual:
        return None
    else:
        # TODO: think about whether this should be normalized
        return (time,lc, xbar - np.nanmean(xbar), ybar - np.nanmean(ybar), regnum, lab)

if __name__=="__main__":
    pass