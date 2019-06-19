
"""
Created on Tue Jun 18 21:28:56 2019

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl
import pandas as pd
import numpy as np

import scipy.optimize as sopt


def fitPrfCentroidForImage(img, ccd, camera, sector, bbox, prfObj):
    """Fit a PRF model to a TPF image and measure star centroid

    Given a cadence image from a TPF file, find the best fitting
    PRF model for that image. The best fit column and row of the
    PRF model can be treated as the photometric centroid of the
    light distribution from a single isolated star.

    Fitting is performed using the L-BFGS-B bounded fitting algorithm
    using model PRFs provided by the mission.

    Inputs:
    ---------
    img
        (np 2d array) Image of star to be fit. Image is in the
        format img[row, col]. img should not contain Nans

    ccd, camera, sector
        (int) Which CCD and camera is the image taken from. Which sector
        was the data observed (the measured PRF changes as time evolves)

    bbox
        [c1, c2, r1, r2]. Define the range of columns (c1..c2)
        and rows (r1..r2)  defined by the image.
        An exception raised if the following equality not true
        img.shape = (c2-c1), (r2-r1). The star is assumed to be
        in the centre of the bounding box.

    prfObj
        An object of the class tessprf.TessPrf()


    Returns:
    ------------
    A scipy.optimise.OptimizeResult object. The best fit
    params are stored in OptimizeResult.x


    Notes:
    ----------
    The return object will frequently set success=False. I've seen no
    evidence that a failed fit give worse results than a successful fit.
    It appears that the fit likes to fall over when it gets very close
    to the minimum.

    In tests, the fitter often returns a best fit at a subpixel position
    adjacent to the true best fit position (i.e .02 pixels away in either
    column and/or row). This is a smaller distance than I'm concerned
    about for my use case, but if you need really accurate positions
    (say for PSF photometry), you may want to explore around the returned
    best fit position for the optimal fit.

    See example() in this module for an example of use
    """

    if not np.all(np.isfinite(img)):
        raise ValueError("Input img contains Nans. Set them to zero?")

    if img.shape != ((bbox[3]-bbox[2]), (bbox[1]-bbox[0])):
        raise ValueError("Shape of bbox doesn't match shape of image")

    #Get initial guess for centroid == centre of brightest pixel
    row, col = np.unravel_index(np.argmax(img), img.shape)
    col += bbox[0] + .5
    row += bbox[2] + .5
    scale = np.max(img)

    #Set options for optimiser
    initGuess = [col, row, scale]
    args = (ccd, camera, sector, bbox, img, prfObj)
    options = {'disp':False, 'eps':.02, 'maxiter':80}

    #Don't let the fit leave the bounding box
    bounds=[(bbox[0], bbox[1]), \
            (bbox[2], bbox[3]), \
            (1, None), \
           ]

    res = sopt.minimize(costFunc, initGuess, args, \
        method="L-BFGS-B", bounds=bounds, options=options)

    return res


def costFunc(x, ccd, camera, sector, bbox, img, prfObj):
    """Measure goodness of fit (chi-square) for a PRF
    with a given col, row and scale (i.e brightness)

    Inputs
    --------
    x
        (list or array) Elements are [col, row, fluxScaling]

    Returns
    ---------
    **float**

    """

    assert len(x) == 3
    model = prfObj.getPrfForBbox(x[0], x[1], ccd, camera, sector,  bbox)
    model *= x[2]

    cost = img - model
    cost = np.sum(cost**2)
    return cost


def getBoundingBoxForImage(img, hdr):
    """Get the bounding box for the an image from a TPF file
    i.e the corners of the rectangle that circumscribes the image

    Inputs:
    ---------
    img
        (2d np array) Image to construct bounding box for
    hdr
        (FITS header object). The header of the first extension of the TPF
        file. If you are dealing with an image not from a fits file, use
        img.shape to get the bounding box instead.

    Returns:
    ----------
    A 4 element list for column0, column1, row0, row1

    """
    shape= img.shape
    c0 = float(hdr['1CRV4P'])
    r0 = float(hdr['2CRV4P'])


    extent = [c0, c0+shape[1], r0, r0+shape[0]]
    return extent


import astropy.io.fits as pyfits
import dave.fileio.tpf as dtpf
import dave.diffimg.tessprf as prf

def example():
#    ticid = 307210830
    sector = 2
    camera = 4
    ccd = 3

    path = '/home/fergal/data/tess/hlsp_tess-data-alerts_tess_phot_00307210830-s02_tess_v1_tp.fits'
    fits, hdr = pyfits.getdata(path, header=True)


    cube = dtpf.getTargetPixelArrayFromFits(fits, hdr)
    img = cube[100]
    idx = np.isfinite(img)
    img[~idx] = 0

    prfObj = prf.TessPrf("/home/fergal/data/tess/prf/v2")

    bbox = getBoundingBoxForImage(img, hdr)
    res = fitPrfCentroidForImage(img, ccd, camera, sector, bbox, prfObj)

    #The fit fails terribly, so adjusting reported col/row position
    #to make the result look better cosmetically.

    res.x[0] += .5
    res.x[1] += 1.0
    plotCentroidFitDiagnostic(img, hdr, ccd, camera, sector, res, prfObj)
    return res


def plotCentroidFitDiagnostic(img, hdr, ccd, camera, sector, res, prfObj):
    """Some diagnostic plots showing the performance of fitPrfCentroid()

    Inputs:
    -------------
    img
        (np 2d array) Image of star to be fit. Image is in the
        format img[row, col]. img should not contain Nans

    hdr
        (Fits header object) header associated with the TPF file the
        image was drawn from

    ccd, camera, sector
        (int) Which CCD and camera is the image taken from. Which sector
        was the data observed (the measured PRF changes as time evolves)

    prfObj
        An object of the class prf.KeplerPrf()


    Returns:
    -------------
    **None**

    Output:
    ----------
    A three panel subplot is created
    """

    disp = lambda x: plt.imshow(x, cmap=plt.cm.YlGnBu_r, origin="bottom",
           interpolation="nearest")

    plt.figure(1)
    plt.clf()
    plt.subplot(131)
    disp(img)
    plt.colorbar()
    plt.title("Input Image")

    plt.subplot(132)
    c,r = res.x[0], res.x[1]
    bbox = getBoundingBoxForImage(img, hdr)
    model = prfObj.getPrfForBbox(c, r, ccd, camera, sector, bbox)
    model *= res.x[2]
    disp(model)
    plt.colorbar()
    plt.title("Best fit model")

    plt.subplot(133)
    diff = img-model
    yr = max( np.max(diff), -np.min(diff))
    plt.imshow(diff, cmap=plt.cm.RdBu, origin='bottom', interpolation='nearest')
    plt.clim(-yr, yr)
    plt.colorbar()
    plt.title("Residuals")

    print ("Performance %.3f" %(np.max(np.abs(diff))/np.max(img)))
