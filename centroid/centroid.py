# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 15:03:35 2015

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib as mp
import numpy as np


import matplotlib.pyplot as mp
import numpy as np

import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.centroid.prf as prf
import plotTpf

import scipy.optimize as sopt


def example():
    kepid = 8554498
    quarter = 16

    ar = mastio.KeplerArchive()
    fits, hdr = ar.getLongTpf(kepid, quarter, header=True)
    hdr0 = ar.getLongTpf(kepid, quarter, ext=0)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    module = hdr0['MODULE']
    output = hdr0['OUTPUT']
    img = cube[100]
    idx = np.isfinite(img)
    img[~idx] = 0

    prfObj = prf.KeplerPrf("/home/fergal/data/keplerprf")

    bbox = getBoundingBoxForImage(img, hdr)
    res = fitPrfCentroidForImage(img, module, output, bbox, prfObj)


    plotCentroidFitDiagnostic(img, hdr, module, output, res, prfObj)


    return res


def fitPrfCentroidForImage(img, ccdMod, ccdOut, bbox, prfObj):
    """Fit a PRF model to a TPF image and measure star centroid

    Given a cadence image from a TPF file, find the best fitting
    PRF model for that image. The best fit column and row of the
    PRF model can be treated as the photometric centroid of the
    light distribution from a single isolated star.

    Inputs:
    ---------
    img
        (np 2d array) Image of star to be fit. Image is in the
        format img[row, col]. img should not contain Nans

    ccdMod, ccdOut
        (int) CCD module and output of image. Needed to
        create the correct PRF model

    bbox
        [c1, c2, r1, r2]. Define the range of columns (c1..c2)
        and rows (r1..r2)  defined by the image.
        An exception raised if the following equality not true
        img.shape = (c2-c1), (r2-r1)

    prfObj
        An object of the class prf.KeplerPrf()


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
    args = (ccdMod, ccdOut, bbox, img, prfObj)
    options = {'disp':False, 'eps':.02, 'maxiter':80}

    #Don't let the fit leave the bounding box
    bounds=[(bbox[0], bbox[1]), \
            (bbox[2], bbox[3]), \
            (1, None), \
           ]

    res = sopt.minimize(costFunc, initGuess, args, \
        method="L-BFGS-B", bounds=bounds, options=options)

    return res


def getBoundingBoxForImage(img, hdr):
    shape= img.shape
    c0 = float(hdr['1CRV4P'])
    r0 = float(hdr['2CRV4P'])


    extent = [c0, c0+shape[1], r0, r0+shape[0]]
    return extent




def plotCentroidFitDiagnostic(img, hdr, ccdMod, ccdOut, res, prfObj):
    """Some diagnostic plots showing the performance of fitPrfCentroid()

    Inputs:
    -------------
    img
        (np 2d array) Image of star to be fit. Image is in the
        format img[row, col]. img should not contain Nans

    hdr
        (Fits header object) header associated with the TPF file the
        image was drawn from

    ccdMod, ccdOut
        (int) CCD module and output of image. Needed to
        create the correct PRF model

    prfObj
        An object of the class prf.KeplerPrf()


    Returns:
    -------------
    **None**

    Output:
    ----------
    A three panel subplot is created
    """
    mp.figure(1)
    mp.clf()
    mp.subplot(131)
    plotTpf.plotCadence(img, hdr)
    mp.colorbar()
    mp.title("Input Image")

    mp.subplot(132)
    c,r = res.x[0], res.x[1]
    bbox = getBoundingBoxForImage(img, hdr)
    model = prfObj.getPrfForBbox(ccdMod, ccdOut, c, r, bbox)
    model *= res.x[2]
    plotTpf.plotCadence(model, hdr)
    mp.colorbar()
    mp.title("Best fit model")

    mp.subplot(133)
    diff = img-model
    plotTpf.plotCadence(diff, hdr)
    mp.colorbar()
    mp.title("Residuals")

    print "Performance %.3f" %(np.max(np.abs(diff))/np.max(img))



def costFunc(x, module, output, bbox, img, prfObj):
    """Measure goodness of fit (chi-square) for a PRF
    with a given col, row and scale (i.e brightness)"""
    model = prfObj.getPrfForBbox(module, output, x[0], x[1], bbox)
    model *= x[2]

    cost = img-model
    cost = np.sum(cost**2)
#    print "%.3e" %(cost)
#    cost = np.log10(cost)
    return cost


def costFunc1(x, module, output, col, row, bbox, img, prfObj):
    """Debugging function.

    Does the same as costFunc, but col and row are constants,
    and only the brightness of the prf can be changed.
    """

    model = prfObj.getPrfForBbox(module, output, col, row, bbox)
    model *= x[0]

    cost = img-model
    cost = np.sum(cost**2)
    return cost
