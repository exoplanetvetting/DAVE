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
import dave.diffimg.diffimg as diffimg
import dave.diffimg.arclen as arclen
import plotTpf

import scipy.optimize as sopt
np.ones

def exampleDiffImgCentroiding():
    k2id =  206103150
    campaign = 3

    ar = mastio.K2Archive()
    fits, hdr = ar.getLongTpf(k2id, campaign, header=True)
    hdr0 = ar.getLongTpf(k2id, campaign, ext=0)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    idx = np.isfinite(cube)
    cube[~idx] = 0  #Remove Nans

    flags = fits['QUALITY']
    ccdMod = hdr0['module']
    ccdOut = hdr0['output']

    #Compute roll phase
    llc = ar.getLongCadence(k2id, campaign)
    time= llc['TIME']
    cent1 = llc['MOM_CENTR1']
    cent2 = llc['MOM_CENTR2']
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = arclen.computeArcLength(centColRow, flags>0)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    prfObj = prf.KeplerPrf("/home/fergal/data/keplerprf")
    bbox = getBoundingBoxForImage(cube[0], hdr)

    period =  	4.1591409
    epoch = fits['time'][491]
    dur = 3.0

    out = measureDiffOffset(period, epoch, dur, time, prfObj, \
        ccdMod, ccdOut, cube, bbox, rollPhase, flags)

    idx = out[:,1] > 0
    mp.clf()
    mp.plot(out[:,3]-out[:,1], out[:,4]- out[:,2], 'ro')
    return out


def exampleFitting():
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


def measureDiffOffset(period_days, epoch_bkjd, duration_hrs, \
    time, prfObj, ccdMod, ccdOut, cube, bbox, rollPhase, flags):
    """Measure Centroid shift between intransit and difference image
    for every in-transit cadence

    Inputs:
    -----------
    period_days, epoch_bkjd, duration_hrs
        (floats) Properties of transit

    time_bkjd
        Array of times per cadence for the given campaign

    prfObj
        An object of the class prf.KeplerPrf()

    ccdMod, ccdOut
        (int) CCD module and output of image. Needed to
        create the correct PRF model

    cube
        (3d np array) A data cube created from a TPF file.
        See fileio.tpf.getTargetPixelArrayFromFits()

    bbox
        [c1, c2, r1, r2]. Define the range of columns (c1..c2)
        and rows (r1..r2)  defined by the image.
        An exception raised if the following equality not true
        img.shape = (c2-c1), (r2-r1)

    rollPhase
        (1d np array) An array of roll phases for each row
        of cube. len(rollPhase) == len(cube). Units of this
        array don't matter, so long as cadences with similar
        roll angles have similar values of rollPhase. Roll phases
        for bad cadences should be set to a bad value

    flags
        (1d array) flag values indicating bad cadences.
        Currently a non-zero value of flags indicates a bad
        cadence.

    Returns:
    -------------
    A array with 5 columns, and as many rows as there are
    in transit cadences. The columns are

    0: Relative cadence number
    1: In transit centroid column
    2: In transit centroid row
    3: Diff img centroid column
    4: Diff img centroid row

    If there is a statisically significant difference between the intransit
    and difference image centroids then the transit is most likely not
    on the target.
    """
    idx = getIndicesInTransit(period_days, epoch_bkjd, duration_hrs, time)
    wh = np.where(idx)[0]
    out = -1 * np.ones((len(wh), 5))
    for i,w in enumerate(wh):
        out[i,0] = w
        try:
            out[i, 1:] = measureInTransitAndDiffCentroidForOneImg(\
                prfObj, ccdMod, ccdOut, cube, w, bbox, rollPhase, flags, \
                hdr=None, plot=False)
        except ValueError:
            pass
        print i, len(wh)

    return out


def getIndicesInTransit(period_days, epoch_bkjd, duration_hrs, time_bkjd):
    """
    Find the cadences affected by a transit.

    Inputs::
    --------
    period_days, epoch_bkjd, duration_hrs
        (floats) Properties of transit

    time_bkjd
        Array of times per cadence for the given campaign


    Returns:
    ------------
    An array of booleans of length equal to length of time_bkjd.
    Cadences in transit are set to true, all other cadences to false
    """

    time = time_bkjd #Mneumonic
    n1 = int(np.floor((time[0] - epoch_bkjd)/period_days))
    n2 = int(np.ceil((time[-1] - epoch_bkjd)/period_days))
    dur_days = duration_hrs/24.

    inTransit = np.zeros_like(time, dtype=bool)
    for n in range(n1, n2+1):
        t0 = epoch_bkjd + n*period_days - .5*dur_days
        t1 = epoch_bkjd + n*period_days + .5*dur_days

        idx = (time >= t0) & (time <= t1)
#        import pdb; pdb.set_trace()
        inTransit |= idx

    return inTransit


def measureInTransitAndDiffCentroidForOneImg(prfObj, ccdMod, ccdOut, cube, rin, bbox, rollPhase, flags, hdr=None, plot=False):
    """Measure image centroid of in-transit and difference images

    Inputs:
    -----------
    prfObj
        An object of the class prf.KeplerPrf()


    ccdMod, ccdOut
        (int) CCD module and output of image. Needed to
        create the correct PRF model


    cube
        (3d np array) A TPF data cube as returned by
        dave.fileio.getTargetPixelArrayFromFits()

    rin
        (int) Which image to process. rin should be in the range 0..len(cube)

    bbox
        [c1, c2, r1, r2]. Define the range of columns (c1..c2)
        and rows (r1..r2)  defined by the image.
        An exception raised if the following equality not true
        img.shape = (c2-c1), (r2-r1)

    rollPhase
        (1d np array) An array of roll phases for each row
        of cube. len(rollPhase) == len(cube). Units of this
        array don't matter, so long as cadences with similar
        roll angles have similar values of rollPhase

    flags
        (1d array) flag values indicating bad cadences.
        Currently a non-zero value of flags indicates a bad
        cadence.

    Optional Inputs:
    ---------------
    hdr
        Fits header object for TPF file. Useful if you want to plot

    plot
        (bool) Request plots.


    Returns:
    -------------
    A 4 element numpy array
    ic  In transit centroid column
    ir  In transit centroid row
    dc  Difference image centroid column
    dr  Difference image centroid row
    """
    inTrans = cube[rin]
    diff, oot= diffimg.constructK2DifferenceImage(cube, rin, rollPhase, flags)

    ootRes = fitPrfCentroidForImage(oot, ccdMod, ccdOut, bbox, prfObj)
    diffRes = fitPrfCentroidForImage(diff, ccdMod, ccdOut, bbox, prfObj)

    if plot:
        mp.subplot(121)
        plotTpf.plotCadence(inTrans, hdr)
        mp.colorbar()
        mp.subplot(122)
        plotTpf.plotCadence(diff, hdr)
        mp.colorbar()

        mp.plot(ootRes.x[0], ootRes.x[1], 'ro', ms=12, alpha=.4)
        mp.plot(diffRes.x[0], diffRes.x[1], 'r^', ms=12, alpha=.4)

    return np.array([ootRes.x[0], ootRes.x[1], diffRes.x[0], diffRes.x[1]])





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
