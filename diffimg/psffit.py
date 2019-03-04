"""
Created on Wed Nov  7 14:26:29 2018

Fit an analytic model of a point spread function (PSF) to an image.

Nomenclature
-------------------------
A point spread function is the analytic model of the distribution of light on the focal planet
from a point source.

A Pixel Response Function 9PRF) is the distribution of light on the pixels of the detector, with one
value for each pixel. In general, a PRF it may include the jitter, or the intrapixel response,
but those aren't included in this model yet.


Usage
---------
Call `fitPrf()` with an image and a PSF model
Example PSF functions are `gaussianPsf()` and `gaussianWithConstantSkyPsf()`


Notes
----------
* The signature of a model is function(col, row, *args)

* The fitting function is extremely slow. To speed it up, we use numba to on-the-fly compile
  the PSF model function to C-code. This speeds up the fitting by roughly a factor of 7.
  To have numba precompile the model function you need to decorate the function appropriately.
  The function `jit_psf5args()` decorates `gaussianWithConstantSkyPsf()`. You need to write
  a separate decorator for each model that contains a different number of arguments.
  (At least until I can figure out a better way to decorate).Each   argument must be a
  double.  If the model can't be described in this fashion it can't be compiled.
  the

* The code uses the following coordinate conventions

Image coordinates


  row
  ^
  |
  |
  |
  |___________
     col -->

The bottom left of a pixel is the origin


  (0,1)
  +-----+ (1,1)
  |     |
  |     |
  |     |
  +-----+ (1,0)

All functions take (col,row), never (row, col), unless copiously documented

"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import numpy as np

import scipy.integrate as spInt
import scipy.optimize as spOpt
from numba import  njit, cfunc
from numba.types import CPointer, float64, intc
from scipy import LowLevelCallable

#https://stackoverflow.com/questions/49683653/how-to-pass-additional-parameters-to-numba-cfunc-passed-as-lowlevelcallable-to-s
def jit_psf5args(func):
    """A complicated piece of code used by numba to compile the PSF.

    This one works for `gaussianWithConstantSkyPsf()`
    """
    jitted_function = njit(func)

    @cfunc(float64(intc, CPointer(float64)))
    def wrapped(n, xx):
        return jitted_function(xx[0], xx[1], xx[2], xx[3], xx[4], xx[5], xx[6])

    return LowLevelCallable(wrapped.ctypes)


@jit_psf5args
def gaussianWithConstantSkyPsf(col, row, col0, row0, sigma, flux0, sky):
    """A model PSF of a 2d symettric gaussian with a constant background

    Inputs
    ---------
    col, row
        (floats) Which location to compute PSF for
    col0, row0
        (floats) Centre of psf
    sigma
        (float) Width of the Gaussian
    flux0
        (float) Height of the gaussian (sort of)
    sky
        (float) Background level

    Returns
    --------
    (float)
    """
    assert sigma > 0

    z_col = .5 * (col - col0) / sigma
    z_row = .5 * (row - row0) / sigma

    return sky + flux0 * np.exp(- z_col**2) * np.exp( - z_row**2)


def fitPrf(img, prfFunc, guess):
    """Fit a PRF to an image

    Inputs
    --------
    img
        (2d numpy array) Image to fit
    prfFunc
        (function) Model to fit. See module level documentation for more details.
    guess
        (tuple or array) Arguments prfFunc

    Returns
    ------------
    A scipy.optiminze.ResultsObject. The .x attribute contains the best fit parameters

    Example
    ----------
    To fit a model with the signature (col, row, flux, width), the guess array would be of
    length 2.
    """

    nr, nc = img.shape
    soln = spOpt.minimize(costFunc, guess,args=(prfFunc, img), method='Nelder-Mead', bounds=None)
    return soln


def gaussianPsf(col, row, col0, row0, sigma, flux0):

    assert sigma > 0

    z_col = .5 * (col - col0) / sigma
    z_row = .5 * (row - row0) / sigma

    return flux0 * np.exp(- z_col**2) * np.exp( - z_row**2)



def costFunc(arglist, func, img, mask=None):
    """Compute difference between image and its model for given model params

    Inputs
    ----------
    arglist
        (tuple or array) Tunable parameters of model
    func
        (function) Model to fit
    img
        (2d np array) Image to fit


    Optional Inputs
    ----------------
    mask
        (2d np array) Zero elements of mask indicate bad data which should not be
        included in the fit


    Returns
    ----------
    float
    """

    nr, nc = img.shape
    model = computeModel(nc, nr, func, arglist)
    diff = img - model

    if mask is not None:
        assert np.all( mask.shape == img.shape)
        diff[~mask] = 0
        img[~mask] = 0  #In case bad values are set to Nan

    cost = np.sqrt( np.sum(diff**2) )
    return cost





#@profile
def computeModel(numCols, numRows, func, arglist):
    """Compute model flux for an image with size (numCols, numRows)

    Inputs
    -------
    numCols, numRows
        (ints) Shape of the image to compute the model PRF for
    func
        (function) Model PRF
    arglist
        (tuple or array) Tunable parameters of the model

    Returns
    ----------
    A 2d numpy array representing the model PRF image.
    """

    model = np.zeros( (numRows, numCols) )

    for i in range(numCols):
        def gfun(x):
            return i

        def hfun(x):
            return i+1

        for j in range(numRows):
            val = spInt.dblquad(func, j, j+1, gfun, hfun, args=arglist)[0]
            model[j, i] = val  #Numpy flips row and column

    return model
