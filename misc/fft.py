# -*- coding: utf-8 -*-
"""
A wrapper around numpy's fft to compute the frequencies associated
with the FFT. This is hard enough that I get it wrong if I try to
rederive it each time."""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

def computeFft(y, expTime_sec, overres=3):
    """Compute the overresolved FFT of input data and scale the amplitude

    Inputs:
    y           (1d np array)   Input data. Treated as equally spaced data
    expTime_Sec (float) Spacing between input points

    Optional Inputs:
    overres     (int) How many times to overresolve the FT. Values below
                3 cause FTs that are difficult to interpret, above 10
                yields no better information.


    Returns:
    A 2d numpy array. The zeroth column is frequency in microHz. The
    second column is fractional amplitude at that frequency

    Notes:
    Calls numpy's rfft function, which is based on fftw
    """

    ySize = len(y)
    pow2YSize = numPointsInFt(ySize, overres)
    ftSize = pow2YSize/2. + 1

    ft = np.abs(np.fft.rfft(y, pow2YSize))

    ft /= .5*ySize

    #Package as 2d array
    out = np.empty( (ftSize, 2) )
    out[:,1] = ft
    out[:,0] = computeFrequencies_Hz(ySize, overres, expTime_sec)

    return out


def numPointsInFt(nPoints, overres):
    """Compute number of points needed to input to an fft algorithm to ensure
    a) The data is overesolved by at least a factor of overres
    b) The number of poitns is an integer power of 2
    """

    nPts = overres*nPoints+1
    log = np.ceil( np.log2(nPts))
    return int(2**log)



def computeFrequencies_Hz(numPoints, overres, expTime_sec):
    """Compute the frequencies corresponding to the output of fft.

    Inputs:
    numPoints:  (int) Number of points in input data.
    overes      (int) How much the data was overresolved.
    expTime_sec (float) Spacing between points in seconds.

    Returns:
    1d numpy array

    Notes
    numPoints is the number of points in the input data before the input
    is padded to account for overresolution and ensuring the length is a power
    of two.
    """

    nPts = numPointsInFt(numPoints, overres)
    timespan = nPts*expTime_sec
    res = 1/float(timespan)  #In Hz

    out = np.arange( nPts/2.+1) * res
    return out


def computeFrequencies_uHz(numPoints, overres, expTime_sec):
    """Compute the frequencies corresponding to the output of fft.

    Inputs:
    numPoints:  (int) Number of points in input data.
    overes      (int) How much the data was overresolved.
    expTime_sec (float) Spacing between points in seconds.

    Returns:
    1d numpy array

    Notes
    numPoints is the number of points in the input data before the input
    is padded to account for overresolution and ensuring the length is a power
    of two.
    """

    return 1e6 * computeFrequencies_Hz(numPoints, overes, expTime_sec)
