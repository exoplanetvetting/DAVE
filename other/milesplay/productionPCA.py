
"""
Created on Tue Aug  2 09:55:43 2016

@author: Miles
"""

import numpy as np
import matplotlib.pyplot as plt
import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.fileio.kplrfits as kplrfits
import dave.misc.noise as noise
import gapfill

def getData(epic, campaign):
    """Obtains the data for the star

    For a particular star, this obtains its pixel time series data cube as well
    as time points which may be bad due to thruster firings.


    Inputs:
    ----------
    epic
        (int) The k2 epic number for the desired star
    campaign
        (int) The k2 campaign for the desired star

    Returns:
    ------------
    cube
        (3d numpy array) The k2 pixel time series data for the star with 
        dimensions (time, number of rows per cadence, number of columns per 
        cadence)
    badIdx
        (1d boolean numpy array) Boolean array with true corresponding to a bad
        time point due to thruster firings
    """
    
    ar = mastio.K2Archive()
    fits, hdr = ar.getLongTpf(epic, campaign, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    
    # get the thruster firing indeces
    q = fits['QUALITY'].astype(np.uint32)
    rollTweakIdx = q & kplrfits.SapQuality['DefiniteRollTweak']    
    momDumpIdx = q & kplrfits.SapQuality['MomDump']
    badIdx = rollTweakIdx | momDumpIdx
    badIdx = badIdx.astype(bool)
    return cube, badIdx

def takeOutNaNs(pixSeries):
    """Takes out NaNs from a pixel time series

    Takes out NaNs in each pixel's time series and replaces them by filling in 
    the gaps using interpolation methods.
    
    If a single pixel's time series is made up of all NaNs, replace all of the
    NaNs with 0

    Notes:
    -----------
    For more information, see gapfill.py


    Inputs:
    ----------
    pixSeries
        (2d numpy array) Size: (number of pixels in cadence, number of time 
        points)
        The 3d data cube reshaped to a 2d array to make it easy to iterate over

    Returns:
    ------------
    newPixSeries
        (2d numpy array) Size: (number of pixels in cadence, number of time 
        points)
        The input array with all of the NaNs taken out
    """
    newPixSeries = []
    for singlePixTimeSeries in pixSeries:
        if np.all(np.isnan(singlePixTimeSeries)):
            singlePixTimeSeries = np.zeros(len(singlePixTimeSeries))
        else:
            singlePixTimeSeries, badIdx = gapfill.gapFill(singlePixTimeSeries)
        newPixSeries.append(singlePixTimeSeries)
                        
    return np.array(newPixSeries)
        
    
def hasNaNs(pixSeries):
    """Checks if any NaNs are present in the pixel time series
    
    Iterates through a list of single pixel time series and takes the sum of 
    the series. If the sum is NaN, we know there is a NaN present in the series
    
    Inputs:
    ----------
    pixSeries
        (2d numpy array) Size: (number of pixels in cadence, number of time 
        points) The 3d data cube reshaped to a 2d array to make it easy to 
        iterate over.
        
    Returns:
    ----------
    A boolean. True if there are NaNs, false if there are no NaNs
    """
    NaNcheckArray = []
    for singlePixSeries in pixSeries:
        if np.isnan(np.sum(pixSeries)):
            print np.sum(pixSeries)
            NaNcheckArray.append(True)
        else:
            NaNcheckArray.append(False)
    return np.any(NaNcheckArray)
    
def normalize(pixSeries):
    """Normalizes the pixel time series by the standard deviation method

    Normalizes by taking the pixel time series, subtracting the mean, and
    dividing by the standard deviation of the series

    Inputs:
    ----------
    pixSeries
        (2d numpy array) Size: (number of pixels in cadence, number of time 
        points) The 3d data cube reshaped to a 2d array to make it easy to 
        iterate over.
        Note: This series should have the NaNs taken out of it first


    Returns:
    ------------
    pixSeriesNorm 
        (2d numpy array) Size: (number of pixels in cadence, number of time 
        points)
        The normalized input array
    """
    assert not hasNaNs(pixSeries)
    pixSeriesNorm = [] 
    for singlePixSeries in pixSeries:
       
        singlePixSeries = (singlePixSeries - np.mean(singlePixSeries))/ \
                               (np.std(singlePixSeries) + 1E-10)
        pixSeriesNorm.append(singlePixSeries)       
        
    return np.array(pixSeriesNorm)

def getRawLightcurve(pixSeries):
    """Gets the total light curve
    
    Sums up the flux for each cadence to get a total light curve
    
    Inputs:
    ----------
    pixSeries
        (2d numpy array) Size: (number of pixels in cadence, number of time 
        points) The 3d data cube reshaped to a 2d array to make it easy to 
        iterate over.
        Note: This series should have the NaNs taken out of it first
        
    Returns:
    ----------
    totLightcurve
        (1d numpy array) Length equal to number of time steps. The sum of each
        cadence at each time step
    
    """
    totLightcurve = np.sum(pixSeries, axis=0)
    return totLightcurve

def performSVD(pixSeriesNorm):
    """Obtains the principal components of the pixel time series

    Uses the numpy.linalg.svd function to compute the principal components of
    the single pixel time series which are non-zero. The svd function returns
    a tuple of values, the second index being a 2d numpy array of the principal
    components with size (number of components, time steps).

    Notes:
    -----------
    More information on the numpy svd function can be found 
    `here <http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html>`_ 

    Inputs:
    ----------
    pixSeriesNorm
        (2d numpy array) Size: (number of pixels in cadence, number of time 
        points)
        The 3d data cube reshaped to a 2d array to make it easy to iterate over
        Note: The time series should be normalized first using the normalize()
        function

    Returns:
    ------------
    prinComps
        (2d numpy array) Size: (number of components, time steps)
        Numpy array of the principal components, ranked from most influential 
        to least (i.e. prinComps[0] is the most dominant component in the light
        curve)
    """
    # don't include the single pixel time series which are all zero
    idx = np.any(pixSeriesNorm, axis=1)
    pixSeriesNormNonzero = pixSeriesNorm[idx, :]

    prinComps = np.linalg.svd(pixSeriesNormNonzero, full_matrices=0,compute_uv = 1)[2]
    return prinComps

def curveFit(prinCompsToFit, totLightcurve):
    """Fits the desired number of principal components to the raw light curve

    Uses the numpy.linalg.lstsq function to get a vector of coefficients which
    is then used to generate the best fit for the desired number of principal
    components by multilying the matrix of principal components by the vector
    of coefficients. The resulting vector is the best fit light curve.

    Notes:
    -----------
    More information on the numpy lstsq function can be found 
    `here <http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html>`_ 



    Inputs:
    ----------
    prinCompsToFit
        (2d numpy array) Size: (number of principal components to fit, time 
        steps)
        This is a 2d matrix of the principal components to use in generating
        the best fit to the total raw light curve.
    totLightcurve
        (1d numpy array) Array of summed raw light curve generated by summing
        each cadence, the length of which should be equal to the number of time
        steps for this particular data set.


    Returns:
    ------------
    fittedCurve
        (1d numpy array) This array is the best fit to the raw light curve 
        using the specified principal components, the length of which is equal 
        to the number of time steps for this particular data set.
    """
    coeffs = np.linalg.lstsq(prinCompsToFit.T, totLightcurve)[0]
    fittedCurve = prinCompsToFit.T.dot(coeffs.T)
    return fittedCurve

def chooseBestSGcdppValue(sgcdppList, thresh=0.5):
    """Chooses the number of principal components to use in the final best fit

    This function expects a list of values generated by using Jeff Van Cleve's 
    Savitzy-Golay technique. Each of these values correspond to using a 
    different number of principal components in the fit and correction of the 
    light curve. This function then chooses the optimal number of principal 
    components by determining at what value in this list the values become 
    close enough to one another to effectively treat it as a flat line. At this
    point, it is futile and even detrimental to continue adding principal 
    components to the fit calculation.
    
    This is achieved by implementing a threshold, i.e. if the difference 
    between two values becomes less than a certain threshold, the number of 
    principal components used in the fit corresponding to that sgcdpp value is 
    the optimal number of principal components to use.

    Notes:
    -----------
    More information on the Savitzy-Golay technique can be found in 
    dave.misc.noise



    Inputs:
    ----------
    sgcdppList
        (1d array) A list of values obtained by the Savitzy-Golay technique. 
        Each value corresponds to the calculation done with a fit using a 
        different number of principal components.


    Optional Inputs:
    -------------------
    thresh
        (float) Threshold to use in determining where in the list of sgcdpp 
        values it starts to become effectively flat

    Returns:
    ------------
    optimalNumPC
        (int) The optimal number of principal components to use in the fit and 
        correction to the raw light curve
    """
    
    # pick where it starts to flatline based on some threshold
    flatLineInd = np.where(np.abs(np.diff(sgcdppList[1])) < thresh)[0][0]

    optimalNumPC = sgcdppList[0][flatLineInd]
    
    return optimalNumPC
    
def generatePrinComps(dataCube, thrusterFireInds):
    """ Prepares the data for SVD calculation and gets principal components
    
    Prepares the data for SVD calculation by reshaping the data cube into a 2d
    numpy array, taking out thruster firing indeces, taking out NaNs, 
    generating the total raw light curve, and normalizing the pixel time series
    
    Ultimately gets the principal components by doing an SVD calculation
    
    Inputs:
    ----------
    dataCube
        (3d numpy array) Size: (timesteps, num rows, num cols) the raw data 
        cube obtained from Kepler
    thrusterFireInds
        (1d boolean numpy array) True where there are known bad indeces due to
        thruster firings
        
    Returns:
    ----------
    prinComps
        (2d numpy array) Size: (number of components, time steps)
        Numpy array of the principal components, ranked from most influential 
        to least (i.e. prinComps[0] is the most dominant component in the light
        curve)
    totLightcurve
        (1d numpy array) Array of summed raw light curve generated by summing
        each cadence, the length of which should be equal to the number of time
        steps for this particular data set.
    t
        (1d numpy array) the time steps
    numPixels
        (int) the number of pixels per cadence
    
    """
    
    nt, nr, nc = dataCube.shape
    
    # flatten the cube 
    pixSeries = dataCube.reshape((nt, nc*nr)).T
    
    # take out the points with known thruster firings 
    pixSeries = pixSeries[:, ~thrusterFireInds]
    
    # replace NaNs with mean of each pixel series
    pixSeries = takeOutNaNs(pixSeries)

    # get the total number of pixels for each cadence
    numPixels = pixSeries.shape[0]

    # get the total light curve by summing all of the cadences
    totLightcurve = getRawLightcurve(pixSeries)
    
    # get the time axis for light curves
    t = np.arange(len(totLightcurve))

    # normalize the time series for each pixel
    pixSeriesNorm = normalize(pixSeries)

    # get the principal components
    prinComps = performSVD(pixSeriesNorm)
    
    return prinComps, totLightcurve, t, numPixels
    
def inputData():
    """Place to specify what star to get data for
    
    Gets the pixel time series data cube and indeces of known thruster firings
    for the star
    
    Inputs:
    ----------
    None
        
    Returns:
    ----------
    dataCube
        (3d numpy array) 
    """
    epic = 206103150
    campaign = 3 
    dataCube, thrusterFireInds = getData(epic, campaign)
    return dataCube, thrusterFireInds
    
def chooseNumPrinComps(prinComps, totLightcurve, numPixels):
    """Chooses the optimal # of principal components to use in the final fit
    
    Generates a list of sgcdpp values by doing an sgcdpp calculation (see
    noise.computeSgCdpp_ppm)
    
    Iterates through the list of sgcdpp values and determines where the slope 
    between two values starts to "flatten out" based on a certain threshold
    
    Inputs:
    ----------
    prinComps
        (2d numpy array) Size: (number of components, time steps)
        Numpy array of the principal components, ranked from most influential 
        to least (i.e. prinComps[0] is the most dominant component in the light
        curve)
    totLightcurve
        (1d numpy array) Array of summed raw light curve generated by summing
        each cadence, the length of which should be equal to the number of time
        steps for this particular data set.

    numPixels
        (int) the number of pixels per cadence
        
    Returns:
    ----------
    optimalNumPC
        (int) the optimal number of principal components to use in the best fit
        and correction to the raw light curve
    
    """
    # in my experience, the optimal number of principal components is always 
    # less than half of the total number of pixels. To save computation time, 
    # I only generate sgcdpp values for the first half of principal components    
    numPrinComps_sgcdppCalc = numPixels/2
    
    
    # threshold for sigma clip
    sigmaClipThresh = 5.    
    
    # get a list of sgcdpp values for fitting different numbers of prin comps
    n = 1
    sgcdppList = []
    while n < numPrinComps_sgcdppCalc:
        
        prinCompsToFit = prinComps[:n]
        
        # fit the light curve to the amount of prin comps selected
        fittedCurve = curveFit(prinCompsToFit, totLightcurve)

        # correct the light curve by subtracting the prin comp fit
        correctedCurve = totLightcurve - fittedCurve
        
        # make the mean of the corrected curve = 0, 
        # necessary for sgcdpp calculation
        correctedCurveMeanZero = correctedCurve / np.mean(correctedCurve) - 1
        
        # get sigma clip true/false values on the corrected curve with mean 0
        sigClip_tf = noise.sigmaClip(correctedCurveMeanZero, sigmaClipThresh)
        
        # perform sigma clip
        correctedMeanZeroSigClip = correctedCurveMeanZero[~sigClip_tf]
        
        # get the sgcdpp value
        sgcdppValue = noise.computeSgCdpp_ppm(correctedMeanZeroSigClip)
        
        sgcdppList.append([n,sgcdppValue])
        
        n += 1
    
    sgcdppList = np.array(sgcdppList).T
  
    # choose the number of principal components to use
    optimalNumPC = chooseBestSGcdppValue(sgcdppList)
    
    return optimalNumPC
    
def main():
    dataCube, thrusterFireInds = inputData()
    
    prinComps, totLightcurve, t, numPixels = generatePrinComps(dataCube, 
                                                              thrusterFireInds)
    
    optimalNumPC = chooseNumPrinComps(prinComps, totLightcurve, numPixels)
    
    bestFit = curveFit(prinComps[:optimalNumPC], totLightcurve)
    
    bestFitCorrectedLightcurve = totLightcurve - bestFit
    
    # plot the raw light curve with best fit in green
    plt.figure()
    plt.plot(t, bestFitCorrectedLightcurve, ".", label="Corrected Lightcurve")
    plt.title("sgcdpp best correction, PC = %i"%optimalNumPC)
    plt.xlabel("time")
    plt.ylabel("flux")
    #plt.legend()
    plt.show()
     

main()