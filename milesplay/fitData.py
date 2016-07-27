#==============================================================================
# 
# Author: Miles Currie
# Created: 06/21/2016
#
# Cleaned up a little
#
#
#==============================================================================



import numpy as np
import matplotlib.pyplot as plt
import dave.fileio.mastio as mastio
import dave.fileio.tpf as tpf
import dave.trapezoidFit.trapfit as tf
import sep
import dave.fileio.kplrfits as kplrfits
import dave.misc.noise as noise


def getData(epic, campaign):
    ar = mastio.K2Archive()
    fits, hdr = ar.getLongTpf(epic, campaign, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    
    # get the thruster firing indeces
    q = fits['QUALITY'].astype(np.uint32)
    idx1 = q & kplrfits.SapQuality['DefiniteRollTweak']    
    idx2 = q & kplrfits.SapQuality['MomDump']
    idx3 = idx1 + idx2
    
    return cube, idx3
    
def normalize(lcvMatrix):
    lcvMatrixNorm = [] 
    for lcv in lcvMatrix:
        if np.std(lcv) == 0.0:
            continue
        else:
            lcv = (lcv - np.mean(lcv))/np.std(lcv)
            lcvMatrixNorm.append(lcv)       
    return np.array(lcvMatrixNorm)
    
def getPrinComps(lcvMatrixNorm):
    return np.linalg.svd(lcvMatrixNorm, full_matrices=0,compute_uv = 1)
    
def make_plot(x, y, new=True, show=True, title=None, marker=".", label=None):
    if new == True:
        plt.figure()
    plt.title(title)
    plt.plot(x,y,marker,label=label)
    plt.legend()
    if show == True:
        plt.show()
         

def curveFit(pcaMatrix, cadenceSum, numPrinComps, plotPC=False):
    pcaMatrixNorm = normalize(pcaMatrix)
    prinComps = getPrinComps(pcaMatrixNorm)[2]
    
    n = 0
    x = []
    while n < numPrinComps:
        x.append(prinComps[n])        
        n += 1
    if plotPC == True:
        n = 1
        for comp in x:
            plt.figure()
            plt.plot(range(3386), comp)
            plt.title("PC %s"%str(n))
            plt.show()
            n += 1
        
    y = cadenceSum
    x = np.array(x).T
    lstsq_return = np.linalg.lstsq(x,y)
    #print lstsq_return
    coeffs = lstsq_return[0]
    #print coeffs
    residuals = lstsq_return[1]
    fittedCurve = x.dot(coeffs.T)
    #print fittedCurve
    return fittedCurve, residuals

def plotCadence(singleImage, hdr=None, axis='Fixed', *args, **kwargs):
    """
    Plot a single cadence from a TPF data cube

    Inputs:
    ---------
    singleImage
        (np 2d array) Image to plot

    Optional Inputs:
    ----------------
    hdr
        Header from fits file
    axis
        ('relative' of 'fixed'). If set to fixed, and hdr is set,
        the axis labels are in units of pixels. If relative,
        the bottom left of the image is labeled as pixel 0,0

    Additional arguments passed to mp.imshow
    """
    #Default interpolation is to plot square pixels
    interpolation = kwargs.pop('interpolation', 'nearest')
    cmap = kwargs.pop('cmap', plt.cm.YlGnBu_r)

    shape = singleImage.shape   #Gives the shape as (row,col)

    #What numbers to put on x and y axes?
    extent=None
    if axis.lower()=="relative":
        nr, nc = singleImage.shape
        extent=[0,nc, 0, nr]
    elif axis.lower() == "fixed":
        if hdr is None:
            raise ValueError("Set axis to 'relative' if header not available")
        #For kepler TPF files, all the columns in the binary
        #table contain images of the same shape, so it
        #doesn't matter which CRPIXs I take.
        c0 = float(hdr['1CRV4P'])
        r0 = float(hdr['2CRV4P'])

        extent = [c0, c0+shape[1], r0, r0+shape[0]]
    else:
        raise ValueError("Unrecognised axis style requested")
    plt.figure()
    plt.imshow(singleImage, origin='lower', interpolation=interpolation, \
        extent=extent, cmap=cmap, *args, **kwargs)

def noMoreNaN(lcvMatrix):
    newMatrix = []
    n = 0
    for lcv in lcvMatrix:
        n += 1
        if np.isnan(np.sum(lcv)):  # quickly checks for a nan in the list
            #print "Removing NaNs from row %s"%str(n-1)            
            where_are_nans = np.isnan(lcv) # find indeces of nans
            nanMean = np.nanmean(lcv) # gets mean of array ignoring nans
            if np.isnan(nanMean): # if the mean itself is nan, set nanMean to 0
                nanMean = 0
            lcv[where_are_nans] = nanMean # replace nans with nanMean
            newMatrix.append(lcv) 
        else:
            newMatrix.append(lcv)
           # print "Row %s does not contain NaNs"%str(n-1)
            continue
        
    return np.array(newMatrix)
            
            
            
def checkNaN(lcvMatrix):
    for lcv in lcvMatrix:
        if np.isnan(np.sum(lcv)):
            print "Found a NaN. Your code doesn't work."
            return False
        else:
            return True

def correctedCurve(rawData, sumCoeffPC):
    #corrected = rawData/sumCoeffPC - 1.
    corrected = rawData - sumCoeffPC
    return corrected

def createMask(data, thresh=100, title="Mask", plotMask=True):
    extraction = sep.extract(data, thresh, segmentation_map=True)
    mask = extraction[1]
    if plotMask:
        plt.figure()
        plt.imshow(mask, aspect='auto', interpolation='nearest', origin='lower')
        plt.title(title)
        plt.show()
    return mask

def reduceAperture(mask, data):
    newData = []
    mask = np.ndarray.flatten(mask)
    for cadence in data[:]:       
        newAperture = []
        cadence = np.ndarray.flatten(cadence)
        inds = np.where(mask == 1)
        for ind in inds[0]:
            newAperture.append(cadence[ind])
            
        newData.append(newAperture)
    lcvMatrix = np.array(newData).T
    return lcvMatrix
    
 
def trapFit(time_days, flux_frac, period_days, phase_bkjd, duration_hrs, depth_ppm):
    trapFit = tf.trapezoid_fit(time_days, flux_frac, np.ones(len(time_days)), \
                   period_days, phase_bkjd, duration_hrs, \
                   depth_ppm, fitTrialN=13, fitRegion=10.0, \
                   errorScale=1.0, debugLevel=3, \
                   sampleN=15)
     #print trapFit
    return trapFit
       
       
def doBeforeAnything(epic, campaign, plotMask):
    # get the data
    getData_return = getData(epic, campaign)
    data = getData_return[0]
    inds_to_eliminate = getData_return[1].astype(bool)
    
    # plot the image
    #plotCadence(data[0], axis='relative')    
    
    # create a mask for the light curve
    mask = createMask(data[0], thresh=100, title="Light Curve Mask", plotMask=plotMask)
    lcvMatrix = reduceAperture(mask, data)
    
    # create a mask for PCA
    pca_mask = createMask(data[0], thresh=100, title="PCA Mask", plotMask=plotMask)
    pcaMatrix = reduceAperture(pca_mask, data)

    # take out data points with thruster firings
    lcvMatrix = lcvMatrix[:,~inds_to_eliminate]
    pcaMatrix = pcaMatrix[:,~inds_to_eliminate]
    
    #get rid of nans in matrices
    lcvMatrix = noMoreNaN(lcvMatrix)
    assert checkNaN(lcvMatrix)
    pcaMatrix = noMoreNaN(pcaMatrix)
    assert checkNaN(pcaMatrix)
    
    # make the raw light curve
    cadenceSum = np.sum(lcvMatrix, axis=0)
    return cadenceSum, lcvMatrix, pcaMatrix, mask, pca_mask

def varylcvMask(epic, campaign, plotMask):
    # get the data
    getData_return = getData(epic, campaign)
    data = getData_return[0]
    inds_to_eliminate = getData_return[1].astype(bool)
    
    numPrinComps = 35
    
    # plot the image
    plotCadence(data[0], axis='relative')    

    thresh_list = np.arange(0, 300, 25)
        
    
    # create a mask for PCA
    pca_mask = createMask(data[0], thresh=100, title="PCA Mask", plotMask=plotMask)
    pcaMatrix = reduceAperture(pca_mask, data)
    

    # take out data points with thruster firings
    pcaMatrix = pcaMatrix[:,~inds_to_eliminate]
    
    #get rid of nans in matrices

    pcaMatrix = noMoreNaN(pcaMatrix)
    assert checkNaN(pcaMatrix)
    


    
    print "\n#pixels_lcv\tdiff_std"
    print "-----------------"
    corrected_list = []
    results_list = []
    for thresh in thresh_list:
        mask = createMask(data[0], thresh=thresh, title="Light Curve Mask, thresh=%s"%str(thresh), plotMask=plotMask)
        lcvMatrix = reduceAperture(mask, data)
        lcvMatrix = lcvMatrix[:,~inds_to_eliminate]
        lcvMatrix = noMoreNaN(lcvMatrix)
        assert checkNaN(lcvMatrix)
        cadenceSum = np.sum(lcvMatrix, axis=0)
        try:
            fittedCurve = curveFit(pcaMatrix, cadenceSum, numPrinComps)[0]
        except IndexError or AssertionError:
            print "Can't use specified amount of prin comps. Continuing with list."
            continue
        corrected = correctedCurve(cadenceSum, fittedCurve)
        results_list.append([np.sum(mask),np.std(np.diff(corrected))])
        print np.sum(mask), np.std(np.diff(corrected))
        corrected_list.append(corrected)
        
   
   #plt.figure()
    
   # plt.title("Varying thresh for light curve mask")
    #for corrected in corrected_list:
    #    plt.plot(dataPoints, corrected + n*offset, ".", markersize=4, label = "%s prin comps"%str(numPrinComps))
    #    n += 1
    #plt.savefig("varylcvMask_epic%s_c%s.png"%(str(epic), str(campaign)), clobber=True, dpi=600)
    #plt.show()
    return np.array(results_list)
    
    
def varypcaMask(epic, campaign, plotMask):
    # get the data
    getData_return = getData(epic, campaign)
    data = getData_return[0]
    inds_to_eliminate = getData_return[1].astype(bool)
    
    numPrinComps = 25
    
    # plot the image
    plotCadence(data[0], axis='relative')    

    cadence_med = np.nanmedian(data[0])
    thresh_list = np.arange(0, 300, 25)
        
    
    # create a mask for PCA
    lcv_mask = createMask(data[0], thresh=100, title="Light Curve Mask", plotMask=plotMask)
    lcvMatrix = reduceAperture(lcv_mask, data)
    

    # take out data points with thruster firings
    lcvMatrix = lcvMatrix[:,~inds_to_eliminate]
    
    #get rid of nans in matrices

    lcvMatrix = noMoreNaN(lcvMatrix)
    assert checkNaN(lcvMatrix)
    


    n = 0
    offset = 5000
    cadenceSum = np.sum(lcvMatrix, axis=0)

    n = 0
    corrected_list = []
    print "\n#pixels_pca\tdiff_std"
    print "-----------------"
    for thresh in thresh_list:
        mask = createMask(data[0], thresh=thresh, title="PCA Mask, thresh=%s"%str(thresh), plotMask=plotMask)
        pcaMatrix = reduceAperture(mask, data)
        pcaMatrix = pcaMatrix[:,~inds_to_eliminate]
        pcaMatrix = noMoreNaN(pcaMatrix)
        assert checkNaN(pcaMatrix)
        fittedCurve = curveFit(pcaMatrix, cadenceSum, numPrinComps)[0]
        corrected = correctedCurve(cadenceSum, fittedCurve)
        print np.sum(mask), np.std(np.diff(corrected))
        corrected_list.append(corrected)
        
    dataPoints = range(len(lcvMatrix[0]))
    plt.figure()
    plt.title("Varying thresh for PCA mask")
    for corrected in corrected_list:
        plt.plot(dataPoints, corrected + n*offset, ".", markersize=4, label = "%s prin comps"%str(numPrinComps))
        n += 1
    plt.savefig("varypcaMask_epic%s_c%s.png"%(str(epic), str(campaign)), clobber=True, dpi=600)
    plt.show()
    
def varyNumPieces(epic, campaign, plotMask):
    
    cadenceSum, lcvMatrix, pcaMatrix, lcvMask, pcaMask = doBeforeAnything(epic, campaign, plotMask)
    
    numPrinComps = 10

    numPieces_list = range(1,11)
    plt.figure()
    plt.title("Varying the number of pieces light curve was broken into")
    n = 0
    offset = 5000
    print "\n#Pieces\tdiff_std"
    print "-------\t--------"
    for numPieces in numPieces_list:
        lcv_pieces = np.array_split(lcvMatrix, numPieces, axis=1)
        pca_pieces = np.array_split(pcaMatrix, numPieces, axis=1)
        fitted_pieces = []
        for lcv_part, pca_part in zip(lcv_pieces, pca_pieces):
            cadenceSum = np.sum(lcv_part, axis=0)
            fittedCurve = curveFit(pca_part, cadenceSum, numPrinComps)[0]
            corrected = correctedCurve(cadenceSum, fittedCurve)
            #corrected = normalize([corrected])[0]
            fitted_pieces.append(corrected)
        corrected_lcv = np.concatenate(fitted_pieces)
        diff_std = np.std(np.diff(corrected_lcv))
        print numPieces, "\t", diff_std
        
        plt.plot(range(len(corrected_lcv)), corrected_lcv + n*offset, ".",markersize=4)
        n += 1
        
    plt.savefig("varyNumPieces_epic%s_c%s.png"%(str(epic), str(campaign)), clobber=True, dpi=600)
    plt.show()
            
def getStarz(txtFile):
    epic_list = []
    campaign_list = []
    with open(txtFile) as file:
        for line in file:
            line = line.split("\t")
            epic_list.append(line[0])
            campaign_list.append(line[1])
    return epic_list, campaign_list

def makeMeanZero(corrected):
    """input: 1d light curve array 
        output: 1d array scaled so the mean is zero"""
    corrected /= np.mean(corrected)
    corrected -= 1
    return corrected

def getNoiseParams(lcv):
    # sigma clip it
    sigmaClip_tf = noise.sigmaClip(lcv, 5.)
    lcv = lcv[~sigmaClip_tf]
    
    # get params
    rta = noise.computeRollTweakAmplitude(lcv)
    sgcdpp = noise.computeSgCdpp_ppm(lcv)
    scatter = noise.estimateScatterWithMarshallMethod(lcv)
    
    return rta, sgcdpp, scatter
    
def paramSubplot(cadenceSum, numPrinCompList, results_list):
    y = cadenceSum/np.mean(cadenceSum) - 1
    raw_rta = noise.computeRollTweakAmplitude(y)  
    raw_sgcdpp = noise.computeSgCdpp_ppm(noise.sigmaClip(y,5.))
    raw_scatter = noise.estimateScatterWithMarshallMethod(y)
    f, axarr = plt.subplots(3, sharex=True)    
    axarr[0].axhline(raw_rta, label="RTA for raw light curve")
    axarr[1].axhline(raw_sgcdpp, label="sgcdpp for raw light curve")
    axarr[2].axhline(raw_scatter, label="scatter for raw light curve")
    axarr[0].plot(numPrinCompList, results_list.T[1], ".", color="green")
    axarr[1].plot(numPrinCompList, results_list.T[2], ".", color="green")
    axarr[2].plot(numPrinCompList, results_list.T[3], ".", color="green")
    plt.xlabel("Number of Principal Components")
    axarr[0].set_ylabel("Roll Tweak Amplitude")
    axarr[1].set_ylabel("sgcdpp")
    axarr[2].set_ylabel("scatter")
    plt.show()
    
def plotOptimalPCforParam(n, results_list, t, cadenceSum, fittedCurve, numPrinCompList, pcaMatrix):
    if n == 1:
        param = "rta"
    elif n == 2:
        param = "sgcdpp"
    elif n == 3:
        param = "MM scatter"
    
    plt.figure()
    # plot lowest value of param on top of raw lcv
    plt.plot(t, cadenceSum, ".")
    optimalPC = numPrinCompList[np.argmin(results_list.T[n])]
    fittedCurve = curveFit(pcaMatrix, cadenceSum, optimalPC)[0]
    optimal_lcv = correctedCurve(cadenceSum, fittedCurve)
    plt.plot(t, optimal_lcv, ".")
    plt.title("lowest %s, PC=%i"%(param, optimalPC))
    plt.show()
    
    # plot lowest slope between params on top of raw lcv
    plt.figure()
    smallest_slope = np.argmin(np.abs(np.diff(results_list.T[n])))
    plt.plot(t, cadenceSum, ".")
    optimalPC = numPrinCompList[smallest_slope]
    fittedCurve = curveFit(pcaMatrix, cadenceSum, optimalPC)[0]
    optimal_lcv = correctedCurve(cadenceSum, fittedCurve)
    plt.plot(t, optimal_lcv, ".")
    plt.title("smallest slope %s, PC=%i"%(param, optimalPC))
    plt.show()

def foldLightcurve(t, period, lcv):
    phase = np.fmod(t, period)
    plt.figure()
    plt.plot(phase, lcv)
    plt.show()
    
    
def varyPrinComp(epic, campaign, plotMask):
    
    cadenceSum, lcvMatrix, pcaMatrix, lcvMask, pcaMask = doBeforeAnything(epic, campaign, plotMask)
    t = np.arange(len(cadenceSum))
    
    # plot the raw light curve
    make_plot(t,  cadenceSum,show=False, title="Raw Light Curve, e: %s, c: %s"%(epic, campaign))
    
    # make a list of number of principal comps to try
    numPrinCompList = np.arange(2,np.sum(lcvMask), 2)


    results_list = []
    n = 0
    
    #plt.figure()
    for numPrinComps in numPrinCompList:
        curveFit_return = curveFit(pcaMatrix, cadenceSum, numPrinComps)
        fittedCurve = curveFit_return[0]
        
        # correct curve for PC's
        corrected = correctedCurve(cadenceSum, fittedCurve)
        
        # make the mean zero
        correctedMeanZero = makeMeanZero(corrected)
        
        # get rta, sgcdpp, MM scatter parameters:
        rta, sgcdpp, scatter = getNoiseParams(correctedMeanZero)

        results_list.append([numPrinComps,rta, sgcdpp, scatter])
        #plt.plot(range(len(corrected)), corrected+ n*offset, ".", markersize=4, label = "%f,%i"%(sgcdpp, int(numPrinComps)))
        n += 1
        
    #plt.legend()
    #plt.show()
    
    results_list = np.array(results_list)
    
    # plot the three params wrt the raw values
    if plotMask==True:
        paramSubplot(cadenceSum, numPrinCompList, results_list)
    
    # plot results
    n = 1
    numParams = 3
    while n < numParams + 1:
        plotOptimalPCforParam(n,results_list, t, cadenceSum, fittedCurve, numPrinCompList, pcaMatrix)
        n += 1

    return np.array(results_list).T
    
    

def plotCumSumPCA(epic, campaign,n):
    dimReturnCutoff_diff = 0.01
    dimReturnCutoff_dec = 0.99
    method = "dec"
    cadenceSum, lcvMatrix, pcaMatrix = doBeforeAnything(epic, campaign, False)[:3]
    t = np.arange(len(cadenceSum))

    lcvMatrixNorm = normalize(lcvMatrix)
    U, s, V = getPrinComps(lcvMatrixNorm)
    prinComps = V
    numPrinComps = np.arange(len(prinComps))
    evalues = np.square(s)
    cumSum = np.cumsum(evalues)
    totSum = np.sum(evalues)
    diff =  np.diff(cumSum/totSum)
    if method == "diff":
        for element in diff:
            if element < dimReturnCutoff_diff:
                idx = np.where(diff == element)[0][0]
                break
    elif method == "dec":
        for element in cumSum/totSum:
            if element > dimReturnCutoff_dec:
                idx = np.where((cumSum/totSum) == element)[0][0]
                break
    plt.plot(numPrinComps,cumSum/totSum, label="e:%i"%epic)
    plt.legend(loc=4)
    plt.show()
    optimalPC= idx + 1
    
    # fit the curve
    curveFit_return = curveFit(pcaMatrix, cadenceSum, optimalPC)
    fittedCurve = curveFit_return[0]
        
    # correct curve for PC's
    corrected = correctedCurve(cadenceSum, fittedCurve)
    plt.figure()
    plt.title("e:%i PC=%i"%(epic, optimalPC))
    plt.plot(t, cadenceSum, ".", label = "raw")
    plt.plot(t, corrected, ".", label = "corrected")
    plt.legend()
    plt.show()
        

def cumSumPCA():
    txtFile = "k2goodTransits.txt"
    #txtFile = "k2_targets_sorted.txt"
    getList = getStarz(txtFile)
    epic_list = getList[0]
    campaign_list = getList[1]
    
    print len(epic_list)
    n = 0
    plt.figure()
    plt.title("Cumulative Sum of Eigenvalues")
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Sum of Eigenvalues / Total Sum")
    plt.ylim(0.5, 1.1)
    plt.xlim(0,20)
    for epic, campaign in zip(epic_list[12:22], campaign_list[12:22]):
        print n
        plotCumSumPCA(int(epic), int(campaign), n)
        n += 1

        
        
        
def main():
    
    # EPICS AND CAMPAIGNS TO TRY
#==============================================================================
    txtFile = "k2goodTransits.txt"
    txtFile = "k2_targets_sorted.txt"
    getList = getStarz(txtFile)

    epic_list = [206103150, 211351816]
    epic_list = np.concatenate((epic_list, getList[0]))
    campaign_list = [3,5]
    campaign_list= np.concatenate((campaign_list, getList[1]))
    # WASP-47 K2 
    #epic = 206103150
    #campaign = 3
    # Red Giant
    #epic = 211351816
    #campaign = 5
    
#==============================================================================
#     epicsGoodTransits = [206103150,205996447,205947214,210789323,
#                          210744674,211418729,211399359,212351405,
#                          201182911,201488265,202843107,]
#     campaignsGoodTransits = [3,3,3,4,4,5,5,6,1,]
#==============================================================================
#==============================================================================    
    plotMask=False
    
    results_list = []
    n = 1
    #for epic, campaign in zip(epic_list, campaign_list):
    for epic, campaign in zip(epic_list, campaign_list):
        print "\nCOUNTER: %s"%str(n)
        print "epic, campaign = ",epic, campaign
        results = varyPrinComp(int(epic), int(campaign), plotMask)
        results_list.append(results)
        answer = raw_input("Finished with epic %s. Want to move on to the next one (y or n)? "%epic)
        n += 1
        if answer == "n":
            break
        else:
            continue

#main()
cumSumPCA()
    
